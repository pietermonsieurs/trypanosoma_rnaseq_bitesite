library('DESeq2')
library('ggplot2')
library("RColorBrewer")
library('pheatmap')
library('vsn')
library('VennDiagram')
library('reshape2')
library('stringr')
# library("countToFPKM")
# library("apeglm")

setwd('/Users/pmonsieurs/programming/trypanosoma_rnaseq/results/deseq2/')

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
# BiocManager::install("ComplexHeatmap")
# BiocManager::install("apeglm")

#### Section 1: Input Data  ####-------------------------

# ---- 1.1: set file names ----
# read in the file names from the directory based on a pattern 
# match. Only the ones -----ending with .csv should be kept. 
src_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq/results/star/'
out_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq/results/deseq2/'
gff_file = '/Users/pmonsieurs/programming/trypanosoma_rnaseq/data/TriTrypDB-48_TbruceiTREU927.gff'
setwd(src_dir)

# select the ouput files of STAR
sample_files <- list.files( path=src_dir, pattern = "*ReadsPerGene.out.tab$", full.names = TRUE )
sample_files

# read in the data, and store them in a list, one entry 
# per sample file. The column number is an important parameters
# as it reflects the overall count, sense count or antisense
# count. 
counts.files <- lapply(sample_files, read.table, skip = 4 )
col_number = 3
counts <- as.data.frame( sapply( counts.files, function(x) x[ , col_number ] ) )

# set rownames of the count matrix
row.names(counts) <- counts.files[[1]]$V1

# create column names
sample_names = gsub( "[.]ReadsPerGene[.]out[.]tab", "", sample_files )
sample_names = gsub( paste0(src_dir, "/"), "", sample_names)
sample_names
colnames(counts) = sample_names
head(counts)

# find the samples split over multiple sequencing runs, and sum the 
# corresponding columns. 001 = 1 + 4 // 005 = 2 + 8 // 009 = 3 + 12 
grep("R_", colnames(counts))
colnames(counts)
counts[,4] = counts[,4] + counts[,1]
counts[,8] = counts[,8] + counts[,2]
counts[,12] = counts[,12] + counts[,3]

# remove the first 3 columsn, as they are now summed up to the other
# corresponding replicates
counts = counts[,-seq(1,3)]
column_names = c(paste0("Metac", seq(1:5)), paste0("Blood", seq(1:5)))
colnames(counts) = column_names
head(counts)


# ---- 1.2: create sampleTable ----
# define the sample name and the biological condition and use this 
# information to build a data frame as input for the DESeq2 object
coldata = as.data.frame(colnames(counts))
conditions = c(rep("Metac", 5), rep("Blood", 5))
coldata$condition = as.factor(conditions)
colnames(coldata) = c('sample', 'condition')
coldata

# ---- 1.3: create DESeq object with design ----
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design= ~ condition)
dds$condition = relevel(dds$condition, ref = "Blood")
dds = DESeq(dds)

# the following selection is done to reduce the number of genes
# with a low expression level. Has *no* effect in the case of quiescence
# experiment as no genes have less than 100 reads over all samples. 
keep <- rowSums(counts(dds)) >= 500
sum(keep)
dds <- dds[keep,]

# do differential expression between both conditions
res = results(dds)




cutoff_fc = 1
cutoff_pval = 0.05
### up ###
sum(res$log2FoldChange > cutoff_fc & res$padj < cutoff_pval, na.rm=TRUE)

### down ###
sum(res$log2FoldChange < -cutoff_fc & res$padj < cutoff_pval, na.rm=TRUE)


# ---- 1.4: read and parse GFF file ----
# add gene lengths to calculate FPKM values 
gff = read.csv(gff_file, sep="\t", comment.char="#", header=FALSE)
gff = gff[gff$V3 == "gene",]
field9 = str_split_fixed(gff[,9], ";",n=Inf)
gene_id = gsub("ID=", "", field9[,1])
description = gsub("description=", "", field9[,2])

gff = gff[,c(1,3,4,5,7)]
colnames(gff) = c("chrom", "feature", "start", "stop", "strand")
gff$length = gff$stop - gff$start + 1

gff$gene_id = gene_id
gff$description = description

rownames(gff) = gff$gene_id
head(gff)






# match between geneid and dds features, to add lengths to the
# DESeq object to calculate FPKM
index = match(rownames(dds), gff$geneid)
mcols(dds)$basepairs = gff$length[index]


# first run section 4.1 --> create the gff_ribo object only 
# containing ribosomal RNA

dds = dds[! rownames(dds) %in% gff_ribo$gene,]




#### Section 2: Quality Control #### 

setwd('/Users/pmonsieurs/programming/trypanosoma_rnaseq/results/deseq2/')


# ---- 2.1 conversion for QC ---- 
# first make some conversion that make it easier to do some 
# quality control plots
ntd <- normTransform(dds)
vsd <- vst(dds, blind=TRUE)
rld <- rlog(dds, blind=TRUE)
head(assay(vsd), 3)



# ---- 2.2 PCA / correlation matrix ----
sampleDists <- dist(t(assay(vsd)))

# correlation plot shown as heatmap
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(vsd))
colnames(sampleDistMatrix) <- NULL
#colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
p = pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             fontsize=12) #,
#         col=colors)
p
png_file_cor = paste0(out_dir, 'correlation_replicates.png')
ggsave(png_file_cor, p)


# make a PCA plot with all 10 samples
p = plotPCA(vsd, intgroup="condition")
p = p + geom_point(size=10, alpha=0.50)
p = p + theme_bw()
p = p + theme(text = element_text(size = 18)) 
p = p + geom_text(aes(label=colnames(vsd)), hjust=-0.25, vjust=-0.10, size=4)
p = p + coord_cartesian(xlim=c(-25,30))
p

png_file_pca = paste0(out_dir, 'PCA_replicates.png')
ggsave(png_file_pca, p)


# heatmap of count matrix - first do a selection of the genes
# with the highest expression
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]
df <- as.data.frame(colData(dds)$condition)
rownames(df) = colnames(dds)
#pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#         cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(ntd)[select,], 
         cluster_cols=TRUE, 
         cluster_rows=FALSE,
         show_rownames=FALSE,
         annotation_col = df
)



# ---- 2.3 MA plot ----
DESeq2::plotMA(res)





# ---- 2.4: density plots (including FPKM) ----
# density plots of the counts (raw and normalized)
density_raw = melt(log10(counts(dds)))
colnames(density_raw) = c('genename', 'condition', 'count')
p = ggplot(density_raw, aes(x=count, color=condition))
# p = p + geom_density()
p = p + stat_density(geom="line",position="identity", size=1)
p = p + xlim(c(0,5))
p = p + ggtitle("Raw read count per gene")
p = p + labs(x="log 10 read count per gene")
p = p + theme_bw()
p = p + theme(axis.text=element_text(size=12),
              axis.title=element_text(size=16,face="bold"),
              title=element_text(size=16))
p

dens_plot_raw_file = paste0(out_dir, 'density_plot_raw.png')
ggsave(dens_plot_raw_file, p)

density_raw = melt(log10(counts(dds, normalized=TRUE)))
colnames(density_raw) = c('genename', 'condition', 'count')
p = ggplot(density_raw, aes(x=count, color=condition))
p = p + stat_density(geom="line",position="identity", size=1)
p = p + xlim(c(0,5))
p = p + theme_bw()
p = p + ggtitle("Normalized read count per gene")
p = p + labs(x="log 10 read count per gene")
p = p + theme(axis.text=element_text(size=12),
              axis.title=element_text(size=16,face="bold"),
              title=element_text(size=16))
p

dens_plot_norm_file = paste0(out_dir, 'density_plot_normalized.png')
ggsave(dens_plot_norm_file, p)

# FPKM values and density values
dds_fpkm = fpkm(dds)
fpkm_melted = melt(dds_fpkm)
colnames(fpkm_melted) = c('genename', 'condition', 'fpkm')

p = ggplot(fpkm_melted, aes(x=fpkm, color=condition))
p = p + geom_density()
p = p + xlim(c(0,400))
p = p + theme_bw()
p

dens_plot_norm_file = paste0(out_dir, 'density_plot_fpkm.png')
ggsave(dens_plot_norm_file, p)




#### Section 3: Analyse DE genes ####

# ---- 3.1 export the data to a .csv file ----
res_ann = as.data.frame(res)
head(res_ann)
gff_index = match(rownames(res_ann), rownames(gff))
res_ann <- cbind(res_ann, gff[gff_index,])
# write.csv(as.data.frame(res), file='blood_vs_metacyclic.csv')
write.csv(res_ann, file='blood_vs_metacyclic.csv')

# ---- 3.2. analyse DE genes ----
number_of_DE_genes = 20
# topGene <- rownames(res)[which.min(res$padj)]

# based on selection / sorting on FDR value. However, this only results
# in overexpressed genes
topGenes = rownames(res)[order(res$padj)][1:number_of_DE_genes]

# same approach, but first filtering for downregulated genes, and 
# then sort on p-value / FDR
res_down = res[res$log2FoldChange < -1,]
topGenes = rownames(res_down)[order(res_down$padj)][1:number_of_DE_genes]

figure_list = c()
i = 0
for (topGene in topGenes) {
  i = i+1
  data <- plotCounts(dds, gene=topGene, intgroup=c("condition"), returnData=TRUE)

  p = ggplot(data, aes(x=condition, y=count, fill=condition))
  p = p + scale_y_log10() 
  p = p + geom_dotplot(binaxis="y", stackdir="center", dotsize = 2)
  p = p + ggtitle(topGene)
  p = p + theme(legend.position = "none")
  p = p + theme(plot.title = element_text(size=10))
  
  
  figure_list[[i]] = local({
      print(p)  
    })  
}

source("http://peterhaschke.com/Code/multiplot.R") #load multiplot function
multiplot(plotlist = figure_list, cols = 5)

library('ggpubr')
figure <- ggarrange(figure_list,
                    labels = c("A", "B", "C"),
                    ncol = 2, nrow = 2)



# ---- 3.3. export DE genes for GO analysis ----

# remove the results with a NA value in the adjusted p-value. 
res_filter = res[-which(is.na(res$padj)),]

# do filter with a log2 fold change of 1
cutoff_fc = 2
cutoff_pval = 0.05

### up ###
de_up = rownames(res_filter[res_filter$log2FoldChange > cutoff_fc & res_filter$padj < cutoff_pval,])
out_file = paste0('DE_up_fc', cutoff_fc, '.csv')
write.table(de_up, file=out_file, row.names = FALSE, quote=FALSE, col.names=FALSE)

### down ###
de_down= rownames(res_filter[res_filter$log2FoldChange < -cutoff_fc & res_filter$padj < cutoff_pval,])
out_file = paste0('DE_down_fc', cutoff_fc, '.csv')
write.table(de_down, file=out_file, row.names = FALSE, quote=FALSE, col.names=FALSE)



#### 4. specialized analyses on demand ####

# ---- 4.1 count plot per replicate ----

# specify the output directory where the figures should be 
# stored
figure_dir = paste0(out_dir, 'figures/')
# gene_list = 'fam10'
gene_list = 'pad'

select_genes_file = paste0(out_dir, gene_list, '_genes.csv')
selectGenes_df = read.csv(select_genes_file, header=FALSE)
selectGenes = selectGenes_df[,1]


figure_list = c()
i = 0
for (selectGene in selectGenes) {
  i = i+1
  data <- plotCounts(dds, gene=selectGene, intgroup=c("condition"), returnData=TRUE)
  
  p = ggplot(data, aes(x=condition, y=count, fill=condition))
  p = p + scale_y_log10() 
  p = p + geom_dotplot(binaxis="y", stackdir="center", dotsize = 2)
  p = p + ggtitle(selectGene)
  p = p + theme(legend.position = "none")
  p = p + theme(plot.title = element_text(size=10))
  p = p + geom_text(aes(label=rownames(data)), hjust=-0.25, vjust=-0.10, size=4)
  

  figure_file = paste0(figure_dir, gene_list, '.', selectGene, '.countplot.png')  
  ggsave(figure_file, plot=p)
  
  figure_list[[i]] = local({
    print(p)  
  })  
}

source("http://peterhaschke.com/Code/multiplot.R") #load multiplot function
multiplot(plotlist = figure_list, cols = 3)


##### OLD CODE ######





# ---- Section 2.3: Sd plots ----

# show dispersion estimates
plotDispEsts(dds)

# make a plot of the variance stabilized values 
meanSdPlot(assay(dds), ranks=FALSE) # variance depending on value
meanSdPlot(assay(ntd), ranks=TRUE) # variance not depending anymore on value
meanSdPlot(assay(vsd))








#### Section 3: Differential expression ####

# ---- Section 3.1: make contrasts ----
# the default command to do differential expression. This will calculate
# some random constrast between two conditions (first versus last?). 
# obsolet... will now be run with default parameters above
# dds = DESeq(dds, betaPrior = FALSE, full=design_matrix)
# dds = DESeq(dds, betaPrior = FALSE)
# res = results(dds)

# find contrasts between two individual factor e.g. prosta and prolog. 
# such a contrast can be found without having to set this condition
# in the initial design matrix above. 
res_prosta_prolog = results(dds, contrast=c("condition_", "prosta", "prolog"))
res_amasta_amalog = results(dds, contrast=c("condition_", "amasta", "amalog"))
res_amaint_amalog = results(dds, contrast=c("condition_", "amaint", "amalog"))

DESeq2::plotMA(res_amasta_amalog, main="amasta vs amalog - housekeeping genes")

head(res_prosta_prolog)

write.csv(res_prosta_prolog, file="res_prosta_prolog.csv")
write.csv(res_amasta_amalog, file="res_amasta_amalog.csv")
write.csv(res_amaint_amalog, file="res_amaint_amalog.csv")


# output of prosta_vs_prolog above should be the same (almost) as the
# prosta_vs_prolog comparision which is integrated into the design 
# matrix. prosta_vs_prolog is the first contrast (of 4) in the design matrix
# (first column), so contrast vector should be c(1,0,0,0). However, this
# approach is **NOT WORKING**
### res_prosta_prolog_design = results(dds, contrast=c(1,0,0,0))
### head(res_prosta_prolog_design)

# with other approach using a design like design = ~ 0 + condition, this 
# gives the same results. For prosta-prolog:
res_prosta_prolog_design = results(dds, contrast=list("condition_prosta", "condition_prolog"))
head(res_prosta_prolog_design)
res_prosta_prolog = res_prosta_prolog_design
hist(res_prosta_prolog$log2FoldChange, breaks=100)
DESeq2::plotMA(res_prosta_prolog_design)

# contrast amasta - amalog
res_amasta_amalog = results(dds, contrast=list("condition_amasta", "condition_amalog"))
head(res_amasta_amalog)
hist(res_amasta_amalog$log2FoldChange, breaks=100)
DESeq2::plotMA(res_amasta_amalog)


# contrast amaint - amalog
res_amaint_amalog = results(dds, contrast=list("condition_amaint", "condition_amalog"))
head(res_amaint_amalog)
hist(res_amaint_amalog$log2FoldChange, breaks=100)
DESeq2::plotMA(res_amaint_amalog)


# test with complex design, where all quiescent states are compared
# with the prolog condition. Does not seem to work very well..
res_quisc_prolog = results(dds, 
                                  contrast=list(c("condition_prosta", "condition_amalog", "condition_amasta"),
                                  c("condition_prolog")),
                                  listValues=c(1/3,-1)
                                 )
head(res_quisc_prolog)
hist(res_quisc_prolog$log2FoldChange, n=100)


# do some filtering on the data, and get some feeling on the
# amount of genes up or down regulated. Produce Venn diagrams
# to show the overlap between different comparisons

# ---- Section 3.2: overlap of DE genes (Venn diagrams) ----

# cutoffs_pval = c(1, 0.01, 1e-5)
# cutoffs_fc = c(1,2)
cutoffs_pval = c(0.05, 0.01)
cutoffs_fc = c(1)


cutoff_fc = 1
cutoff_pval = 0.01
### up ###
sum(res_quisc_prolog$log2FoldChange > cutoff_fc & res_quisc_prolog$padj < cutoff_pval, na.rm=TRUE)
sum(res_amasta_amalog$log2FoldChange > cutoff_fc  & res_amasta_amalog$padj < cutoff_pval, na.rm=TRUE)
sum(res_amaint_amalog$log2FoldChange > cutoff_fc & res_amaint_amalog$padj < cutoff_pval, na.rm=TRUE)
sum(res_prosta_prolog$log2FoldChange > cutoff_fc & res_prosta_prolog$padj < cutoff_pval, na.rm=TRUE)

### down ###
cutoff_fc = -1
sum(res_quisc_prolog$log2FoldChange < cutoff_fc & res_quisc_prolog$padj < cutoff_pval, na.rm=TRUE)
sum(res_amasta_amalog$log2FoldChange < cutoff_fc  & res_amasta_amalog$padj < cutoff_pval, na.rm=TRUE)
sum(res_amaint_amalog$log2FoldChange < cutoff_fc & res_amaint_amalog$padj < cutoff_pval, na.rm=TRUE)
sum(res_prosta_prolog$log2FoldChange < cutoff_fc & res_prosta_prolog$padj < cutoff_pval, na.rm=TRUE)


for (cutoff_pval in cutoffs_pval) {
  for (cutoff_fc in cutoffs_fc) {  
    # cutoff_fc = 1
    # cutoff_pval = 1
    # index_quisc_prolog = which(res_quisc_prolog$log2FoldChange > cutoff_fc & res_quisc_prolog$padj < cutoff_pval)
    index_amasta_amalog = which(res_amasta_amalog$log2FoldChange > cutoff_fc & res_amasta_amalog$padj < cutoff_pval)
    index_amaint_amalog = which(res_amaint_amalog$log2FoldChange > cutoff_fc & res_amaint_amalog$padj < cutoff_pval)
    index_prosta_prolog = which(res_prosta_prolog$log2FoldChange > cutoff_fc & res_prosta_prolog$padj < cutoff_pval)
    
    # VENN diagrams #
    # Chart
    myCol <- brewer.pal(4, "Set1")
    myCol <- brewer.pal(3, "Set1")
    
    venn.diagram(
      x = list(# index_quisc_prolog, 
               index_amasta_amalog,
               index_amaint_amalog,
               index_prosta_prolog),
      category.names = c(# "Quiscence vs prolog",
                         "amasta vs amalog",
                         "amaint vs amalog",
                         "prosta vs amalog"),
      filename = paste0(out_dir, 
                        'venn_diagramm',
                        #'.two_conditions',
                        '.v2',
                        '.pvalue_', cutoff_pval, 
                        '.fc_', cutoff_fc, 
                        '.png'),
      
      # add general title 
      main = paste0('Venn p-value ', cutoff_pval, ' -- logFC ', cutoff_fc),
      main.cex = 0.4,
      main.fontface = "bold",
      # filename = NULL,
      output=TRUE,
      
      # # Output features
      imagetype="png" ,
      height = 480 , 
      width = 480 , 
      resolution = 300,
      compression = "lzw",
      # 
      # # Circles
      lwd = 2,
      lty = 'blank',
      # lty = 1,
      fill = myCol,
      # fill=myCol[1:2],
      # 
      # # Numbers
      cex = .3,
      #fontface = "regular",
      fontfamily = "sans",
      # 
      # # Set names
      cat.cex = 0.4,
      cat.fontface = "bold",
      # cat.default.pos = "outer",
      # cat.pos = c(-27, 27, 135),
      # cat.dist = c(0.055, 0.055, 0.085),
      # cat.fontfamily = "sans",
      # rotation = 1
      )
  }
}


# ---- Section 3.3: detailed ananlysis of subset ----

# create heatmap for overlap of 8 genes, or select the overlapping gene
# between all genes
cutoff_pval = 0.05
cutoff_fc = 1

index_quisc_prolog = which(res_quisc_prolog$log2FoldChange > cutoff_fc & res_quisc_prolog$padj < cutoff_pval)
index_amasta_amalog = which(res_amasta_amalog$log2FoldChange > cutoff_fc & res_amasta_amalog$padj < cutoff_pval)
index_amaint_amalog = which(res_amaint_amalog$log2FoldChange > cutoff_fc & res_amaint_amalog$padj < cutoff_pval)
index_prosta_prolog = which(res_prosta_prolog$log2FoldChange > cutoff_fc & res_prosta_prolog$padj < cutoff_pval)

index_intersection = intersect(index_amasta_amalog, index_prosta_prolog)
index_intersection = intersect(index_intersection, index_amaint_amalog)



png(paste0(out_dir, "heatmap_cutoff", cutoff_pval, '_cutofflogfc', cutoff_fc, ".png"))
data_plot = counts(dds)[index_intersection,]
data_plot = log10(counts(dds)[index_intersection,])
data_plot[is.na(data_plot)] = 1
data_plot[is.infinite(data_plot)] = 1
pheatmap(data_plot, 
         cluster_cols=TRUE, 
         cluster_rows=TRUE,
         show_rownames=TRUE
         #col=rev(heat.colors(n=100)),
         #breaks=seq(1.5,3.5,by=0.02)
         # annotation_col = counts(dds)[index_intersection,]
)
dev.off()


### zoom in on single gene
ii = index_intersection
raw_counts_ii = as.data.frame(counts(dds)[ii,])
norm_counts_ii = as.data.frame(counts(dds, normalized=TRUE)[ii,])
fpkm_ii = as.data.frame(dds_fpkm[ii,])


df_selected = raw_counts_ii
colnames(df_selected) = c('raw_counts')
df_selected$sample = rownames(df_selected)
df_selected$condition = c(rep('amaint',3), rep('amalog', 3), rep('amasta', 3), rep('prolog',3), rep('prosta', 3))
p = ggplot(data=df_selected, aes(x=sample, y=raw_counts, fill=condition))
p = p + geom_bar(stat="identity")
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

bar_plot = paste0(out_dir, 'barplot_rawcounts.png')
ggsave(bar_plot, p)


df_selected = norm_counts_ii
colnames(df_selected) = c('norm_counts')
df_selected$sample = rownames(df_selected)
df_selected$condition = c(rep('amaint',3), rep('amalog', 3), rep('amasta', 3), rep('prolog',3), rep('prosta', 3))
p = ggplot(data=df_selected, aes(x=sample, y=norm_counts, fill=condition))
p = p + geom_bar(stat="identity")
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

bar_plot = paste0(out_dir, 'barplot_norm_counts.png')
ggsave(bar_plot, p)
                  



# ---- Section 3.4: MA- and volcano plots ----


volcano_amasta_amalog = paste0(out_dir, 'volcano.amasta_vs_amalog.png')
make_volcano(res_amasta_amalog, volcano_amasta_amalog)

volcano_prosta_prolog = paste0(out_dir, 'volcano.prosta_vs_prolog.png')
make_volcano(res_prosta_prolog, volcano_prosta_prolog)


make_volcano = function (res, png_file) {
  png(png_file, width=600, height=600)
  alpha <- 0.05 # Threshold on the adjusted p-value
  cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
  plot(res$log2FoldChange, -log10(res$padj), col=cols, panel.first=grid(),
       main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
       pch=20, cex=0.6)
  abline(v=0)
  abline(v=c(-1,1), col="brown")
  abline(h=-log10(alpha), col="brown")
  
  gn.selected <- abs(res$log2FoldChange) > 2.5 & res$padj < alpha 
  text(res$log2FoldChange[gn.selected],
       -log10(res$padj)[gn.selected],
       lab=rownames(res)[gn.selected ], cex=0.4)
  dev.off()
}


#### Section 4: RNA mapping ####

# ---- Section 4.1: heatmap of ribosomal genes ----

# select ribosomal genes based on the description in the gff-file. This
# might include some artefacts, but most of them are correct. 
# First parse out the geneID, as we need this to link to the read 
# count matrix
gff = read.csv(gff_file, sep="\t", comment.char="#", header=FALSE)
gene_id = str_split_fixed(gff[,9], ";", n=Inf)[,2]
gene_id = gsub(" gene_id ", "", gene_id)
gff$gene_id = gene_id 

# add an extra column to indicate whether ribo or not
gff$ribo = 0
gff[grep("ribosomal", gff[,9]),]$ribo = 1
gff[grep("Ribosomal", gff[,9]),]$ribo = 1
ribo_genes = gff[gff$ribo==1, ]$gene_id

# create gff_ribo for backward compatibility with previous code below
gff_ribo = gff[gff$ribo==1,]

# heatmap based on raw counts
counts_raw = counts(dds, normalized=FALSE)
counts_raw = counts_raw[rownames(counts_raw) %in% gff_ribo$gene,]
ribo_heatmap_raw = paste0(out_dir, 'heatmap_ribo_rawcounts.png')
make_pheatmap(counts_raw, ribo_heatmap_raw)

# normalized counts
counts_norm = counts(dds, normalized=TRUE)
counts_norm = counts_norm[rownames(counts_norm) %in% gff_ribo$gene,]
ribo_heatmap_norm = paste0(out_dir, 'heatmap_ribo_normcounts.png')
make_pheatmap(counts_norm, ribo_heatmap_norm)

# FPKM values
# this next line is broken. still needs to be fixed
dds_fpkm = fpkm(dds)
dds_fpkm_ribo = dds_fpkm[rownames(dds_fpkm) %in% gff_ribo$gene,]
ribo_heatmap_fpkm = paste0(out_dir, 'heatmap_ribo_fpkm.png')
make_pheatmap(dds_fpkm_ribo, ribo_heatmap_fpkm)

    
# venn diagrams based on the FPKM values where the count should
# be low in prolog condition and high in amasta and prosta (and
# maybe also )

# ---- Section 4.2: percentages ----
counts_all = counts(dds)
ribo_counts = colSums(counts_all[rownames(counts_all) %in% ribo_genes,])
total_counts = colSums(counts_all)
percentages = data.frame(100*ribo_counts/total_counts)

mrna_levels = 100 - percentages
percentages$condition = rownames(percentages)
mrna_levels$condition = rownames(mrna_levels)
mrna_levels$type = 'mRNA'
percentages$type = 'rRNA'
df_barplot = rbind(mrna_levels, percentages)
colnames(df_barplot)[1] = 'percentage'

p = ggplot(df_barplot, aes(fill=type, y=percentage, x=condition)) 
p = p + geom_bar(position="stack", stat="identity", alpha=0.80)
p = p + theme_bw()
p = p + theme(text = element_text(size=18), axis.text.x = element_text(angle = 90, hjust = 1))
p

# ---- Section 4.3: density plots ----
density_raw = melt(counts_all)
colnames(density_raw) = c('gene_id', 'condition', 'count')
density_raw$ribo = 0
density_raw$ribo = factor(density_raw$ribo)
density_raw[density_raw$gene_id %in% ribo_genes, ]$ribo = 1

# p = ggplot(density_raw, aes(x=count, color=condition, linetype=ribo))
p = ggplot(density_raw, aes(x=count, color=ribo, fill = ribo))
p = p + geom_density(alpha=0.3)
p = p + xlim(c(0,25000))
p = p + theme_bw()
p

# ---- Section 4.4: find top expressed rRNA genes ----
head(counts_all)
counts_ribo = counts_all[rownames(counts_all) %in% ribo_genes,]
counts_ribo_average_sort = as.data.frame(rowMeans(counts_ribo[order(rowMeans(counts_ribo)), ]))
for (i in seq(1,20)) {
  gene_id = rownames(counts_ribo_average_sort)[i]
  description = as.character(gff[gff$gene_id == gene_id,]$V9)
  data = strsplit(description, ";")
  print(data[[1]][3])
}


#### Section X: functions ###

make_pheatmap = function (counts, png_file) {
  png(png_file)
  dev.set(which=2)
  pheatmap(counts, 
           cluster_cols=TRUE, 
           cluster_rows=TRUE,
           show_rownames=FALSE
  )
  dev.copy(which=4)
  dev.off()
}








