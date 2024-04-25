library('DESeq2')
library('ggplot2')
library('pheatmap')
library('reshape2')
library('stringr')
library('openxlsx')


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
install.packages("ggplot2")

#### Section 1: Input Data  ####-------------------------

# ---- 1.1: set file names ----
# read in the file names from the directory based on a pattern 
# match. Only the ones -----ending with .csv should be kept. 
src_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/star/'
out_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/deseq2/'
metadata_file = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/metadata.xlsx'

# select the ouput files of STAR
sample_files <- list.files( path=src_dir, pattern = "*ReadsPerGene.out.tab$", full.names = TRUE )
sample_files

# read in the data, and store them in a list, one entry 
# per sample file. The column number is an important parameters
# as it reflects the overall count, sense count or antisense
# count. 
counts.files <- lapply(sample_files, read.table, skip = 4 )
col_number = 2
counts <- as.data.frame( sapply( counts.files, function(x) x[ , col_number ] ) )

# set rownames of the count matrix
row.names(counts) <- counts.files[[1]]$V1

# create column names
sample_names = gsub( "_ReadsPerGene[.]out[.]tab", "", sample_files )
sample_names = gsub( paste0(src_dir, "/"), "", sample_names)
sample_names
colnames(counts) = sample_names
head(counts)


# ---- 1.2: create sampleTable ----
## define the sample name and the biological condition as extracted 
## from the meta data excel file
coldata = read.xlsx(metadata_file)
coldata$new_column = paste0(coldata$tissue, "_", coldata$Condition)
coldata$new_column2 = paste0(coldata$tissue, "_", coldata$Condition, "_", coldata$Timepoint)

coldata

# ---- 1.3: create DESeq object with design ----

## preprocessing the counts
counts = counts[,-grep("BS", colnames(counts))]
counts = counts[,-grep("MC", colnames(counts))]
coldata = coldata[-grep("BS", coldata$sample_name),]
coldata = coldata[-grep("MC", coldata$sample_name),]

colnames(counts) == coldata$sample_name

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design=  ~ tissue)
dds$condition = relevel(dds$tissue, ref = "ear")
dds = DESeq(dds)

# the following selection is done to reduce the number of genes
# with a low expression level.  
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


## make correlation heatmap without the outlier
index_to_remove = which(rownames(sampleDistMatrix) == "E96I2")
sampleDistMatrix = sampleDistMatrix[-index_to_remove, -index_to_remove]
rownames(sampleDistMatrix)  = colnames(vsd)[-index_to_remove]
colnames(sampleDistMatrix)  = colnames(vsd)[-index_to_remove]
sampleDists <- dist(t(assay(vsd)[,-index_to_remove]))


p = pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             fontsize=12) #,
#         col=colors)
p

# make a PCA plot with all 10 samples
p = plotPCA(vsd, intgroup="condition")
# p = p + geom_point(size=10, alpha=0.50)
p = p + theme_bw()
p = p + theme(text = element_text(size = 18)) 
p = p + geom_text(aes(label=colnames(vsd)), hjust=-0.25, vjust=-0.10, size=4)
# p = p + coord_cartesian(xlim=c(-25,30))
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
         show_rownames=TRUE,
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


# ---- 2.5 densityplots per sample ----
sample_counts = as.data.frame(colSums(counts))
colnames(sample_counts)[1] = 'counts'
sample_counts$tissue = coldata[match(rownames(sample_counts), coldata$sample_name),]$tissue

ggplot(data = sample_counts, aes(x=counts)) +
  geom_density(aes(color=tissue), alpha=.50) + 
  theme_bw()

ggplot(data = sample_counts, aes(y=counts)) +
  geom_boxplot(aes(group=tissue, fill=tissue), alpha=.50) + 
  theme_bw()

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

