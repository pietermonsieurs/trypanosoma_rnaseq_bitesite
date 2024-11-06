library('Mfuzz')
library('reshape')
library('pheatmap')
library('openxlsx')
library('ggplot2')
library(gridExtra)
library(DESeq2)

## GSEA enrichment libraries
library(clusterProfiler)
library(org.Mm.eg.db)  # Replace with the appropriate organism database
library(AnnotationDbi)
library(enrichplot)

#### 1. process input data ####

setwd('/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/star/')

## specify whether you run the analysis for ears or for lymph nodes, in order
## to do analysis separately for both approaches
tissue = 'ear'
tissue = 'lymphnode'

# ---- ⊢ 1.1: set file names ----
# read in the file names from the directory based on a pattern 
# match. Only the ones -----ending with .csv should be kept. 
src_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/star/'
out_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/deseq2/'
metadata_file = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/metadata.xlsx'
# excel_out = paste0(out_dir, 'rnaseq_immune_response.xlsx')
setwd(src_dir)

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

## save the source counts data frame
counts_src = counts



# ---- ⊢ 1.2 create sampleTable ----
## define the sample name and the biological condition as extracted 
## from the meta data excel file
coldata = read.xlsx(metadata_file)
coldata$Timepoint = as.factor(coldata$Timepoint)
coldata

# ---- ⊢ 1.3 create DESeq object with design ----

## preprocessing the counts. Remove the pure Tb samples
## and only keep the mouse genes
counts = counts_src
counts = counts[,-grep("BS", colnames(counts_src))]
counts = counts[,-grep("MC", colnames(counts))]
coldata = coldata[-grep("BS", coldata$sample_name),]
coldata = coldata[-grep("MC", coldata$sample_name),]

## optional: select only the ear samples or the lymph
## node samples
if (tissue == "ear") {
  counts = counts[,grep("^E", colnames(counts))]
  coldata = coldata[grep("^E", coldata$sample_name),]
}else if (tissue == "lymphnode") {
  counts = counts[,grep("^L", colnames(counts))]
  coldata = coldata[grep("^L", coldata$sample_name),]
}

## only select the infected samples as they have the most time points
coldata = coldata[coldata$Condition == "I" | coldata$Condition == "control",]
counts = counts[,colnames(counts) %in% coldata$sample_name]
dim(coldata)
dim(counts)

## create a new column that merges the information from the time point and 
## whether sample is infected or not. 
coldata$condition_per_timepoint = paste0(coldata$Timepoint, "_", coldata$Condition)

## create DESeq2 object and run it
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design= ~ condition_per_timepoint)

## remove outlier sample
if (tissue == 'ear') {
  dds = dds[, -which(colnames(dds) == "E96I2")]
}



# ---- ⊢ 1.4 select most significant genes ----
## first remove those genes with a too low expression value counted over 
## all samples
keep <- rowSums(counts(dds)) >= 1000
sum(keep)
dds <- dds[keep,]

## remove the Tb genes
dds = dds[- grep("Tb", rownames(dds)),]

## in this step a subselection can be made that consists of those genes that
## are most varying over time. This cna
dds <- DESeq(dds, test="LRT", reduced = ~ 1)
res <- results(dds)
sig_genes <- rownames(res[which(res$padj < 0.01), ])
normalized_counts <- counts(dds, normalized = TRUE)
filtered_counts <- normalized_counts[sig_genes, ]


## additional selection to limit the number of input genes
## based on the standard deviation
gene_variances <- apply(normalized_counts, 1, var)
var_threshold <- quantile(gene_variances, 0.90)
variable_genes <- normalized_counts[gene_variances > var_threshold, ]


# ---- ⊢ 1.5 average expression ----
## calculate the average expression over all the biological replicates
## and return one value per time point
if (tissue == "ear") {
  variable_genes_average = data.frame(
    "0h" = rowMeans(variable_genes[,c("E01", "E02", "E03")]),
    "12h" = rowMeans(variable_genes[,c("E12I1", "E12I2", "E12I3")]),
    "24h" = rowMeans(variable_genes[,c("E24I1", "E24I2", "E24I3")]),
    "48h" = rowMeans(variable_genes[,c("E48I1", "E48I2", "E48I3")]),
    "72h" = rowMeans(variable_genes[,c("E72I1", "E72I2", "E72I3")]),
    "96h" = rowMeans(variable_genes[,c("E96I1", "E96I3")])
  )
}


## do normalisation and extract the normalised read

#### 2. run Mfuzz #####

# ---- ⊢ 2.1 preprocessing & fuzzifier ----
# preprocessing and removal of NA values and genes with too many NA values
data <- ExpressionSet(assayData = as.matrix(variable_genes_average))
data = filter.NA(data, thres=0.25) # removes 15 genes with all NA values = 75 NA values
data = fill.NA(data, mode="mean") # no NA values present anymore so this step obsolete
data.s = standardise(data)

# a fuzifier gives an estimate on how difficult it is to cluster the 
# transcriptomics data into cluster. Closer to 1 = hard to cluster. 
m1 <- mestimate(data.s)
m1

# ---- ⊢ 2.2 filtering ----
# filter out those genes which are expressed at low levels, or only
# varying to a limited extend. 


# ---- ⊢ 2.3 find number of clusters ---- 
# find the optimal number of clusters using the Dmin funtion
Dmin(data.s, m=m1, crange=seq(2,22,1), repeats=3, visu=TRUE)
ncluster = 9

# ---- ⊢ 2.4 do clustering ----
## set seed is need to make reproducible results. If you do not set this to the 
## same number, every run of mfuzz will results in a different kind of clustering.
## overall clustering profiles will be the same, however different numbering
set.seed(5)
cl = mfuzz(data.s, c=ncluster,m=m1) 

## mfuzz.plot is the integrated plotting function of the 
# mfuzz.plot(data.s, mfrow=c(3,4),
#            cl=cl, 
#            time.labels=sampleNames(data.s), 
#            # new.window=TRUE,
#            min.mem=0.50)

length(which(cl$cluster == 9))


#### 3. Mfuzz Visualization ####

# ---- ⊢ 3.1 make a visualation of all the clusters ----

## ? also additional cutoff on membership?

figure_list=c()
i = 0

for (cluster_id in seq(1, length(unique(cl$cluster)))) {
  print(cluster_id)
  class(cl$cluster)
  genes_in_cluster = names(cl$cluster)[cl$cluster == cluster_id]
  expr_index = which(rownames(exprs(data.s)) %in% genes_in_cluster)
  expr_data = exprs(data.s)[expr_index,]
  
  expr_data_plot = melt(expr_data)
  
  average_expr = as.data.frame(colMeans(expr_data))
  average_expr$condition = rownames(average_expr)
  colnames(average_expr) = c('expression', 'condition')
  
  colnames(expr_data_plot) = c('geneID', 'condition', 'expression')
  
  nr_of_genes = dim(expr_data)[1]
  title = paste0("cluster #", cluster_id, "(n=", nr_of_genes, ")")

  p = ggplot(data = expr_data_plot)
  p = p + geom_line( aes(x=condition, y=expression, group=geneID), col="blue", alpha=0.2)
  p = p + stat_summary(aes(x=condition, y=expression,group = 1), fun = mean, color="red", geom = 'line', size=2, alpha=0.75)
  p = p + theme_bw()
  p = p + ggtitle(title) + theme(plot.title = element_text(size = 10, face = "bold"))
  p = p + theme(axis.title.x = element_blank())
  # p = p + theme(axis.text.x = element_text(angle = 30))
  p
  
  i = i + 1
  figure_list[[i]] = local({
    print(p)
  })
  

}  

# source("http://peterhaschke.com/Code/multiplot.R") #load multiplot function
# multiplot(plotlist = figure_list, cols = 3)

grid.arrange(grobs = figure_list)

## DEGreport is not working as it shoudld. resulting in difficult to solve errors.
# library('DEGreport')
# ma = assays(rlog(dds[top_genes_index,]))
# res <- degPatterns(ma, design(dds), time = "time")

# ---- ⊢ 3.2 GSEA enrichment mfuzz ----
out_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/deseq2/gsea_mfuzz/'
go_plot_files = list()
## extract ensembl IDs and convert them to entrez IDs
for (cluster_id in seq(1, length(unique(cl$cluster)))) {
  print(cluster_id)
  class(cl$cluster)
  
  ## get ensemblIDs and convert to entrezIDs
  ensembl_ids = names(cl$cluster)[cl$cluster == cluster_id]
  entrez_ids <- mapIds(org.Mm.eg.db, keys = ensembl_ids, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")

  
  ## GO enrichment analysis based on EntrezIDs
  go_enrich <- enrichGO(gene = entrez_ids, 
                        OrgDb = org.Mm.eg.db,  
                        ont = "BP",  ## vary over BP, CC, or MF
                        pAdjustMethod = "BH", 
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.2, 
                        readable = TRUE,
                        minGSSize = 20)
  
  ## filter the go_enrich output
  go_enrich_simplified <- simplify(go_enrich, cutoff = 0.5, by = "p.adjust", select_fun = min)
  dim(go_enrich_simplified)

  # Visualize GO enrichment results
  go_plot_enriched = dotplot(go_enrich_simplified, 
                             showCategory = 25,
                             font.size = 6,
  )

  ## save output to png file
  go_plot_file_enriched = paste0('mfuzz_time_series_cluster', cluster_id, '_enricheGO.png')
  ggsave(go_plot_file_enriched, plot = go_plot_enriched)
  go_plot_files[[cluster_id]] = go_plot_enriched

}

grid.arrange(grobs = go_plot_files)

grid.arrange(grobs = go_plot_files[1:6], ncol=3)
grid.arrange(grobs = go_plot_files[7:9], ncol=3, nrow=2)



#### 4. TCseq ####

# ---- ⊢ 4.1 
library(TCseq)

time_points = unique(coldata$Timepoint)
tca <- TCA(design = data.frame(time = unique(time_points)))





##### OLD CODE #####
##### OLD CODE #####
##### OLD CODE #####




# ---- 3.2 make boxplots of selected clusters -----
cluster_select = 11
membership_cutoff = 0.50
membership_taken_into_account = 1

genes_in_cluster = names(cl$cluster)[cl$cluster == cluster_select]
expr_index = which(rownames(exprs(data.s)) %in% genes_in_cluster)
expr_data = exprs(data.s)[expr_index,]

## optionally: also take membership into account. Now, a gene can be 
## assigned to a cluster, however, the membership might be low (all 
## genes are "forced" to be at least in one cluster)
if (membership_taken_into_account) {
  gene_ids = rownames(expr_data)
  membership = cl$membership[gene_ids, cluster_select]
  core_genes = names(membership[membership > membership_cutoff])
  expr_data = expr_data[core_genes,]
}

expr_data_plot = melt(expr_data)
colnames(expr_data_plot) = c('geneID', 'condition', 'expression')



## rename with name starting by number to automate the ordering. 
expr_data_plot$condition = as.character(expr_data_plot$condition)
expr_data_plot[expr_data_plot$condition == "prolog",'condition'] = "0_prolog"
expr_data_plot[expr_data_plot$condition == "prosta",'condition'] = "1_prosta"
expr_data_plot[expr_data_plot$condition == "amalog",'condition'] = "2_amalog"
expr_data_plot[expr_data_plot$condition == "amasta",'condition'] = "3_amasta"
expr_data_plot[expr_data_plot$condition == "amaint",'condition'] = "4_amaint"

plot_title = paste0("cluster #", cluster_select, " (n=", dim(expr_data)[1], ")")
p = ggplot(data = expr_data_plot, aes(y=expression, x=condition)) 
p = p + geom_boxplot(aes(fill=condition))
p = p + geom_jitter(width=0.15, alpha=0.40)
p = p + theme(axis.text.x = element_text(angle = 0))
p = p + ggtitle(plot_title)
p = p + theme_bw()
p



# ---- 3.3 heatmap of clusters -----
cluster_select = 9
membership_cutoff = 0.50
membership_taken_into_account = 1

genes_in_cluster = names(cl$cluster)[cl$cluster == cluster_select]
expr_index = which(rownames(exprs(data.s)) %in% genes_in_cluster)
expr_data = exprs(data.s)[expr_index,]

## optionally: also take membership into account. Now, a gene can be 
## assigned to a cluster, however, the membership might be low (all 
## genes are "forced" to be at least in one cluster)
if (membership_taken_into_account) {
  gene_ids = rownames(expr_data)
  membership = cl$membership[gene_ids, cluster_select]
  core_genes = names(membership[membership > membership_cutoff])
  expr_data = expr_data[core_genes,]
}

column_order = c("prolog_1","prolog_2","prolog_3", "prosta_1", "prosta_2", "prosta_3", 
                 "amalog_1", "amalog_2","amalog_3", "amasta_1", "amasta_2", "amasta_3",
                 "amaint_1", "amaint_2", "amaint_3")
rpkm_core = rpkm_all[core_genes,column_order]
pheatmap(rpkm_core, cluster_cols=FALSE, scale="row")

annotation_table = gff[gff$geneid %in% core_genes, c('geneid', 'chrom', 'start', 'stop', 'annotation')]
summary_annotation = as.data.frame(table(annotation_table$annotation))
summary_annotation[rev(order(summary_annotation$Freq)),]

# ---- 3.4 make dotplot per gene with good contribution ---- 
## this dotplot directly uses the dds object, which is created in another 
## R script, i.e. deseq2_timeseries.R
cluster_select = 9
membership_cutoff = 0.70
membership_taken_into_account = 1

genes_in_cluster = names(cl$cluster)[cl$cluster == cluster_select]
expr_index = which(rownames(exprs(data.s)) %in% genes_in_cluster)
expr_data = exprs(data.s)[expr_index,]

## optionally: also take membership into account. Now, a gene can be 
## assigned to a cluster, however, the membership might be low (all 
## genes are "forced" to be at least in one cluster)
if (membership_taken_into_account) {
  gene_ids = rownames(expr_data)
  membership = cl$membership[gene_ids, cluster_select]
  core_genes = names(membership[membership > membership_cutoff])
  expr_data = expr_data[core_genes,]
}


## create empty figure list, that will contain all pictures that should
## be shown in a multi-plot
figure_list = c()
i = 0

nr_of_plot_columns = 4
start = 1
end = 12
core_genes_select = core_genes[start:end]

for (core_gene in core_genes_select) {
  
  ## get the geneID to be used as title of the plot
  plot_title = core_gene
  index = which(rownames(dds) == core_gene)
  
  plot_data <- plotCounts(dds, index, 
                          #intgroup = c("condition", "time"), returnData = TRUE)
                          intgroup = c("condition"), returnData = TRUE)
  
  plot_data$condition = gsub("prolog", "0_prolog", plot_data$condition)
  plot_data$condition = gsub("prosta", "1_prosta", plot_data$condition)
  plot_data$condition = gsub("amalog", "2_amalog", plot_data$condition)
  plot_data$condition = gsub("amasta", "3_amasta", plot_data$condition)
  plot_data$condition = gsub("amaint", "4_amaint", plot_data$condition)
  
  # p = ggplot(plot_data, aes(x = time, y = count, group=1))
  # p = p + geom_point(aes(fill=time))
  # p = p + stat_summary(fun=mean, geom="line")
  # # p = p + geom_dotplot(binaxis="y", stackdir="center", dotsize = 2)
  # p = p + scale_y_log10()
  # p = p + theme(legend.position = "none")
  # p = p + ggtitle(plot_title)
  # p = p + theme(axis.text.x = element_text(angle = 45))
  # p
  
  # p = ggplot(plot_data, aes(x=time, y=count, fill=time))
  p = ggplot(plot_data, aes(x=condition, y=count, fill=condition))
  p = p + scale_y_log10()
  p = p + geom_dotplot(binaxis="y", stackdir="center", dotsize = 2)
  p = p + ggtitle(plot_title)
  p = p + theme(legend.position = "none")
  # p = p + stat_summary(fun=mean, geom="line")
  p = p + theme(plot.title = element_text(size=10))
  p = p + scale_x_discrete(labels=c("0_prolog" = "0_PL", 
                                "1_prosta" = "1_PS",
                                "2_amalog" = "2_AL",
                                "3_amasta" = "3_AS",
                                "4_amaint" = "4_AI")
                                )
  p

  i = i + 1
  figure_list[[i]] = local({
    print(p)
  })
  
}

source("http://peterhaschke.com/Code/multiplot.R") #load multiplot function
multiplot(plotlist = figure_list, cols = nr_of_plot_columns)






#### 4. export data ####

# ---- 4.0 simple export per cluster ----
cluster_of_interest = 11
membership_cutoff = 0.40
membership = cl$membership[, cluster_of_interest]
core_genes = names(membership[membership > membership_cutoff])
length(core_genes)
out_dir = '/Users/pmonsieurs/programming/leishmania_q/results/mfuzz/'
out_file = paste0(out_dir, 'export_mfuzz_cluster', cluster_of_interest, '.20210617.csv')
write.table(core_genes, out_file, quote=FALSE, row.names=FALSE, col.names=FALSE)

## the output of this file can be run on to extract the GO classes and the 
## KEGG pathways using the following scripts:
## 1) GO_genes_to_GOenrichment.sh & GO_genes_to_GOenrichment.py
## 2) GO_enrichment_analysis.R

# ---- 4.1 extra memberships -----
cl_memberships = as.data.frame(cl$membership)
head(cl_memberships)

# ----- 4.2 add pvalues from DESeq2 time series analysis -----
top_res = as.data.frame(res[top_genes_index,c('padj')])
rownames(top_res) = rownames(res[top_genes_index,])
top_res = res

# ----- 4.3 add normalized read counts -----
df_out = cbind(top_res[rownames(top_res) %in% rownames(cl_memberships),], cl_memberships)
colnames(df_out)[1] = "padj"
head(df_out)
# rownames(top_res)[500]
# rownames(cl_memberships)[499]

counts_top = counts(dds, normalized=TRUE)[rownames(counts(dds)) %in% rownames(df_out),]
counts_top_ordered = counts_top[match(rownames(df_out),rownames(counts_top)),]
head(counts_top_ordered)
head(df_out)

df_out = cbind(df_out, counts_top_ordered)


# ----- 4.4 add gff annotation ---- 
gff_file = '/Users/pmonsieurs/programming/leishmania_q/data/gene_annotations.gtf'
gff = read.csv(gff_file, sep="\t", comment.char="#", header=FALSE)

# extract the gene_id 
geneid = str_split_fixed(gff[,9], ";",n=Inf)[,2]
geneid = str_split_fixed(geneid, " ", n=Inf)[,3]

# extract the gene annotation
gene_annotation = str_split_fixed(gff[,9], ";",n=Inf)[,3]
gene_annotation = str_split_fixed(gene_annotation, " gene_name ",n=Inf)[,2]

gff = gff[,1:5]
colnames(gff) = c("chrom", "type", "feature", "start", "stop")
gff$geneid = geneid
gff$length = gff$stop - gff$start + 1
gff$annotation = gene_annotation
head(gff)

gff[match(rownames(df_out), gff$geneid),]$geneid[1:6]
head(df_out)
gff_out = gff[match(rownames(df_out), gff$geneid),]

df_out = cbind(df_out, gff_out)

# ---- 4.5 add the old gene ID ----
mapping_file = '/Users/pmonsieurs/programming/leishmania_10X/data/blast/Mapping_Sanger_vs_TriTrypDB.csv'
mapping_data = read.csv(mapping_file, sep="\t")
gene_ids_v1 = mapping_data[match(paste0(rownames(df_out), ".1"),mapping_data[,1]),2]
head(cbind(df_out, gene_ids_v1))
df_out = cbind(df_out, gene_ids_v1)

# ---- 4.6 write to excel ----w
write.xlsx(df_out, 'mfuzz_clustering.xlsx', row.names=TRUE)



# ---- 4.7 check the genes selected by Gosia
## check the number of genes that Gosia selected, but were not present in the 
## top 500 list:
setwd('/Users/pmonsieurs/programming/leishmania_q/results/deseq2_custom/')
missing_genes = read.csv('missing_genes_in_timeseries.csv', strip.white = TRUE, header=FALSE)
colnames(missing_genes) = 'gene_id'

index_genes_gosia = which(df_out$geneid %in% missing_genes[,1])
df_out_gosia = df_out[index_genes_gosia,]
write.xlsx(df_out_gosia, 'mfuzz_genes_gosia.xlsx', row.names=TRUE) 

