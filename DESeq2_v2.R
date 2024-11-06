library('DESeq2')
library('ggplot2')
library("RColorBrewer")
library('pheatmap')
# library('vsn')
library('VennDiagram')
library('reshape2')
library('stringr')
library(openxlsx)


setwd('/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/star/')

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
# BiocManager::install("ComplexHeatmap")
# BiocManager::install("apeglm")

#### Section 1: Input Data  ####-------------------------

## specify whether you run the analysis for ears or for lymph nodes, in order
## to do analysis separately for both approaches
tissue = 'ear'
tissue = 'lymphnode'

# ---- 1.1: set file names ----
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



#### 2. infected versus non-infected ####

# ---- ⊢ 2.1 create sampleTable ----
## define the sample name and the biological condition as extracted 
## from the meta data excel file
coldata = read.xlsx(metadata_file)
coldata$Timepoint = as.factor(coldata$Timepoint)
coldata

# ---- ⊢ 2.2 create DESeq object with design ----

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
dds = DESeq(dds)


## the following selection is done to reduce the number of genes
## with a low expression level. Optional, as this will also have the 
## consequence that the number of genes is different for different runs
keep <- rowSums(counts(dds)) >= 500
sum(keep)
dds <- dds[keep,]

## only select those genes that are coming from Mouse, and omit the ones that
## are originating from Tb
dds <- dds[-grep("Tb", rownames(dds)),]

# do differential expression between different infected and non-infected
## conditions and this for the different time points
res = results(dds, contrast = c("condition_per_timepoint", "4_I", "4_NI"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = c('4h_fc', '4h_pval', '4h_fdr')
df_out = res

res = results(dds, contrast = c("condition_per_timepoint", "12_I", "12_NI"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = c('12h_fc', '12h_pval', '12h_fdr')
df_out = cbind.data.frame(df_out, res)

res = results(dds, contrast = c("condition_per_timepoint", "24_I", "24_NI"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = c('24h_fc', '24h_pval', '24h_fdr')
df_out = cbind.data.frame(df_out, res)

res = results(dds, contrast = c("condition_per_timepoint", "96_I", "96_NI"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = c('96h_fc', '96h_pval', '96h_fdr')
df_out = cbind.data.frame(df_out, res)


# ---- ⊢ 2.3 get annotaition information and gene lengths / Biomart ----
library(biomaRt)

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
gene_lengths <- getBM(attributes=c("ensembl_gene_id", 
                                   "transcript_length", 
                                   "external_gene_name",
                                   "description"
                                   ),
                      filters="ensembl_gene_id", 
                      values=rownames(counts(dds, normalized=FALSE)), 
                      mart=mart)

gene_metadata <- getBM(attributes=c("ensembl_gene_id", 
                                   "external_gene_name",
                                   "description"
                                    ),
                      filters="ensembl_gene_id", 
                      values=rownames(counts(dds, normalized=FALSE)), 
                      mart=mart)

# ---- ⊢ 2.4 calculate FPKM values ----
## calculate total number of mapped reads per sample, and convert 
## gene lengths from base pairs to kilobases
total_mapped_reads <- colSums(counts(dds, normalized=FALSE))
gene_lengths_kb <- gene_lengths$transcript_length / 1000

library(dplyr)  
# Group by gene and select the longest transcript for each gene
longest_transcripts <- gene_lengths %>%
  group_by(ensembl_gene_id) %>%
  filter(transcript_length == max(transcript_length)) %>%
  ungroup()

# If there are ties (i.e., multiple transcripts with the same longest length), this will return all.
# To keep just one transcript per gene, you can add:
longest_transcripts <- longest_transcripts %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# Check the result
longest_transcripts$transcript_length_kb = longest_transcripts$transcript_length/1000

## calculate FPKM in different steps
sample_fpkm <- (counts(dds, normalized=FALSE) / longest_transcripts$transcript_length_kb)
sample_fpkm <- sample_fpkm/(total_mapped_reads/1e6)
dim(sample_fpkm)
sum(rownames(sample_fpkm) == rownames(df_out))

gene_metadata_sorted = gene_metadata[match(rownames(df_out), gene_metadata$ensembl_gene_id),]

# ---- ⊢ 2.5 write to output ----

rownames_df_out = rownames(df_out)
df_out_full = cbind.data.frame(gene_metadata_sorted, df_out, sample_fpkm)
rownames(df_out_full) = rownames_df_out
head(df_out_full)

excel_out = paste0(out_dir, 'rnaseq_immune_response.', tissue,'.I_vs_NI.xlsx')


write.xlsx(df_out_full, excel_out, 
           colNames = TRUE,
           rowNames = TRUE)



#### 3. Time series experiment ####

## specify whether you run the analysis for ears or for lymph nodes, in order
## to do analysis separately for both approaches
tissue = 'ear'
tissue = 'lymphnode'

# ---- ⊢ 3.1 create sampleTable ----
## define the sample name and the biological condition as extracted 
## from the meta data excel file
coldata = read.xlsx(metadata_file)
coldata$Timepoint = as.factor(coldata$Timepoint)
coldata

# ---- ⊢ 3.2 create DESeq object with design ----

## preprocessing the counts. Remove the pure Tb samples
## and only keep the mouse genes
counts = counts_src
counts = counts[,-grep("BS", colnames(counts_src))]
counts = counts[,-grep("MC", colnames(counts))]
coldata = coldata[-grep("BS", coldata$sample_name),]
coldata = coldata[-grep("MC", coldata$sample_name),]

# dds <- DESeqDataSetFromMatrix(countData = counts,
#                               colData = coldata,
#                               design= ~ tissue)
# dds$condition = relevel(dds$tissue, ref = "ear")
# dds = DESeq(dds)


## optional: select only the ear samples or the lymph
## node samples
if (tissue == "ear") {
  counts = counts[,grep("^E", colnames(counts))]
  coldata = coldata[grep("^E", coldata$sample_name),]
}else if (tissue == "lymphnode") {
  counts = counts[,grep("^L", colnames(counts))]
  coldata = coldata[grep("^L", coldata$sample_name),]
}

coldata$condition_per_timepoint = paste0(coldata$Timepoint, "_", coldata$Condition)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design= ~ condition_per_timepoint)

if (tissue == 'ear') {
  dds = dds[, -which(colnames(dds) == "E96I2")]
}

# dds$condition = relevel(dds$Condition, ref = "control")
dds = DESeq(dds)


## the following selection is done to reduce the number of genes
## with a low expression level. Optional, as this will also have the 
## consequence that the number of genes is different for different runs
keep <- rowSums(counts(dds)) >= 500
sum(keep)
dds <- dds[keep,]

## only select those genes that are coming from Mouse, and omit the ones that
## are originating from Tb
dds <- dds[-grep("Tb", rownames(dds)),]

# do differential expression between both conditions
res = results(dds, contrast = c("condition_per_timepoint", "4_I", "0_control"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = c('4h_fc', '4h_pval', '4h_fdr')
df_out = res

res = results(dds, contrast = c("condition_per_timepoint", "12_I", "0_control"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = c('12h_fc', '12h_pval', '12h_fdr')
df_out = cbind.data.frame(df_out, res)

res = results(dds, contrast = c("condition_per_timepoint", "24_I", "0_control"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = c('24h_fc', '24h_pval', '24h_fdr')
df_out = cbind.data.frame(df_out, res)

res = results(dds, contrast = c("condition_per_timepoint", "48_I", "0_control"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = c('48h_fc', '48h_pval', '48h_fdr')
df_out = cbind.data.frame(df_out, res)

res = results(dds, contrast = c("condition_per_timepoint", "72_I", "0_control"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = c('72h_fc', '72h_pval', '72h_fdr')
df_out = cbind.data.frame(df_out, res)

res = results(dds, contrast = c("condition_per_timepoint", "96_I", "0_control"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = c('96h_fc', '96h_pval', '96h_fdr')
df_out = cbind.data.frame(df_out, res)



## extract the gene lengths and the gene annotation from BioMart to 
## add to output file or use length to calculate FPKM values

library(biomaRt)

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
gene_lengths <- getBM(attributes=c("ensembl_gene_id", 
                                   "transcript_length", 
                                   "external_gene_name",
                                   "description"
                                   ),
                      filters="ensembl_gene_id", 
                      values=rownames(counts(dds, normalized=FALSE)), 
                      mart=mart)

gene_metadata <- getBM(attributes=c("ensembl_gene_id", 
                                   "external_gene_name",
                                   "description"
                                    ),
                      filters="ensembl_gene_id", 
                      values=rownames(counts(dds, normalized=FALSE)), 
                      mart=mart)

## calculate total number of mapped reads per sample, and convert 
## gene lengths from base pairs to kilobases
total_mapped_reads <- colSums(counts(dds, normalized=FALSE))
gene_lengths_kb <- gene_lengths$transcript_length / 1000

library(dplyr)  
# Group by gene and select the longest transcript for each gene
longest_transcripts <- gene_lengths %>%
  group_by(ensembl_gene_id) %>%
  filter(transcript_length == max(transcript_length)) %>%
  ungroup()

## If there are ties (i.e., multiple transcripts with the same longest length), 
## this will return all. To keep just one transcript per gene, you can add:
longest_transcripts <- longest_transcripts %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# Check the result
longest_transcripts$transcript_length_kb = longest_transcripts$transcript_length/1000


sample_fpkm <- (counts(dds, normalized=FALSE) / longest_transcripts$transcript_length_kb)
sample_fpkm <- sample_fpkm/(total_mapped_reads/1e6)
dim(sample_fpkm)
sum(rownames(sample_fpkm) == rownames(df_out))

gene_metadata_sorted = gene_metadata[match(rownames(df_out), gene_metadata$ensembl_gene_id),]

rownames_df_out = rownames(df_out)
df_out_full = cbind.data.frame(gene_metadata_sorted, df_out, sample_fpkm)
rownames(df_out_full) = rownames_df_out


excel_out = paste0(out_dir, 'rnaseq_immune_response.', tissue,'.timeseries.xlsx')


write.xlsx(df_out_full, excel_out, 
           colNames = TRUE,
           rowNames = TRUE)


#### 4. Time series experiment - NON infected ####

## specify whether you run the analysis for ears or for lymph nodes, in order
## to do analysis separately for both approaches
tissue = 'ear'
tissue = 'lymphnode'

# ---- ⊢ 4.1 create sampleTable ----
## define the sample name and the biological condition as extracted 
## from the meta data excel file
coldata = read.xlsx(metadata_file)
coldata$Timepoint = as.factor(coldata$Timepoint)
coldata

# ---- ⊢ 4.2 create DESeq object with design ----

## preprocessing the counts. Remove the pure Tb samples
## and only keep the mouse genes
counts = counts_src
counts = counts[,-grep("BS", colnames(counts_src))]
counts = counts[,-grep("MC", colnames(counts))]
coldata = coldata[-grep("BS", coldata$sample_name),]
coldata = coldata[-grep("MC", coldata$sample_name),]

# dds <- DESeqDataSetFromMatrix(countData = counts,
#                               colData = coldata,
#                               design= ~ tissue)
# dds$condition = relevel(dds$tissue, ref = "ear")
# dds = DESeq(dds)


## optional: select only the ear samples or the lymph
## node samples
if (tissue == "ear") {
  counts = counts[,grep("^E", colnames(counts))]
  coldata = coldata[grep("^E", coldata$sample_name),]
}else if (tissue == "lymphnode") {
  counts = counts[,grep("^L", colnames(counts))]
  coldata = coldata[grep("^L", coldata$sample_name),]
}

coldata$condition_per_timepoint = paste0(coldata$Timepoint, "_", coldata$Condition)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design= ~ condition_per_timepoint)

if (tissue == 'ear') {
  dds = dds[, -which(colnames(dds) == "E96I2")]
}

# dds$condition = relevel(dds$Condition, ref = "control")
dds = DESeq(dds)


## the following selection is done to reduce the number of genes
## with a low expression level. Optional, as this will also have the 
## consequence that the number of genes is different for different runs
keep <- rowSums(counts(dds)) >= 500
sum(keep)
dds <- dds[keep,]

## only select those genes that are coming from Mouse, and omit the ones that
## are originating from Tb
dds <- dds[-grep("Tb", rownames(dds)),]

# do differential expression between both conditions
res = results(dds, contrast = c("condition_per_timepoint", "4_NI", "0_control"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = c('4h_fc', '4h_pval', '4h_fdr')
df_out = res

res = results(dds, contrast = c("condition_per_timepoint", "12_NI", "0_control"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = c('12h_fc', '12h_pval', '12h_fdr')
df_out = cbind.data.frame(df_out, res)

res = results(dds, contrast = c("condition_per_timepoint", "24_NI", "0_control"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = c('24h_fc', '24h_pval', '24h_fdr')
df_out = cbind.data.frame(df_out, res)

res = results(dds, contrast = c("condition_per_timepoint", "96_NI", "0_control"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = c('96h_fc', '96h_pval', '96h_fdr')
df_out = cbind.data.frame(df_out, res)



## extract the gene lengths and the gene annotation from BioMart to 
## add to output file or use length to calculate FPKM values

library(biomaRt)

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
gene_lengths <- getBM(attributes=c("ensembl_gene_id", 
                                   "transcript_length", 
                                   "external_gene_name",
                                   "description"
                                  ),
                      filters="ensembl_gene_id", 
                      values=rownames(counts(dds, normalized=FALSE)), 
                      mart=mart)

gene_metadata <- getBM(attributes=c("ensembl_gene_id", 
                                    "external_gene_name",
                                    "description"
                      ),
                      filters="ensembl_gene_id", 
                      values=rownames(counts(dds, normalized=FALSE)), 
                      mart=mart)

## calculate total number of mapped reads per sample, and convert 
## gene lengths from base pairs to kilobases
total_mapped_reads <- colSums(counts(dds, normalized=FALSE))
gene_lengths_kb <- gene_lengths$transcript_length / 1000

library(dplyr)  
# Group by gene and select the longest transcript for each gene
longest_transcripts <- gene_lengths %>%
  group_by(ensembl_gene_id) %>%
  filter(transcript_length == max(transcript_length)) %>%
  ungroup()

## If there are ties (i.e., multiple transcripts with the same longest length), 
## this will return all. To keep just one transcript per gene, you can add:
longest_transcripts <- longest_transcripts %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# Check the result
longest_transcripts$transcript_length_kb = longest_transcripts$transcript_length/1000


sample_fpkm <- (counts(dds, normalized=FALSE) / longest_transcripts$transcript_length_kb)
sample_fpkm <- sample_fpkm/(total_mapped_reads/1e6)
dim(sample_fpkm)
sum(rownames(sample_fpkm) == rownames(df_out))

gene_metadata_sorted = gene_metadata[match(rownames(df_out), gene_metadata$ensembl_gene_id),]

rownames_df_out = rownames(df_out)
df_out_full = cbind.data.frame(gene_metadata_sorted, df_out, sample_fpkm)
rownames(df_out_full) = rownames_df_out


excel_out = paste0(out_dir, 'rnaseq_immune_response.', tissue,'.timeseries_noninfected.xlsx')


write.xlsx(df_out_full, excel_out, 
           colNames = TRUE,
           rowNames = TRUE)





#### 5. T. brucei MC vs BS ####

# ---- ⊢ 5.1 create sampleTable ----
## define the sample name and the biological condition as extracted 
## from the meta data excel file
coldata = read.xlsx(metadata_file)
coldata$Timepoint = as.factor(coldata$Timepoint)
coldata

# ---- ⊢ 5.2 create DESeq object with design ----

## preprocessing the counts. Remove the pure Tb samples
## and only keep the mouse genes
counts = counts_src
counts = counts[,grep("BS|MC", colnames(counts_src))]
coldata = coldata[grep("BS|MC", coldata$sample_name),]

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design= ~ tissue)


## only select those genes that are coming from Tb, and omit the ones coming
## from mice
dds <- dds[grep("Tb", rownames(dds)),]

# dds$condition = relevel(dds$Condition, ref = "control")
dds = DESeq(dds)


## the following selection is done to reduce the number of genes
## with a low expression level. Optional, as this will also have the 
## consequence that the number of genes is different for different runs
# keep <- rowSums(counts(dds)) >= 500
# sum(keep)
# dds <- dds[keep,]



# do differential expression between both conditions
res = results(dds, contrast = c("tissue", "bloodstream", "metacyclics"))
res = as.data.frame(res[,c(2,5,6)])
colnames(res) = c('BS_vs_MC_fc', 'BS_vs_MC_pval', 'BS_vs_MC_fdr')
df_out = res
df_out

# ----  ⊢ 5.3 read and parse GFF file ----
# add gene lengths to calculate FPKM values 
gff_file = '/Users/pmonsieurs/programming/trypanosoma_rnaseq/data/TriTrypDB-48_TbruceiTREU927.gff'

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
dds_fpkm = fpkm(dds)



# ---- ⊢  5.4 export the data to a .csv file ----

## add the annotation to the output file
res_ann = as.data.frame(res)
head(res_ann)
gff_index = match(rownames(res_ann), rownames(gff))
res_ann <- cbind(res_ann, gff[gff_index,])
# write.csv(as.data.frame(res), file='blood_vs_metacyclic.csv')

dds_fpkm_sorted = dds_fpkm[match(rownames(res_ann), rownames(dds_fpkm)),]
res_ann = cbind.data.frame(res_ann, dds_fpkm_sorted)
head(res_ann)

tb_file_out = paste0(out_dir, 'Bloodstream_vs_metayclics.xlsx')
write.xlsx(res_ann, file=tb_file_out, rowNames = TRUE)



############## @@@@@ OLD CODE@@@@ ###########

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
index = match(rownames(dds), gff$gene_id)
mcols(dds)$basepairs = gff$length[index]
dds_fpkm = fpkm(dds)




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
colnames(sampleDistMatrix) <- paste(colnames(vsd))


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


## add annotation data
coldata = coldata[-which(coldata$sample_name== "E96I2"),]
annotation_data = coldata[, c('Condition', 'Timepoint')]
annotation_data$Timepoint = as.factor(annotation_data$Timepoint)
rownames(annotation_data) = colnames(sampleDistMatrix)

p = pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             fontsize=12, 
             annotation_row = annotation_data) #,
#         col=colors)
p

# make a PCA plot with all samples
p = plotPCA(vsd, intgroup="condition")
p = plotPCA(vsd, intgroup="Timepoint")

# p = p + geom_point(size=10, alpha=0.50)
p = p + theme_bw()
p = p + theme(text = element_text(size = 18)) 
p = p + geom_text(aes(label=colnames(vsd)), hjust=-0.25, vjust=-0.10, size=4)
p = p + coord_cartesian(xlim=c(-20,20), ylim=c(-15,15))
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

