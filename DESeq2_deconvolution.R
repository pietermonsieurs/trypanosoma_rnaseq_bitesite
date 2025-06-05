library(DESeq2)
library(openxlsx)
library(pheatmap)

## most tools for deconvolution are difficult to install, so strategy now is to
## convert (part of) the mouse expression data to the immunology data set in 
## humans as used in CybersortX 

src_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/star/'
cibersort_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/cibersortx/'


#### 1. Input data ####

## specify whether you run the analysis for ears or for lymph nodes, in order
## to do analysis separately for both approaches
tissue = 'ear'
# tissue = 'lymphnode'

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



#### 2. subselect the data ####

## only select the relevant conditions from the expression matrix + focus only 
## on either ears or lymph nodes
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


#### 3. merge / conversion using the LM22 dataset ####

## LM22 is a dataset on CybersortX containing immune cell signatures for 
## 22 different immune cell types. This 

## get LM22 genes as download from CyberSortX. The full signature matrix has 
## been downloaded and stored under 
## /Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/cibersortx/LM22.txt
## but gene list (LM22_genes.txt) only contains the first column
M22_genes = as.vector(unlist(read.table('/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/cibersortx/LM22_genes.txt')))

## for the M22 genes, look them up with gprofiler and use g:orth to convert to 
## mouse (including Ensembl ID ENSMUSxxx)
M22_genes_with_orthologs = read.csv("/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/cibersortx/gProfiler_hsapiens_mmusculus_orthologs.csv")

## select the ENSMUS counter part in mouse for each of the genes in the M22 
## biomarker dataset
M22_in_mouse_ensmus = rownames(counts)[rownames(counts) %in% M22_genes_with_orthologs$ortholog_ensg]
counts_M22 = counts[M22_in_mouse_ensmus,]


## some gene names of M22 match with multiples ENSMUS id, so the rownames will
## not be unique anymore. Only select the first occurence of this ENSMUS row
counts_M22$new_rowname = M22_genes_with_orthologs[match(rownames(counts_M22), M22_genes_with_orthologs$ortholog_ensg),]$initial_alias
counts_M22_unique <- counts_M22[!duplicated(counts_M22$new_rowname), ]
rownames(counts_M22_unique) = counts_M22_unique$new_rowname
counts_M22_unique = counts_M22_unique[, -ncol(counts_M22_unique)]

## make sure that "Gene" is printed as column name in the output

counts_M22_unique_out =cbind(Gene = rownames(counts_M22_unique), counts_M22_unique)
write.table(counts_M22_unique_out, file = paste0(cibersort_dir, "mixture_file_only_matches_", tissue, ".txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)


## cibersortX might work better if you do not only use those genes which match
## with the signature matrix M22, but also include the other genes that do not 
## match with the M22 genes (it will confuse the internal normalization / 
## regression when only matched genes are used --> no difference with the other
## setup, but this one used to be sure... (so full datasets imported in cibersortX)
counts_non_match <- counts[! rownames(counts) %in% M22_in_mouse_ensmus, ]
counts_full = rbind.data.frame(counts_M22_unique, counts_non_match)
counts_full_out <- cbind(Gene = rownames(counts_full), counts_full)
write.table(counts_full_out, file = paste0(cibersort_dir, "mixture_file_full_", tissue, ".txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)



#### 4. visualization ####

## running CibersortX: 1) login, 2) go to Menu, 3) run CibersortX, 4) click on 
## number 2. Impute Cell Fractions, 5) click on "Custom", 6) select as 
## reference sample file "LM22", and as mixture "mixture_file_full_ear". Run and 
## download .txt

## after running the CiberSortX online, do visualisation of the txt-file output
## of the heatmap
tissue = 'ears'
tissue = 'lymphnodes'

cibersortx_file = paste0(cibersort_dir, tissue, "_full_out.txt")
cs_data = read.table(cibersortx_file, 
                     sep="\t", 
                     header = TRUE)
rownames(cs_data) = cs_data$Mixture
cs_data = cs_data[, -1]
cs_data = cs_data[,1:(ncol(cs_data) - 3)]
head(cs_data)

## create heatmap
heatmap = pheatmap(cs_data, 
         cluster_cols = FALSE,
         cluster_rows = FALSE)
heatmap_file = paste0(cibersort_dir, 'heatmap_', tissue, '.png')
ggsave(filename = heatmap_file, plot = heatmap, width=16, height=9)


## melt and make barplot per type
plot_data = melt(as.matrix(cs_data))
colnames(plot_data) = c('sample', 'immune_type', 'fraction')

barplot = ggplot(plot_data, aes(x=sample, y=fraction)) + 
  geom_col() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_wrap(~ immune_type, scales="free")

barplot_file = paste0(cibersort_dir, 'barplot_', tissue, '.png')
ggsave(filename = barplot_file, plot = barplot, width=16, height=9)


