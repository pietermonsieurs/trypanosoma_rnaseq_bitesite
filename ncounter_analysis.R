## install tools required for doing the analysis of the nCounter 
## independt of the nSolver software
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("NanoStringQCPro")
# BiocManager::install("NanoStringDiff")


# Load the packages
library(NanoStringQCPro)
library(NanoStringDiff)
library(limma)
library(pheatmap)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(openxlsx)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(EnhancedVolcano)

 
#### 1. input data and metadata ####

# ---- ├ 1.1 parameter settings ----
out_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/ncounter/'

## check the meta data. In contrast to what the manual says, the data should
## be stored in a tab-delimited file instead of a comma separated file, so 
## use the .txt metadata and not the .csv version
meta_data_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/'
meta_data_file = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/ncounter_metadata.txt'
meta_data_file_xlsx = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/ncounter_metadata.xlsx'
meta_data = read.xlsx(meta_data_file_xlsx)
head(meta_data)


# ---- ├ 1.2. nanostring data -----

## specity the rlf file, which is downloaded from the nanostring website
## https://nanostring.com/products/ncounter-assays-panels/immunology/host-response/
rlf_file = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/ncounter/NS_Mm_HostResponse_v1.0.rlf'

rcc_directory <- "/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/ncounter/combined/"
rcc_data <- newRccSet(rccFiles = dir(rcc_directory, full.names = TRUE),
                      extraPdata = meta_data_file,
                      rlf = rlf_file,
                      blankLabel = "blank"
                      )
checkRccSet(rcc_data)

## inspect the resulting data set
str(max.level=2, example_rccSet)
rcc_data@phenoData@data



#### 2. normalisation and QC ####

# ---- ├ 2.1 normaliation ----

rcc_data_norm <- preprocRccSet(rccSet = rcc_data, 
                               normMethod = "housekeeping",
                               bgReference="negatives")


# ---- ├ 2.2 QC function ----

rcc_data_norm_qc <- makeQCReport(rcc_data_norm, 
                                 paste0(out_dir, "nCounter_QC_report.html"),
                                 sampleNameCol = "sample_id")



# ---- ├ 2.3 correlation analysis -----

## correlation between biological replicates
meta_data$group = paste0(meta_data$timepoint, meta_data$infection)
expr_matrix = rcc_data_norm@assayData$normData
colnames(expr_matrix) = unname(sapply(colnames(expr_matrix), function(x) unlist(strsplit(x, split = "\\_"))[3]))
colnames(expr_matrix)[24] = "Standard2"

# filter per "biological condition", e.g. per time point, or timepoint x infection
for (group in unique(meta_data$group)) {
  samples = colnames(expr_matrix)[grep(group, colnames(expr_matrix))]
  print(samples)
}



## calculate the correlation between all samples, and plot the correlations
## in a heatmap. Remove the samples containing the "standard" and optionally
## remove 96hI2 as this seems to be an outlier. Also 24hI2 is present in both
## batches, so one should be renamed (hard coded) to 24hI2.2 (this is also the
## name in the metadata file)
expr_matrix = rcc_data_norm@assayData$normData
colnames(expr_matrix) = unname(sapply(colnames(expr_matrix), function(x) unlist(strsplit(x, split = "\\_"))[3]))
colnames(expr_matrix)[24] = "Standard2"
colnames(expr_matrix)[grep("24hI2", colnames(expr_matrix))[2]] = "24hI2.2"
expr_matrix = expr_matrix[, -grep("Standard", colnames(expr_matrix))]
expr_matrix = expr_matrix[, - which(colnames(expr_matrix) == "96hI2")]
cor_matrix <- cor(expr_matrix, method = "pearson")

annotation_data = meta_data[match(colnames(cor_matrix), meta_data$sample_id),c('timepoint', 'infection', 'batch')]
rownames(annotation_data) = colnames(cor_matrix)

pheatmap(cor_matrix,
         clustering_distance_rows = "euclidean",  
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         main = "Sample Correlation Heatmap",
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Choose a color palette
         show_rownames = TRUE,  
         show_colnames = TRUE,
         annotation_col = annotation_data
)

grep("Standard", colnames(cor_matrix))


## calculate the R^2 value 

rsq <- function (x, y) cor(x, y) ^ 2
samples = colnames(expr_matrix)
all_plots = list()

## Create an empty grid to store the plots
grid <- expand.grid(i = 1:24, j = 1:24)
count = 0
for (i in 1:(length(samples))) {
  ## first method
  sample1 = samples[i]

  
  for (j in 1:length(samples)) {
    ## second method
    sample2 = samples[j]
    
    # if (i == j) {next}
    
    r2 = rsq(expr_matrix[, sample1], expr_matrix[,sample2])
    
    count = count + 1
    print(paste0("[[", count, "]] ", i, " -- ", j, " = ", r2))
    
    plot_data = as.data.frame(expr_matrix[,c(sample1, sample2)])
    colnames(plot_data) = c('sample1', 'sample2')
    
    p = ggplot(plot_data, aes(x = sample1, y = sample2))    
    p = p + geom_point(alpha=0.3, size=1.5, color="black")
    p = p + geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) 
    p = p + theme_bw()
    # p = p + stat_cor(aes(label = ..rr.label..), color = "red", geom = "label")
    p = p + stat_cor(method = "pearson", label.x = 0, label.y = 5, color="red")
    # p =    
    p
    
    all_plots[[paste0(i, "_", j)]] = p
  }

  
}

p = wrap_plots(all_plots, ncol=24)
correlation_file = paste0(out_dir, 'correlation_samples.png')
ggsave(p, file=correlation_file, width=16, height=9)



# ---- ├  2.4 PCA ----

## do PCA analysis
expr_matrix_scaled <- scale(expr_matrix)
pca_result <- prcomp(t(expr_matrix_scaled))  # Transpose to treat samples as rows
summary(pca_result)  

## add metadata for colouring in PCA plot 
meta_data_ordered <- meta_data[match(rownames(pca_data), meta_data$sample_id), c('timepoint', 'infection', 'batch')]
pca_data <- as.data.frame(pca_result$x)  
pca_data <- cbind.data.frame(pca_data, meta_data_ordered)
pca_data$batch = as.factor(pca_data$batch)

ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = timepoint), size = 3) +
  # geom_point(aes(color = batch), size = 3) +
  # geom_point(aes(color = infection), size = 3) +
  # geom_text(aes(label = rownames(pca_data)), vjust = -0.5, label.size=0, check_overlap = TRUE) +
  # geom_text(aes(label = rownames(pca_data)), vjust = -0.5) +
  geom_text_repel(aes(label = rownames(pca_data)),                # Add labels with lines
                  box.padding = 0.35,                             # Extra space around labels
                  point.padding = 0.3,                            # Space between points and labels
                  max.overlaps = Inf) +
  theme_minimal() +
  labs(title = "PCA of Samples", x = "PC1", y = "PC2")


#### 3. differential expression ####

# ---- ├  3.1 using NanoStringDiff ---- 

## obsolete: glm.LRT gets stuck in an endless loop
## do differential expression analysis using NanoStringDiff. For this 
## you need to split up between endogenous probes, positve, negative 
## and housekeeping probbes. 
ncounter_data = rcc_data@assayData[["exprs"]]
endogenous = ncounter_data[grep("Endogenous", rownames(ncounter_data)),]
negative = ncounter_data[grep("Negative", rownames(ncounter_data)),]
positive = ncounter_data[grep("Positive", rownames(ncounter_data)),]
housekeeping = ncounter_data[grep("Housekeeping", rownames(ncounter_data)),]

designs = data.frame(time_point = meta_data$timepoint,
                     infection = meta_data$infection, 
                     group = paste0(meta_data$timepoint, "_", meta_data$infection))
                     
                       
ns_data = createNanoStringSet(endogenous,
                    positive,
                    negative,
                    housekeeping,
                    designs)


pData(ns_data)
head(exprs(ns_data))

## do normalisation on the nanostring data
ns_data = estNormalizationFactors(ns_data)

## create design matrix to perform differential expression analysis
pheno = pData(ns_data)
group = pheno$group
design.full = model.matrix(~0+group)

## do differential expression. This takes forever, somewhere stuck in 
## an endless loop? 
result_12hI =glm.LRT(ns_data, design.full, 
                     Beta=ncol(design.full), 
                     contrast=c(-1,1,0,0,0,0,0,0,0))


# ---- ├ 3.2 edgeR ----
## RUVseq
## edgeR on normalised data, but edgeR does not except negative
## values
group <- designs$group # Adjust to your conditions
dge <- DGEList(counts = rcc_data_norm@assayData$normData, group = designs$group)

## Skip TMM normalization for already normalized data. normalised using the 
## NanoStringQCPro package
dge$samples$norm.factors <- 1

design <- model.matrix(~ 0 + group)
dge <- estimateDisp(dge, design)


# ---- ├ 3.3 limma ----

## configure the experimental design
designs = data.frame(time_point = meta_data$timepoint,
                     infection = meta_data$infection, 
                     group = paste0(meta_data$timepoint, "_", meta_data$infection))

g_ = factor(designs$group)
g_ = relevel(g_, ref= "0h_I")
design <- model.matrix(~ 0 + g_)


## fit the linear model
fit <- lmFit(rcc_data_norm@assayData$normData, design)
fit <- eBayes(fit)
fit

## make the contrasts, each time comparing with the control at 0h
contrast_matrix <- makeContrasts(g_4h_I-g_0h_I,
                                g_4h_NI-g_0h_I,
                                g_12h_I-g_0h_I,
                                g_12h_NI-g_0h_I,
                                g_24h_I-g_0h_I,
                                g_24h_NI-g_0h_I,
                                g_96h_I-g_0h_I,
                                g_96h_NI-g_0h_I,
                                levels = design)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

fit2$coefficients
fit2$p.value
dim(fit2$coefficients)
dim(fit2$p.value)


## create heatmaps on all genes
pheatmap(fit2$coefficients)

## create heatmap on diff expressed genes
diff_expressed = topTable(fit2, p.value=0.01, lfc=1, n=200) 
rownames(diff_expressed) = gsub("Endogenous_", "", rownames(diff_expressed))
pheatmap(diff_expressed[,1:8],
         fontsize_row = 7, 
         treeheight_col = 0)


## check with output Emma using nSolver. She compared for time point 12h and 
## time point 4h the infected with the non-infected samples
contrast_matrix_emma <- makeContrasts(g_4h_I-g_4h_NI,
                                     g_12h_I-g_12h_NI,
                                     levels = design)
fit_emma <- contrasts.fit(fit, contrast_matrix_emma)
fit_emma <- eBayes(fit_emma)
topTable(fit_emma, coef=2, n=100)
topTable(fit_emma, coef=1)

## volcano plots per comparison
for (coeff in 1:8) {
  condition = colnames(contrast_matrix)[coeff]
  condition = gsub("g_", "", condition)
  condition = gsub(" ", "", condition)
  
  diff_expressed = topTable(fit2, coef=coeff, n=1000, p.value=1,lfc=0)
  diff_expressed = diff_expressed[grep("Endogenous_", rownames(diff_expressed)),]
  rownames(diff_expressed) = sapply(strsplit(rownames(diff_expressed), "_"), `[`, 2)

  x = EnhancedVolcano(diff_expressed,
                  rownames(diff_expressed),
                  x= "logFC",
                  y="adj.P.Val", 
                  xlim = c(min(diff_expressed$logFC, na.rm = TRUE) - 0.5, max(diff_expressed$logFC, na.rm = TRUE) + 0.5),
                  ylim = c(0, max(-log10(diff_expressed$adj.P.Val), na.rm = TRUE) + 0.2),
                  pCutoff = 5e-02, 
                  title = condition,
                  subtitle = "",
                  titleLabSize = 14,
                  labSize = 4,
                  subtitleLabSize = 2,
                  captionLabSize = 10,
                  legendPosition = "right",
                  legendLabels = c("NS", "FC", "Pval", "FC+pval"),
                  drawConnectors = TRUE,
                  arrowheads = FALSE

  )
  
  png_out_file = paste0(out_dir, "volcano_", condition, ".png")
  ggsave(png_out_file, plot = p_volcano, width= 10, height=8)
  
  
}

## correlation with the RNA-seq data sets. First read the DESeq2 output, which
## is derived for infected (different time points) versus 0h
tissue = 'ear'
src_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/deseq2/'
excel_file = paste0(src_dir, 'rnaseq_immune_response.',tissue,'.timeseries.xlsx')
deseq_data = read.xlsx(excel_file)

count = 0
plot_list = list()
for (timepoint in c("4h", "12h", "24h", "96h")) {
  column_name_deseq = paste0(timepoint, "_fc")
  column_name_limma = paste0("g_", timepoint,"_I - g_0h_I")
  
  ## get RNA-seq data
  rnaseq_data = deseq_data[, c('external_gene_name', column_name_deseq)]
  
  ## get nCounter data
  ncounter_list = fit2$coefficients[, column_name_limma]
  ncounter_list = ncounter_list[grep("Endogenous", names(ncounter_list))]
  ncounter_data = as.data.frame(ncounter_list)
  rownames(ncounter_data) = sapply(strsplit(rownames(ncounter_data), "_"), `[`, 2)
  rnaseq_sorted = rnaseq_data[match(rownames(ncounter_data), rnaseq_data$external_gene_name),]
  dim(rnaseq_sorted)
  dim(ncounter_data)
  plot_data = cbind.data.frame(rnaseq_sorted, ncounter_data)
  colnames(plot_data) = c('gene_id', 'rnaseq_fc', 'ncounter_fc')
  
  ## calculate and plot correlation
  cor_value = cor(plot_data[,2], plot_data[,3], use="complete.obs")
  
  p = ggplot(plot_data, aes(x = rnaseq_fc, y = ncounter_fc))    
  p = p + geom_point(alpha=0.3, size=1.5, color="black")
  p = p + geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) 
  p = p + theme_bw()
  # p = p + stat_cor(aes(label = ..rr.label..), color = "red", geom = "label")
  p = p + stat_cor(method = "pearson", label.x = 0, label.y = 6, color="red")
  p = p + ggtitle(timepoint)
  p
  
  plot_list[[timepoint]] = p
}

wrap_plots(plot_list, ncol=2)

#### 4. GSEA analysis ####

# ---- ├ 4.1 using publicly available databases ----

## Extract gene symbols and convert to EntrezIDS

# for (coeff in c(1,3,5,7)) {
for (coeff in c(2,4,6,8)) {
    
  sample_name = colnames(fit2$contrasts)[coeff]
  sample_name = gsub("g_", "", sample_name)
  sample_name = gsub(" ", "", sample_name)
  print(sample_name)
  
  diff_expressed = topTable(fit2, coef = coeff, number = 10000,  p.value=0.05, lfc=1, sort.by = "logFC")
  gene_symbols <- sapply(strsplit(rownames(diff_expressed), "_"), `[`, 2)
  entrez_ids <- mapIds(org.Mm.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  print(entrez_ids)
  
  ## GO enrichment analysis
  go_enrich <- enrichGO(gene = entrez_ids,
                        OrgDb = org.Mm.eg.db,
                        ont = "ALL",  ## vary over BP, CC, or MF
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2,
                        readable = TRUE,
                        minGSSize = 20)
  
  ## you can usse entrez_ids here, as the genes are already sorted when you 
  ## export them via top
  gene_list = diff_expressed[,1]
  gene_list = setNames(gene_list, entrez_ids)
  
  # go_gse <- gseGO(gene = gene_list,
  #                 OrgDb = org.Mm.eg.db,
  #                 ont = "ALL",  ## vary over BP, CC, or MF
  #                 pAdjustMethod = "BH",
  #                 pvalueCutoff = 0.05)
  # dim(go_gse)
  
  ## filter the go_enrich output
  go_enrich_simplified <- simplify(go_enrich, cutoff = 0.5, by = "p.adjust", select_fun = min)
  dim(go_enrich_simplified)
  # 
  # go_gse_simplified <- simplify(go_gse, cutoff = 0.5, by = "p.adjust", select_fun = min)
  # dim(go_gse_simplified)
  
  
  
  
  # Visualize GO enrichment results
  go_plot_enriched = dotplot(go_enrich_simplified,
                             showCategory = 25,
                             font.size = 7,
  )
  # 
  # 
  # go_plot_gse = dotplot(go_gse_simplified,
  #                       showCategory = 25,
  #                       font.size = 10,
  #)
  
  go_plot_file_enriched = paste0(out_dir, 'gsea/', sample_name, '_enricheGO.png')
  # go_plot_file_gse = paste0(out_dir, 'gsea/', sample_name, '_gsego.png')
  
  ggsave(go_plot_enriched, file = go_plot_file_enriched)
  # ggsave(go_plot_gse, file = go_plot_file_gse)

}


  
  
# ---- ├ 4.2 using categories of Nanostring using hypergeometric distribution ----

# ---- ├ ├ 4.2.1 read in data and parse ----

## pathway data
class_file = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/ncounter/LBL-10835-01_Mouse-Host-Response-Gene-List.xlsx'
sheet = 'Annotations'
class_data_df = read.xlsx(class_file, sheet = sheet, colNames=TRUE, startRow=2)
head(class_data_df) 

class_data = melt(class_data_df, id.vars=c('Gene', 'Human.Gene', 'Cell.Type'))
colnames(class_data) = c('Gene', 'Human_Gene', 'Cell_type', 'pathway', 'value')
class_data = class_data[class_data$value == "+", ]
dim(class_data)
head(class_data)

## get the counts of gene per pathway
pathway_count = as.data.frame(table(class_data$pathway))
colnames(pathway_count) = c('pathway', 'class_count')



# ---- ├ ├ 4.2.2 run GSEA using phyper ----

plot_data = as.data.frame(pathway_count$pathway)
for (coeff in 1:8) {
  condition = colnames(contrast_matrix)[coeff]
  condition = gsub("g_", "", condition)
  condition = gsub(" ", "", condition)
  
  diff_expressed = topTable(fit2, coef=coeff, n=1000, p.value=0.05,lfc=0)
  diff_expressed = diff_expressed[grep("Endogenous_", rownames(diff_expressed)),]
  rownames(diff_expressed) = sapply(strsplit(rownames(diff_expressed), "_"), `[`, 2)
  genes_diff = rownames(diff_expressed)
  
  class_data_sub = class_data[class_data$Gene %in% genes_diff,]
  pathway_count_cond = as.data.frame(class_data_sub$pathway)
  dim(pathway_count_cond)
  
  ## count the number of genes that are diff expressed in the ncounter
  ## dataset and merge with the overall classes. Also add the total 
  ## number of genes assigned to pathways and the number of genes
  ## that are diff expressed
  pathway_count_sub = as.data.frame(table(class_data_sub$pathway))
  colnames(pathway_count_sub) = c('pathway', 'diff_count')
  gsea_data = merge(pathway_count, pathway_count_sub, by="pathway")
  
  gsea_data$class_count_total = sum(gsea_data$class_count)
  gsea_data$diff_count_total = sum(gsea_data$diff_count)
  head(gsea_data)
  
  if (sum(gsea_data$diff_count) > 0) {
    gsea_data$pval = 1-phyper(gsea_data$diff_count, 
                              gsea_data$class_count,
                              gsea_data$class_count_total-gsea_data$class_count, 
                              gsea_data$diff_count_total)
  }else{
    gsea_data$pval = 1
  }
  
  # print(gsea_data[gsea_data$pval < 0.05,])
  plot_data = cbind.data.frame(plot_data, gsea_data$pval)

}


## make heatmap
rownames(plot_data) = plot_data[,1]
plot_data = plot_data[,-1]
new_colnames = gsub(" - g_0h_I", "", colnames(contrast_matrix))
colnames(plot_data) = new_colnames
head(plot_data)

plot_data_log10 = -log10(plot_data)
breaks = seq(0,3,0.03)
pheatmap(plot_data_log10,
         cluster_cols = FALSE, 
         cluster_rows = FALSE,
         fontsize_row = 10, 
         breaks = breaks)


plot_data_log10_sign = plot_data_log10[rowSums(plot_data_log10 > -log10(0.05)) > 0,]
breaks = seq(0,3,0.03)
pheatmap(plot_data_log10_sign,
         cluster_cols = FALSE, 
         cluster_rows = FALSE,
         fontsize_row = 12, 
         breaks = breaks)

## To Do ##
## correlation between normalised values of the biological replicates
## enhanced volcanoplot per condition






## load normalised data into 



#### example data & code snippets ####
exampleDataDir <- system.file("extdata", package="NanoStringQCPro")
rccDir <- file.path(exampleDataDir, "RCC")
example_rccSet <- newRccSet(
  rccFiles = dir(rccDir, full.names=TRUE)
  ,rlf = file.path(exampleDataDir, "RLF", "NQCP_example.rlf")
  ,cdrDesignData = file.path(exampleDataDir, "CDR", "CDR-DesignData.csv")
  ,extraPdata = file.path(exampleDataDir, "extraPdata", "SampleType.txt")
  ,blankLabel = "blank"
)



