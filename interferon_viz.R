library(ggplot2)
library(openxlsx)
library(pheatmap)
library(reshape)


library(NanoStringQCPro)
library(NanoStringDiff)
library(limma)



#### 1. RNA-seq data ####

# ----  ├ 1.1. input data ----
tissue = 'ear'
# tissue = 'lymphnode'
src_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/deseq2/'
excel_file = paste0(src_dir, 'rnaseq_immune_response.',tissue,'.timeseries.xlsx')
deseq_data = read.xlsx(excel_file)
head(deseq_data)

## interferon genes obtained from Eva - 24/11/2024
inf_gene_file = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/data_eva_AllCellTypesFC3.txt'
inf_gene_file = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/data_eva_AtLeastOneCellTypeFC3.txt'
inf_genes = read.csv(inf_gene_file, header=FALSE)
colnames(inf_genes) = c('inf_gene')

time_points = c('4h', '12h', '24h', '48h', '72h', '96h')
cutoff_fdr = 0.05
cutoff_fc = 1

match_index = match(inf_genes$inf_gene, deseq_data$external_gene_name)
match_index = match_index[- which(is.na(match_index))]
plot_data = deseq_data[match_index,]


# ----  ├ 1.2 plot FC data ----

plot_data = plot_data[,grep("fc", colnames(plot_data))]
rownames(plot_data) = deseq_data$external_gene_name[match_index]
pheatmap(plot_data,
         breaks = seq(0,5, 0.05),
         fontsize_row = 8,
         cluster_cols = FALSE)

## plot only those data where there is a fold change > 1 or < -1
plot_data_fc = plot_data[rowSums(abs(plot_data) > 2) > 0,]
pheatmap(t(plot_data_fc),
         breaks = seq(-5,5, 0.10),
         fontsize_col = 9,
         cluster_rows = FALSE)

## melt the data to make a line plot
line_data = melt(as.matrix(plot_data_fc))
colnames(line_data) = c('gene', 'time', 'fold_change')
head(line_data)

## calculate the average line
average_data <- aggregate(fold_change ~ time, data = line_data, FUN = mean)


ggplot(data = line_data, aes(x=time, y=fold_change)) + 
  geom_line(aes(group = gene), alpha = 0.50) + 
  geom_line(data = average_data, aes(x = time, y = fold_change, group = 1), color = "red", size = 2) + 
  coord_cartesian(ylim = c(-2.5,10)) + 
  theme_bw()
  



# ----  ├ 1.3 plot RPKM data ----
match_index = match(inf_genes$inf_gene, deseq_data$external_gene_name)
match_index = match_index[- which(is.na(match_index))]
plot_data = deseq_data[match_index,]

plot_data = plot_data[,grep("E", colnames(plot_data))]
rownames(plot_data) = deseq_data$external_gene_name[match_index]
pheatmap(t(log10(plot_data+0.01)),
         breaks = seq(-2,2, 0.04),
         fontsize_row = 8,
         fontsize_col = 8,
         cluster_cols = TRUE)


## melt the data to make a line plot
line_data = melt(as.matrix(plot_data))
colnames(line_data) = c('gene', 'time', 'fold_change')
head(line_data)

## calculate the average line
average_data <- aggregate(fold_change ~ time, data = line_data, FUN = mean)

## NOT VERY USEFUL PLOT
ggplot(data = line_data, aes(x=time, y=fold_change)) + 
  geom_line(aes(group = gene), alpha = 0.50) + 
  geom_line(data = average_data, aes(x = time, y = fold_change, group = 1), color = "red", size = 2) + 
  # coord_cartesian(ylim = c(-2.5,10)) + 
  theme_bw()



#### 2. ncounter data ####



#----  ├ 2.1 input data and metadata ####

out_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/ncounter/'

## check the meta data. In contrast to what the manual says, the data should
## be stored in a tab-delimited file instead of a comma separated file, so 
## use the .txt metadata and not the .csv version
meta_data_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/'
meta_data_file = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/ncounter_metadata.txt'
meta_data_file_xlsx = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/ncounter_metadata.xlsx'
meta_data = read.xlsx(meta_data_file_xlsx)
head(meta_data)


# ---- ├ 2.2 process and normalize data -----

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

rcc_data_norm <- preprocRccSet(rccSet = rcc_data, 
                               normMethod = "housekeeping",
                               bgReference="negatives")

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



# ----  ├ 2.3 plot fold change data ----
plot_data = fit2$coefficients
plot_data = plot_data[grep("Endogenous_", rownames(plot_data)),]
rownames(plot_data) = sapply(strsplit(rownames(plot_data), "_"), function(x) x[2])
head(plot_data)


match_index = match(inf_genes$inf_gene, rownames(plot_data))
match_index = match_index[- which(is.na(match_index))]
plot_data = plot_data[match_index,]


pheatmap(plot_data,
         breaks = seq(0,5, 0.05),
         fontsize_row = 8,
         cluster_cols = FALSE)

## plot only those data where there is a fold change > 1 or < -1
# plot_data_fc = plot_data[rowSums(abs(plot_data) > 2) > 0,]
pheatmap(t(plot_data_fc),
         breaks = seq(-4,4, 0.08),
         fontsize_col = 9,
         cluster_rows = FALSE)

## melt the data to make a line plot
line_data = melt(as.matrix(plot_data_fc))
colnames(line_data) = c('gene', 'time', 'fold_change')
head(line_data)

## calculate the average line
average_data <- aggregate(fold_change ~ time, data = line_data, FUN = mean)


ggplot(data = line_data, aes(x=time, y=fold_change)) + 
  geom_line(aes(group = gene), alpha = 0.50) + 
  geom_line(data = average_data, aes(x = time, y = fold_change, group = 1), color = "red", size = 2) + 
  coord_cartesian(ylim = c(-2.5,10)) + 
  theme_bw()




# ----  ├ 2.4 plot raw expression data data ----
plot_data = rcc_data_norm@assayData$normData
colnames(plot_data) = sapply(strsplit(colnames(plot_data), "_"), function(x) x[3])
index_standard = grep("Standard", colnames(plot_data))
plot_data = plot_data[, -index_standard]

plot_data = plot_data[grep("Endogenous_", rownames(plot_data)),]
rownames(plot_data) = sapply(strsplit(rownames(plot_data), "_"), function(x) x[2])
head(plot_data)


match_index = match(inf_genes$inf_gene, rownames(plot_data))
match_index = match_index[- which(is.na(match_index))]
plot_data = plot_data[match_index,]

pheatmap(plot_data,
         breaks = seq(-2,2, 0.04),
         #breaks = seq(-5,5, 0.10),
         fontsize_row = 8,
         fontsize_col = 8,
         cluster_cols = TRUE)


