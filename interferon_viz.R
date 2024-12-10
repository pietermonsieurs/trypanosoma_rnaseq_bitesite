library(ggplot2)
library(openxlsx)
library(pheatmap)
library(reshape)


#### 1. read in fold changes ####
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


## plot all the data
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
  

