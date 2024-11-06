library(ggplot2)
library(openxlsx)
library(reshape2)
library(scales) 


## STAR Log out parsing using grep and cut
mapping_data_file= '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/mapping_stats.xlsx'


mapping_data = read.xlsx(mapping_data_file)

## remove last 2 columns (with sum + OK)
mapping_data = mapping_data[, -ncol(mapping_data)]
mapping_data = mapping_data[, -ncol(mapping_data)]


mapping_data_melt = melt(mapping_data)
mapping_data_melt = mapping_data_melt[- which(mapping_data_melt$variable == "total"),]
head(mapping_data_melt)


## create stacked bar plot
ggplot(mapping_data_melt, aes(y=value, fill=variable, x=sample)) + 
  geom_bar(position="stack", stat="identity", alpha=1) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=20)) + 
  scale_y_continuous(labels = comma) 

## stacked barplot - percentage 
ggplot(mapping_data_melt, aes(y=value, fill=variable, x=sample)) + 
  geom_bar(position="fill", stat="identity", alpha=1) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=20)) + 
  scale_y_continuous(labels = comma) 

