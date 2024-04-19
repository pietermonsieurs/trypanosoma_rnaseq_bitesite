library(ggplot2)
library(openxlsx)

qc_file = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/data/NVNL2024031201_SampleQCresult.xlsx'
qc_data = read.xlsx(qc_file, sheet="Data QC")
qc_data = read.xlsx(qc_file, sheet="Data QC2")
head(qc_data)

ggplot(data = qc_data, aes(x=PE150, group=QC)) + 
  geom_density(aes(color=QC)) + 
  theme_bw()


ggplot(data = qc_data, aes(y=PE150, x=QC)) + 
  geom_boxplot(aes(fill=QC)) + 
  theme_bw()

