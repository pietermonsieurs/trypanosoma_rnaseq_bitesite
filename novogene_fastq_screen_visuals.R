library(ggplot2)
library(openxlsx)

data_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/fastq_screen/'
data_file = paste0(data_dir, 'fastq_screen_summary.csv')

plot_data = read.csv(data_file)


## plot the data per species
org = 'Mouse'
org = 'Tbrucei_TREU927'
plot_data = plot_data[plot_data$organism == org,]

## define the mapping colors
mapping_colors = c("01_unmapped_pc" = "gray", 
                   "05_one_hit_one_genome" = "lightskyblue1", 
                   "04_multiple_hits_one_genome" = "cornflowerblue",
                   "03_one_hit_multiple_genomes" = "darksalmon",
                   "02_multiple_hits_multiple_genomes" = "red")

## create as factor to revert the order of coord_flip
plot_data$sample = as.factor(plot_data$sample)

# plot_data_scr = plot_data
# plot_data = plot_data[grep("Ldon", plot_data$organism),]

p = ggplot(plot_data, aes(x=sample,y=percentage,fill=type)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual("legend", values = mapping_colors) + 
  scale_x_discrete(limits=rev(levels(plot_data$sample))) + ## change order to put sample 001 on top
  coord_flip() + 
  theme(text = element_text(size=12)) + 
  theme(legend.text=element_text(size=12)) + 
  ggtitle(paste0("FastQ screen -- ", org)) 

p


#### L. aethiopica ####
## plot only mixes
plot_data = plot_data[plot_data$organism == "Laethiopica", ]
plot_data = plot_data[grep("Mix", plot_data$sample),]
plot_data$sample_short = gsub("SuSL.Laeth.Laeth.", "", plot_data$sample)
plot_data$sample_short = gsub(".R1", "", plot_data$sample_short)

p = ggplot(plot_data, aes(x=sample_short,y=percentage,fill=type))
p = p + geom_bar(position="stack", stat="identity")
p = p + scale_fill_manual("legend", values = mapping_colors)
p = p + scale_x_discrete(limits=rev(levels(plot_data$sample))) ## change order to put sample 001 on top
p = p + coord_flip()
p = p + theme(text = element_text(size=16))
p = p + theme(legend.text=element_text(size=16))
p = p + ggtitle(paste0("FastQ screen -- ", org))
p

#### L.aethiopica clinical and hyrax ####
plot_data$sample = gsub("_R1", "", plot_data$sample)
meta_data_file = '/Users/pmonsieurs/programming/leishmania_susl/data/sample_submission_form_SuSL.20230116.xlsx'
meta_data = read.xlsx(meta_data_file)
meta_data$GS_ID = gsub("105328-001-", "", meta_data$GS_ID)
meta_data$GS_ID



#### L. braziliensis ####
org = 'Lbraziliensis'
plot_data = plot_data[plot_data$organism == org,]
plot_data$sample_short = gsub("SuSL.Lbraz.Lbraz.", "", plot_data$sample)
plot_data$sample_short = gsub(".R1", "", plot_data$sample_short)

p = ggplot(plot_data, aes(x=sample_short,y=percentage,fill=type))
p = p + geom_bar(position="stack", stat="identity")
p = p + scale_fill_manual("legend", values = mapping_colors)
# p = p + scale_x_discrete(limits=rev(levels(plot_data$sample))) ## change order to put sample 001 on top
p = p + coord_flip()
p = p + theme(text = element_text(size=12))
p = p + theme(legend.text=element_text(size=12))
p = p + ggtitle(paste0("FastQ screen -- ", org))
p




#### L. tropica ####

## only take the Ltropica with mixes
## compare the two different designs
plot_Ltrop = plot_data[grep("Ltrop", plot_data$sample),]
plot_Ltrop = plot_Ltrop[grep("ix", plot_Ltrop$sample),]
plot_Ltrop = plot_Ltrop[-grep("Ampure", plot_Ltrop$sample),]
plot_Ltrop = plot_Ltrop[-grep("FPPE", plot_Ltrop$sample),]

## compare the ampure with FPPE
plot_Ltrop = plot_data[grep("Ltrop", plot_data$sample),]
# plot_Ltrop = plot_Ltrop[grep("ix", plot_Ltrop$sample),]
plot_Ltrop = plot_Ltrop[grep("Ampure|FPPE", plot_Ltrop$sample),]
# plot_Ltrop = plot_Ltrop[-grep("FPPE", plot_Ltrop$sample),]


plot_Ltrop

p = ggplot(plot_Ltrop, aes(x=sample,y=percentage,fill=type))
p = p + geom_bar(position="stack", stat="identity")
p = p + scale_fill_manual("legend", values = mapping_colors)
p = p + scale_x_discrete(limits=rev(levels(plot_data$sample))) ## change order to put sample 001 on top
p = p + coord_flip()
p = p + theme(text = element_text(size=12))
p = p + theme(legend.text=element_text(size=15))
p = p + ggtitle(paste0("FastQ screen -- ", org))
p



## only take the Ltropica clinical samples
plot_Ltrop = plot_data[grep("Ltrop", plot_data$sample),]
plot_Ltrop = plot_Ltrop[grep("Laeth.F|Laeth.I", plot_Ltrop$sample),]

plot_Ltrop$sample_short = gsub("SuSL.Ltrop.Laeth.", "", plot_Ltrop$sample)
plot_Ltrop$sample_short = gsub(".R1", "", plot_Ltrop$sample_short)

plot_Ltrop

p = ggplot(plot_Ltrop, aes(x=sample_short,y=percentage,fill=type))
p = p + geom_bar(position="stack", stat="identity")
p = p + scale_fill_manual("legend", values = mapping_colors)
p = p + scale_x_discrete(limits=rev(levels(plot_data$sample))) ## change order to put sample 001 on top
p = p + coord_flip()
p = p + theme(text = element_text(size=15))
p = p + theme(legend.text=element_text(size=15))
#p = p + ggtitle(paste0("FastQ screen -- ", org))
p = p + ggtitle("Percentage of reads mapping to L. tropica genome")
p = p + xlab("percentage of reads") + ylab("clinical sample")
p


#### L. donovania - sWGA ####
plot_swga = plot_data[grep("sWGA", plot_data$sample),]

plot_swga$sample = gsub("sWGA.Ldon.", "", plot_swga$sample)
plot_swga$sample = gsub(".R1$", "", plot_swga$sample)

## remove erroneous samples - those are duplicates
plot_swga = plot_swga[- which(plot_swga$sample == "Mix01_BPK282"),]
plot_swga = plot_swga[- which(plot_swga$sample == "Mix006_BPK282"),]


p = ggplot(plot_swga, aes(x=sample,y=percentage,fill=type))
p = p + geom_bar(position="stack", stat="identity")
p = p + scale_fill_manual("legend", values = mapping_colors)
p = p + scale_x_discrete(limits=rev(levels(plot_data$sample))) ## change order to put sample 001 on top
p = p + coord_flip()
p = p + theme(text = element_text(size=12))
p = p + theme(legend.text=element_text(size=15))
p = p + ggtitle(paste0("FastQ screen -- ", org))
p




#### 16.01.2023 ####
## analysis for the samples of 16.01.2023, with combination of different
## species (Ldon, Laethiopica, Lbraz)
plot_data = read.csv(data_file)

## selection if you only want the three donovani samples
org = 'Ldonovani'
org = 'Laethiopica'
org = 'Lbraziliensis'
plot_data = plot_data[plot_data$organism == org,]

## selection if you want all the 30 samples combined, as they are 
## coming from different species
plot_data = plot_data[plot_data$organism != 'Human',]

## metadata to add to the plot_data so that we can print meaningfull
## sample names instead of numbers
meta_data_file = '/Users/pmonsieurs/programming/leishmania_susl/data/sample_submission_form_SuSL.20230116.xlsx'
meta_data = read.xlsx(meta_data_file)
meta_data$GS_ID_short = gsub("105328-001-", "", meta_data$GS_ID)
meta_data$GS_ID_short = paste0(meta_data$GS_ID_short, "_R1")
plot_data$sample_name = meta_data[match(plot_data$sample, meta_data$GS_ID_short), ]$Customer.ID
plot_data$sample_name_full = paste0(plot_data$sample_name, "::", plot_data$organism, "::", plot_data$sample)
head(plot_data)

## create as factor to revert the order of coord_flip
plot_data$sample_name = as.factor(plot_data$sample_name)

# plot_data_scr = plot_data
# plot_data = plot_data[grep("Ldon", plot_data$organism),]

p = ggplot(plot_data, aes(x=sample_name_full,y=percentage,fill=type))
p = p + geom_bar(position="stack", stat="identity")
p = p + scale_fill_manual("legend", values = mapping_colors)
p = p + scale_x_discrete(limits=rev(levels(plot_data$sample))) ## change order to put sample 001 on top
p = p + coord_flip()
p = p + theme(text = element_text(size=12))
p = p + theme(legend.text=element_text(size=12))
p = p + ggtitle(paste0("FastQ screen"))
p


#### 05.06.2023 - SpatialCL slide ####
## create Lbraziliensis figure with SureSelect samples coming
## from two different runs

## first batch of data - coming from a batch with a mixture
## of different species. Also containing mixed species samples
## which should be removed in this context
data_file = '/Users/pmonsieurs/programming/leishmania_susl/results/fastq_screen_2022/fastq_screen_summary.csv'
plot_data = read.csv(data_file)
org = 'Lbraziliensis'
plot_data = plot_data[plot_data$organism == org,]
plot_data$sample_short = gsub("SuSL.Lbraz.Lbraz.", "", plot_data$sample)
plot_data$sample_short = gsub(".R1", "", plot_data$sample_short)
# plot_data$sample = plot_data$sample_short

## remove the mixed samples, only keep the clinical 
## samples of the first batch
plot_data = plot_data[-grep("Mix", plot_data$sample_short),]


## second batch with clinical Lbraziliensis samples. Also containing
## data from other species. Use the excel file of the sample submission 
## form to extract the sample name
data_file2 = '/Users/pmonsieurs/programming/leishmania_susl/results/fastq_screen_20230116/fastq_screen_summary.csv'
plot_data2 = read.csv(data_file2)
org = 'Lbraziliensis'
plot_data2 = plot_data2[plot_data2$organism == org,]

## metadata to add to the plot_data so that we can print meaningfull
## sample names instead of numbers
meta_data_file = '/Users/pmonsieurs/programming/leishmania_susl/data/sample_submission_form_SuSL.20230116.xlsx'
meta_data = read.xlsx(meta_data_file)
meta_data$GS_ID_short = gsub("105328-001-", "", meta_data$GS_ID)
meta_data$GS_ID_short = paste0(meta_data$GS_ID_short, "_R1")
plot_data2$sample_short = meta_data[match(plot_data2$sample, meta_data$GS_ID_short), ]$Customer.ID
# plot_data$sample_name_full = paste0(plot_data$sample_name, "::", plot_data$organism, "::", plot_data$sample)
head(plot_data2)

## create as factor to revert the order of coord_flip
plot_data_all = rbind.data.frame(plot_data, plot_data2)
plot_data_all$sample = as.factor(plot_data_all$sample)
plot_data_all$sample_short = as.factor(plot_data_all$sample_short)

## plot the data
p = ggplot(plot_data_all, aes(x=sample_short,y=percentage,fill=type))
p = p + geom_bar(position="stack", stat="identity")
p = p + scale_fill_manual("legend", values = mapping_colors)
p = p + scale_x_discrete(limits=rev(levels(plot_data$sample_short))) ## change order to put sample 001 on top
p = p + coord_flip()
p = p + theme(text = element_text(size=12))
p = p + theme(legend.text=element_text(size=12))
p = p + ggtitle(paste0("FastQ screen"))
p




#### 30.08.2023: L. aethiopica samples with skinslit and biopsies

data_dir = '/Users/pmonsieurs/programming/leishmania_susl/results/fastq_screen_laethiopica_skinslit/'
data_file = paste0(data_dir, 'fastq_screen_summary.csv')

plot_data = read.csv(data_file)
# org = 'Human'
# org = 'Ldonovani'
# org = 'Lbraziliensis'
# org = 'Laethiopica'
# org = 'Ltropica'

# plot_data = plot_data[plot_data$organism == org,]
plot_data = plot_data[plot_data$organism != 'Human',]

meta_data_file = '/Users/pmonsieurs/programming/leishmania_susl/data/Laethiopica_skinslit_metadata.xlsx'
meta_data = read.xlsx(meta_data_file)
meta_data$Initial_Id = gsub("_R1.fastq.gz", "", meta_data$Initial_Id)
meta_data$Full_name

## convert to meaning full sample name
new_cols_full = strsplit(meta_data$Full_name, "/")
meta_data$sample_short = sapply(new_cols_full, function(x) x[3])

## match with plot_data
plot_data$sample = gsub("_R1", "", plot_data$sample)
plot_data$sample_short = meta_data[match(plot_data$sample, meta_data$Initial_Id),]$Workname

## optional: select the good samples
good_samples = c("07_biopsy", "08_biopsy", "07_skinslit", "08_skinslit")
plot_data = plot_data[plot_data$sample_short %in% good_samples,]

p = ggplot(plot_data, aes(x=sample_short,y=percentage,fill=type))
p = p + geom_bar(position="stack", stat="identity")
p = p + scale_fill_manual("legend", values = mapping_colors)
p = p + scale_x_discrete(limits=rev(levels(plot_data$sample))) ## change order to put sample 001 on top
p = p + coord_flip()
p = p + theme(text = element_text(size=12))
p = p + theme(legend.text=element_text(size=12))
p = p + ggtitle(paste0("FastQ screen -- ", org))
p = p + facet_wrap(~ organism)
p

