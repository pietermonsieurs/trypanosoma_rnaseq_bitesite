#### 0. installation of the libraries ####

# Install necessary packages (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
# BiocManager::install("clusterProfiler")
## BiocManager::install("org.Hs.eg.db")  # For human data
# BiocManager::install("org.Mm.eg.db")  # For mouse data 
# BiocManager::install("AnnotationDbi")
# BiocManager::install("enrichplot")  
# if (!requireNamespace("rrvgo", quietly = TRUE)) {
#   remotes::install_github("ssayols/rrvgo")
# }
if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
  install.packages("ggVennDiagram")
}

# Load libraries
library(clusterProfiler)
library(org.Mm.eg.db)  # Replace with the appropriate organism database
library(AnnotationDbi)
library(enrichplot)
library(openxlsx)
library(ggplot2)
library(rrvgo)
library(GOSemSim)
library(ggVennDiagram)

#### 1. GSEA for time series approach ####
tissue = 'ear'
# tissue = 'lymphnode'
src_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/deseq2/'
excel_file = paste0(src_dir, 'rnaseq_immune_response.',tissue,'.timeseries.xlsx')
deseq_data = read.xlsx(excel_file)

time_points = c('4h', '12h', '24h', '48h', '72h', '96h')
cutoff_fdr = 0.05
cutoff_fc = 1

enriched_GO_list = list()
gsego_list = list()
enriched_GO_plot_list = list()
gsego_plot_list = list()

GO_types = c('BP', 'MF', 'CC')
# GO_type = 'ALL'

for (time_point in time_points) {
  print(time_point)
  
  col_name_fc = paste0(time_point, '_fc')
  col_name_fdr = paste0(time_point, '_fdr')
  
  ## select the differentially expressed genes
  diff_indices = which(deseq_data[,col_name_fc] >= cutoff_fc & deseq_data[,col_name_fdr] <= cutoff_fdr)
  ensembl_ids = deseq_data[diff_indices,]$ensembl_gene_id
  
  ## select all ensembl IDs for the gene set enrichment (for gseGO)
  ensembl_ids_fc = deseq_data[, c('ensembl_gene_id', col_name_fc)]
  ensembl_ids_fc_sorted = ensembl_ids_fc[order(ensembl_ids_fc[, col_name_fc], decreasing = TRUE),]
  ensembl_ids_sorted = ensembl_ids_fc_sorted$ensembl_gene_id
  
  ## extract ensembl IDs and convert them to entrez IDs
  entrez_ids <- mapIds(org.Mm.eg.db, keys = ensembl_ids, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
  entrez_ids_sorted <- mapIds(org.Mm.eg.db, keys = ensembl_ids_sorted, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
  
  geneList = ensembl_ids_fc_sorted[, col_name_fc]
  names(geneList) = entrez_ids_sorted
  geneList <- geneList[!is.na(names(geneList))]
  
  
  for (GO_type in GO_types) {
  #   print(GO_type)
    
    ## GO enrichment analysis
    go_enrich <- enrichGO(gene = entrez_ids,
                          OrgDb = org.Mm.eg.db,
                          ont = "ALL",  ## vary over BP, CC, or MF
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2,
                          readable = TRUE,
                          minGSSize = 20)

    go_gse <- gseGO(gene = geneList,
                          OrgDb = org.Mm.eg.db,
                          ont = "ALL",  ## vary over BP, CC, or MF
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05)
    dim(go_gse)
    
    ## filter the go_enrich output
    go_enrich_simplified <- simplify(go_enrich, cutoff = 0.5, by = "p.adjust", select_fun = min)
    dim(go_enrich_simplified)
  
    go_gse_simplified <- simplify(go_gse, cutoff = 0.5, by = "p.adjust", select_fun = min)
    dim(go_gse_simplified)
    
      
    
    
    # Visualize GO enrichment results
    go_plot_enriched = dotplot(go_enrich_simplified,
            showCategory = 25,
            font.size = 7,
    )
    
    
    go_plot_gse = dotplot(go_gse_simplified,
            showCategory = 25,
            font.size = 10,
    )
    
  
    go_plot_file_enriched = paste0(src_dir, 'gsea/', 'time_series_', tissue, '_', time_point, '_', GO_type, '_enricheGO.png')
    go_plot_file_gse = paste0(src_dir, 'gsea/', 'time_series_', tissue, '_', time_point, '_', GO_type, '_gsego.png')
    
    ggsave(go_plot_enriched, file = go_plot_file_enriched)
    ggsave(go_plot_gse, file = go_plot_file_gse)
    
    
    enriched_GO_list[[time_point]][[GO_type]] = go_enrich_simplified
    gsego_list[[time_point]][[GO_type]] = go_gse_simplified
    enriched_GO_plot_list[[time_point]][[GO_type]] = go_plot_enriched
    gsego_plot_list[[time_point]][[GO_type]] = go_plot_gse
    
  }
  
}

for (GO_type in GO_types) {
  plot_list = list()
  for (time_point in time_points) {
    plot_list[[time_point]] = enriched_GO_plot_list[[time_point]][[GO_type]]
  }
  grid.arrange(grobs = plot_list[1:3], ncol=3)
  
}


## ===> viz for GO enrichment approach <====
## extract all the biological processes out of the GO output for the GO
## enrichment approach
GO_type = 'BP'
GO_type = 'ALL'
GO_list_all = c()
GO_list_all_desc = c()
for (time_point in time_points) {
  GO_list = enriched_GO_list[[time_point]][[GO_type]]$ID
  GO_list_desc = enriched_GO_list[[time_point]][[GO_type]]$Description
  
  GO_list_all = c(GO_list_all, GO_list)
  GO_list_all_desc = c(GO_list_all_desc, GO_list_desc)
}

GO_list_all_unique = unique(GO_list_all)
GO_list_all_desc_unique = unique(GO_list_all_desc)

GO_matrix = matrix(1, nrow=length(GO_list_all_unique), ncol=length(time_points))
GO_df = as.data.frame(GO_matrix, row.names = GO_list_all_unique)
colnames(GO_df) = time_points
head(GO_df)

for (time_point in time_points) {
  index_values = match(rownames(GO_df), enriched_GO_list[[time_point]][[GO_type]]$ID)
  p_values = enriched_GO_list[[time_point]][[GO_type]][index_values]$p.adjust
  p_values[is.na(p_values)] = 1
  GO_df[,time_point] = p_values
  
}

rownames(GO_df) = GO_list_all_desc_unique
GO_df_log10 = -log10(GO_df)
hist(as.matrix(GO_df_log10))
GO_df_log10[GO_df_log10 > 10] = 10
GO_df_log10 = GO_df_log10[rowSums(GO_df_log10 < 2) < 6,]
pheatmap(GO_df_log10,
         cluster_cols = FALSE, 
         fontsize_row = 6)


## zoom in on one pathway
pathway = "cytokine-mediated signaling pathway"
pathway = "leukocyte degranulation"
pathway = "wound healing"

GO_id = GO_list_all_unique[grep(pathway, GO_list_all_desc_unique)]

genes_all = c()
for (time_point in time_points) {
  enrich_id = which(enriched_GO_list[[time_point]][[GO_type]]$ID == GO_id)
  print(time_point)
  if (length(enrich_id) > 0) {
    print(enrich_id)
    genes_cat = enriched_GO_list[[time_point]][[GO_type]][enrich_id]$geneID
    genes = strsplit(genes_cat, "/")[[1]]
    genes_all = c(genes_all, genes)
  }
}

genes_all_unique = unique(genes_all)
plot_data = deseq_data[match(genes_all_unique, deseq_data$external_gene_name),]
plot_data = plot_data[,grep("fc", colnames(plot_data))]
rownames(plot_data) = genes_all_unique
pheatmap(plot_data,
         breaks = seq(0,5, 0.05),
         fontsize_row = 8,
         cluster_cols = FALSE)




## ===> viz GSEA enrichment approach <====
## extract all the biological processes out of the GO output for the GO
## enrichment approach
GO_type = 'BP'
GO_type = 'ALL'
GO_list_all = c()
GO_list_all_desc = c()
for (time_point in time_points) {
  GO_list = gsego_list[[time_point]][[GO_type]]$ID
  GO_list_desc = gsego_list[[time_point]][[GO_type]]$Description
  
  GO_list_all = c(GO_list_all, GO_list)
  GO_list_all_desc = c(GO_list_all_desc, GO_list_desc)
}

GO_list_all_unique = unique(GO_list_all)
GO_list_all_desc_unique = unique(GO_list_all_desc)

GO_matrix = matrix(1, nrow=length(GO_list_all_unique), ncol=length(time_points))
GO_df = as.data.frame(GO_matrix, row.names = GO_list_all_unique)
colnames(GO_df) = time_points
head(GO_df)

for (time_point in time_points) {
  index_values = match(rownames(GO_df), gsego_list[[time_point]][[GO_type]]$ID)
  p_values = gsego_list[[time_point]][[GO_type]][index_values]$p.adjust
  p_values[is.na(p_values)] = 1
  GO_df[,time_point] = p_values
  
}

rownames(GO_df) = GO_list_all_desc_unique
GO_df_log10 = -log10(GO_df)
hist(as.matrix(GO_df_log10))
GO_df_log10[GO_df_log10 > 10] = 10
GO_df_log10 = GO_df_log10[rowSums(GO_df_log10 < 2) < 6,]
pheatmap(GO_df_log10,
         breaks = seq(0,2,0.02),
         cluster_cols = FALSE, 
         fontsize_row = 6)


# ---- >  zoom in on one pathway < ----
pathway = "cytokine-mediated signaling pathway"
pathway = "leukocyte degranulation"
pathway = "wound healing"

GO_id = GO_list_all_unique[grep(pathway, GO_list_all_desc_unique)]

genes_all = c()
for (time_point in time_points) {
  enrich_id = which(gsego_list[[time_point]][[GO_type]]$ID == GO_id)
  print(time_point)
  if (length(enrich_id) > 0) {
    print(enrich_id)
    genes_cat = gsego_list[[time_point]][[GO_type]][enrich_id]$geneID
    genes = strsplit(genes_cat, "/")[[1]]
    genes_all = c(genes_all, genes)
  }
}

genes_all_unique = unique(genes_all)
plot_data = deseq_data[match(genes_all_unique, deseq_data$external_gene_name),]
plot_data = plot_data[,grep("fc", colnames(plot_data))]
rownames(plot_data) = genes_all_unique
pheatmap(plot_data,
         breaks = seq(0,5, 0.05),
         fontsize_row = 8,
         cluster_cols = FALSE)


# ---- > REVIGO < -----
GO_types = c('BP')
for (time_point in time_points) {
  for (GO_type in GO_types) {
    # ego = enriched_GO_list[[time_point]][[GO_type]]
    ego = gsego_list[[time_point]][[GO_type]]
    go_results <- ego@result
    go_results = go_results[go_results$ONTOLOGY == GO_type,]
    
    ## extract GO IDs and their adjusted p-values
    go_terms <- go_results$ID
    p_values <- go_results$p.adjust
    
    # Ensure GO terms and p-values are in a data frame
    go_data <- data.frame(GO_ID = go_terms, p_value = p_values)
    head(go_data)  # Inspect the data

    ## set up semantic similarity computation
    simMatrix <- calculateSimMatrix(go_data$GO_ID,
                                    orgdb = "org.Mm.eg.db",  
                                    ont = GO_type,  
                                    method = "Wang")
    
    
    # Subset the scores vector to only include terms present in the similarity matrix
    sim_go_terms <- rownames(simMatrix)  # GO terms in the similarity matrix
    go_data_filtered <- go_data[go_data$GO_ID %in% sim_go_terms, ]
    
    scores = go_data_filtered$p_value
    names(scores) = go_data_filtered$GO_ID
    
    ## Reduce GO terms by clustering based on similarity
    reducedTerms <- reduceSimMatrix(simMatrix,
                                    scores = scores,  # P-values as scores for importance
                                    threshold = 0.7,  # Adjust this threshold for more/less reduction
                                    orgdb = "org.Mm.eg.db")
    
    # Inspect the reduced terms
    head(reducedTerms)
    dim(reducedTerms)
    
    # Scatterplot of the reduced GO terms
    scatter_p = scatterPlot(simMatrix, reducedTerms, size="size")
    scatter_file = paste0(src_dir, '/gsea/', 'scatter_', time_point, '_', GO_type, '.png')
    ggsave(filename = scatter_file, plot = scatter_p)
        
    # Treemap visualization of reduced GO terms
    # tree_p = treemapPlot(reducedTerms)
    # tree_file = paste0(src_dir, '/gsea/', 'treemap_', time_point, '_', GO_type, '.png')
    # ggsave(filename = tree_file, plot = tree_p)
    # 
    
  }
}  


# ----> Venn Diagram diff expressed genes <----
venn_list = list()
for (time_point in time_points) {
  diff_indices = which(deseq_data[,col_name_fc] >= cutoff_fc & deseq_data[,col_name_fdr] <= cutoff_fdr)
  ensembl_ids = deseq_data[diff_indices,]$ensembl_gene_id
  venn_list[[time_point]] = ensembl_ids
}


# Create the Venn diagram
ggVennDiagram(my_lists) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")

  
#### 2. GSEA for I vs NI approach ####

tissue = 'ear'
tissue = 'lymphnode'
src_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/deseq2/'
excel_file = paste0(src_dir, 'rnaseq_immune_response.',tissue,'.I_vs_NI.xlsx')
deseq_data = read.xlsx(excel_file)

time_points = c('4h', '12h', '24h', '48h', '72h', '96h')
cutoff_fdr = 0.05
cutoff_fc = 1

enriched_GO_list = list()
gsego_list = list()
enriched_GO_plot_list = list()
gsego_plot_list = list()

GO_types = c('BP', 'MF', 'CC')

for (time_point in time_points) {
  print(time_point)
  
  col_name_fc = paste0(time_point, '_fc')
  col_name_fdr = paste0(time_point, '_fdr')
  
  ## select the differentially expressed genes
  diff_indices = which(deseq_data[,col_name_fc] >= cutoff_fc & deseq_data[,col_name_fdr] <= cutoff_fdr)
  ensembl_ids = deseq_data[diff_indices,]$ensembl_gene_id
  
  ## select all ensembl IDs for the gene set enrichment (for gseGO)
  ensembl_ids_fc = deseq_data[, c('ensembl_gene_id', col_name_fc)]
  ensembl_ids_fc_sorted = ensembl_ids_fc[order(ensembl_ids_fc[, col_name_fc], decreasing = TRUE),]
  ensembl_ids_sorted = ensembl_ids_fc_sorted$ensembl_gene_id
  
  ## extract ensembl IDs and convert them to entrez IDs
  entrez_ids <- mapIds(org.Mm.eg.db, keys = ensembl_ids, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
  entrez_ids_sorted <- mapIds(org.Mm.eg.db, keys = ensembl_ids_sorted, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
  
  geneList = ensembl_ids_fc_sorted[, col_name_fc]
  names(geneList) = entrez_ids_sorted
  geneList <- geneList[!is.na(names(geneList))]
  
  
  for (GO_type in GO_types) {
    print(GO_type)
    
    ## GO enrichment analysis
    go_enrich <- enrichGO(gene = entrez_ids, 
                          OrgDb = org.Mm.eg.db,  
                          ont = GO_type,  ## vary over BP, CC, or MF
                          pAdjustMethod = "BH", 
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.2, 
                          readable = TRUE,
                          minGSSize = 20)
    
    # go_gse <- gseGO(gene = geneList, 
    #                       OrgDb = org.Mm.eg.db,  
    #                       ont = "BP",  ## vary over BP, CC, or MF
    #                       pAdjustMethod = "BH", 
    #                       pvalueCutoff = 0.05) 
    # dim(go_gse)
    
    ## filter the go_enrich output
    go_enrich_simplified <- simplify(go_enrich, cutoff = 0.5, by = "p.adjust", select_fun = min)
    dim(go_enrich_simplified)
    
    # go_gse_simplified <- simplify(go_gse, cutoff = 0.5, by = "p.adjust", select_fun = min)
    # dim(go_gse_simplified)
    
    
    
    
    # Visualize GO enrichment results
    go_plot_enriched = dotplot(go_enrich_simplified, 
                               showCategory = 25,
                               font.size = 10,
    )
    
    
    # go_plot_gse = dotplot(go_gse_simplified, 
    #         showCategory = 25,
    #         font.size = 10,
    # )
    
    
    go_plot_file_enriched = paste0(src_dir, 'gsea/', 'I_vs_NI_', tissue, '_', time_point, '_', GO_type, '_enricheGO.png')
    # go_plot_file_gse = paste0(src_dir, 'gsea/', 'I_vs_NI_', tissue, '_', time_point, '_', GO_type, '_gsego.png')
    
    ggsave(go_plot_enriched, file = go_plot_file_enriched)
    # ggsave(go_plot_gse, file = go_plot_file_gse)
    
    
    enriched_GO_list[[time_point]][[GO_type]] = go_enrich_simplified
    # gsego_list[[time_point]][[GO_type]] = go_gse_simplified
    enriched_GO_plot_list[[time_point]][[GO_type]] = go_plot_enriched
    # gsego_plot_list[[time_point]][[GO_type]] = go_plot_gse
    
  }
  
}




