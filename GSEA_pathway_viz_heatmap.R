library(GO.db)
library(org.Mm.eg.db)
library(pheatmap)

# Your list of GO description keywords
pathways <- c("extracellular matrix",
              "neutrophil migration",
              # "ERK1/2 cascade",
              "ERK1 and ERK2 cascade",
              "response to interferon-beta",
              "myeloid leukocyte migration",
              "defense response to virus",
              "adaptive immune response",
              "antigen processing and presentation")

src_dir = '/Users/pmonsieurs/programming/trypanosoma_rnaseq_bitesite/results/deseq2/'


# Get all GO IDs with their terms
go_terms <- as.list(GOTERM)
go_id_to_term <- sapply(go_terms, Term)

for (pathway in pathways) {
  cat("\n➡️ Processing pathway:", pathway, "\n")
  
  # Match GO descriptions
  matched_go_ids <- names(go_id_to_term[grepl(pathway, go_id_to_term, ignore.case = TRUE)])
  
  if (length(matched_go_ids) == 0) {
    cat("❌ No GO term match for:", pathway, "\n")
    next
  }
  
  all_genes <- c()
  
  # Only keep valid GO IDs that exist in org.Mm.eg.db
  valid_go_ids <- intersect(matched_go_ids, keys(org.Mm.eg.db, keytype = "GOALL"))
  
  for (go_id in valid_go_ids) {
    entrez_df <- tryCatch(
      AnnotationDbi::select(org.Mm.eg.db,
                            keys = go_id,
                            columns = "ENTREZID",
                            keytype = "GOALL"),
      error = function(e) NULL
    )
    if (is.null(entrez_df) || nrow(entrez_df) == 0) next
    
    gene_df <- AnnotationDbi::select(org.Mm.eg.db,
                                     keys = unique(entrez_df$ENTREZID),
                                     columns = "SYMBOL",
                                     keytype = "ENTREZID")
    all_genes <- c(all_genes, gene_df$SYMBOL)
  }
  
  genes_of_interest <- unique(na.omit(all_genes))
  
  # Remove duplicated gene names first
  deseq_unique <- deseq_data[!duplicated(deseq_data$external_gene_name), ]
  
  # Subset only genes of interest
  plot_data <- deseq_unique[deseq_unique$external_gene_name %in% genes_of_interest, ]
  
  # Filter FC columns
  plot_data_fc <- plot_data[, grep("_fc$", colnames(plot_data))]
  
  # Set unique rownames
  rownames(plot_data_fc) <- plot_data$external_gene_name
  
  # Skip empty
  if (nrow(plot_data_fc) == 0) {
    cat("⚠️ No matching genes in expression data for:", pathway, "\n")
    next
  }
  
  # Draw heatmap
  p_heat = pheatmap(plot_data_fc,
           main = pathway,
           fontsize_row = 8,
           fontsize_col = 24,
           fontsize = 28,   
           cluster_cols = FALSE,
           breaks = seq(0, 5, 0.05))
  
  png_file = paste0(src_dir, 'pathway_heatmaps/', pathway, '.png')
  ggsave(filename = png_file, plot = p_heat, width=32, height=18)
}
