perform_enrichment_analysis <- function(gene_list, universe_genes, subset_name, output_prefix,
                                        significant_threshold = 0.05, m_t2g = NULL,
                                        do_GO = TRUE, species = "Homo sapiens", height = 5, width = 7) {
  # Load required OrgDb based on species
  if (species == "Homo sapiens") {
    OrgDb <- org.Hs.eg.db
  } else if (species == "Mus musculus") {
    OrgDb <- org.Mm.eg.db
  } else {
    stop("Unsupported species. Use 'Homo sapiens' or 'Mus musculus'.")
  }
  
  # Load default Hallmark gene sets if m_t2g not provided
  if (is.null(m_t2g)) {
    m_t2g <- msigdbr(species = species, collection = "H") %>%
      mutate(gs_name = gsub("HALLMARK_", "", gs_name)) %>%
      dplyr::select(gs_name, gene_symbol)
  }
  
  # GO:BP Enrichment
  if (do_GO) {
    enrich_go <- enrichGO(
      gene          = gene_list,
      OrgDb         = OrgDb,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      universe      = universe_genes,
      pvalueCutoff  = 1,
      qvalueCutoff  = 1
    )
  } else {
    enrich_go <- NULL
  }
  
  # Hallmark (or custom) enrichment
  enrich_hallmark <- enricher(
    gene          = gene_list,
    TERM2GENE     = m_t2g,
    pAdjustMethod = "BH",
    universe      = universe_genes,
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  )
  
  # Format function
  format_descriptions <- function(enrich_result, threshold) {
    enrich_result@result$Description <- ifelse(
      enrich_result@result$p.adjust < threshold,
      paste0(enrich_result@result$Description, " *"),
      enrich_result@result$Description
    )
    return(enrich_result)
  }
  
  if (!is.null(enrich_go)) enrich_go <- format_descriptions(enrich_go, significant_threshold)
  enrich_hallmark <- format_descriptions(enrich_hallmark, significant_threshold)
  
  # Plot GO if done
  if (!is.null(enrich_go)) {
    go_plot <- dotplot(enrich_go) + ggtitle(paste(subset_name, "GO:BP")) + theme_pubr() +
      theme(legend.position = "right",
            axis.text.x = element_text(size = 9),
            axis.text.y = element_text(size = 9),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            plot.caption = element_text(size = 9),
            legend.text = element_text(size = 9),
            plot.title = element_text(size = 11))
    
    ggsave(filename = paste0(output_prefix, "_GO_BP.svg"), height = 5, width = 7, bg = "transparent",
           plot = go_plot)
  }
  
  hallmark_plot <- dotplot(enrich_hallmark) + ggtitle(paste(subset_name, "Hallmark")) + theme_pubr() +
    theme(legend.position = "right",
          axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 9),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          plot.caption = element_text(size = 9),
          legend.text = element_text(size = 9),
          plot.title = element_text(size = 11))
  
  ggsave(filename = paste0(output_prefix, "_Hallmark.svg"), height = height, width = width, bg = "transparent",
         plot = hallmark_plot)
  
  return(list(GO_BP = enrich_go, Hallmark = enrich_hallmark))
}
