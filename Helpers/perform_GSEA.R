gsea_barplot_multi <- function(
    df,                         # limma topTable with Symbol, logFC, P.Value
    m_t2g = NULL,               # Hallmark TERM2GENE-like df (gs_name, gene_symbol OR term,gene)
    use = c("H", "GO:BP"),      # which sources to run
    per_source = TRUE,          # pick top/bottom within each source
    n_each = 20,                # # of top & bottom terms by NES
    add_source_prefix = FALSE,  # keep FALSE to avoid adding "Hallmark: " / "GO:BP: "
    plot_title = "GSEA - Ranked by NES",
    seed = 123,
    p_label = 0.05,
    palette = c("#00468BFF", "#ED0000FF"),
    # gseGO controls:
    keyType = "SYMBOL",
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pvalueCutoff = 1,
    # NEW: wrapping control
    wrap_width = 35
) {
  suppressPackageStartupMessages({
    library(dplyr); library(ggplot2); library(ggpubr)
    library(clusterProfiler); library(org.Hs.eg.db)
    library(stringr)
  })
  
  # ---- checks ----
  needed <- c("Symbol", "logFC", "P.Value")
  if (!all(needed %in% names(df)))
    stop("Input `df` must contain: ", paste(needed, collapse = ", "))
  
  use <- unique(use)
  
  # ---- helpers ----
  std_t2g <- function(x, source_name) {
    if (all(c("gs_name","gene_symbol") %in% names(x))) {
      out <- x %>% transmute(term = gs_name, gene = gene_symbol)
    } else if (all(c("term","gene") %in% names(x))) {
      out <- x %>% dplyr::select(term, gene)
    } else {
      stop("TERM2GENE for ", source_name, " must have (gs_name, gene_symbol) or (term, gene).")
    }
    if (add_source_prefix) out <- out %>% mutate(term = paste0(source_name, ": ", term))
    out
  }
  
  rank_gene_list <- function(df) {
    df %>%
      filter(!is.na(Symbol), !is.na(logFC), !is.na(P.Value), P.Value > 0) %>%
      mutate(score = sign(logFC) * -log10(P.Value)) %>%
      group_by(Symbol) %>% summarize(score = mean(score), .groups = "drop") %>%
      arrange(desc(score)) %>%
      { setNames(.$score, .$Symbol) }
  }
  
  # (kept for completeness; now wraps labels too if used)
  build_plot_tbl <- function(gsea_res, plot_title, n_each, p_label, palette) {
    dfres <- as.data.frame(gsea_res)
    if (nrow(dfres) == 0) {
      return(list(
        plot  = ggplot() + theme_void() + labs(title = paste0(plot_title, " (no enriched terms)")),
        table = tibble::tibble()
      ))
    }
    label_col <- if ("Description" %in% names(dfres)) "Description" else "ID"
    
    nes <- dfres %>%
      filter(!is.na(NES)) %>%
      dplyr::select(ID, NES, p.adjust, !!label_col) %>%
      mutate(
        Label = .data[[label_col]],
        Label = gsub("_", " ", Label),
        Label = stringr::str_to_title(Label),
        Label = ifelse(p.adjust < p_label, paste0(Label, " *"), Label)
      ) %>%
      arrange(desc(NES))
    
    top_n <- nes %>% slice_head(n = n_each)
    bot_n <- nes %>% slice_tail(n = n_each)
    combined <- bind_rows(top_n, bot_n) %>% distinct(ID, .keep_all = TRUE)
    
    plt <- ggplot(combined, aes(x = reorder(Label, NES), y = NES)) +
      geom_col(aes(fill = p.adjust), width = 0.8) +
      scale_fill_gradient(low = palette[1], high = palette[2], guide = "colorbar") +
      coord_flip() +
      scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = wrap_width)) +
      labs(x = "Pathway", y = "Normalized Enrichment Score",
           title = plot_title, fill = "Adjusted p-value") +
      theme_pubr() +
      theme(panel.spacing = grid::unit(0.1, "lines"),
            legend.position = "right",
            axis.text.y = element_text(lineheight = 0.9))
    
    list(plot = plt, table = combined)
  }
  
  # ---- rank genes ----
  gene_list <- rank_gene_list(df)
  
  # ---- run per source ----
  sources <- character(0)
  if ("H" %in% use) {
    if (is.null(m_t2g)) stop("Hallmark selected in `use`, but `m_t2g` is NULL.")
    sources <- c(sources, "Hallmark")
  }
  if ("GO:BP" %in% use) {
    sources <- c(sources, "GO_BP")
  }
  if (!length(sources)) stop("Nothing to run. Set `use` to include 'H' and/or 'GO:BP'.")
  
  set.seed(seed)
  gsea_by_src <- lapply(sources, function(src) {
    if (src == "Hallmark") {
      TERM2GENE <- std_t2g(m_t2g, "Hallmark")
      res <- GSEA(geneList = gene_list, TERM2GENE = TERM2GENE,
                  pvalueCutoff = 1, seed = TRUE)
    } else { # GO_BP via gseGO
      res <- gseGO(
        geneList     = gene_list,
        ont          = ont,
        keyType      = keyType,
        OrgDb        = OrgDb,
        pvalueCutoff = pvalueCutoff,
        seed         = TRUE
      )
    }
    list(source = src, res = res)
  })
  names(gsea_by_src) <- sources
  
  # ---- combine results ----
  df_all <- bind_rows(lapply(gsea_by_src, function(x) {
    if (is.null(x$res)) return(NULL)
    as.data.frame(x$res) %>% mutate(Source = x$source)
  }))
  
  if (is.null(df_all) || nrow(df_all) == 0) {
    empty_plot <- ggplot() + theme_void() +
      labs(title = paste0(plot_title, " (no enriched terms)"))
    return(list(plot = empty_plot, table = tibble::tibble(), gsea = gsea_by_src))
  }
  
  # ---- labels (no source prefixes added to terms) ----
  nes_tbl <- df_all %>%
    dplyr::select(Source, ID, NES, p.adjust,
                  Description = dplyr::any_of("Description")) %>%
    mutate(
      Label = if (!all(is.na(Description))) Description else ID,
      Label = gsub("_", " ", Label),
      Label = stringr::str_to_title(Label),
      Label = ifelse(p.adjust < p_label, paste0(Label, " *"), Label)
    ) %>%
    dplyr::select(-Description)
  
  # ---- choose top/bottom by NES ----
  if (per_source) {
    pick_tb <- function(d) {
      bind_rows(
        d %>% arrange(desc(NES)) %>% slice_head(n = n_each),
        d %>% arrange(NES)       %>% slice_head(n = n_each)
      )
    }
    combined <- nes_tbl %>%
      group_by(Source) %>%
      do(pick_tb(.)) %>%
      ungroup() %>%
      distinct(Source, ID, .keep_all = TRUE)
  } else {
    combined <- bind_rows(
      nes_tbl %>% arrange(desc(NES)) %>% slice_head(n = n_each),
      nes_tbl %>% arrange(NES)       %>% slice_head(n = n_each)
    ) %>% distinct(Source, ID, .keep_all = TRUE)
  }
  
  # ---- plot (with wrapped labels) ----
  facet_needed <- (length(unique(combined$Source)) > 1)
  plt <- ggplot(combined, aes(x = reorder(Label, NES), y = NES)) +
    geom_col(aes(fill = p.adjust), width = 0.8) +
    scale_fill_gradient(low = palette[1], high = palette[2], guide = "colorbar") +
    coord_flip() +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = wrap_width)) +
    labs(x = "Pathway", y = "Normalized Enrichment Score",
         title = plot_title, fill = "Adjusted p-value") +
    theme_pubr() +
    theme(panel.spacing = grid::unit(0.1, "lines"),
          legend.position = "right",
          #axis.text.y = element_text(lineheight = 0.9),
          axis.text.y     = element_text(size = 9, lineheight = 0.8),  # pathway labels
          axis.text.x     = element_text(size = 9))
  
  if (facet_needed) {
    plt <- plt + facet_grid(Source ~ ., scales = "free_y", space = "free_y")
  }
  
  return(list(plot = plt, table = combined, gsea = gsea_by_src))
}
