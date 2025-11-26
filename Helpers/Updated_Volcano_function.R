plot_volcano <- function(data, plot_title,
                         output_file = NULL, height = 6, width = 9,
                         p_cut = 0.05, fc_cut = 1.2) {
  
  # ---- deps ----
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(EnhancedVolcano)
    library(ggsci)
  })
  
  # ---- prepare data safely ----
  df <- data
  
  # ensure rownames exist and are unique (needed by colCustom)
  if (is.null(rownames(df))) {
    if (!"Symbol" %in% names(df)) {
      df$Symbol <- paste0("gene_", seq_len(nrow(df)))
    }
    rn <- make.unique(ifelse(is.na(df$Symbol) | df$Symbol == "",
                             paste0("gene_", seq_len(nrow(df))),
                             df$Symbol))
    rownames(df) <- rn
  }
  
  # ---- colors: grey for NS, red for up, blue for down (only when both cutoffs pass) ----
  lancet <- ggsci::pal_lancet()(3)
  blue   <- lancet[1]  # "#00468BFF"
  red    <- lancet[2]  # "#ED0000FF"
  
  col_vec <- rep("grey70", nrow(df))
  both_sig <- (df$P.Value < p_cut) & (abs(df$logFC) >= fc_cut)
  col_vec[both_sig & df$logFC > 0] <- red
  col_vec[both_sig & df$logFC < 0] <- blue
  names(col_vec) <- rownames(df)
  
  # ---- labels (top by raw P; guard against NAs/empty) ----
  top_genes <- df %>%
    filter(!is.na(P.Value), P.Value < p_cut) %>%
    arrange(P.Value) %>%
    slice_head(n = 30) %>%
    pull(Symbol) %>%
    { .[!is.na(.) & nzchar(.)] }
  
  lab_italics       <- paste0("italic('", df$Symbol, "')")
  selectLab_italics <- if (length(top_genes)) paste0("italic('", top_genes, "')") else NULL
  draw_connect      <- !is.null(selectLab_italics) && length(selectLab_italics) > 0
  
  # ---- plot ----
  p <- EnhancedVolcano(
    df,
    subtitle = "",
    lab = lab_italics,
    parseLabels = TRUE,
    selectLab = selectLab_italics,
    x = "logFC",
    y = "P.Value",
    title = plot_title,
    xlab = bquote(~Log[2]~"fold change"),
    ylab = expression(-log[10](P[val])),
    caption = paste(
      "Total significant genes (Padj): ",
      sum(data$adj.P.Val < 0.05), ", Up: ",
      sum(data$adj.P.Val < 0.05 & data$logFC > 0), " Down: ",
      sum(data$adj.P.Val < 0.05 & data$logFC < -0)
    ),
    pCutoff = p_cut,
    FCcutoff = fc_cut,
    cutoffLineType = "twodash",
    cutoffLineWidth = 0.8,
    pointSize = 1.0,
    labSize = 3.0,
    colAlpha = 0.50,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    legendPosition = "none",     
    drawConnectors = draw_connect,
    widthConnectors = 0.5,
    max.overlaps = 10,
    colCustom = col_vec          
  ) +
    theme(
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.caption = element_text(size = 9),
      plot.title = element_text(size = 11)
    )
  
  
  
  if (!is.null(output_file)) {
    ggplot2::ggsave(filename = output_file, plot = p,
                    height = height, width = width,
                    dpi = 300,
                    bg = "transparent")
  }
  
  p
}
