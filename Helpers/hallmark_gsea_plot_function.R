library(tidyverse)
library(clusterProfiler)
library(ggpubr)
hallmark_gsea_plot <- function(gsea_object, top_n = 10, low_color = "#00468BFF", high_color = "#ED0000FF", title = "GSEA Plot") {
  # Get the name of the input gsea_object as a string
  gsea_name <- deparse(substitute(gsea_object))
  
  # Extract relevant columns
  nes <- gsea_object[, c("ID","NES", "p.adjust")]
  
  # Format and clean pathway IDs
  nes <- nes %>%
    mutate(
      ID = gsub("HALLMARK_", "", ID),       # Remove "HALLMARK_" prefix
      ID = gsub("_", " ", ID),             # Replace underscores with spaces
      ID = str_to_sentence(ID),
      ID = ifelse(p.adjust < 0.05, paste0(ID, " *"), ID) # Mark significant pathways
    )
  
  # Sort by NES
  nes_sorted <- nes %>% arrange(desc(NES))
  
  # Select top and bottom pathways based on NES
  top_n_nes <- nes_sorted %>% slice_head(n = top_n)
  bottom_n_nes <- nes_sorted %>% slice_tail(n = top_n)
  
  # Combine top and bottom pathways, ensuring no duplicates
  combined_nes <- bind_rows(
    top_n_nes %>% mutate(Category = "Top"),
    bottom_n_nes %>% mutate(Category = "Bottom")
  ) %>% distinct(ID, .keep_all = TRUE)
  
  # Dynamically assign combined_nes to a global variable using the GSEA object name
  assign(paste0(gsea_name, "_combined_nes"), combined_nes, envir = .GlobalEnv)
  
  # Create the plot
  plot <- ggplot(combined_nes, aes(reorder(ID, NES), NES)) +
    geom_col(aes(fill = p.adjust), width = 0.8) +   # Fill by p.adjust
    scale_fill_gradient(low = low_color, high = high_color, guide = "colorbar") +
    coord_flip() +
    labs(
      x = "Pathway",
      y = "Normalized Enrichment Score (NES)",
      title = title,
      fill = "Adjusted p-value"
    ) +
    theme_pubr() +
    theme(
      panel.spacing = unit(0.1, "lines"),
      legend.position = "right",
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.caption = element_text(size = 9),
      legend.text = element_text(size = 9),
      plot.title = element_text(size = 11)
    )
  
  return(plot)
}

go_gsea_plot <- function(gsea_object, top_n = 10, low_color = "#00468BFF", high_color = "#ED0000FF", title = "GSEA Plot", line_break_after = 20) {
  # Get the name of the input gsea_object as a string
  gsea_name <- deparse(substitute(gsea_object))
  
  # Extract relevant columns
  nes <- gsea_object@result
  nes$ID <- nes$Description
  
  nes <- nes[, c("ID","NES", "p.adjust")]
  # Format and clean pathway IDs
  nes <- nes %>%
    mutate(
      ID = gsub("_", " ", ID),
      ID = str_to_sentence(ID), 
      ID = ifelse(p.adjust < 0.05, paste0(ID, " *"), ID) # Mark significant pathways
    )
  
  nes$ID <- sapply(nes$ID, function(x) {
    paste(strwrap(x, width = line_break_after), collapse = "\n")
  })
  # Sort by NES
  nes_sorted <- nes %>% arrange(desc(NES))
  
  # Select top and bottom pathways based on NES
  top_n_nes <- nes_sorted %>% slice_head(n = top_n)
  bottom_n_nes <- nes_sorted %>% slice_tail(n = top_n)
  
  # Combine top and bottom pathways, ensuring no duplicates
  combined_nes <- bind_rows(
    top_n_nes %>% mutate(Category = "Top"),
    bottom_n_nes %>% mutate(Category = "Bottom")
  ) %>% distinct(ID, .keep_all = TRUE)
  
  # Dynamically assign combined_nes to a global variable using the GSEA object name
  assign(paste0(gsea_name, "_combined_nes"), combined_nes, envir = .GlobalEnv)
  
  # Create the plot
  plot <- ggplot(combined_nes, aes(reorder(ID, NES), NES)) +
    geom_col(aes(fill = p.adjust), width = 0.8) +   # Fill by p.adjust
    scale_fill_gradient(low = low_color, high = high_color, guide = "colorbar") +
    coord_flip() +
    labs(
      x = "Pathway",
      y = "Normalized Enrichment Score (NES)",
      title = title,
      fill = "Adjusted p-value"
    ) +
    theme_pubr() +
    theme(
      panel.spacing = unit(0.1, "lines"),
      legend.position = "right",
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.caption = element_text(size = 9),
      legend.text = element_text(size = 9),
      plot.title = element_text(size = 11)
    )
  
  return(plot)
}
