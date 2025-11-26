library(tidyverse)
library(msigdbr)
library(ggpubr)

thp1_top.treat <- ko_wt_m1_top.treat
mgi <- read.csv2("C:/Users/ssa214/UiT Office 365/O365-PhD Saikat - General/Proteomics/Mouse_LPS/Mouse_Human_MGI.csv")


create_mapped_dataframe <- function(genes_vector, mgi) {
  # Filter the `mgi` data frame to keep only matching mouse gene symbols
  mapped_df <- mgi[mgi$gene_symbol %in% genes_vector, ]
  return(mapped_df)
}


mouse_upregulated <- mouse_top.treat %>% filter(adj.P.Val < 0.05 & logFC > 0) %>% pull(Symbol)
mouse_downregulated <- mouse_top.treat %>% filter(adj.P.Val < 0.05 & logFC < 0) %>% pull(Symbol)

mouse_upregulated_df <- create_mapped_dataframe(mouse_upregulated, mgi)
mouse_downregulated_df <- create_mapped_dataframe(mouse_downregulated, mgi)
mouse_merged <- rbind(mouse_upregulated_df, mouse_downregulated_df)
universe_mouse_df <- create_mapped_dataframe(rownames(lcpm), mgi)

common <- unique(c(rownames(lcpm), universe_mouse_df$human_gene_symbol))

THP1_upregulated <- thp1_top.treat %>% filter(adj.P.Val < 0.05 & logFC > 0) %>% pull(Symbol)
THP1_downregulated <- thp1_top.treat %>% filter(adj.P.Val < 0.05 & logFC < 0) %>% pull(Symbol)


Up_list <- list(Mouse_Up = mouse_upregulated_df$human_gene_symbol,
                THP1_Up = THP1_upregulated)

Down_list <- list(Mouse_Down = mouse_downregulated_df$human_gene_symbol,
                  THP1_Down = THP1_downregulated)

library(VennDetail)
library(ggpubr)
mouse_merged$X <- NULL
mouse_merged$DB.Class.Key <- NULL
mouse_df <- merge(mouse_top.treat, mouse_merged, by.x = "Symbol",by.y = "gene_symbol")
#mouse_df_duplicates <- mouse_df[mouse_df$Symbol %in% mouse_df$Symbol[duplicated(mouse_df$Symbol)], ]

mouse_df_subset <- mouse_df[,c(7,2,6)]
names(mouse_df_subset) <- c("Symbol", "Mouse_logFC", "Mouse_adj.P.Val")
mouse_df_subset_unique <- mouse_df_subset[!duplicated(mouse_df_subset$Symbol), ]

thp1_df_subset <- thp1_top.treat[,c(6,1,5)]
names(thp1_df_subset) <- c("Symbol", "THP1_logFC", "THP1_adj.P.Val")

Up_plot <- venndetail(
  Up_list)
Up_df <- result(Up_plot)
Up_mouse <- mouse_df_subset %>% filter(Symbol %in% Up_df$Detail)
Up_df_added <- merge(Up_df, Up_mouse, by.x = "Detail", by.y = "Symbol", all.x = T) 
Up_df_added <- merge(Up_df_added, thp1_df_subset, by.x = "Detail", by.y = "Symbol", all.x = T)
Up_genes_merged <- Up_df_added %>% distinct(Detail, .keep_all = T)

Down_plot <- venndetail(
  Down_list)
Down_df <- result(Down_plot)

Down_mouse <- mouse_df_subset %>% filter(Symbol %in% Down_df$Detail)
Down_df_added <- merge(Down_df, Down_mouse, by.x = "Detail", by.y = "Symbol", all.x = T) 
Down_df_added <- merge(Down_df_added, thp1_df_subset, by.x = "Detail", by.y = "Symbol", all.x = T)
Down_genes_merged <- Down_df_added %>% distinct(Detail, .keep_all = T)




library(svglite)
ggsave(filename = "Mouse 4Hours vs THP1 M1 Upregulated Venn_Treat.svg", height = 5, width = 7, bg = "transparent",
       plot = plot(Up_plot, percentage = TRUE, mycol = c("#00468BFF", 
                                                         "#ED0000FF"),order = FALSE,
                   alpha = 0.3,
                   gp = gpar(fontsize = 9)))

ggsave(filename = "Mouse 4Hours vs THP1 M1 Downregulated Venn_Treat.svg", height = 5, width = 7, bg = "transparent",
       plot = plot(Down_plot, percentage = TRUE, mycol = c("#00468BFF", 
                                                           "#ED0000FF"),order = FALSE,
                   alpha = 0.3,
                   gp = gpar(fontsize = 10)))





m_t2g_tnfa <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  mutate(gs_name=gsub("HALLMARK_","",gs_name)) %>%
  dplyr::select(gs_name,gene_symbol) %>% filter(gs_name == "TNFA_SIGNALING_VIA_NFKB") %>%
  pull(gene_symbol)




mouse_upregulated <- mouse_top.treat %>% filter(adj.P.Val < 0.05 & logFC > 0) %>% pull(Symbol)
mouse_downregulated <- mouse_top.treat %>% filter(adj.P.Val < 0.05 & logFC < 0) %>% pull(Symbol)

mouse_upregulated_df <- create_mapped_dataframe(mouse_upregulated, mgi)
mouse_downregulated_df <- create_mapped_dataframe(mouse_downregulated, mgi)

universe_mouse_df <- create_mapped_dataframe(rownames(lcpm), mgi)

common <- unique(c(rownames(lcpm), universe_mouse_df$human_gene_symbol))

THP1_upregulated <- thp1_top.treat %>% filter(adj.P.Val < 0.05 & logFC > 0) %>% pull(Symbol)
THP1_downregulated <- thp1_top.treat %>% filter(adj.P.Val < 0.05 & logFC < 0) %>% pull(Symbol)


Up_list <- list(Mouse_Up = intersect(m_t2g_tnfa, mouse_upregulated_df$human_gene_symbol),
                THP1_Up = intersect(m_t2g_tnfa, THP1_upregulated))

Down_list <- list(Mouse_Down = intersect(m_t2g_tnfa, mouse_downregulated_df$human_gene_symbol),
                  THP1_Down = intersect(m_t2g_tnfa, THP1_downregulated))

library(VennDetail)
library(ggpubr)
Up_plot <- venndetail(
  Up_list)
Up_df <- result(Up_plot)
Up_mouse <- mouse_df_subset %>% filter(Symbol %in% Up_df$Detail)
Up_df_added <- merge(Up_df, Up_mouse, by.x = "Detail", by.y = "Symbol", all.x = T) 
Up_df_added <- merge(Up_df_added, thp1_df_subset, by.x = "Detail", by.y = "Symbol", all.x = T)
Up_genes_merged_tnfa <- Up_df_added %>% distinct(Detail, .keep_all = T)


Down_plot <- venndetail(
  Down_list)

Down_df <- result(Down_plot)
Down_mouse <- mouse_df_subset %>% filter(Symbol %in% Down_df$Detail)
Down_df_added <- merge(Down_df, Down_mouse, by.x = "Detail", by.y = "Symbol", all.x = T) 
Down_df_added <- merge(Down_df_added, thp1_df_subset, by.x = "Detail", by.y = "Symbol", all.x = T)
Down_genes_merged_tnfa <- Down_df_added %>% distinct(Detail, .keep_all = T)


library(svglite)
ggsave(filename = "Mouse 4Hours vs THP1 M1 TNFA Upregulated Venn_Treat.svg", height = 5, width = 7, bg = "transparent",
       plot = plot(Up_plot, percentage = TRUE, mycol = c("#00468BFF", 
                                                         "#ED0000FF"),order = FALSE,
                   alpha = 0.3,
                   gp = gpar(fontsize = 9)))

ggsave(filename = "Mouse 4Hours vs THP1 M1 TNFA Downregulated Venn_Treat.svg", height = 5, width = 7, bg = "transparent",
       plot = plot(Down_plot, percentage = TRUE, mycol = c("#00468BFF", 
                                                           "#ED0000FF"),order = FALSE,
                   alpha = 0.3,
                   gp = gpar(fontsize = 10)))


up_tnfa <- Up_df
down_tnfa <- Down_df

tnfa_lcpm <- as.data.frame(lcpm)[intersect(m_t2g_tnfa, rownames(lcpm)),]
tnfa_lcpm <- rownames_to_column(tnfa_lcpm, var = "Symbol")
tnfa_thp1_de <- ko_wt_m1_top.treat %>% filter(Symbol %in% m_t2g_tnfa)

up_tnfa
dataframes <- list("Upregulated_All" = Up_genes_merged,
                   "Downregulated_All" = Down_genes_merged,
                   "Upregulated_TNFA" = Up_genes_merged_tnfa,
                   "Downregulated_TNFA" = Down_genes_merged_tnfa)
write.xlsx(dataframes, file= "Venn Diagrams Mouse and THP1_Treat.xlsx")


dataframes2 <- list("TNFA_DE_Treat_THP1" = tnfa_thp1_de,
                    "TNFA_LogCPM" = tnfa_lcpm)

write.xlsx(dataframes2, file= "TNFA Genes in THP1 M1.xlsx")

####Common ORA
Up_df_common <- Up_df %>% filter(Subset == "Shared") %>% pull(Detail)
Down_df_common <- Down_df %>% filter(Subset == "Shared") %>% pull(Detail)


m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  mutate(gs_name=gsub("HALLMARK_","",gs_name)) %>%
  dplyr::select(gs_name,gene_symbol)

up_enrichment <- perform_enrichment_analysis(gene_list = Up_df_common, universe_genes = common,species = "Homo sapiens",
                                             subset_name = "Common Upregulated - Mouse and THP1",
                                             output_prefix = "Common Upregulated from Venn_Mouse and THP1")

down_enrichment <- perform_enrichment_analysis(gene_list = Down_df_common, universe_genes = common,species = "Homo sapiens",
                                               subset_name = "Common Downregulated - Mouse and THP1",
                                               output_prefix = "Common Downregulated from Venn_Mouse and THP1")






library(nichenetr)


mouse_de <- mouse_top.treat
mouse_df <- merge(mouse, mouse_de, by.x = "gene_symbol", by.y = "Symbol")

load(file="C:/Users/ssa214/UiT Office 365/O365-PhD Saikat - General/Proteomics/20250504_THP1/DE_KO_WT_M1_Treat.RData")

thp1 <- ko_wt_m1_top.treat



m_t2g_chol <- msigdbr(species = "Mus musculus", collection = "H") %>%
  mutate(gs_name = gsub("HALLMARK_", "", gs_name)) %>%
  dplyr::select(gs_name, gene_symbol)


tnfa_thp1 <- m_t2g %>% 
  filter(gs_name == "TNFA_SIGNALING_VIA_NFKB") %>%
  pull(gene_symbol)

tnfa_thp1_df <- thp1 %>% filter(Symbol %in% tnfa_thp1)

tnfa_thp1_df_sig <- tnfa_thp1_df %>% filter(adj.P.Val < 0.05)


tnfa_mouse <- m_t2g %>% 
  filter(gs_name == "TNFA_SIGNALING_VIA_NFKB") %>%
  pull(gene_symbol)

tnfa_mouse_df <- mouse %>% filter(Symbol %in% tnfa_mouse)

tnfa_mouse_df_sig <- tnfa_mouse_df %>% filter(adj.P.Val < 0.05)

mouse_f <- mouse %>% filter(gene_symbol %in% tnfa_mouse_df_sig$Symbol) %>% select(gene_symbol, human_gene_symbol)
mouse_df <- merge(tnfa_mouse_df_sig, mouse_f, by.x = "Symbol", by.y = "gene_symbol")
mouse_df$mouse_symbol <- mouse_df$Symbol
mouse_df$Symbol <- mouse_df$human_gene_symbol



mac_merged <- inner_join(tnfa_thp1_df_sig, mouse_df, by = "Symbol", suffix = c("_THP1", "_Mouse"))
mac_merged <- mac_merged %>%
  mutate(direction_match = sign(logFC_THP1) == sign(logFC_Mouse))


features_list <- list(
  TNFA_Genes = sort(mac_merged$Symbol)
)



long_df <- mac_merged %>%
  select(Symbol, logFC_THP1, logFC_Mouse) %>%
  pivot_longer(
    cols = starts_with("logFC"),
    names_to = "Origin",
    values_to = "logFC"
  ) %>%
  mutate(
    Origin = recode(Origin,
                       "logFC_THP1" = "THP1",
                       "logFC_Mouse" = "Mouse")
  )

long_df$Origin <- factor(long_df$Origin, levels = c("Mouse", "THP1"))
category_colors <- c(
  "THP1" = "sienna1",            # Red
  "Mouse" = "skyblue"#,            # Blue
  #"Promoter/enhancer" = "#fdae61"    # Orange
)

#head(meta)
#head(expr)

#write.csv2(result, "LogFC1_Four_and_eight hours Treated_Significant_Clusters.csv")


pdf("20250718_Overlap_TNFA_Mouse_and_THP1.pdf", width = 14, height = 12)
#pdf("20250710_KO_OE_Proteoglycans_Boxplots.pdf", width = 14, height = 12)

for (category in names(features_list)) {
  genes <- features_list[[category]]
  
  # Subset to relevant genes
  sub_df <- long_df %>%
    filter(Symbol %in% genes)
  
  # Chunk into groups of 20
  gene_chunks <- split(genes, ceiling(seq_along(genes) / 18))
  
  for (chunk in gene_chunks) {
    plot_df <- sub_df %>% filter(Symbol %in% chunk)
    
    if (nrow(plot_df) == 0) {
      message("Empty chunk for genes: ", paste(chunk, collapse = ", "))
      next
    }
    
    # Skip empty chunks
    if (nrow(plot_df) == 0) next
    
    plot_df$Symbol <- factor(plot_df$Symbol, levels = unique(chunk))
    
    p <- ggplot(plot_df, aes(x = Origin, y = logFC)) +
      geom_bar(stat = "identity", width = 0.8, alpha = 0.7, aes(fill = Origin)) +
      geom_text(
        aes(label = sprintf("%.2f", logFC)),
        vjust = -0.3,
        size = 3
      ) +
      scale_y_continuous(breaks = seq(-6, 6, by = 1)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
      ggforce::facet_wrap_paginate(
        ~ Symbol, 
        ncol = 6, 
        nrow = 3,
        scales = "fixed"
      ) +
      scale_fill_manual(values = category_colors) +
      theme_pubr(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "top",
        legend.title = element_blank()
      ) +
      labs(title = category, x = NULL, y = "log2 Fold Change")
    
    print(p)
  }
}

dev.off()

###Venn
mouse_up <- mouse_df %>% filter(logFC > 0) %>% pull(Symbol)
mouse_down <- mouse_df %>% filter(logFC < 0) %>% pull(Symbol)

thp1_up <- tnfa_thp1_df_sig %>% filter(logFC > 0) %>% pull(Symbol)
thp1_down <- tnfa_thp1_df_sig %>% filter(logFC < 0) %>% pull(Symbol)


Up_m1 <- list(Mouse_Up = mouse_up,
                THP1_Up = thp1_up)

Down_m1 <- list(Mouse_Down = mouse_down,
              THP1_Down = thp1_down)
library(VennDetail)
Up_plot <- venndetail(Up_m1)
Up_df <- result(Up_plot)

Down_plot <- venndetail(Down_m1)
Down_df <- result(Down_plot)

pdf("Mouse vs THP1 TNFA Genes_Venn.pdf")
plot(Up_plot, percentage = TRUE, mycol = c("#0073C2FF", 
                                                "#CD534CFF"),order = TRUE,
     alpha = 0.3)
grid::grid.newpage()
plot(Down_plot, percentage = TRUE, mycol = c("#0073C2FF", 
                                                "#CD534CFF"),order = TRUE,
     alpha = 0.3)

dev.off()


sheets <- list("Upregulated" = Up_df,
               "Downregulated" = Down_df)

writexl::write_xlsx(sheets, path= "Mouse vs THP1 TNFA Genes_Venn.xlsx")






