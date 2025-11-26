library(dplyr)
library(limma)
library(edgeR)
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)

source("MDS_plot.R")
source("Generate_barplot.R")
source("perform_enrichment_function_correlated.R")
source("Updated_Volcano_function.R")
source("hallmark_gsea_plot_function.R")

ctrl <- readxl::read_xlsx("Srgn_correlations_control.xlsx")
pos_ctrl <- ctrl %>% filter(bicor > 0 & pvalue < 0.05) %>% dplyr::select(gene_symbol_2,bicor, pvalue)
neg_ctrl <- ctrl %>% filter(bicor < 0 & pvalue < 0.05) %>% dplyr::select(gene_symbol_2,bicor, pvalue)
lps <- readxl::read_xlsx("Srgn_correlations_LPS.xlsx")
pos_lps <- lps %>% filter(bicor > 0 & pvalue < 0.05) %>% dplyr::select(gene_symbol_2,bicor, pvalue)
neg_lps <- lps %>% filter(bicor < 0 & pvalue < 0.05) %>% dplyr::select(gene_symbol_2,bicor, pvalue)

pos_ctrl_enrichment <- perform_enrichment_analysis(gene_list = pos_ctrl$gene_symbol_2, 
                                             subset_name = "Positively Correlated in Control",
                                             output_prefix = "Positively Correlated in Control_ORA",
                                             height = 4, width = 6)

neg_ctrl_enrichment <- perform_enrichment_analysis(gene_list = neg_ctrl$gene_symbol_2, 
                                                   subset_name = "Negatively Correlated in Control",
                                                   output_prefix = "Negatively Correlated in Control_ORA",
                                                   height = 4, width = 6)

pos_lps_enrichment <- perform_enrichment_analysis(gene_list = pos_lps$gene_symbol_2, 
                                                   subset_name = "Positively Correlated in LPS",
                                                   output_prefix = "Positively Correlated in LPS_ORA",
                                                   height = 4, width = 6)

neg_lps_enrichment <- perform_enrichment_analysis(gene_list = neg_lps$gene_symbol_2, 
                                                   subset_name = "Negatively Correlated in LPS",
                                                   output_prefix = "Negatively Correlated in LPS_ORA",
                                                   height = 4, width = 6)


