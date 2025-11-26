library(dplyr)
library(limma)
library(edgeR)
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)

source("MDS_plot.R")
source("Generate_barplot.R")
source("perform_enrichment_function.R")
source("Updated_Volcano_function.R")
source("hallmark_gsea_plot_function.R")

metadata <- read.csv2("THP1_Metadata.csv",as.is = T, check.names = F)
expr <- read.csv2("THP1_RNAseq_raw_counts.csv",as.is = T, check.names = F)

expr <- subset_df

expr$KO_M0_Rep3 <- NULL ##MDS revealed KO_M0_Rep3 is an outlier

expr <- column_to_rownames(expr, var = "Name")
expr$Identifier <- NULL
names(expr)
metadata <- read_xlsx("Frode_Metadata_Combined.xlsx")
metadata <- metadata %>% filter(Modified_Names %in% names(expr)) %>% arrange(match(Modified_Names, names(expr)))
metadata <- column_to_rownames(metadata, var = "Modified_Names")
all(rownames(metadata) == colnames(expr))

group <- as.factor(metadata$Group)

y <- DGEList(expr, group = group) #Converting to DGE list for edgeR

y$samples$Sample_Names <- rownames(y$samples)

keep.genes <- filterByExpr(y, group=y$samples$group)
table(keep.genes)

y <- y[keep.genes, , keep=FALSE]
y <- calcNormFactors(y, method = "TMM")

y$samples$replicate <- metadata$Replicate
replicates <- as.factor(y$samples$replicate)

y$samples$Lane <- metadata$Lane
lanes <- as.factor(y$samples$Lane)

lcpm <-edgeR::cpm(y, log = T)
p1 <- plotMDS_ggpubr(
  lcpm = lcpm, 
  sample_data = y$samples, 
  color_var = "replicate", # Use 'Condition' column for colors
  label_var = "Sample_Names", # Use 'Replicate' column for labels
  palette = "Set2",         # Change color palette
  title = "Original MDS Plot",
  dims=c(1,2)
)
p1



write.csv2(lcpm, "Original M0_CPM.csv", row.names = T)

library(pheatmap)
raw_cor <- cor(lcpm)
p1 <- pheatmap(raw_cor, main = "Correlation - M0 Samples")
p1

design <- model.matrix(~0+group+replicates, y$samples)
design
colnames(design) <- gsub("group","", colnames(design))
design

v <- voom(y, design, plot=TRUE)
vfit <- lmFit(v, design)


contr <- makeContrasts(KO_m1vsm0=KO_M1 - KO_M0,
                       WT_m1vsm0=WT_M1 - WT_M0,
                       KO_m1vsWT_M1=KO_M1 - WT_M1,
                       KO_m0vsWT_M0=KO_M0 - WT_M0,
                       levels = colnames(coef(vfit)))
contr
cfit <- contrasts.fit(vfit, contrasts = contr)

efit <- eBayes(cfit, robust = T)

wt_top.table <- topTable(efit, sort.by = "P", n = Inf, coef = "WT_m1vsm0", adjust.method = "BH")
wt_top.table$Symbol <- rownames(wt_top.table)

ko_top.table <- topTable(efit, sort.by = "P", n = Inf, coef = "KO_m1vsm0", adjust.method = "BH")
ko_top.table$Symbol <- rownames(ko_top.table)

ko_wt_m1_top.table <- topTable(efit, sort.by = "P", n = Inf, coef = "KO_m1vsWT_M1", adjust.method = "BH")
ko_wt_m1_top.table$Symbol <- rownames(ko_wt_m1_top.table)

ko_wt_m0_top.table <- topTable(efit, sort.by = "P", n = Inf, coef = "KO_m0vsWT_M0", adjust.method = "BH")
ko_wt_m0_top.table$Symbol <- rownames(ko_wt_m0_top.table)

#Saving the table method DE data for GSEA
save(lcpm, wt_top.table, metadata,ko_top.table, ko_wt_m1_top.table, ko_wt_m0_top.table, file="C:/Users/ssa214/UiT Office 365/O365-PhD_Saikat - General/Proteomics/20250218/DE_KO_WT_M1.RData")

#Treat Data
tfit <- treat(cfit, robust = T)

wt_top.treat <- topTreat(tfit, sort.by = "P", n = Inf, coef = "WT_m1vsm0", adjust.method = "BH")
wt_top.treat$Symbol <- rownames(wt_top.treat)

ko_top.treat <- topTreat(tfit, sort.by = "P", n = Inf, coef = "KO_m1vsm0", adjust.method = "BH")
ko_top.treat$Symbol <- rownames(ko_top.treat)

ko_wt_m1_top.treat <- topTreat(tfit, sort.by = "P", n = Inf, coef = "KO_m1vsWT_M1", adjust.method = "BH")
ko_wt_m1_top.treat$Symbol <- rownames(ko_wt_m1_top.treat)

ko_wt_m0_top.treat <- topTreat(tfit, sort.by = "P", n = Inf, coef = "KO_m0vsWT_M0", adjust.method = "BH")
ko_wt_m0_top.treat$Symbol <- rownames(ko_wt_m0_top.treat)
save(lcpm, wt_top.treat, metadata,ko_top.treat, ko_wt_m1_top.treat, ko_wt_m0_top.treat, file="DE_KO_WT_M1_Treat.RData")

plot_volcano(data=ko_wt_m1_top.treat, plot_title = "THP1 Knockout vs Wildtype - M1", 
             output_file = "Volcano THP1 Wildtype vs Knockout - M1.svg", 
             height = 4, width = 4, fc_cut = 1.2)



##ORA on up and downregulated
up_genes <- ko_wt_m1_top.treat %>% filter(adj.P.Val < 0.05 & logFC > 0) %>% pull(Symbol)
down_genes <- ko_wt_m1_top.treat %>% filter(adj.P.Val < 0.05 & logFC < 0) %>% pull(Symbol)

library(org.Hs.eg.db)
library(msigdbr)
up_enrichment <- perform_enrichment_analysis(gene_list = up_genes, universe_genes = ko_wt_m1_top.treat$Symbol,
                                             species = "Homo sapiens",
                                             subset_name = "Upregulated - Knockout vs Wildtype THP1 M1",
                                             output_prefix = "Upregulated - THP1 - KO vs WT M1 ORA",
                                             height = 4, width = 6)


down_enrichment <- perform_enrichment_analysis(gene_list = down_genes, universe_genes = ko_wt_m1_top.treat$Symbol,
                                             species = "Homo sapiens",
                                             subset_name = "Downregulated - Knockout vs Wildtype THP1 M1",
                                             output_prefix = "Downregulated - THP1 - KO vs WT M1 ORA",
                                             height = 4, width = 6)


#GSEA
#Hallmark and GO BP Datasets
m_t2g <- msigdbr(species = "Homo sapiens", collection = "H") %>% 
  mutate(gs_name=gsub("HALLMARK_","",gs_name)) %>%
  dplyr::select(gs_name,gene_symbol)

m_t2go_bp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") %>% 
  mutate(gs_name=gsub("GOBP_","",gs_name)) %>%
  dplyr::select(gs_name,gene_symbol)

wt_m1_gsea_go <- gsea_barplot_multi(
  df = wt_top.table, plot_title = "M1 vs M0 Wildtype THP1 - GO:BP",
  use = "GO:BP",
  add_source_prefix = FALSE,
  n_each = 10
)
ggsave("THP1 M1 vs M0 Wildtype - GO_BP.svg",wt_m1_gsea_go[["plot"]], dpi = 300, height = 6, width = 8, bg = "transparent")

wt_m1_gsea <- gsea_barplot_multi(
  df = wt_top.table, plot_title = "M1 vs M0 Wildtype THP1 - Hallmark",
  m_t2g = m_t2g,
  #m_t2go_bp = m_t2go_bp,
  use = c("H"),
  add_source_prefix = FALSE,
  n_each = 10
)
ggsave("THP1 M1 vs M0 Wildtype - Hallmark.svg",wt_m1_gsea[["plot"]], dpi = 300, height = 6, width = 8, bg = "transparent")

ko_m1_gsea_go <- gsea_barplot_multi(
  df = ko_top.table, plot_title = "M1 vs M0 Knockout THP1",
  #m_t2go_bp = m_t2go_bp,
  use = "GO:BP",
  add_source_prefix = FALSE,
  n_each = 10
)
ggsave("THP1 M1 vs M0 Knockout - GO_BP.svg",ko_m1_gsea_go[["plot"]], dpi = 300, height = 6, width = 8, bg = "transparent")

ko_m1_gsea <- gsea_barplot_multi(
  df = ko_top.table, plot_title = "M1 vs M0 Wildtype THP1 - Hallmark",
  m_t2g = m_t2g,
  #m_t2go_bp = m_t2go_bp,
  use = c("H"),
  add_source_prefix = FALSE,
  n_each = 10
)
ggsave("THP1 M1 vs M0 Knockout - Hallmark.svg",ko_m1_gsea[["plot"]], dpi = 300, height = 6, width = 8, bg = "transparent")

save(lcpm, ko_wt_m0_top.table, ko_wt_m0_top.treat,
     ko_wt_m1_top.table, ko_wt_m1_top.treat,
     wt_top.table, wt_top.treat, ko_top.table, ko_top.treat, file = "THP1_DE_KO.RData")


##Comparisons between GSEA Objects

m_go_list <- list(Wildtype = wt_m1_gsea_go[["gsea"]][["GO_BP"]][["res"]],
                  Knockout = ko_m1_gsea_go[["gsea"]][["GO_BP"]][["res"]])


library(GseaVis)
p1 <- GseaVis::GSEAmultiGP(gsea_list = m_go_list,
                           geneSetID = "GO:0042116",
                           exp_name = c("Wildtype","Knockout"),
                           curve.col = ggsci::pal_lancet()(3))

p1

ggsave("WT_KO_GSEA_Macrophage_activation.svg", p1, width = 5, height = 4, dpi = 300)
m_gsea_list <- list(Wildtype = wt_m1_gsea[["gsea"]][["Hallmark"]][["res"]],
                    Knockout = ko_m1_gsea[["gsea"]][["Hallmark"]][["res"]])


p2 <- GseaVis::GSEAmultiGP(gsea_list = m_gsea_list,
                           geneSetID = "TNFA_SIGNALING_VIA_NFKB",
                           exp_name = c("Wildtype","Knockout"),
                           curve.col = ggsci::pal_lancet()(3))


p2


#Pulling out the gene lists associated with macrophages
macrophage_dif <- m_t2go_bp %>% 
  filter(gs_name == "MACROPHAGE_DIFFERENTIATION") %>%
  pull(gene_symbol)

macrophage_activation <- m_t2go_bp %>% 
  filter(gs_name == "MACROPHAGE_ACTIVATION") %>%
  pull(gene_symbol)

#Extracting the gene lists for later
save(macrophage_dif, macrophage_activation, m_t2go_bp, file="Macrophage_Genes.RData")

