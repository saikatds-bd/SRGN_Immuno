library(dplyr)
library(limma)
library(magrittr)
library(edgeR)
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)

#Loading the helper function
source("MDS_plot.R")
source("Generate_barplot.R")
source("perform_enrichment_function.R")
source("Updated_Volcano_function.R")
source("hallmark_gsea_plot_function.R")


##We will be doing this analysis separately for four and eight hours data because they will be considered separately (i,e Four Hours KO vs WT, and Eight Hours KO vs WT)
expr <- as.data.frame(read_xlsx("C:/Users/ssa214/UiT Office 365/O365-PhD Saikat - General/Proteomics/Mouse_LPS/BMDM_Raw Data LPS.xlsx", sheet = "4hours Raw"))
expr <- column_to_rownames(expr, var = "Gene")

expr$WT_4h_Rep_4 <- NULL #Outlier based on the MDS, removed after confirming
expr$KO_4h_Rep_1 <- NULL #Outlier, removed after confirming

#Metadata for these samples
metadata <- read.csv2("info_4hours.csv", as.is = T, check.names = F)
metadata <- metadata %>% filter(new_name %in% names(expr)) %>% arrange(match(new_name, names(expr)))
metadata <- column_to_rownames(metadata, var = "new_name")

all(rownames(metadata) == colnames(expr))

metadata$Group <- metadata$Condition

group <- as.factor(metadata$Condition)

y <- DGEList(expr, group = group) #Converting to DGE list for edgeR

y$samples$Sample_Names <- rownames(y$samples)

keep.genes <- filterByExpr(y, group=y$samples$group)
table(keep.genes)

y <- y[keep.genes, , keep=FALSE]
y <- calcNormFactors(y, method = "TMM")
lcpm <-edgeR::cpm(y, log = T)

metadata$replicate <- paste0("Rep_", metadata$replicate)
y$samples$replicate <- metadata$replicate
replicates <- as.factor(y$samples$replicate)

y$samples$Lane <- metadata$Run
lanes <- as.factor(y$samples$Lane)


p1 <- plotMDS_ggpubr(
  lcpm = lcpm, 
  sample_data = y$samples, 
  color_var = "replicate", 
  label_var = "Sample_Names", 
  palette = "Set2",      
  title = "Original MDS Plot",
  dims=c(1,2)
)
p1

#Saving the MDS plot
ggsave("Mouse_Four_Hours_MDS.svg", plot = p1,height = 5, width = 6, dpi = 300)

#Exporting the log CPM values
write.csv2(lcpm, "Original Four_Hours_CPM.csv", row.names = T)

#Further checking it with the pearson correlation between samples
library(pheatmap)
raw_cor <- cor(lcpm)
p1 <- pheatmap(raw_cor, main = "Correlation - Four Samples")
p1

ggsave("Mouse_Four_Hours_Correlation.svg", plot = p1,height = 5, width = 6, dpi = 300)

#No apparent effect found for the replicate or lane
design <- model.matrix(~0+group, y$samples)

#Making the design matrix more readable
colnames(design) <- gsub("group","", colnames(design))
design

##Still the samples are not showing very good overlap. So using voomwithqualityweights to reduce the effect of outliers
vwts <- voomWithQualityWeights(y, design=design,normalize.method="none", plot=TRUE)
vfit <- lmFit(vwts, design)
contr <- makeContrasts(KO - WT, levels = colnames(coef(vfit)))
contr
cfit <- contrasts.fit(vfit, contrasts = contr)


##Performing toptable so that it can be used for GSEA
efit <- eBayes(cfit, robust = T)
plotSA(efit, main="Final model: Mean-variance trend")
summary(decideTests(efit))
top.table <- topTable(efit, sort.by = "P", n = Inf, coef = "KO - WT", adjust.method = "BH") #Since we want to see the DE genes between KO and WT
top.table$Symbol <- rownames(top.table) 
save(lcpm, top.table, metadata, file="DE_Mouse_Four_Table.RData")

write.csv2(top.table,"DE Analysis Four_toptable.csv")


##Treat will be used for calling a gene significant or non-significant
tfit <- treat(cfit, robust = T) #Deafult foldchange is 1.2
dt <- decideTests(tfit)
summary(dt)
top.treat<- topTreat(tfit, coef="KO - WT", n=Inf)
top.treat$Symbol <- rownames(top.treat)
write.csv2(top.treat,"DE Analysis Four_toptreat.csv")

save(lcpm, top.treat, metadata, file="DE_Mouse_Four_Treat.RData")


####Making the Volcano Plot
plot_volcano(data=mouse_top.treat, plot_title = "Mouse Wildtype vs Knockout - Four Hours", 
             output_file = "Volcano Mouse Wildtype vs Knockout - Four Hours.svg", 
             height = 4, width = 5)



###Over-representation analysis
mouse_top.treat <- top.treat #Saving it in a new object to save as RData later and use for other analysis
mouse_upregulated <- mouse_top.treat %>% filter(adj.P.Val < 0.05 & logFC > 0) %>% pull(Symbol) #Since we have already incorporated the FC threshold in treat, not using any additional thresholding
mouse_downregulated <- mouse_top.treat %>% filter(adj.P.Val < 0.05 & logFC < 0) %>% pull(Symbol)

library(org.Mm.eg.db) #As it is a mouse data
up_genes <- mouse_upregulated

up_enrichment <- perform_enrichment_analysis(gene_list = up_genes, universe_genes = mouse_top.treat$Symbol,species = "Mus musculus", #Using all genes as background
                                             subset_name = "Upregulated in KO vs WT 4Hours - Mouse",
                                             output_prefix = "Upregulated in KO vs WT_Mouse_ORA",
                                             height = 4, width = 6)
down_genes <- mouse_downregulated
down_enrichment <- perform_enrichment_analysis(gene_list = down_genes, universe_genes = mouse_top.treat$Symbol,species = "Mus musculus",
                                               subset_name = "Downregulated in KO vs WT 4Hours - Mouse",
                                               output_prefix = "Downregulated in KO vs WT_Mouse_ORA",
                                               height = 4, width = 6)






###GSEA
library(clusterProfiler)
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)
library(dplyr)
library(stats)
library(readxl)
library(ggplot2)
library(viridis)

back <-  as.data.frame(read_xlsx("BMDM_Raw Data LPS.xlsx", sheet = "4hours Raw")) %>% pull(Gene)
#df <- read.csv2("limma_Outlier and batch Removed_4Hours.csv", as.is = T, check.names = F, row.names = 1)
df <- top.treat
eg <- bitr(df$Symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
df2 <- left_join(df, eg, by = c("Symbol" = "SYMBOL"))
m_t2g_r <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% 
  mutate(gs_name=gsub("REACTOME_","",gs_name)) %>%
  dplyr::select(gs_name,gene_symbol)

m_t2g_reactome <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME") %>% 
  mutate(gs_name=gsub("REACTOME_","",gs_name)) %>%
  dplyr::select(gs_name,gene_symbol)

m_t2g_hallmark <- msigdbr(species = "Mus musculus", category = "H") %>% 
  mutate(gs_name=gsub("HALLMARK_","",gs_name)) %>%
  dplyr::select(gs_name,gene_symbol)
term2gene_list <- list(
  Hallmark = m_t2g_hallmark,
  Reactome = m_t2g_reactome
)

df2_unique <- df2 %>% 
  distinct(ENTREZID, .keep_all = TRUE)  %>%  filter(!is.na(ENTREZID))
back_2 <- bitr(back, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")


##GSEA
gs <- df %>% 
  mutate(signed_rank_stats = sign(logFC) * -log10(P.Value)) %>%
  #left_join(background_genes, by= c("gene" = "SYMBOL")) %>%
  arrange(desc(signed_rank_stats))

gene_list<- gs$signed_rank_stats
names(gene_list)<- gs$Symbol

set.seed(123)
em2 <- GSEA(gene_list, TERM2GENE=m_t2g_hallmark, seed = 123,
            pvalueCutoff = 1)

four_hours_gsea <- as.data.frame(em2@result)
four_hours_top_treat <- top.treat
four_hours_lcpm <- lcpm
four_hours_counts <- expr

save(four_hours_gsea, metadata, four_hours_top_treat, 
     four_hours_lcpm, four_hours_counts, file = "Four_Hourse_DE.RData")

library(GseaVis)
nes <- em2[,c(1,5,7)]
nes <- nes %>% mutate(ID = gsub("HALLMARK_", "", ID),      # Remove "HALLMARK_" prefix
                      ID = gsub("_", " ", ID))
nes$ID <- ifelse(
  nes$p.adjust < 0.05,
  paste0(nes$ID, " *"),
  nes$ID
)

nes_sorted <- nes %>% arrange(desc(NES))

# Select the top 10 IDs with the highest NES values
top_10_nes <- nes_sorted %>% slice_head(n = 20)

# Select the bottom 10 IDs with the lowest NES values
bottom_10_nes <- nes_sorted %>% slice_tail(n = 20)

# We can combine these into a single dataframe if needed:
combined_top_bottom_nes <- bind_rows(
  top_10_nes %>% mutate(Category = "Top 10"),
  bottom_10_nes %>% mutate(Category = "Bottom 10")
)

plot <- ggplot(combined_top_bottom_nes, aes(reorder(ID, NES), NES)) +
  geom_col(aes(fill= `p.adjust`), width = 0.8) +           # Fill by actual p.adjust values
  scale_fill_gradient(low = "#347eba", high="#df6664",
                      guide = "colorbar") +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="KO vs WT Four Hours - Hallmark", fill="Adjusted p-value") + 
  theme_bw()+theme(panel.spacing = unit(0.1, "lines")) #+ coord_fixed(ratio=0.4)

plot
pdf("GSEA Four_Hours_Hallmakr.pdf",  width=10,height=8)
print(plot)
dev.off()

##Reactome
gs2 <- df2_unique %>% 
  mutate(signed_rank_stats = sign(logFC) * -log10(P.Value)) %>%
  #left_join(background_genes, by= c("gene" = "SYMBOL")) %>%
  arrange(desc(signed_rank_stats))

gene_list2<- gs2$signed_rank_stats
names(gene_list2)<- gs2$ENTREZID
set.seed(123)
library(ReactomePA)
em3 <- gsePathway(gene_list2, organism = "human", seed = 123,
                  pvalueCutoff = 1)

library(GseaVis)
nes <- em3[,c(2,5,7)]
nes$Description  <- ifelse(
  nes$p.adjust < 0.05,
  paste0(nes$Description, " *"),
  nes$Description
)

nes_sorted <- nes %>% arrange(desc(NES))

# Select the top 10 IDs with the highest NES values
top_10_nes <- nes_sorted %>% slice_head(n = 20)

# Select the bottom 10 IDs with the lowest NES values
bottom_10_nes <- nes_sorted %>% slice_tail(n = 20)

# Optionally, you can combine these into a single dataframe if needed:
combined_top_bottom_nes <- bind_rows(
  top_10_nes %>% mutate(Category = "Top 10"),
  bottom_10_nes %>% mutate(Category = "Bottom 10")
)

combined_top_bottom_nes <- distinct(combined_top_bottom_nes, Description, .keep_all = T)

plot2 <- ggplot(combined_top_bottom_nes, aes(reorder(Description, NES), NES)) +
  geom_col(aes(fill= `p.adjust`), width = 0.8) +           # Fill by actual p.adjust values
  scale_fill_gradient(low = "#347eba", high="#df6664",
                      guide = "colorbar") +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="KO vs WT M0 - Reactome", fill="Adjusted p-value") + 
  theme_bw()+theme(panel.spacing = unit(0.1, "lines")) #+ coord_fixed(ratio=0.4)


plot2
library(patchwork)
gs2 <- df2_unique %>% 
  mutate(signed_rank_stats = sign(logFC) * -log10(P.Value)) %>%
  #left_join(background_genes, by= c("gene" = "SYMBOL")) %>%
  arrange(desc(signed_rank_stats))

gene_list2<- gs2$signed_rank_stats
names(gene_list2)<- gs2$ENTREZID
set.seed(123)
library(ReactomePA)
em3 <- gseKEGG(gene_list2, organism = "hsa", seed = 123,
               pvalueCutoff = 1)

library(GseaVis)
nes <- em3[,c(2,5,7)]
nes$Description  <- ifelse(
  nes$p.adjust < 0.05,
  paste0(nes$Description, " *"),
  nes$Description
)

nes_sorted <- nes %>% arrange(desc(NES))

# Select the top 10 IDs with the highest NES values
top_10_nes <- nes_sorted %>% slice_head(n = 20)

# Select the bottom 10 IDs with the lowest NES values
bottom_10_nes <- nes_sorted %>% slice_tail(n = 20)

# Optionally, you can combine these into a single dataframe if needed:
combined_top_bottom_nes <- bind_rows(
  top_10_nes %>% mutate(Category = "Top 10"),
  bottom_10_nes %>% mutate(Category = "Bottom 10")
)

combined_top_bottom_nes <- distinct(combined_top_bottom_nes, Description, .keep_all = T)

plot3 <- ggplot(combined_top_bottom_nes, aes(reorder(Description, NES), NES)) +
  geom_col(aes(fill= `p.adjust`), width = 0.8) +           # Fill by actual p.adjust values
  scale_fill_gradient(low = "#347eba", high="#df6664",
                      guide = "colorbar") +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="KO vs WT M0 - KEGG", fill="Adjusted p-value") + 
  theme_bw()+theme(panel.spacing = unit(0.1, "lines")) #+ coord_fixed(ratio=0.4)


plot3


pdf("GSEA M0_toptable.pdf",  width=10,height=8)
print(plot)
print(plot2)
print(plot3)
dev.off()

# Overrepresentation Analysis

m_t2g_r <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% 
  mutate(gs_name=gsub("REACTOME_","",gs_name)) %>%
  dplyr::select(gs_name,gene_symbol)


m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  mutate(gs_name=gsub("HALLMARK_","",gs_name)) %>%
  dplyr::select(gs_name,gene_symbol)
significant_threshold <- 0.05
up <-  df %>% filter(adj.P.Val < 0.05 & logFC > 0.5) %>% pull(Symbol)
up_ent <- df2 %>% filter(adj.P.Val < 0.05 & logFC > 0.5) %>% pull(ENTREZID)

enricher_up <- enricher(up, TERM2GENE=m_t2g, 
                        universe = back,
                        pvalueCutoff = 1,
                        qvalueCutoff = 1
)


enricher_up@result$Description <- ifelse(
  enricher_up@result$p.adjust < significant_threshold,
  paste0(enricher_up@result$Description, " *"),
  enricher_up@result$Description
)

up_plot <- dotplot(enricher_up) + ggtitle("Upregulated - 4 Hours")
up_plot


library(ReactomePA)
enricher_up2 <- enrichPathway(up_ent, #TERM2GENE=m_t2g,
                                organism = "mouse",
                                universe = back_2$ENTREZID,
                                pvalueCutoff = 1,
                                qvalueCutoff = 1
)

enricher_up2@result$Description <- ifelse(
  enricher_up2@result$p.adjust < significant_threshold,
  paste0(enricher_up2@result$Description, " *"),
  enricher_up2@result$Description
)

up_plot2 <- dotplot(enricher_up2) + ggtitle("Upregulated - 4 Hours Reactome")
up_plot2


down <-  df %>% filter(adj.P.Val < 0.05 & logFC < -0.5) %>% pull(Symbol)

enricher_down <- enricher(down, TERM2GENE=m_t2g, 
                          universe = df$Symbol,
                          pvalueCutoff = 1,
                          qvalueCutoff = 1
)

enricher_down@result$Description <- ifelse(
  enricher_down@result$p.adjust < significant_threshold,
  paste0(enricher_down@result$Description, " *"),
  enricher_down@result$Description
)

down_plot <- dotplot(enricher_down) + ggtitle("Downregulated - 4 Hours")
down_plot

down_ent <- df2 %>% filter(adj.P.Val < 0.05 & logFC < -0.5) %>% pull(ENTREZID)
library(ReactomePA)
enricher_down2 <- enrichPathway(down_ent, #TERM2GENE=m_t2g,
                              organism = "mouse",
                              universe = back_2$ENTREZID,
                              pvalueCutoff = 1,
                              qvalueCutoff = 1
)

enricher_down2@result$Description <- ifelse(
  enricher_down2@result$p.adjust < significant_threshold,
  paste0(enricher_down2@result$Description, " *"),
  enricher_down2@result$Description
)

down_plot2 <- dotplot(enricher_down2, showCategory=10) + ggtitle("Downregulated - 4 Hours Reactome")
down_plot2

pdf(paste0("4hours_Enrich.pdf"),  width=7,height=6)
print(up_plot)
print(down_plot)
print(up_plot2)
print(down_plot2)
dev.off()

