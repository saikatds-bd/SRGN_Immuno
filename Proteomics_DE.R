library(tidyverse)
library(dplyr)
library(ggpubr)
library(limma)
library(magrittr)
library(edgeR)
library(readxl)
library(DEP)
library(RColorBrewer)
source("MDS_plot.R")
source("Generate_barplot.R")
source("perform_enrichment_function.R")
source("Updated_Volcano_function.R")
source("hallmark_gsea_plot_function.R")


solid_colors <- c("KO_M0" = "#5094c2",  
                  "KO_M1" = "#fca050",  
                  "KO_M2" = "#40a940",  
                  "WT_M0" = "#db4647",  
                  "WT_M1" = "#a885c8", 
                  "WO_M2" = "#9a6c62")

fill_colors <- c("KO_M0" = "#88b4d3",  
                 "KO_M1" = "#f9b880",  
                 "KO_M2" = "#8fc98f",  
                 "WT_M0" = "#e48c8d",  
                 "WT_M1" = "#c3acd8", 
                 "WT_M2" = "#bfa49e")

# Define colors for the three replicates
replicate_colors <- c("Replicate1" = "#A6CEE3",
                      "Replicate2" = "#1F78B4",
                      "Replicate3" = "#B2DF8A")

raw <- read.delim("proteinGroups.txt", as.is = T, check.names = F)#Maxquant output file
design <- read.csv2("metadata.csv", as.is = T, check.names = F)

x <- raw[,c(150:167)]
names(x) <- gsub("LFQ intensity ","", names(x))

design <- design[match(names(x), design$label), ] #making sure the orders are same


#Filtering based on Only identified by site, reverse and potential contaminant
sub_raw <- raw %>%
  filter(is.na(Reverse) & 
           is.na(`Potential contaminant`) & 
           is.na(`Only identified by site`))

#Checking if there is any duplicated gene names
table(sub_raw$`Gene names` %>% duplicated()) #Only 2 duplicated gene names

#Renaming the duplicated gene names
raw_unique <- make_unique(sub_raw, "Gene names", "Protein IDs", delim = ";") #Will create two new columns

#Rechecking if there is any more duplicated gene names
raw_unique$name %>% duplicated() %>% any()

LFQ_columns <- grep("LFQ ", colnames(raw_unique)) # get LFQ column numbers
x <- raw_unique[,c(150:167)]
all(names(x)==design$rename)
names(x) <- gsub("LFQ intensity ","", names(x))
all(names(x)==design$rename)


experimental_design <- design #It has column names: label, condition and replicate
experimental_design$label <- experimental_design$rename
experimental_design <- experimental_design[,-c(4:7)] #It has column names: label, condition and replicate
all(names(x)==experimental_design$label)


raw_unique[LFQ_columns] <- lapply(raw_unique[LFQ_columns], function(x) {
  x <- gsub("eTRUE", "e+", x, fixed = TRUE)
  
  x <- suppressWarnings(as.numeric(x))
  #x[x == 0] <- NA
  x
})


#names(experimental_design) <- c("label", "condition", "replicate")
data_se <- DEP::make_se(raw_unique, LFQ_columns, experimental_design) #Here the assay data is log2 transformed
data_se
dim(data_se)

data_filt <- filter_missval(data_se, thr = 1)
y <- get_df_wide(data_filt)
data_filt_man <- DEP::impute(data_filt, fun = "man", shift = 1.8, scale = 0.3)
y <- get_df_wide(data_filt_man)

write.csv2(y, "Imputed_unnormalized_mq.csv")

###Limma DE
KD <- read.csv2("C:/Users/ssa214/UiT Office 365/O365-PhD Saikat - General/Proteomics/20240130/Imputed_unnormalized_mq.csv", as.is = T,check.names = F, row.names = 1)
gc()

rownames(KD) <- KD$name
KD$name <- NULL

expr <- KD[,c(1:18)]

hist(as.matrix(expr), col="blue", border="white", breaks=50, main = "Histogram of expression values")

#Get the metadata
metadata <- read.csv2("C:/Users/ssa214/UiT Office 365/O365-PhD Saikat - General/Proteomics/20240130/metadata.csv", as.is = T, check.names = F) #Sample metadatarmation.
metadata <- metadata %>% filter(rename %in% names(expr)) %>% arrange(match(rename, names(expr)))
metadata <- column_to_rownames(metadata, var = "rename")
all(rownames(metadata) == colnames(expr))

metadata$Sample_Names <- rownames(metadata)
metadata$group <- metadata$condition
metadata$replicate <- paste0("Rep_",metadata$replicate)

replicates <- as.factor(metadata$replicate)
phase <- as.factor(metadata$Phase)
condition <- as.factor(metadata$condition)

p1 <- plotMDS_ggpubr(
  lcpm = expr, 
  sample_data = metadata, 
  color_var = "replicate", # Use 'Condition' column for colors
  label_var = "Sample_Names", # Use 'Replicate' column for labels
  palette = "Set2",         # Change color palette
  title = "Original MDS Plot",
  dims=c(1,2)
)
p1

group <- as.factor(metadata$condition)
design <- model.matrix(~0+condition+replicates)
colnames(design) <- gsub("condition","", colnames(design))

vfit <- lmFit(expr, design, method = "robust")

contr <- makeContrasts(M0=KO_M0 - WT_M0,
                       M1=KO_M1 - WT_M1,
                       M2=KO_M2 - WT_M2,
                       KO_M1 = KO_M1 - KO_M0,
                       KO_M2 = KO_M2 - KO_M0,
                       WT_M1 = WT_M1 - WT_M0,
                       WT_M2 = WT_M2 - WT_M0,
                       levels = colnames(coef(vfit)))
contr


vfit <- contrasts.fit(vfit, contr)

efit <- eBayes(vfit)
summary(decideTests(efit))

KO_vs_WT_M0_top.table <- topTable(efit, sort.by = "P", n = Inf, coef = "M0")
KO_vs_WT_M0_top.table$Symbol <- rownames(KO_vs_WT_M0_top.table)

KO_vs_WT_M1_top.table <- topTable(efit, sort.by = "P", n = Inf, coef = "M1")
KO_vs_WT_M1_top.table$Symbol <- rownames(KO_vs_WT_M1_top.table)

plot_volcano(data=KO_vs_WT_M1_top.table, plot_title = "Secretome Wildtype vs Knockout - M1", 
             output_file = "Volcano Secretome Wildtype vs Knockout - M1.svg", 
             height = 4, width = 4, fc_cut = 0.5)

##ORA on up and downregulated
up_genes <- ko_wt_m1_top.treat %>% filter(adj.P.Val < 0.05 & logFC > 0) %>% pull(Symbol)
down_genes <- ko_wt_m1_top.treat %>% filter(adj.P.Val < 0.05 & logFC < 0) %>% pull(Symbol)

library(org.Hs.eg.db)
library(msigdbr)
up_enrichment <- perform_enrichment_analysis(gene_list = up_genes, universe_genes = ko_wt_m1_top.treat$Symbol,
                                             species = "Homo sapiens",
                                             subset_name = "Upregulated - Knockout vs Wildtype Secretome",
                                             output_prefix = "Upregulated - Secretome - KO vs WT M1 ORA",
                                             height = 4, width = 6)


down_enrichment <- perform_enrichment_analysis(gene_list = down_genes, universe_genes = ko_wt_m1_top.treat$Symbol,
                                               species = "Homo sapiens",
                                               subset_name = "Downregulated - Knockout vs Wildtype Secretome",
                                               output_prefix = "Downregulated - Secretome - KO vs WT M1 ORA",
                                               height = 4, width = 6)
