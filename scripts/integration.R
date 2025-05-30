# Load libraries and data ----
library(AnnotationDbi)
library(affy)
library(biomaRt)
library(Biobase)
library(car)
library(ComplexHeatmap)
library(colorRamp2)
library(forcats)
library(gprofiler2)
library(ggpubr)
library(gplots)
library(ggrepel)
library(gridExtra)
library(ggvenn)
library(GOplot) # Probar si puedo usar esto
library(hgu133plus2.db)
library(limma)
library(openxlsx)
library(pheatmap)
library(plyr)
library(paletteer)
library(plotly)
library(pheatmap)
library(sva)
library(tidyverse)
library(treemap)
library(rstatix)

load("data/rdata/HtoM.RData")
load("data/rdata/glyco_terms.RData")
load("~/projects/tesis/output/integration/glycogenes2.RData")
load("~/projects/tesis/data/pancolitis/pan_uc_toptable.RData")
load("~/projects/tesis/data/pancolitis/cacrc_uc_toptable.RData")

firma <- c("DEFA4", "DEFA5", "DEFA6", "SLC9A3", "LGALS4", "DEGS2", "SLC35G1",
           "LGALS12", "HNF4A", "SLC26A2", "ICAM1", "SELL", "ST3GAL1", "SELP",
           "CHI3L1", "CHIA", "GCNT2")
info <- data.frame(author = c("Digby-Bell", "Watanabe", "Pekow", "Bjerrum", "Arijs", "Zhang", "Tang"),
                   gse = c("GSE224758", "GSE3629", "GSE37283", "GSE47908", "GSE73661", "GSE87473, GSE92415, and GSE206285", "GSE31106"),
                   phenotypes = c("h-uc", "uc-cac-sporadic", "h-uc-neo", "h-uc-d", "h-uc", "h-uc", "h-uc-dys-aden"))
process_data <- function(data, id) {
  c1 <- paste0(id, ".lfc")
  c2 <- paste0(id, ".adjP")
  data %>% dplyr::select(c("logFC", "adj.P.Val")) %>%
    dplyr::rename(!!c1 := "logFC") %>%
    dplyr::mutate(adj.P.Val = ifelse(adj.P.Val < 0.05, "*", "")) %>%
    dplyr::rename(!!c2 := "adj.P.Val")
}
process_data_nes <- function(data, id) {
  c1 <- paste0(id, ".nes")
  c2 <- paste0(id, ".adjP")
  data %>% dplyr::select(c("NES", "p.adjust")) %>%
    dplyr::rename(!!c1 := "NES") %>%
    dplyr::rename(!!c2 := "p.adjust") %>%
    dplyr::mutate(
      !!c2 := case_when(
        !!sym(c2) <= 0.00000000001 ~ "****",
        !!sym(c2) <= 0.00000001 ~ "***",
        !!sym(c2) <= 0.00001  ~ "**",
        !!sym(c2) <= 0.01  ~ "*",
        TRUE               ~ ""))
}
load("data/rdata/HtoM.RData")
colnames(HtoM) <- c("human", "mouse")

# Import gene data----
# GSE224758 Digby-Bell
digby.uc_c <- read.xlsx(xlsxFile = "output/GSE224758/tables/uc_c_tables.xlsx", sheet = "all_genes", rowNames = T)

# GSE3629 Watanabe
watanabe.cac_sporadic <- read.xlsx(xlsxFile = "output/GSE3629/tables/cac_sporadic_tables.xlsx", sheet = "all_genes", rowNames = T)
watanabe.ca_pan <- read.xlsx(xlsxFile = "output/GSE3629/tables/cac_uc_tables.xlsx", sheet = "all_genes", rowNames = T)
watanabe.ca_pan$comparison <- "CACRC - Pancolitis"

# GSE37283 Pekow
pekow.n_c <- read.xlsx(xlsxFile = "output/GSE37283/tables/n_c_tables.xlsx", sheet = "all_genes", rowNames = T)
pekow.n_uc <- read.xlsx(xlsxFile = "output/GSE37283/tables/n_uc_tables.xlsx", sheet = "all_genes", rowNames = T)
pekow.uc_c <- read.xlsx(xlsxFile = "output/GSE37283/tables/uc_c_tables.xlsx", sheet = "all_genes", rowNames = T)

# GSE47908 Bjerrum
bjerrum.d_c <- read.xlsx(xlsxFile = "output/GSE47908/tables/d_c_tables.xlsx", sheet = "all_genes", rowNames = T)
bjerrum.d_uc <- read.xlsx(xlsxFile = "output/GSE47908/tables/d_uc_tables.xlsx", sheet = "all_genes", rowNames = T)
bjerrum.uc_c <- read.xlsx(xlsxFile = "output/GSE47908/tables/uc_c_tables.xlsx", sheet = "all_genes", rowNames = T)
load("~/projects/tesis/data/pancolitis/pan_uc_toptable.RData")
bjerrum.pan_uc <- pan_uc_toptable
load("~/projects/tesis/data/pancolitis/dys_pan_toptable.RData")
bjerrum.dys_pan <- dys_pan_toptable

# GSE73661 roman 
arijs.uc_c <- read.xlsx(xlsxFile = "output/roman/uc_c_tables.xlsx", sheet = "all", rowNames = T)

# zhang
zhang.uc_c <- read.xlsx(xlsxFile = "output/zhang/uc_c_tables.xlsx", sheet = "all", rowNames = T)

# GSE31106 - Tang
tang.a_c <- read.xlsx(xlsxFile = "output/GSE31106/tables/a_c_tables.xlsx", sheet = "all_genes", rowNames = T)
tang.a_i <- read.xlsx(xlsxFile = "output/GSE31106/tables/a_i_tables.xlsx", sheet = "all_genes", rowNames = T)
tang.d_i <- read.xlsx(xlsxFile = "output/GSE31106/tables/d_i_tables.xlsx", sheet = "all_genes", rowNames = T)
tang.a_d <- read.xlsx(xlsxFile = "output/GSE31106/tables/a_d_tables.xlsx", sheet = "all_genes", rowNames = T)
tang.d_c <- read.xlsx(xlsxFile = "output/GSE31106/tables/d_c_tables.xlsx", sheet = "all_genes", rowNames = T)
tang.i_c <- read.xlsx(xlsxFile = "output/GSE31106/tables/i_c_tables.xlsx", sheet = "all_genes", rowNames = T)

# Import NES data----
# GSE224758 Digby-Bell
digby.uc_c <- read.xlsx(xlsxFile = "output/GSE224758/tables/uc_c_tables.xlsx", sheet = "full_gsea", rowNames = T)

# GSE3629 Watanabe
watanabe.cac_sporadic <- read.xlsx(xlsxFile = "output/GSE3629/tables/cac_sporadic_tables.xlsx", sheet = "full_gsea", rowNames = T)
watanabe.ca_pan <- read.xlsx(xlsxFile = "output/GSE3629/tables/cac_uc_tables.xlsx", sheet = "full_gsea", rowNames = T)
watanabe.ca_pan$comparison <- "CACRC - Pancolitis"

# GSE37283 Pekow
pekow.n_c <- read.xlsx(xlsxFile = "output/GSE37283/tables/n_c_tables.xlsx", sheet = "full_gsea", rowNames = T)
pekow.n_uc <- read.xlsx(xlsxFile = "output/GSE37283/tables/n_uc_tables.xlsx", sheet = "full_gsea", rowNames = T)
pekow.uc_c <- read.xlsx(xlsxFile = "output/GSE37283/tables/uc_c_tables.xlsx", sheet = "full_gsea", rowNames = T)

# GSE47908 Bjerrum
bjerrum.d_c <- read.xlsx(xlsxFile = "output/GSE47908/tables/d_c_tables.xlsx", sheet = "full_gsea", rowNames = T)
bjerrum.d_uc <- read.xlsx(xlsxFile = "output/GSE47908/tables/d_uc_tables.xlsx", sheet = "full_gsea", rowNames = T)
bjerrum.uc_c <- read.xlsx(xlsxFile = "output/GSE47908/tables/uc_c_tables.xlsx", sheet = "full_gsea", rowNames = T)
bjerrum.dys_pan <- read.xlsx(xlsxFile = "output/GSE47908/tables/dys_pan_tables.xlsx", sheet = "full_gsea", rowNames = T)

# GSE31106 - Tang
tang.a_c <- read.xlsx(xlsxFile = "output/GSE31106/tables/a_c_tables.xlsx", sheet = "full_gsea", rowNames = T)
tang.a_i <- read.xlsx(xlsxFile = "output/GSE31106/tables/a_i_tables.xlsx", sheet = "full_gsea", rowNames = T)
tang.d_i <- read.xlsx(xlsxFile = "output/GSE31106/tables/d_i_tables.xlsx", sheet = "full_gsea", rowNames = T)
tang.a_d <- read.xlsx(xlsxFile = "output/GSE31106/tables/a_d_tables.xlsx", sheet = "full_gsea", rowNames = T)
tang.d_c <- read.xlsx(xlsxFile = "output/GSE31106/tables/d_c_tables.xlsx", sheet = "full_gsea", rowNames = T)
tang.i_c <- read.xlsx(xlsxFile = "output/GSE31106/tables/i_c_tables.xlsx", sheet = "full_gsea", rowNames = T)

# Make orthologize function ----
# Get mouse orthologs
listEnsembl()
datasets <- listDatasets(useEnsembl(biomart = "genes", mirror = "www"))
ensembl_connection <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
attributes <- listAttributes(ensembl_connection)
filters <- listFilters(ensembl_connection)

orthologs <- getBM(
  attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"),
  filters = "with_hsapiens_homolog",
  values = TRUE,
  mart = ensembl_connection)
colnames(orthologs) <- c("mouse", "human")

orthologize <- function(input) {
  input$mouse <- rownames(input)
  df <- left_join(input, orthologs, by = "mouse")
  df <- drop_na(df)
  df <- df[!df$human == "",]
  # Note: when there are multiple possible human orthologs for the mouse gene
  # (e.g. Defa5), left_join() performed a many-to-one join, which resulted in
  # multiple rows for each matching pair
  df <- df %>%
    group_by(human) %>%
    dplyr::summarise(across(where(is.numeric), mean, na.rm = TRUE),  # Calculate mean for numeric columns
                     across(where(is.character), ~ unique(.)[1]))  # Keep the first unique entry for text columns
  df <- column_to_rownames(df, "human")
  return(df)
}
# UC/C DEGs ----
di <- process_data(digby.uc_c, id = "digby")
pe <- process_data(pekow.uc_c, id = "pekow")
bj <- process_data(bjerrum.uc_c, id = "bjerrum")
ar <- process_data(arijs.uc_c, id = "arijs")
zh <- process_data(zhang.uc_c, id = "zhang")
degs <- unique(c(rownames(pe[pe$pekow.adjP == "*", ]), rownames(bj[bj$bjerrum.adjP == "*", ]),
                 rownames(ar[ar$arijs.adjP == "*", ]), rownames(zh[zh$zhang.adjP == "*", ])))
# degs <- unique(c(rownames(bj[bj$bjerrum.adjP == "*", ]), rownames(ar[ar$arijs.adjP == "*", ])))

dif <- di[rownames(di) %in% degs, ] %>% rownames_to_column(var = "Gene")
pef <- pe[rownames(pe) %in% degs, ] %>% rownames_to_column(var = "Gene")
bjf <- bj[rownames(bj) %in% degs, ] %>% rownames_to_column(var = "Gene")
arf <- ar[rownames(ar) %in% degs, ] %>% rownames_to_column(var = "Gene")
zhf <- zh[rownames(zh) %in% degs, ] %>% rownames_to_column(var = "Gene")

all <- plyr::join_all(list(bjf, arf), by = "Gene", type = "inner")
rownames(all) <- all$Gene
all <- all[all$Gene %in% glycogenes2,]
all <- all[ , -1]
all_lfc <- all[,grep(x = colnames(all), pattern = "*.lfc", value = T)]
all_lfc[is.na(all_lfc)] <- 0
all_pval <- all[,grep(x = colnames(all), pattern = "*.adjP", value = T)]
all_pval[is.na(all_pval)] <- ""

heatmap_colors <- as.character(paletteer_c("pals::coolwarm", 100, direction = 1))
hm <- pheatmap::pheatmap(as.matrix(all_lfc),
               show_rownames = T,
               cutree_rows = 2,
               treeheight_row = 0,
               color = heatmap_colors,
               breaks = seq(-4, 4, length.out = 101),
               display_numbers = as.matrix(all_pval))

row_tree <- hm$tree_row
clusters <- cutree(row_tree, k = 2)
df <- data.frame(cluster = clusters, gene = row_tree$labels) %>% arrange(cluster)
# cluster 1 rojo 277, cluster 2 azul 144

ora1 <- gprofiler2::gost(df[df$cluster == 1, "gene"], organism = "hsapiens", sources = c("GO"), significant = TRUE)
ora2 <- gprofiler2::gost(df[df$cluster == 2, "gene"], organism = "hsapiens", sources = c("GO"), significant = TRUE)
ora1.res <- ora1$result[ora1$result$term_size < 500,]
ora2.res <- ora2$result[ora2$result$term_size < 500,]
gprofiler2::gostplot(ora1, interactive = T, capped = T)
gprofiler2::gostplot(ora2, interactive = T, capped = T)

View(ora1.res) #ora1 es up
ora1.res$query <- "up"
ora2.res$query <- "down"
ora.res <- rbind(ora1.res, ora2.res)
openxlsx::write.xlsx(ora.res,
                     file = "~/projects/tesis/final_outputs/uc_control/ora_results.xlsx", rowNames = TRUE)

# make a venn (all)
sig_arijs <- rownames(ar[ar$arijs.adjP == "*", ])
sig_bj <- rownames(bj[bj$bjerrum.adjP == "*", ])
sig_pe <- rownames(pe[pe$pekow.adjP == "*", ])
sig_zh <- rownames(zh[zh$zhang.adjP == "*", ])
gene_lists <- list(GSE47908 = sig_bj, GSE73661 = sig_arijs, GSE37283 = sig_pe, Zhang = sig_zh)
ggvenn(gene_lists, fill_color = c("skyblue", "darkgoldenrod2", "pink", "purple3"), text_size = 7, show_percentage = F)

# make a venn (glyco)
sig_arijs <- sig_arijs[sig_arijs %in% glycogenes2]
sig_bj <- sig_bj[sig_bj %in% glycogenes2]
gene_lists <- list(GSE47908 = sig_bj, GSE73661 = sig_arijs)
ggvenn(gene_lists, fill_color = c("skyblue", "darkgoldenrod2"), text_size = 7)

# ANORMAL/UC DEGs ----
pe <- process_data(pekow.n_uc, "pekow.n_uc")
bj <- process_data(bjerrum.d_uc, "bjerrum.d_uc")
wa <- process_data(watanabe.ca_pan, "watanabe.ca_pan")
bj2 <- process_data(bjerrum.dys_pan, "bjerrum.d_pan")
ta <- process_data(tang.a_i, "tang_a_i")
ta <- orthologize(ta)

degs <- unique(c(rownames(pe[pe$pekow.n_uc.adjP == "*", ]),
                 rownames(bj[bj$bjerrum.d_uc.adjP == "*", ]),
                 rownames(wa[wa$watanabe.ca_pan.adjP == "*", ]),
                 rownames(bj2[bj2$bjerrum.dys_pan.adjP == "*", ]),
                 rownames(ta[ta$tang_a_i.adjP == "*",])))

pef <- pe[rownames(pe) %in% degs, ] %>% rownames_to_column(var = "Gene")
bjf <- bj[rownames(bj) %in% degs, ] %>% rownames_to_column(var = "Gene")
waf <- wa[rownames(wa) %in% degs, ] %>% rownames_to_column(var = "Gene")
bj2f <- bj2[rownames(bj2) %in% degs, ] %>% rownames_to_column(var = "Gene")
taf <- ta[rownames(ta) %in% degs, ] %>% rownames_to_column(var = "Gene")

all.df <- plyr::join_all(list(waf, pef, bjf, bj2f, taf), by = "Gene", type = "full")
rownames(all.df) <- all.df$Gene
all.df <- all.df[ , -1]

all_lfc <- all.df[,grep(x = colnames(all.df), pattern = "*.lfc", value = T)]
all_lfc[is.na(all_lfc)] <- 0
all_pval <- as.matrix(all.df[,grep(x = colnames(all.df), pattern = "*.adjP", value = T)])
all_pval[is.na(all_pval)] <- ""

heatmap_colors <- as.character(paletteer_c("pals::coolwarm", 100, direction = 1))
hm <- pheatmap::pheatmap(as.matrix(all_lfc),
                         show_rownames = F,
                         cutree_rows = 4,
                         treeheight_row = 0,
                         color = heatmap_colors,
                         breaks = seq(-1.5, 1.5, length.out = 101),
                         display_numbers = as.matrix(all_pval))

# ANORMAL/UC NES ----
pe <- process_data_nes(pekow.n_uc, "pekow.n_uc") %>% rownames_to_column(var = "ID")
bj <- process_data_nes(bjerrum.d_uc, "bjerrum.d_uc") %>% rownames_to_column(var = "ID")
wa <- process_data_nes(watanabe.ca_pan, "watanabe.ca_pan") %>% rownames_to_column(var = "ID")
bj2 <- process_data_nes(bjerrum.dys_pan, "bjerrum.d_pan") %>% rownames_to_column(var = "ID")
#ta <- process_data_nes(tang.a_i, "tang_a_i") %>% rownames_to_column(var = "ID")

# selection <- c("GOBP_CELL_KILLING",
#                "REACTOME_CELL_CYCLE_CHECKPOINTS",
#                "REACTOME_SYNTHESIS_OF_DNA",
#                "REACTOME_S_PHASE",
#                "GOBP_OSSIFICATION",
#                "GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT",
#                "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
#                "GOBP_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN",
#                "GOBP_SISTER_CHROMATID_SEGREGATION",
#                "GOBP_LEUKOCYTE_MIGRATION",
#                "GOBP_LEUKOCYTE_CHEMOTAXIS",
#                "GOBP_POSITIVE_REGULATION_OF_CELL_ADHESION",
#                "REACTOME_DNA_REPAIR",
#                "GOBP_NUCLEAR_CHROMOSOME_SEGREGATION",
#                "GOBP_DNA_REPLICATION",
#                "GOBP_CELLULAR_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN",
#                "GOBP_CELLULAR_RESPONSE_TO_BIOTIC_STIMULUS",
#                "GOBP_ADAPTIVE_IMMUNE_RESPONSE",
#                "GOBP_LEUKOCYTE_MEDIATED_IMMUNITY",
#                "GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY")

#pval <- 0.00005
#watanabe <- watanabe.ca_pan[watanabe.ca_pan$p.adjust < pval, ]
#bjerrum <- bjerrum.dys_pan[bjerrum.dys_pan$p.adjust < pval, ]
#wa <- process_data_nes(watanabe, "CCRAC/Pancolitis") %>% rownames_to_column(var = "ID")
#bj <- process_data_nes(bjerrum, "Displasia/Pancolitis") %>% rownames_to_column(var = "ID")

all.df <- plyr::join_all(list(wa, bj, pe, bj2), by = "ID", type = "inner")
rownames(all.df) <- all.df$ID
#all.df <- all.df[!grepl("^HP", rownames(all.df)),]
#all.df <- all.df[!grepl("^REACTOME", rownames(all.df)),]
all.df <- all.df[grepl("^GOBP", rownames(all.df)),]
all.df <- all.df[ , -1]

selection <- c("GOBP_ADAPTIVE_IMMUNE_RESPONSE", "GOBP_HUMORAL_IMMUNE_RESPONSE",
               "GOBP_COMPLEMENT_ACTIVATION", "GOBP_LEUKOCYTE_MEDIATED_IMMUNITY",
               "GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY", "GOBP_IMMUNOGLOBULIN_PRODUCTION",
               "GOBP_POSITIVE_REGULATION_OF_CELL_ACTIVATION", "GOBP_B_CELL_ACTIVATION",
               "GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION", "GOBP_GRANULOCYTE_MIGRATION",
               "GOBP_REGULATION_OF_HEMOPOIESIS", "GOBP_LEUKOCYTE_CHEMOTAXIS",
               "GOBP_MYELOID_LEUKOCYTE_MIGRATION", "GOBP_GRANULOCYTE_CHEMOTAXIS",
               "GOBP_PHAGOCYTOSIS","GOBP_NEUTROPHIL_MIGRATION", "GOBP_RESPONSE_TO_VIRUS",
               "GOBP_DEFENSE_RESPONSE_TO_SYMBIONT", "GOBP_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY",
               "GOBP_MONOCARBOXYLIC_ACID_CATABOLIC_PROCESS", "GOBP_FATTY_ACID_CATABOLIC_PROCESS")
all.df <- all.df[rownames(all.df) %in% selection,]

all_nes <- all.df[,grep(x = colnames(all.df), pattern = "*.nes", value = T)]
all_nes[is.na(all_nes)] <- 0
all_pval <- as.matrix(all.df[,grep(x = colnames(all.df), pattern = "*.adjP", value = T)])
all_pval[is.na(all_pval)] <- ""

heatmap_colors <- as.character(paletteer_c("pals::coolwarm", 100, direction = 1))
hm <- pheatmap::pheatmap(as.matrix(all_nes),
                         show_rownames = F,
                         show_colnames = F,
                         cutree_rows = 2,
                         treeheight_row = 0,
                         border_color = "white",
                         color = heatmap_colors,
                         breaks = seq(-3, 3, length.out = 101),
                         display_numbers = as.matrix(all_pval))

row_tree <- hm$tree_row
clusters <- cutree(row_tree, k = 2)
df <- data.frame(cluster = clusters, gene = row_tree$labels) %>% arrange(cluster)
table(df$cluster)
filter <- df[df$cluster == 2, 2]

# NEO/CTRL y NEO/qUC GENES ----
vs_control <- process_data(pekow.n_c, "Neoplasia/control")
vs_quc <- process_data(pekow.n_uc, "Neoplasia/qUC")

degs <- unique(c(rownames(vs_control[vs_control$`Neoplasia/control.adjP` == "*", ]),
                 rownames(vs_quc[vs_quc$`Neoplasia/qUC.adjP` == "*", ])))
# degs <- intersect(rownames(vs_control[vs_control$`Neoplasia/control.adjP` == "*", ]),
#                   rownames(vs_quc[vs_quc$`Neoplasia/qUC.adjP` == "*", ]))
# degs <- c("REG1A", "CXCL8", "REG1B", "CLC", "SLC6A14", "MMP1", "IL1B", "CXCR4",
#           "SFRP2", "MNDA", "CXCL13", "MMP10", "IL13RA2", "FDCSP", "S100A8", "MMP3",
#           "FCGR3B", "SERPINA3", "IDO1", "SELL", "CXCL1", "S100A9", "MMP12", "BCL2A1",
#           "S100A12", "REG3A", "GPR171", "TCN1", "CD38", "AQP8", "SLC26A2", "VIPR1",
#           "MUC12", "SELENBP1", "PIGZ", "CDHR1", "DHRS11", "CHP2")

vs_control.f <- vs_control[rownames(vs_control) %in% degs,] %>% rownames_to_column(var = "gen")
vs_quc.f <- vs_quc[rownames(vs_quc) %in% degs,] %>% rownames_to_column(var = "gen")

all.df <- plyr::join_all(list(vs_control.f, vs_quc.f), by = "gen", type = "full")
rownames(all.df) <- all.df$gen
all.df <- all.df[ , -1]

all_lfc <- all.df[,grep(x = colnames(all.df), pattern = "*.lfc", value = T)]
all_lfc[is.na(all_lfc)] <- 0
all_pval <- as.matrix(all.df[,grep(x = colnames(all.df), pattern = "*.adjP", value = T)])
all_pval[is.na(all_pval)] <- ""

heatmap_colors <- as.character(paletteer_c("pals::coolwarm", 100, direction = 1))
hm <- pheatmap::pheatmap(as.matrix(all_lfc),
                         show_rownames = F,
                         cutree_rows = 2,
                         treeheight_row = 0,
                         treeheight_col = 0,
                         color = heatmap_colors,
                         breaks = seq(-2, 2, length.out = 101),
                         display_numbers = as.matrix(all_pval))

# ORA
row_tree <- hm$tree_row
clusters <- cutree(row_tree, k = 2)
df <- data.frame(cluster = clusters, gene = row_tree$labels) %>% arrange(cluster)
table(df$cluster)

ora1 <- gprofiler2::gost(df[df$cluster == 1, "gene"], organism = "hsapiens", sources = c("GO"), significant = TRUE)
ora2 <- gprofiler2::gost(df[df$cluster == 2, "gene"], organism = "hsapiens", sources = c("GO"), significant = TRUE)
ora1.res <- ora1$result[ora1$result$term_size < 500,]
ora2.res <- ora2$result[ora2$result$term_size < 500,]
ora1.filtered <- list(result = ora1.res, meta = ora1$meta)
ora2.filtered <- list(result = ora2.res, meta = ora2$meta)
class(ora1.filtered) <- "gostres"
class(ora2.filtered) <- "gostres"
gprofiler2::gostplot(ora1.filtered, interactive = T, capped = T)
gprofiler2::gostplot(ora2.filtered, interactive = T, capped = T)

View(ora1.res)
View(ora2.res)

genes_test <- c("DUOX2", "PI3", "LCN2", "TIMP1", "MMP9", "CXCL9", "DMBT1",
                "S100A9", "OLFM4", "CCL2")
ora_test <- gprofiler2::gost(genes_test, organism = "hsapiens", sources = c("GO"), significant = TRUE)
ora_test.res <- ora_test$result[ora_test$result$term_size < 500,]


degs <- unique(c(rownames(vs_control[vs_control$`Neoplasia/control.adjP` == "*", ]),
                 rownames(vs_quc[vs_quc$`Neoplasia/qUC.adjP` == "*", ])))


neo_control_significant <- rownames(vs_control[vs_control$`Neoplasia/control.adjP` == "*", ])
neo_quc_significant <- rownames(vs_quc[vs_quc$`Neoplasia/qUC.adjP` == "*", ])
gene_lists <- list(Neoplasia_control = neo_control_significant,
                   Neoplasia_qUC = neo_quc_significant)
ggvenn(gene_lists, fill_color = c("skyblue", "darkgoldenrod2"), text_size = 7, show_percentage = F)

selected_sigs <- neo_quc_significant[!neo_quc_significant %in% neo_control_significant]
sixty_quc_c <- pekow.n_uc[selected_sigs, c(1,3,5,6)]

openxlsx::write.xlsx(sixty_quc_c, file = "~/projects/tesis/final_outputs/uc_activa_vs_quiescente/60_desregulados_nq_no_nc.xlsx", rowNames = TRUE)


# Anormal/control DEGs ----
pe <- process_data(pekow.n_c, id = "Neoplasia/Control (H)")
bj <- process_data(bjerrum.d_c, id = "Displasia/Control (H)")
ta.a <- process_data(tang.a_c, id = "Adenocarcinoma/Control (R)")
ta.a <- orthologize(ta.a)
ta.d <- process_data(tang.d_c, id = "Displasia/Control (R)")
ta.d <- orthologize(ta.d)

degs <- unique(c(rownames(pe[pe$`Neoplasia/Control (H).adjP` == "*", ]),
                 rownames(bj[bj$`Displasia/Control (H).adjP` == "*", ]),
                 rownames(ta.a[ta.a$`Adenocarcinoma/Control (R).adjP` == "*", ]),
                 rownames(ta.d[ta.d$`Displasia/Control (R).adjP` == "*", ])))

pef <- pe[rownames(pe) %in% degs, ] %>% rownames_to_column(var = "Gene")
bjf <- bj[rownames(bj) %in% degs, ] %>% rownames_to_column(var = "Gene")
ta.af <- ta.a[rownames(ta.a) %in% degs, ] %>% rownames_to_column(var = "Gene")
ta.df <- ta.d[rownames(ta.d) %in% degs, ] %>% rownames_to_column(var = "Gene")

all.df <- plyr::join_all(list(bjf, pef, ta.af, ta.df), by = "Gene", type = "inner")
rownames(all.df) <- all.df$Gene
all.df <- all.df[ , -1]

all_lfc <- all.df[,grep(x = colnames(all.df), pattern = "*.lfc", value = T)]
all_lfc[is.na(all_lfc)] <- 0
all_pval <- all.df[,grep(x = colnames(all.df), pattern = "*.adjP", value = T)]
all_pval[is.na(all_pval)] <- ""

heatmap_colors <- as.character(paletteer_c("pals::coolwarm", 100, direction = 1))
pheatmap::pheatmap(as.matrix(all_lfc),
                   show_rownames = F,
                   cutree_rows = 2,
                   treeheight_row = 0,
                   treeheight_col = 6,
                   color = heatmap_colors,
                   breaks = seq(-2, 2, length.out = 101),
                   display_numbers = as.matrix(all_pval))

# Anormal/control NES ----
pe <- process_data_nes(pekow.n_c, id = "Neoplasia/Control (H)") %>% rownames_to_column(var = "ID")
bj <- process_data_nes(bjerrum.d_c, id = "Displasia/Control (H)") %>% rownames_to_column(var = "ID")
ta.a <- process_data_nes(tang.a_c, id = "Adenocarcinoma/Control (R)") %>% rownames_to_column(var = "ID")
ta.d <- process_data_nes(tang.d_c, id = "Displasia/Control (R)") %>% rownames_to_column(var = "ID")

all.df <- plyr::join_all(list(bj, pe, ta.a, ta.d), by = "ID", type = "inner")
rownames(all.df) <- all.df$ID
all.df <- all.df[ , -1]
#selection <- rownames(all.df[grepl("^REACTOME_", rownames(all.df)),])

all_lfc <- all.df[,grep(x = colnames(all.df), pattern = "*.nes", value = T)]
all_lfc[is.na(all_lfc)] <- 0
all_pval <- all.df[,grep(x = colnames(all.df), pattern = "*.adjP", value = T)]
all_pval[is.na(all_pval)] <- ""

heatmap_colors <- as.character(paletteer_c("pals::coolwarm", 100, direction = 1))

# Hm general
pheatmap::pheatmap(as.matrix(all_lfc),
                   show_rownames = T,
                   cutree_rows = 1,
                   treeheight_row = 0,
                   treeheight_col = 5,
                   color = heatmap_colors,
                   breaks = seq(-3, 3, length.out = 101),
                   display_numbers = as.matrix(all_pval))

# Bjerrum DEGs ----
d_c <- process_data(bjerrum.d_c, id = "Displasia/Control")
d_uc <- process_data(bjerrum.d_uc, id = "Displasia/CU")
uc_c <- process_data(bjerrum.uc_c, id = "CU/Control")
dis_pan <- process_data(bjerrum.dys_pan, id = "Displasia/Pancolitis")

degs <- unique(c(rownames(d_c[d_c$`Displasia/Control.adjP` == "*", ]),
                 rownames(d_uc[d_uc$`Displasia/CU.adjP` == "*", ]),
                 rownames(uc_c[uc_c$`CU/Control.adjP` == "*", ]),
                 rownames(dis_pan[dis_pan$`Displasia/Pancolitis.adjP` == "*", ])))

d_c.f <- d_c[rownames(d_c) %in% degs, ] %>% rownames_to_column(var = "Gene")
d_uc.f <- d_uc[rownames(d_uc) %in% degs, ] %>% rownames_to_column(var = "Gene")
uc_c.f <- uc_c[rownames(uc_c) %in% degs, ] %>% rownames_to_column(var = "Gene")
dis_pan.f <- dis_pan[rownames(dis_pan) %in% degs, ] %>% rownames_to_column(var = "Gene")

all.df <- plyr::join_all(list(d_c.f, d_uc.f, uc_c.f, dis_pan.f), by = "Gene", type = "inner")
rownames(all.df) <- all.df$Gene
all.df <- all.df[ , -1]

all_lfc <- all.df[,grep(x = colnames(all.df), pattern = "*.lfc", value = T)]
all_lfc[is.na(all_lfc)] <- 0
all_pval <- all.df[,grep(x = colnames(all.df), pattern = "*.adjP", value = T)]
all_pval[is.na(all_pval)] <- ""

heatmap_colors <- as.character(paletteer_c("pals::coolwarm", 100, direction = 1))
pheatmap::pheatmap(as.matrix(all_lfc),
                   show_rownames = F,
                   cutree_rows = 2,
                   treeheight_row = 0,
                   treeheight_col = 6,
                   color = heatmap_colors,
                   breaks = seq(-1, 1, length.out = 101),
                   display_numbers = as.matrix(all_pval))

# Bjerrum NES ----
#d_c <- process_data_nes(bjerrum.d_c, id = "Displasia/Control") %>% rownames_to_column(var = "ID")
d_uc <- process_data_nes(bjerrum.d_uc, id = "Displasia/CU") %>% rownames_to_column(var = "ID")
#uc_c <- process_data_nes(bjerrum.uc_c, id = "CU/Control") %>% rownames_to_column(var = "ID")
dis_pan <- process_data_nes(bjerrum.dys_pan, id = "Displasia/Pancolitis") %>% rownames_to_column(var = "ID")

# choose important nes
my_vector <- bjerrum.d_uc$p.adjust
names(my_vector) <- rownames(bjerrum.d_uc)
my_vector <- sort(my_vector)
my_vector <- my_vector[grepl("^GOBP_", names(my_vector))]
selected_terms <- names(my_vector[1:30])

all.df <- plyr::join_all(list(d_uc, dis_pan), by = "ID", type = "full")
rownames(all.df) <- all.df$ID
#all.df <- all.df[!grepl("^HP_", rownames(all.df)),]
#all.df <- all.df[grepl("^GOBP_", rownames(all.df)),]
#terms <- rownames(all.df)
# ix <- c(which(terms %in% c("GOBP_AEROBIC_RESPIRATION",
#                            "GOBP_SENSORY_PERCEPTION_OF_CHEMICAL_STIMULUS",
#                            "GOBP_LEUKOCYTE_MEDIATED_IMMUNITY",
#                            "GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION",
#                            "GOBP_REGULATION_OF_ADAPTIVE_IMMUNE_RESPONSE",
#                            "GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY",
#                            "GOBP_ADAPTIVE_IMMUNE_RESPONSE",
#                            "GOBP_POSITIVE_REGULATION_OF_RESPONSE_TO_BIOTIC_STIMULUS",
#                            "GOBP_INTERLEUKIN_2_PRODUCTION",
#                            "GOBP_REGULATION_OF_INNATE_IMMUNE_RESPONSE",
#                            "GOBP_PROTEIN_EXIT_FROM_ENDOPLASMIC_RETICULUM",
#                            "GOBP_REGULATION_OF_RESPONSE_TO_BIOTIC_STIMULUS",
#                            "GOBP_CELLULAR_RESPONSE_TO_INTERFERON_GAMMA",
#                            "GOBP_MYD88_DEPENDENT_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
#                            "GOBP_PROTEIN_FOLDING_IN_ENDOPLASMIC_RETICULUM",
#                            "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
#                            "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
#                            "GOBP_RESPONSE_TO_VIRUS",
#                            "GOBP_DEFENSE_RESPONSE_TO_SYMBIONT",
#                            "GOBP_CELL_KILLING",
#                            "GOBP_ESTABLISHMENT_OF_PROTEIN_LOCALIZATION_TO_ENDOPLASMIC_RETICULUM",
#                            "GOBP_PHAGOCYTOSIS",
#                            "GOBP_REGULATION_OF_INFLAMMATORY_RESPONSE",
#                            "GOBP_PROTEIN_LOCALIZATION_TO_CHROMOSOME",
#                            "GOBP_POSITIVE_REGULATION_OF_PROTEIN_LOCALIZATION_TO_NUCLEUS",
#                            "GOBP_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN",
#                            "GOBP_POSITIVE_REGULATION_OF_PROTEOLYSIS",
#                            "GOBP_PROTEIN_N_LINKED_GLYCOSYLATION",
#                            "GOBP_I_KAPPAB_KINASE_NF_KAPPAB_SIGNALING")))
#selected_terms <- terms[ix]
all.df <- all.df[rownames(all.df) %in% selected_terms,]
all.df <- all.df[ , -1]

all_lfc <- all.df[,grep(x = colnames(all.df), pattern = "*.nes", value = T)]
all_lfc[is.na(all_lfc)] <- 0
all_pval <- all.df[,grep(x = colnames(all.df), pattern = "*.adjP", value = T)]
all_pval[is.na(all_pval)] <- ""

heatmap_colors <- as.character(paletteer_c("pals::coolwarm", 100, direction = 1))
pheatmap::pheatmap(as.matrix(all_lfc),
                   show_rownames = T,
                   show_colnames = F,
                   cutree_rows = 1,
                   treeheight_row = 0,
                   treeheight_col = 0,
                   #fontsize_row = 5, 
                   color = heatmap_colors,
                   breaks = seq(-3, 3, length.out = 101),
                   display_numbers = as.matrix(all_pval))

# General NES ----
d_uc <- process_data_nes(bjerrum.d_uc, "Displasia/L-UC") %>% rownames_to_column(var = "ID")
d_c <- process_data_nes(bjerrum.d_c, "Displasia/Control") %>% rownames_to_column(var = "ID")
d_pan <- process_data_nes(bjerrum.dys_pan, "Displasia/Pancolitis") %>% rownames_to_column(var = "ID")
luc_c <- process_data_nes(bjerrum.uc_c, "L-UC/Control") %>% rownames_to_column(var = "ID")
CACRC_pan <- process_data_nes(watanabe.ca_pan, "CACRC/Pancolitis") %>% rownames_to_column(var = "ID")
n_c <- process_data_nes(pekow.n_c, "Neoplasia/Control") %>% rownames_to_column(var = "ID")
n_uc <- process_data_nes(pekow.n_uc, "Neoplasia/qUC") %>% rownames_to_column(var = "ID")
#uc_c <- process_data_nes(pekow.uc_c, "qUC/Control") %>% rownames_to_column(var = "ID")
#a_c <- process_data_nes(tang.a_c, "Adenocarcinoma/Control R") %>% rownames_to_column(var = "ID")
#a_i <- process_data_nes(tang.a_i, "Adenocarcinoma/Inflamado R") %>% rownames_to_column(var = "ID")
#d_i <- process_data_nes(tang.d_i, "Displasia/Inflamado R") %>% rownames_to_column(var = "ID")
#d_c_m <- process_data_nes(tang.d_c, "Displasia/Control R") %>% rownames_to_column(var = "ID")
#a_d <- process_data_nes(tang.a_d, "Adenocarcinoma/Displasia R") %>% rownames_to_column(var = "ID")
#i_c <- process_data_nes(tang.i_c, "Inflamado/Control R") %>% rownames_to_column(var = "ID")

# ALL
all.df <- plyr::join_all(list(d_uc, d_c, d_pan, luc_c, CACRC_pan, n_c, n_uc), by = "ID", type = "full")
all.df <- plyr::join_all(list(d_uc, d_c, d_pan, luc_c, CACRC_pan, n_c, n_uc,
                              uc_c, a_c, a_i, d_i, d_c_m, a_d, i_c),
                         by = "ID", type = "inner")
# HUMAN – reportado en tesis
all.df <- plyr::join_all(list(d_uc, d_pan, CACRC_pan, n_uc), by = "ID", type = "inner")

# MOUSE
all.df <- plyr::join_all(list(a_c, a_i, d_i, d_c_m, a_d, i_c), by = "ID", type = "inner")

# DYS HUMANO
all.df <- plyr::join_all(list(d_uc, d_pan), by = "ID", type = "full")

rownames(all.df) <- all.df$ID
all.df <- all.df[ , -1]
# all.df <- all.df[!grepl("^HP_", rownames(all.df)),]
all.df <- all.df[grepl("^GOBP_", rownames(all.df)),]
#all.df <- all.df[grepl("REACTOME_", rownames(all.df)),]
selection <- c("GOBP_ADAPTIVE_IMMUNE_RESPONSE", "GOBP_HUMORAL_IMMUNE_RESPONSE",
               "GOBP_COMPLEMENT_ACTIVATION", "GOBP_LEUKOCYTE_MEDIATED_IMMUNITY",
               "GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY", "GOBP_IMMUNOGLOBULIN_PRODUCTION",
               "GOBP_POSITIVE_REGULATION_OF_CELL_ACTIVATION", "GOBP_B_CELL_ACTIVATION",
               "GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION", "GOBP_GRANULOCYTE_MIGRATION",
               "GOBP_REGULATION_OF_HEMOPOIESIS", "GOBP_LEUKOCYTE_CHEMOTAXIS",
               "GOBP_MYELOID_LEUKOCYTE_MIGRATION", "GOBP_GRANULOCYTE_CHEMOTAXIS",
               "GOBP_PHAGOCYTOSIS","GOBP_NEUTROPHIL_MIGRATION", "GOBP_RESPONSE_TO_VIRUS",
               "GOBP_DEFENSE_RESPONSE_TO_SYMBIONT", "GOBP_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY",
               "GOBP_MONOCARBOXYLIC_ACID_CATABOLIC_PROCESS", "GOBP_FATTY_ACID_CATABOLIC_PROCESS")
all.df <- all.df[rownames(all.df) %in% selection,]

all_nes <- all.df[,grep(x = colnames(all.df), pattern = "*.nes", value = T)]
all_nes[is.na(all_nes)] <- 0
all_pval <- as.matrix(all.df[,grep(x = colnames(all.df), pattern = "*.adjP", value = T)])
all_pval[is.na(all_pval)] <- ""

heatmap_colors <- as.character(paletteer_c("pals::coolwarm", 100, direction = 1))
hm <- pheatmap::pheatmap(as.matrix(all_nes),
                         show_rownames = T,
                         cutree_rows = 2,
                         treeheight_row = 0,
                         treeheight_col = 0,
                         color = heatmap_colors,
                         breaks = seq(-3, 3, length.out = 101),
                         display_numbers = as.matrix(all_pval))

d_uc_h <- process_data_nes(bjerrum.d_uc, "Displasia/L-UC (H)") %>% rownames_to_column(var = "ID")
d_c_h <- process_data_nes(bjerrum.d_c, "Displasia/Control (H)") %>% rownames_to_column(var = "ID")
d_i_m <- process_data_nes(tang.d_i, "Displasia/Inflamado (R)") %>% rownames_to_column(var = "ID")
d_c_m <- process_data_nes(tang.d_c, "Displasia/Control (R)") %>% rownames_to_column(var = "ID")

all.df <- plyr::join_all(list(d_uc_h, d_c_h, d_i_m, d_c_m), by = "ID", type = "inner")
rownames(all.df) <- all.df$ID
all.df <- all.df[ , -1]
all.df <- all.df[!grepl("^HP_", rownames(all.df)),]

all_nes <- all.df[,grep(x = colnames(all.df), pattern = "*.nes", value = T)]
all_nes[is.na(all_nes)] <- 0
all_pval <- as.matrix(all.df[,grep(x = colnames(all.df), pattern = "*.adjP", value = T)])
all_pval[is.na(all_pval)] <- ""

heatmap_colors <- as.character(paletteer_c("pals::coolwarm", 100, direction = 1))
hm <- pheatmap::pheatmap(as.matrix(all_nes),
                         show_rownames = F,
                         cutree_rows = 1,
                         treeheight_row = 0,
                         treeheight_col = 15,
                         color = heatmap_colors,
                         breaks = seq(-3, 3, length.out = 101),
                         display_numbers = as.matrix(all_pval))

# UC to Pancolitis to CACRC significant ----
load("data/rdata/glycogenes.RData")
load("~/projects/tesis/data/pancolitis/pan_uc_toptable.RData")
load("~/projects/tesis/data/pancolitis/cacrc_uc_toptable.RData")
cacrc_pan <- cac_uc_toptable[,c(1,3,5)]
pan_luc <- pan_uc_toptable[,c(1,3,5)]

cacrc_pan.up <- cacrc_pan[cacrc_pan$direction == "up", -3]
# cacrc_pan.up <- cacrc_pan.up[cacrc_pan.up$adj.P.Val < 0.05,]
cacrc_pan.down <- cacrc_pan[cacrc_pan$direction == "down", -3]
# cacrc_pan.down <- cacrc_pan.down[cacrc_pan.down$adj.P.Val < 0.05,]

pan_luc.up <- pan_luc[pan_luc$direction == "up",-3]
# pan_luc.up <- pan_luc.up[pan_luc.up$adj.P.Val < 0.05,]
pan_luc.down <- pan_luc[pan_luc$direction == "down",-3]
# pan_luc.down <- pan_luc.down[pan_luc.down$adj.P.Val < 0.05,]

# Activados en CACRC vs. L-UC
upregulated <- intersect(rownames(cacrc_pan.up), rownames(pan_luc.up))
ggvenn(list("+ CCRAC/Pancolitis" = rownames(cacrc_pan.up), "+ Pancolitis/L-UC" = rownames(pan_luc.up)), text_size = 10,fill_color = c("red", "red4"),show_percentage = F)
# upregulated_ora <- gprofiler2::gost(upregulated, organism = "hsapiens", sources = c("GO:BP", "REACTOME"), significant = TRUE)
# upregulated_ora.res <- upregulated_ora$result[upregulated_ora$result$term_size < 500,]
# gprofiler2::gostplot(upregulated_ora)

# Downregulados en CACRC vs. L-UC
downregulated <- intersect(rownames(cacrc_pan.down), rownames(pan_luc.down))
ggvenn(list("- CCRAC/Pancolitis" = rownames(cacrc_pan.down), "- Pancolitis/L-UC" = rownames(pan_luc.down)), text_size = 10,fill_color = c("blue", "blue4"),show_percentage = F)
# downregulated_ora <- gprofiler2::gost(downregulated, organism = "hsapiens", sources = c("GO:BP", "REACTOME"), significant = TRUE)
# downregulated_ora.res <- downregulated_ora$result[downregulated_ora$result$term_size < 500,]
# gprofiler2::gostplot(downregulated_ora)

upregulated.g <- upregulated[upregulated %in% glycogenes2]
downregulated.g <- downregulated[downregulated %in% glycogenes2]
length(upregulated.g)
length(downregulated.g)
length(upregulated)
length(downregulated)

# Save genelist
cacrc_luc_genelist <- list("up" = upregulated, "down" = downregulated)
save(cacrc_luc_genelist, file = "~/projects/tesis/output/integration/cacrc_luc_genelist.RData")
options(max.print = 50)
capture.output(print(list("upregulated" = paste0(upregulated, collapse = " "),
                          "downregulated" = paste0(downregulated, collapse = " "))),
               file = "~/projects/tesis/output/integration/cacrc_luc_genes.txt")
genes <- list(up = upregulated, down = downregulated)
genes.g <- list(up = upregulated.g, down = downregulated.g)

ora.up <- gprofiler2::gost(upregulated, organism = "hsapiens", sources = c("GO", significant = TRUE))
ora.up.res <- ora.up$result[ora.up$result$term_size < 500,]
ora.up.res$query <- "up"
ora.up.filtered <- list(result = ora.up.res, meta = ora.up$meta)
class(ora.up.filtered) <- "gostres"
gprofiler2::gostplot(ora.up.filtered, interactive = T, capped = T)

View(ora.up.res)

ora.down <- gprofiler2::gost(downregulated, organism = "hsapiens", sources = c("GO", significant = TRUE))
ora.down.res <- ora.down$result[ora.down$result$term_size < 500,]
ora.down.filtered <- list(result = ora.down.res, meta = ora.down$meta)
class(ora.down.filtered) <- "gostres"
gprofiler2::gostplot(ora.down.filtered, interactive = T, capped = T)
ora.down.res$query <- "down"
View(ora.down.res)

# save sheets
openxlsx::write.xlsx(list("UP" = ora.up.res, "DOWN" = ora.down.res),
                     file = "~/projects/tesis/final_outputs/ccrac_vs_uc/ORA.xlsx",
                     rowNames = FALSE)

# Select mice GSE31106 - Tang
tang.a_i <- read.xlsx(xlsxFile = "output/GSE31106/tables/a_i_tables.xlsx", sheet = "all_significant", rowNames = T)

orthologs <- getBM(
  attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"),
  filters = "with_hsapiens_homolog",
  values = TRUE,
  mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl"))
colnames(orthologs) <- c("mouse", "human")
orthologs.up <- orthologs[orthologs$human %in% genes$up,]$mouse
orthologs.down <- orthologs[orthologs$human %in% genes$down,]$mouse

# volcano
mouse_up.df <- rownames(tang.a_i[tang.a_i$direction == "up",])
mouse_down.df <- rownames(tang.a_i[tang.a_i$direction == "down",])
coregulated.up <- mouse_up.df[mouse_up.df %in% orthologs.up]
coregulated.down <- mouse_down.df[mouse_down.df %in% orthologs.down]
coregulated.df <- tang.a_i[rownames(tang.a_i) %in% c(coregulated.up, coregulated.down),]
coregulated.df <- rownames_to_column(coregulated.df, var = "geneID")

glycogenes2.g <- orthologs[orthologs$human %in% glycogenes2,]$mouse
glycogenes2.g <- glycogenes2.g[!is.na(glycogenes2.g)]
coregulated.df.g <- coregulated.df[coregulated.df$geneID %in% glycogenes2.g,]
table(coregulated.df[coregulated.df$adj.P.Val < 0.05,]$direction)
table(coregulated.df.g[coregulated.df.g$adj.P.Val < 0.05,]$direction)

my_data <- coregulated.df.g
lfc <- 2
pval <- 0.02
to_show <- my_data[my_data$adj.P.Val < pval & abs(my_data$logFC) > lfc,]$geneID
ggplot(my_data) +
  aes(y = -log10(adj.P.Val), x = logFC, text = paste("Symbol:", geneID), customdata = geneID) +
  geom_point(size = 5, color = ifelse(!my_data$geneID %in% to_show, "grey", ifelse(my_data$logFC > 0, "red3", "blue3"))) +
  geom_hline(yintercept = -log10(pval), linetype = "longdash", colour = "grey", linewidth = 1) +
  geom_vline(xintercept = lfc, linetype = "longdash", colour = "#BE684D", linewidth = 1) +
  geom_vline(xintercept = -lfc, linetype = "longdash", colour = "#2C467A", linewidth = 1) +
  geom_text_repel(data = subset(my_data, geneID %in% to_show), aes(label = geneID), size = 8) +
  theme_bw(base_size = 25)

example.df <- my_data[,c(1,2,4,6)]
colnames(example.df) <- c("gen", "logFC", "p", "sentido")

openxlsx::write.xlsx(coregulated.df.g[,!colnames(coregulated.df.g) %in% "rank"],
                     file = "~/projects/tesis/final_outputs/firma/firma.xlsx", rowNames = TRUE)

save(example.df, file = "~/projects/tesis/shiny_app/degs/example.df")

# Sporadic vs CACRC (Watanabe) ----
# Load data
cac_uc <- read.xlsx(xlsxFile = "output/GSE3629/tables/cac_uc_tables.xlsx", sheet = "all_significant", rowNames = T)
sporadic_uc <- read.xlsx(xlsxFile = "output/GSE3629/tables/sporadic_uc_tables.xlsx", sheet = "all_significant", rowNames = T)
cac_sporadic <- read.xlsx(xlsxFile = "output/GSE3629/tables/cac_sporadic_tables.xlsx", sheet = "all_significant", rowNames = T)
cac_genes <- rownames(cac_uc)
sporadic_genes <- rownames(sporadic_uc)
differential_genes <- rownames(cac_sporadic)
gplots::venn(list(cac_genes, sporadic_genes))

# Genes únicos de CAC
cac_unique <- setdiff(cac_genes, sporadic_genes)
length(cac_unique)
length(cac_unique[cac_unique %in% glycogenes])

# Genes únicos de Sporadic
sporadic_unique <- setdiff(sporadic_genes, cac_genes)
length(sporadic_unique)
length(sporadic_unique[sporadic_unique %in% glycogenes2])

# Genes de ambos tipos de cáncer
intersection_genes <- intersect(cac_genes, sporadic_genes)
length(intersection_genes)
length(intersection_genes[intersection_genes %in% glycogenes])

# Genes de ambos tipos de cáncer que se expresan diferencialmente
opposite_genes <- intersect(intersection_genes, differential_genes)
length(opposite_genes)
shared.df <- cac_sporadic[rownames(cac_sporadic) %in% opposite_genes,]
table(shared.df$direction) # 34 down, 20 up
opposite_glyco <- opposite_genes[opposite_genes %in% glycogenes]
length(opposite_glyco)
shared.dfg <- cac_sporadic[rownames(cac_sporadic) %in% opposite_glyco,]
length(opposite_glyco)
table(shared.dfg$direction) # 1 down, 1 up

# Genes de ambos tipos de cáncer que no se expresan diferencialmente
not_opposite_genes <- setdiff(intersection_genes, differential_genes)
length(not_opposite_genes)
length(not_opposite_genes[not_opposite_genes %in% glycogenes])

length(firma)

firma_extendida <- scan("genes.csv", sep=',', what = "", quiet = TRUE)

firma_extendida[firma_extendida %in% cac_unique]
firma_extendida[firma_extendida %in% sporadic_unique]
firma_extendida[firma_extendida %in% intersection_genes]
firma_extendida[firma_extendida %in% opposite_genes]
firma_extendida[firma_extendida %in% not_opposite_genes]

# Otro intento
cac_uc.g <- cac_uc[rownames(cac_uc) %in% glycogenes2,]
sporadic_uc.g <- sporadic_uc[rownames(sporadic_uc) %in% glycogenes2,]
cac_sporadic.g <- cac_sporadic[rownames(cac_sporadic) %in% glycogenes2,]
cac_genes.g <- rownames(cac_uc.g)
sporadic_genes.g <- rownames(sporadic_uc.g)
differential_genes.g <- rownames(cac_sporadic.g)
gplots::venn(list(firma_extendida, cac_genes.g, sporadic_genes.g))

setdiff(firma_extendida, c(cac_genes.g, sporadic_genes.g))
intersect(sporadic_genes.g, cac_genes.g)

# filtro tabla cac/esporádico
df <- cac_sporadic[rownames(cac_sporadic) %in% firma_extendida, c(1,3,5)]
to_save <- paste0(glycogenes2, collapse = ", ")

# Mouse ----
tang.a_i <- read.xlsx(xlsxFile = "output/GSE31106/tables/a_i_tables.xlsx", sheet = "all_significant", rowNames = T)

orthologs <- getBM(
  attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"),
  filters = "with_hsapiens_homolog",
  values = TRUE,
  mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl"))
colnames(orthologs) <- c("mouse", "human")
orthologs <- orthologs[orthologs$human %in% c(genelist$up, genelist$down),]

tang_selected <- tang.a_i[orthologs$mouse, c(1,3)]
tang_selected$genes <- orthologs$human
tang_selected <- tang_selected[!is.na(tang_selected$logFC),]

tang_selected.up <- tang_selected[tang_selected$genes %in% genelist$up,]
tang_selected.up <- tang_selected.up[tang_selected.up$logFC > 0, ]
tang_selected.down <- tang_selected[tang_selected$genes %in% genelist$down,]
tang_selected.down <- tang_selected.down[tang_selected.down$logFC < 0, ]
tang_selected.up.significant <- tang_selected.up[tang_selected.up$adj.P.Val < 0.05, ]
tang_selected.down.significant <- tang_selected.down[tang_selected.down$adj.P.Val < 0.05, ]

pre_signature <- c(rownames(tang_selected.up), rownames(tang_selected.down))
#to_remove <- pre_signature[grepl("\\.", pre_signature)]
#pre_signature <- pre_signature[!pre_signature %in% to_remove]

filtered_mouse <- tang.a_i[pre_signature,]
potential_genes <- rownames(filtered_mouse)

# Volcano
up_n <- sum(filtered_mouse$direction == "up")
down_n <- sum(filtered_mouse$direction == "down")
filtered_mouse <- rownames_to_column(filtered_mouse, var = "geneID")
lfc <- 0.9
p <- 0.05
to_show <- filtered_mouse %>% dplyr::filter(adj.P.Val <= p & abs(logFC) >= lfc)
to_show <- to_show$geneID
#to_show <- c(to_show, "Sell", "Slc26a2", "Lgals4")

plot <- ggplot(filtered_mouse) +
  aes(y = -log10(adj.P.Val), x = logFC, text = paste("Symbol:", geneID), customdata = geneID) +
  geom_hline(yintercept = -log10(p), linetype = "longdash", colour = "grey", linewidth = 1) +
  geom_vline(xintercept = lfc, linetype = "longdash", colour = "#BE684D", linewidth = 1) +
  geom_vline(xintercept = -lfc, linetype = "longdash", colour = "#2C467A", linewidth = 1) +
  geom_point(size = 2, color = ifelse(!filtered_mouse$geneID %in% to_show, "grey", ifelse(filtered_mouse$logFC > 0, "red3", "blue3"))) +
  geom_text_repel(data = subset(filtered_mouse, geneID %in% to_show), aes(label = geneID), vjust = 0) +
  labs(subtitle = paste0("p-valor: ", p, ", logFC = ", lfc),
       caption = paste0("DEGs significativos: ", up_n, " (upregulated), ", down_n, " (downregulated).")) +
  theme_bw()

print(plot)

# Dónde está PANX1?----
"PANX1" %in% uc_c_genes$down
"PANX1" %in% uc_c_genes$up
"PANX1" %in% cac_sporadic_genes$down #TRUE
"PANX1" %in% cac_sporadic_genes$up
"PANX1" %in% anormal_uc_genes$down
"PANX1" %in% anormal_uc_genes$up

digby.uc_c.panx1 <- as_tibble(digby.uc_c[rownames(digby.uc_c) == "PANX1",])
watanabe.cac_sporadic.panx1 <- as_tibble(watanabe.cac_sporadic[rownames(watanabe.cac_sporadic) == "PANX1",])
watanabe.cac_uc.panx1 <- as_tibble(watanabe.cac_uc[rownames(watanabe.cac_uc) == "PANX1",])
pekow.n_c.panx1 <- as_tibble(pekow.n_c[rownames(pekow.n_c) == "PANX1",])
pekow.n_uc.panx1 <- as_tibble(pekow.n_uc[rownames(pekow.n_uc) == "PANX1",])
pekow.uc_c.panx1 <- as_tibble(pekow.uc_c[rownames(pekow.uc_c) == "PANX1",])
bjerrum.d_c.panx1 <- as_tibble(bjerrum.d_c[rownames(bjerrum.d_c) == "PANX1",])
bjerrum.d_uc.panx1 <- as_tibble(bjerrum.d_uc[rownames(bjerrum.d_uc) == "PANX1",])
bjerrum.uc_c.panx1 <- as_tibble(bjerrum.uc_c[rownames(bjerrum.uc_c) == "PANX1",])
arijs.uc_c.panx1 <- as_tibble(arijs.uc_c[rownames(arijs.uc_c) == "PANX1",])
zhang.uc_c.panx1 <- as_tibble(zhang.uc_c[rownames(zhang.uc_c) == "PANX1",])
tang.a_c.panx1 <- as_tibble(tang.a_c[rownames(tang.a_c) == "Panx1",])
tang.a_i.panx1 <- as_tibble(tang.a_i[rownames(tang.a_i) == "Panx1",])
tang.d_i.panx1 <- as_tibble(tang.d_i[rownames(tang.d_i) == "Panx1",])
tang.a_d.panx1 <- as_tibble(tang.a_d[rownames(tang.a_d) == "Panx1",])
tang.d_c.panx1 <- as_tibble(tang.d_c[rownames(tang.d_c) == "Panx1",])
tang.i_c.panx1 <- as_tibble(tang.i_c[rownames(tang.i_c) == "Panx1",])

# No está en zhang ni arijs - los quito
panx1_human <- rbind(digby.uc_c.panx1,
                     pekow.uc_c.panx1,
                     bjerrum.uc_c.panx1,
                     watanabe.cac_sporadic.panx1,
                     watanabe.cac_uc.panx1,
                     bjerrum.d_uc.panx1,
                     pekow.n_uc.panx1,
                     pekow.n_c.panx1,
                     bjerrum.d_c.panx1)
panx1_human$category <- c(rep("UC_C", 3), "CAC_sporadic", rep("anormal_uc", 3), rep("anormal_control", 2))
panx1_human <- panx1_human %>% mutate(comparison = fct_reorder(comparison, logFC, .desc = TRUE))

panx1_mouse <- rbind(rbind(tang.i_c.panx1,
                           tang.a_c.panx1,
                           tang.a_i.panx1,
                           tang.a_d.panx1,
                           tang.d_c.panx1,
                           tang.d_i.panx1))
panx1_mouse <- panx1_mouse %>% mutate(comparison = fct_reorder(comparison, logFC, .desc = TRUE))

p1 <- ggplot(panx1_human, aes(x = logFC, y = comparison, col = category, label = paste("p =", round(adj.P.Val, digits = 4)))) +
  geom_pointrange(aes(xmin = min(logFC) - 1, xmax = max(logFC) + 1)) +
  geom_vline(xintercept = 0, linetype = "longdash", colour = "grey", linewidth = 1) +
  geom_text_repel() + 
  theme_minimal()

p2 <- ggplot(panx1_mouse, aes(x = logFC, y = comparison, label = paste("p =", round(adj.P.Val, digits = 4)))) +
  geom_pointrange(aes(xmin = min(logFC) - 1, xmax = max(logFC) + 1)) +
  geom_vline(xintercept = 0, linetype = "longdash", colour = "grey", linewidth = 1) +
  geom_text_repel() +
  theme_minimal()

panx1_plot <- grid.arrange(p1, p2, top = "Differential expression of Panx1 across all cohorts (top: human, bottom: mouse)")
ggsave(plot = panx1_plot, filename = "panx1.pdf", device = "pdf", dpi = 300,
       height = 10, width = 8, path = "~/projects/tesis/output")

# Donde está GCNT2? ----
# UC / C
GSE47908.uc_c <- read.xlsx(xlsxFile = "output/GSE47908/tables/uc_c_tables.xlsx", sheet = "all_genes")
GSE37283.uc_c <- read.xlsx(xlsxFile = "output/GSE37283/tables/uc_c_tables.xlsx", sheet = "all_genes")
roman.uc_c <- read.xlsx(xlsxFile = "output/roman/uc_c_tables.xlsx", sheet = "all")
colnames(GSE47908.uc_c) <- c("genes", "logFC", "P", "adjP", "rank", "direction", "comparison")
colnames(GSE37283.uc_c) <- c("genes", "logFC", "P", "adjP", "rank", "direction", "comparison")
colnames(roman.uc_c) <- c("genes", "logFC", "adjP")
GSE47908.uc_c.gcnt2 <- GSE47908.uc_c[GSE47908.uc_c$genes == "GCNT2",]
GSE37283.uc_c.gcnt2 <- GSE37283.uc_c[GSE37283.uc_c$genes == "GCNT2",]
roman.uc_c.gcnt2 <- roman.uc_c[roman.uc_c$genes == "GCNT2",]

# CAC/C
GSE37283.n_c <- read.xlsx(xlsxFile = "output/GSE37283/tables/n_c_tables.xlsx", sheet = "all_genes")
GSE47908.d_c <- read.xlsx(xlsxFile = "output/GSE47908/tables/d_c_tables.xlsx", sheet = "all_genes")
colnames(GSE37283.n_c) <- c("genes", "logFC", "P", "adjP", "rank", "direction", "comparison")
colnames(GSE47908.d_c) <- c("genes", "logFC", "P", "adjP", "rank", "direction", "comparison")
GSE37283.n_c.gcnt2 <- GSE37283.n_c[GSE37283.n_c$genes == "GCNT2",]
GSE47908.d_c.gcnt2 <- GSE47908.d_c[GSE47908.d_c$genes == "GCNT2",]

# CAC/UC
GSE3629.cac_uc <- read.xlsx(xlsxFile = "output/GSE3629/tables/cac_uc_tables.xlsx", sheet = "all_genes")
GSE37283.n_uc <- read.xlsx(xlsxFile = "output/GSE37283/tables/n_uc_tables.xlsx", sheet = "all_genes")
GSE47908.d_uc <- read.xlsx(xlsxFile = "output/GSE47908/tables/d_uc_tables.xlsx", sheet = "all_genes")
colnames(GSE3629.cac_uc) <- c("genes", "logFC", "P", "adjP", "rank", "direction", "comparison")
colnames(GSE37283.n_uc) <- c("genes", "logFC", "P", "adjP", "rank", "direction", "comparison")
colnames(GSE47908.d_uc) <- c("genes", "logFC", "P", "adjP", "rank", "direction", "comparison")
GSE3629.cac_uc.gcnt2 <- GSE3629.cac_uc[GSE3629.cac_uc$genes == "GCNT2",]
GSE37283.n_uc.gcnt2 <- GSE37283.n_uc[GSE37283.n_uc$genes == "GCNT2",]
GSE47908.d_uc.gcnt2 <- GSE47908.d_uc[GSE47908.d_uc$genes == "GCNT2",]

# unir
GSE47908.uc_c.gcnt2$gse <- "GSE47908"
GSE37283.uc_c.gcnt2$gse <- "GSE37283"
GSE37283.n_c.gcnt2$gse <- "GSE37283"
GSE47908.d_c.gcnt2$gse <- "GSE47908"
GSE3629.cac_uc.gcnt2$gse <- "GSE3629"
GSE37283.n_uc.gcnt2$gse <- "GSE37283"
GSE47908.d_uc.gcnt2$gse <- "GSE47908"

gcnt2_uc_c.df <- rbind(GSE47908.uc_c.gcnt2, GSE37283.uc_c.gcnt2)
gcnt2_ca_c.df <- rbind(GSE37283.n_c.gcnt2, GSE47908.d_c.gcnt2)
gcnt2_ca_uc.df <- rbind(GSE37283.n_uc.gcnt2, GSE47908.d_uc.gcnt2, GSE3629.cac_uc.gcnt2)
gcnt2_uc_c.df <- gcnt2_uc_c.df[, c(2,4,6,7,8)]
gcnt2_ca_c.df <- gcnt2_ca_c.df[, c(2,4,6,7,8)]
gcnt2_ca_uc.df <- gcnt2_ca_uc.df[, c(2,4,6,7,8)]

gcnt2_uc_c.df[1,4] <- "L-UC/HC"
gcnt2_uc_c.df[2,4] <- "qUC/HC"

g1 <- ggplot(gcnt2_uc_c.df, aes(x = logFC, y = comparison, col = gse, label = paste("p =", round(adjP, digits = 4)))) +
  geom_pointrange(aes(xmin = min(logFC) - 1, xmax = max(logFC) + 1)) +
  geom_vline(xintercept = 0, linetype = "longdash", colour = "grey", linewidth = 1) +
  scale_color_manual(values = c("GSE47908" = "cornflowerblue", "GSE37283" = "coral1")) +
  geom_text_repel() + 
  theme_minimal()

gcnt2_ca_c.df[1,4] <- "Neo/HC"
gcnt2_ca_c.df[2,4] <- "Dis/HC"

g2 <- ggplot(gcnt2_ca_c.df, aes(x = logFC, y = comparison, col = gse, label = paste("p =", round(adjP, digits = 4)))) +
  geom_pointrange(aes(xmin = min(logFC) - 1, xmax = max(logFC) + 1)) +
  geom_vline(xintercept = 0, linetype = "longdash", colour = "grey", linewidth = 1) +
  scale_color_manual(values = c("GSE47908" = "cornflowerblue", "GSE37283" = "coral1")) +
  geom_text_repel() + 
  theme_minimal()

gcnt2_ca_uc.df[1,4] <- "Neo/qUC"
gcnt2_ca_uc.df[2,4] <- "Dis/L-UC"
gcnt2_ca_uc.df[3,4] <- "CACRC/Pan"

g3 <- ggplot(gcnt2_ca_uc.df, aes(x = logFC, y = comparison, col = gse, label = paste("p =", round(adjP, digits = 4)))) +
  geom_pointrange(aes(xmin = min(logFC) - 1, xmax = max(logFC) + 1)) +
  geom_vline(xintercept = 0, linetype = "longdash", colour = "grey", linewidth = 1) +
  scale_color_manual(values = c("GSE47908" = "cornflowerblue", "GSE37283" = "coral1", "GSE3629" = "aquamarine4")) +
  geom_text_repel() + 
  theme_minimal()

gcnt2_plot <- grid.arrange(g1, g2, g3, ncol = 1, top = "Differential expression of GCNT2 across all human cohorts")
plot(gcnt2_plot)


# keep: L-UC/HC, Dis/HC, Neo/qUC (significant)
significant <- rbind(gcnt2_uc_c.df[1,],gcnt2_ca_c.df[2,],gcnt2_ca_uc.df[1,])
ggplot(significant, aes(x = logFC, y = comparison, col = gse, label = paste("p =", round(adjP, digits = 4)))) +
  geom_pointrange(aes(xmin = min(logFC) - 1, xmax = max(logFC) + 1)) +
  geom_vline(xintercept = 0, linetype = "longdash", colour = "grey", linewidth = 1) +
  scale_color_manual(values = c("GSE47908" = "cornflowerblue", "GSE37283" = "coral1", "GSE3629" = "aquamarine4")) +
  geom_text_repel() + 
  theme_minimal()

# Donde está GAL4? ----
# UC / C
GSE47908.uc_c <- read.xlsx(xlsxFile = "output/GSE47908/tables/uc_c_tables.xlsx", sheet = "all_genes")
GSE37283.uc_c <- read.xlsx(xlsxFile = "output/GSE37283/tables/uc_c_tables.xlsx", sheet = "all_genes")
roman.uc_c <- read.xlsx(xlsxFile = "output/roman/uc_c_tables.xlsx", sheet = "all")
colnames(GSE47908.uc_c) <- c("genes", "logFC", "P", "adjP", "rank", "direction", "comparison")
colnames(GSE37283.uc_c) <- c("genes", "logFC", "P", "adjP", "rank", "direction", "comparison")
colnames(roman.uc_c) <- c("genes", "logFC", "adjP")
GSE47908.uc_c.gal4 <- GSE47908.uc_c[GSE47908.uc_c$genes == "LGALS4",]
GSE37283.uc_c.gal4 <- GSE37283.uc_c[GSE37283.uc_c$genes == "LGALS4",]
roman.uc_c.gal4 <- roman.uc_c[roman.uc_c$genes == "LGALS4",]

# CAC/C
GSE37283.n_c <- read.xlsx(xlsxFile = "output/GSE37283/tables/n_c_tables.xlsx", sheet = "all_genes")
GSE47908.d_c <- read.xlsx(xlsxFile = "output/GSE47908/tables/d_c_tables.xlsx", sheet = "all_genes")
colnames(GSE37283.n_c) <- c("genes", "logFC", "P", "adjP", "rank", "direction", "comparison")
colnames(GSE47908.d_c) <- c("genes", "logFC", "P", "adjP", "rank", "direction", "comparison")
GSE37283.n_c.gal4 <- GSE37283.n_c[GSE37283.n_c$genes == "LGALS4",]
GSE47908.d_c.gal4 <- GSE47908.d_c[GSE47908.d_c$genes == "LGALS4",]

# CAC/UC
GSE3629.cac_uc <- read.xlsx(xlsxFile = "output/GSE3629/tables/cac_uc_tables.xlsx", sheet = "all_genes")
GSE37283.n_uc <- read.xlsx(xlsxFile = "output/GSE37283/tables/n_uc_tables.xlsx", sheet = "all_genes")
GSE47908.d_uc <- read.xlsx(xlsxFile = "output/GSE47908/tables/d_uc_tables.xlsx", sheet = "all_genes")
colnames(GSE3629.cac_uc) <- c("genes", "logFC", "P", "adjP", "rank", "direction", "comparison")
colnames(GSE37283.n_uc) <- c("genes", "logFC", "P", "adjP", "rank", "direction", "comparison")
colnames(GSE47908.d_uc) <- c("genes", "logFC", "P", "adjP", "rank", "direction", "comparison")
GSE3629.cac_uc.gal4 <- GSE3629.cac_uc[GSE3629.cac_uc$genes == "LGALS4",]
GSE37283.n_uc.gal4 <- GSE37283.n_uc[GSE37283.n_uc$genes == "LGALS4",]
GSE47908.d_uc.gal4 <- GSE47908.d_uc[GSE47908.d_uc$genes == "LGALS4",]

# unir
# roman.uc_c.gal4
GSE47908.uc_c.gal4$gse <- "GSE47908"
GSE37283.uc_c.gal4$gse <- "GSE37283"
GSE37283.n_c.gal4$gse <- "GSE37283"
GSE47908.d_c.gal4$gse <- "GSE47908"
GSE3629.cac_uc.gal4$gse <- "GSE3629"
GSE37283.n_uc.gal4$gse <- "GSE37283"
GSE47908.d_uc.gal4$gse <- "GSE47908"

gal4_uc_c.df <- rbind(GSE47908.uc_c.gal4, GSE37283.uc_c.gal4)
gal4_ca_c.df <- rbind(GSE37283.n_c.gal4, GSE47908.d_c.gal4)
gal4_ca_uc.df <- rbind(GSE37283.n_uc.gal4, GSE47908.d_uc.gal4, GSE3629.cac_uc.gal4)
gal4_uc_c.df <- gal4_uc_c.df[, c(2,4,6,7,8)]
gal4_ca_c.df <- gal4_ca_c.df[, c(2,4,6,7,8)]
gal4_ca_uc.df <- gal4_ca_uc.df[, c(2,4,6,7,8)]

g1 <- ggplot(gal4_uc_c.df, aes(x = logFC, y = comparison, col = gse, label = paste("p =", round(adjP, digits = 4)))) +
  geom_pointrange(aes(xmin = min(logFC) - 1, xmax = max(logFC) + 1)) +
  geom_vline(xintercept = 0, linetype = "longdash", colour = "grey", linewidth = 1) +
  scale_color_manual(values = c("GSE47908" = "cornflowerblue", "GSE37283" = "coral1")) +
  geom_text_repel() + 
  theme_minimal()

g2 <- ggplot(gal4_ca_c.df, aes(x = logFC, y = comparison, col = gse, label = paste("p =", round(adjP, digits = 4)))) +
  geom_pointrange(aes(xmin = min(logFC) - 1, xmax = max(logFC) + 1)) +
  geom_vline(xintercept = 0, linetype = "longdash", colour = "grey", linewidth = 1) +
  scale_color_manual(values = c("GSE47908" = "cornflowerblue", "GSE37283" = "coral1")) +
  geom_text_repel() + 
  theme_minimal()

g3 <- ggplot(gal4_ca_uc.df, aes(x = logFC, y = comparison, col = gse, label = paste("p =", round(adjP, digits = 4)))) +
  geom_pointrange(aes(xmin = min(logFC) - 1, xmax = max(logFC) + 1)) +
  geom_vline(xintercept = 0, linetype = "longdash", colour = "grey", linewidth = 1) +
  scale_color_manual(values = c("GSE47908" = "cornflowerblue", "GSE37283" = "coral1", "GSE3629" = "aquamarine4")) +
  geom_text_repel() + 
  theme_minimal()

plot(gal4_plot)
gal4_plot <- grid.arrange(g1, g2, g3, ncol = 1, top = "Differential expression of Gal4 across all cohorts")
ggsave(plot = plot.r, filename = "quality.pdf", device = "pdf", dpi = 300,
       height = 7, width = 12, path = "output/GSE37283/quality_plots")

GSE3629.cac_sporadic <- read.xlsx(xlsxFile = "output/GSE3629/tables/cac_sporadic_tables.xlsx", sheet = "all_genes")
colnames(GSE3629.cac_sporadic) <- c("genes", "logFC", "P", "adjP", "rank", "direction", "comparison")
GSE3629.cac_sporadic.gal4 <- GSE3629.cac_sporadic[GSE3629.cac_sporadic$genes == "LGALS4",]

# Analyze qpcr results ----
qpcr_data <- read_csv("~/projects/tesis/data/raw/qpcr.csv") %>%
  filter(target != "B2M") %>%
  mutate(Group = case_when(
    Sample %in% c("A1", "A2", "A3") ~ "AOM/DSS",
    Sample %in% c("C1", "C2", "C3") ~ "Control",
    Sample %in% c("D1", "D2", "D3") ~ "DSS"))

# Gen 1: Defa5
defa5 <- qpcr_data %>% filter(target == "Defa5")
ggboxplot(defa5, x = "Group", y = "Relative Quantity", title = "Expresión de Defa5", xlab = "Grupo", ylab = "Cantidad relativa")
kruskal.test(`Relative Quantity` ~ Group, data = defa5) # No paramétrico (n bajo)
dunn_test(defa5, `Relative Quantity` ~ Group, p.adjust.method = "BH", detailed = F)

# Gen 2: Lgals4
lgals4 <- qpcr_data %>% filter(target == "Lgals4")
ggboxplot(lgals4, x = "Group", y = "Relative Quantity", title = "Expresión de Lgals4", xlab = "Grupo", ylab = "Cantidad relativa")
kruskal.test(`Relative Quantity` ~ Group, data = lgals4)

# Movimiento de genes ----
df <- data.frame(
  group = c("SELL, CHI3L1, MMP1, CXCL1, REG1A, DEFA5, SLC6A14, DUOX, LCN2, LAMP3, TNIP3, LMAN1, C1S, WARS1, etc.",
            "DEFB1, SLC26A2, MT1F, CHP2, HMGCS2, AQP8, CLDN8, NXPE4, ENTREP2, CCDC183, CTSZ, etc."),
  control = c(0, 1),
  CU = c(1, 0),
  displasia = c(0, 1))

df_long <- df %>%
  pivot_longer(cols = -group, names_to = "condition", values_to = "expression")
df_long$condition <- factor(df_long$condition, levels = c("control", "CU", "displasia"))

ggplot(df_long, aes(x = condition, y = expression, group = group, color = group)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  labs(
    title = "Tendencias de expresión génica en la progresión control-CU-displasia",
    x = "Condición",
    y = "Niveles de expresión") +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),axis.text.y = element_blank())

