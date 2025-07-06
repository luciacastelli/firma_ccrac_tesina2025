# GSE37283 "Pekow" cohort - bulk microarray
# Of the patients with neoplasia, four had low-grade dysplasia (LGD), three had high-grade dysplasia (HGD), and four had adenocarcinomas.
# Load packages -----
library(affy)
library(AnnotationDbi)
library(biomaRt)
library(Biobase)
library(clusterProfiler)
library(cowplot)
library(DT)
library(enrichplot)
library(ggrepel)
library(gridExtra)
library(gplots)
library(GSEABase)
library(GSVA)
library(gt)
library(hthgu133pluspm.db)
library(limma)
library(matrixStats)
library(msigdbr)
library(pheatmap)
library(paletteer)
library(plotly)
library(RColorBrewer)
library(shiny)
library(tidyverse)
library(openxlsx)
library(qpdf)
source("scripts/functions.R")

# Set up data ----
logfc <- 2
pval <- 0.05
load("data/rdata/glycogenes.RData")
load("data/rdata/glycogenes2.RData")

# Import counts
expression <- read_tsv("data/raw/GSE37283/metadata", skip = 83, col_names = FALSE)
expression <- as.data.frame(expression[-54615,])
colnames(expression) <- expression[1,]
expression <- expression[-1,]

# Map probes to genenames
expression$ID_REF <- mapIds(hthgu133pluspm.db,
                            keys = expression$ID_REF,
                            column = "SYMBOL",
                            keytype = "PROBEID",
                            multiVals = "first")

expression <- expression[!is.na(expression$ID_REF),]
expression[, 2:ncol(expression)] <- lapply(expression[, 2:ncol(expression)],
                                           function(x) as.numeric(as.character(x)))

expression <- expression %>% group_by(ID_REF) %>% summarise_all(mean) %>% ungroup()
any(duplicated(expression$ID_REF))
geneIDs <- expression$ID_REF
expression <- expression[,-1]
expression <- data.frame(expression)
rownames(expression) <- geneIDs
expression <- expression[, colnames(expression) != "GSM915469"]

expression.t <- expression |> as_tibble(rownames = "geneID")
expression.tp <- expression.t |>
  pivot_longer(cols = colnames(expression.t[, -1]),
               names_to = "samples",
               values_to = "expression")

# Plot raw data
density_plot.r <- ggplot(data = expression.tp, aes(x = expression, color = samples)) +
  geom_density(adjust = 1.5) +
  theme_classic()

eplot.r <- ggplot(expression.tp) +
  aes(x = samples, y = expression, fill = samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", geom = "point", shape = 95, size = 10, color = "black", show.legend = FALSE) +
  labs(y = "log2 CPM", x = "sample") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_bw()

plot.r <- grid.arrange(eplot.r, density_plot.r, ncol = 2, top = "Expression data: quality check and Robust Multichip Normalization using R/Bioconductor affy package")
ggsave(plot = plot.r, filename = "quality.pdf", device = "pdf", dpi = 300,
       height = 7, width = 12, path = "output/GSE37283/quality_plots")

# Import metadata
metadata <- as.data.frame(t(read.csv("data/raw/GSE37283/metadata",sep="\t", skip = 49,nrows = 32)));metadata[1:5,1:5]
study_design <- read_tsv("data/raw/GSE37283/metadata", col_names = FALSE)[c(50, 56), 2]
study_design <- data.frame(sample = str_split(study_design[1,1], '\t'),
                       phenotype = str_split(study_design[2,1], '\t'))
colnames(study_design) <- c("sample", "phenotype")
study_design$phenotype <- factor(study_design$phenotype, labels = c("C", "qUC", "neoplasia"))
colnames(metadata) <- metadata[1,]
metadata <- metadata[-1,]
processing <- metadata[1, 19]
study_design <- study_design[study_design$sample != "GSM915469",]

# Capture samples and phenotypes
sample_names <- study_design$sample
control_samples <- study_design$sample[study_design$phenotype == "C"]
qUC_samples <- study_design$sample[study_design$phenotype == "qUC"]
neoplasia_samples <- study_design$sample[study_design$phenotype == "neoplasia"]
phenotypes <- study_design$phenotype

# PCA test ----
pca_expression <- expression # expression[, colnames(expression) != "GSM915469"]
pca_samples <- colnames(pca_expression)
pca_phenotypes <- study_design$phenotype[study_design$sampleID %in% pca_samples]
pca <- prcomp(t(pca_expression), scale. = T, center = T, retx = T)
eigenvalues <- pca$sdev^2
pc_percentages <- round((eigenvalues) / sum(eigenvalues) * 100, 1)
PCs <- paste0("PC", 1:length(pc_percentages))
pc_percentages <- data.frame(Component = PCs, Percentage = pc_percentages)
pca_data <- as_tibble(pca$x)

pca_plot <- ggplot(pca_data) +
  aes(x = PC1, y = PC2, color = pca_phenotypes) +
  geom_point(size = 3) +
  scale_color_manual(values = c("cornflowerblue", "coral1", "orchid")) +
  xlab(paste0("PC1 (", pc_percentages[1, 2], "%", ")")) +
  ylab(paste0("PC2 (", pc_percentages[2, 2], "%", ")")) +
  #geom_label_repel(aes(label = sample_names), size = 3) +
  theme_bw()

print(pca_plot)
ggsave(path = "output/GSE37283/quality_plots/",
       filename = "pca_after_sample_removal.pdf",
       pca_plot, device = "pdf", width = 6, height = 5, dpi = 300)

# Linear model ----
lm_dmatrix <- model.matrix(~ 0 + phenotypes)
colnames(lm_dmatrix) <- levels(phenotypes)
rownames(lm_dmatrix) <- sample_names
master_lm <- lmFit(expression, lm_dmatrix)
lm_cmatrix <- makeContrasts(UC_control = qUC - C,
                            neoplasia_control = neoplasia - C,
                            neoplasia_UC = neoplasia - qUC,
                            levels = lm_dmatrix)
lm <- contrasts.fit(master_lm, lm_cmatrix)
lm_stats <- eBayes(lm)

# TopTable to view DEGs and filter glycogenes -----
uc_c_toptable <- make_toptable("UC_control", "quiescent UC - control")
n_c_toptable <- make_toptable("neoplasia_control", "neoplasia - control")
n_uc_toptable <- make_toptable("neoplasia_UC", "neoplasia - quiescent UC")

# uc_c_toptable.f <- uc_c_toptable %>% dplyr::filter(adj.P.Val <= pval)
# n_c_toptable.f <- n_c_toptable %>% dplyr::filter(adj.P.Val <= pval)
# n_uc_toptable.f <- n_uc_toptable %>% dplyr::filter(adj.P.Val <= pval)

uc_c_toptable.g <- uc_c_toptable %>% dplyr::filter(row.names(uc_c_toptable) %in% glycogenes2)
n_c_toptable.g <- n_c_toptable %>% dplyr::filter(row.names(n_c_toptable) %in% glycogenes2)
n_uc_toptable.g <- n_uc_toptable %>% dplyr::filter(row.names(n_uc_toptable) %in% glycogenes2)

# uc_c_toptable.fg <- uc_c_toptable.f %>% dplyr::filter(row.names(uc_c_toptable.f) %in% glycogenes)
# n_c_toptable.fg <- n_c_toptable.f %>% dplyr::filter(row.names(n_c_toptable.f) %in% glycogenes)
# n_uc_toptable.fg <- n_uc_toptable.f %>% dplyr::filter(row.names(n_uc_toptable.f) %in% glycogenes)

# Volcano Plot: logFC/p-value ----
table <- as_tibble(uc_c_toptable.g, rownames = "geneID")
pval <- 0.05
lfc <- 2
degs <- table %>% dplyr::filter(adj.P.Val <= pval & abs(logFC) >= lfc)
degs <- degs$geneID

rownames(uc_c_toptable.g[uc_c_toptable.g$adj.P.Val <= pval,])

significant <- table %>% dplyr::filter(adj.P.Val <= pval)
significant_up <- sum(significant$direction == "up")
significant_down <- sum(significant$direction == "down")

ggplot(table) +
  aes(y = -log10(adj.P.Val), x = logFC, text = paste("Symbol:", geneID), customdata = geneID) +
  geom_hline(yintercept = -log10(pval), linetype = "longdash", colour = "grey", linewidth = 1) +
  geom_vline(xintercept = lfc, linetype = "longdash", colour = "#BE684D", linewidth = 1) +
  geom_vline(xintercept = -lfc, linetype = "longdash", colour = "#2C467A", linewidth = 1) +
  geom_point(size = 2, color = ifelse(!table$geneID %in% degs, "grey", ifelse(table$logFC > 0, "red3", "blue3"))) +
  geom_text_repel(data = subset(table, geneID %in% degs), aes(label = geneID), vjust = -1, size = 6) +
  theme_bw(base_size = 20)

volcano(uc_c_toptable, p = pval, lfc = logfc, title = "DEGs in quiescent UC vs. control",
        path = "output/GSE37283/volcano_plots/separate/", filename = "uc_c_volcano.pdf")
volcano(uc_c_toptable.g, p = pval, lfc = logfc, title = "GlycoDEGs in quiescent UC vs. control",
        path = "output/GSE37283/volcano_plots/separate/", filename = "uc_c_volcano_glyco.pdf")
volcano(n_c_toptable, p = pval, lfc = logfc, title = "DEGs in neoplasia vs. control",
        path = "output/GSE37283/volcano_plots/separate/", filename = "n_c_volcano.pdf")
volcano(n_c_toptable.g, p = pval, lfc = logfc, title = "GlycoDEGs in neoplasia vs. control",
        path = "output/GSE37283/volcano_plots/separate/", filename = "n_c_volcano_glyco.pdf")
volcano(n_uc_toptable, p = pval, lfc = logfc, title = "DEGs in neoplasia vs. quiescent UC",
        path = "output/GSE37283/volcano_plots/separate/", filename = "n_uc_volcano.pdf")
volcano(n_uc_toptable.g, p = pval, lfc = logfc, title = "GlycoDEGs in neoplasia vs. quiescent UC",
        path = "output/GSE37283/volcano_plots/separate/", filename = "n_uc_volcano_glyco.pdf")

pdfs <- list.files(path = "output/GSE37283/volcano_plots/separate/", full.names = TRUE)
qpdf::pdf_combine(input = pdfs[grepl("uc_c_", pdfs)],
                  output = "output/GSE37283/volcano_plots/uc_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("n_c_", pdfs)],
                  output = "output/GSE37283/volcano_plots/n_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("n_uc_", pdfs)],
                  output = "output/GSE37283/volcano_plots/n_uc_combined.pdf")

# Perform GSEA (gene set enrichment analysis) ----
load("data/rdata/c5_hs.RData")
load("data/rdata/reactome_hs.RData")
load("data/rdata/glyco_terms.RData")

ranked_uc_c_DEGs <- uc_c_toptable$rank
ranked_n_c_DEGs <- n_c_toptable$rank
ranked_n_uc_DEGs <- n_uc_toptable$rank
names(ranked_uc_c_DEGs) <- rownames(uc_c_toptable)
names(ranked_n_c_DEGs) <- rownames(n_c_toptable)
names(ranked_n_uc_DEGs) <- rownames(n_uc_toptable)

uc_c_gsea <- run_gsea(ranked_uc_c_DEGs, "Quiescent UC vs Control", c5_hs, reactome_hs)
uc_c_gsea.g <- uc_c_gsea %>% dplyr::filter(grepl(paste(glyco_terms, collapse = "|"), as.character(ID), ignore.case = TRUE))
n_c_gsea <- run_gsea(ranked_n_c_DEGs, "Neoplasia vs Control", c5_hs, reactome_hs)
n_c_gsea.g <- n_c_gsea %>% dplyr::filter(grepl(paste(glyco_terms, collapse = "|"), as.character(ID), ignore.case = TRUE))
n_uc_gsea <- run_gsea(ranked_n_uc_DEGs, "Neoplasia vs Quiescent UC", c5_hs, reactome_hs)
n_uc_gsea.g <- n_uc_gsea %>% dplyr::filter(grepl(paste(glyco_terms, collapse = "|"), as.character(ID), ignore.case = TRUE))

gsea_bubble(uc_c_gsea, comp = "uc_c", path = "output/GSE37283/gsea_plots/separate/", n = 30)
gsea_bubble(uc_c_gsea.g, comp = "uc_c_glyco", path = "output/GSE37283/gsea_plots/separate/", n = 30)
gsea_bubble(n_c_gsea, comp = "n_c", path = "output/GSE37283/gsea_plots/separate/", n = 30)
gsea_bubble(n_c_gsea.g, comp = "n_c_glyco", path = "output/GSE37283/gsea_plots/separate/", n = 30)
gsea_bubble(n_uc_gsea, comp = "n_uc", path = "output/GSE37283/gsea_plots/separate/", n = 30)
gsea_bubble(n_uc_gsea.g, comp = "n_uc_glyco", path = "output/GSE37283/gsea_plots/separate/", n = 30)

pdfs <- list.files(path = "output/GSE37283/gsea_plots/separate/", full.names = TRUE)
qpdf::pdf_combine(input = pdfs[grepl("uc_c_", pdfs)],
                  output = "output/GSE37283/gsea_plots/uc_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("n_c_", pdfs)],
                  output = "output/GSE37283/gsea_plots/n_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("n_uc", pdfs)],
                  output = "output/GSE37283/gsea_plots/n_uc_combined.pdf")

#keep <- selection
#selection <- selection[grepl("^GOBP_", selection),]
#my_gsea <- n_c_gsea[rownames(n_c_gsea) %in% selection,]
uc_c_gsea$source <- ifelse(grepl("^GOBP_", rownames(uc_c_gsea)), "GOBP",
                           ifelse(grepl("^GOCC_", rownames(uc_c_gsea)), "GOCC",
                                  ifelse(grepl("^REACTOME_", rownames(uc_c_gsea)), "REACTOME", 
                                         ifelse(grepl("^GOMF_", rownames(uc_c_gsea)), "GOMF", "HP"))))
my_gsea <- uc_c_gsea %>% mutate(NES_category = ifelse(NES > 0, "Positivo", "Negativo"))
my_gsea <- my_gsea[1:20,]
my_gsea$comparison <- ""
ggplot(my_gsea, aes(x = comparison, y = reorder(ID, -log10(p.adjust)))) +
  geom_point(aes(size = setSize, color = NES_category, alpha = -log10(p.adjust))) +
  scale_color_manual(values = c("Positivo" = "red", "Negativo" = "blue")) +
  labs(title = "GSEA: colitis ulcerosa quiescente vs. control", 
       x = element_blank(),
       y = "Vía",
       size = "Tamaño",
       color = "Sentido NES",
       alpha = "-log10(p-ajustado)") +
  facet_wrap(~ my_gsea$source, nrow = 1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 1, size = 14, face = "bold"))

data <- uc_c_gsea %>% dplyr::mutate(NES_category = ifelse(NES > 0, "Positive", "Negative"))
data <- data[1:50,]

ggplot(data, aes(x = comparison, y = reorder(ID, -log10(p.adjust)))) +
  geom_point(aes(size = setSize, color = NES_category, alpha = -log10(p.adjust))) +
  scale_color_manual(values = c("Positive" = "red", "Negative" = "blue")) +
  labs(x = element_blank(),
       y = "Pathway",
       size = "Gene Set Size",
       color = "NES Direction",
       alpha = "-log10(p.adjust)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 1, size = 14, face = "bold"))

# GSVA: pathway enrichment per sample: patient sub-populations ----
GO_BP <- c5_hs |> dplyr::filter(grepl("GOBP", gs_name))
GO_BP <- split(x = GO_BP$gene_symbol, f = GO_BP$gs_name)

GO_MF <- c5_hs |> dplyr::filter(grepl("GOMF", gs_name))
GO_MF <- split(x = GO_MF$gene_symbol, f = GO_MF$gs_name)

GO_CC <- c5_hs |>  dplyr::filter(grepl("GOCC", gs_name))
GO_CC <- split(x = GO_CC$gene_symbol, f = GO_CC$gs_name)

GO_RE <- reactome_hs |>  dplyr::filter(grepl("REACTOME", gs_name))
GO_RE <- split(x = GO_RE$gene_symbol, f = GO_RE$gs_name)

# Run GSVA
gsva_BP <- run_gsva(exp = as.matrix(expression), go = GO_BP, dmatrix = lm_dmatrix, cmatrix = lm_cmatrix, lfc = 0.5)  
gsva_MF <- run_gsva(exp = as.matrix(expression), go = GO_MF, dmatrix = lm_dmatrix, cmatrix = lm_cmatrix, lfc = 0.5)  
gsva_CC <- run_gsva(exp = as.matrix(expression), go = GO_CC, dmatrix = lm_dmatrix, cmatrix = lm_cmatrix, lfc = 0.5)  
gsva_RE <- run_gsva(exp = as.matrix(expression), go = GO_RE, dmatrix = lm_dmatrix, cmatrix = lm_cmatrix, lfc = 0.5)  
summary(gsva_BP$decidetests)
summary(gsva_MF$decidetests)
summary(gsva_CC$decidetests)
summary(gsva_RE$decidetests)

# Pull out enriched GSVA results --> (gsva, decidetests, coeff, samples)
uc_c_GSVA_BP <- pull_gsva(gsva = gsva_BP$res, decidetests = gsva_BP$decidetests,
                          coeff = "UC_control", samples = c(qUC_samples, control_samples))
uc_c_GSVA_MF <- pull_gsva(gsva = gsva_MF$res, decidetests = gsva_MF$decidetests,
                          coeff = "UC_control", samples = c(qUC_samples, control_samples))
uc_c_GSVA_CC <- pull_gsva(gsva = gsva_CC$res, decidetests = gsva_CC$decidetests,
                          coeff = "UC_control", samples = c(qUC_samples, control_samples))
uc_c_GSVA_RE <- pull_gsva(gsva = gsva_RE$res, decidetests = gsva_RE$decidetests,
                          coeff = "UC_control", samples = c(qUC_samples, control_samples))
uc_c_GSVA_TOTAL <- rbind(uc_c_GSVA_BP, uc_c_GSVA_MF, uc_c_GSVA_CC, uc_c_GSVA_RE)
uc_c_GSVA_TOTAL.g <- uc_c_GSVA_TOTAL[grepl(paste(glyco_terms, collapse = "|"), rownames(uc_c_GSVA_TOTAL), ignore.case = TRUE), ]

n_c_GSVA_BP <- pull_gsva(gsva = gsva_BP$res, decidetests = gsva_BP$decidetests,
                        coeff = "neoplasia_control", samples = c(neoplasia_samples, control_samples))
n_c_GSVA_MF <- pull_gsva(gsva = gsva_MF$res, decidetests = gsva_MF$decidetests,
                        coeff = "neoplasia_control", samples = c(neoplasia_samples, control_samples))
n_c_GSVA_CC <- pull_gsva(gsva = gsva_CC$res, decidetests = gsva_CC$decidetests,
                        coeff = "neoplasia_control", samples = c(neoplasia_samples, control_samples))
n_c_GSVA_RE <- pull_gsva(gsva = gsva_RE$res, decidetests = gsva_RE$decidetests,
                        coeff = "neoplasia_control", samples = c(neoplasia_samples, control_samples))
n_c_GSVA_TOTAL <- rbind(n_c_GSVA_BP, n_c_GSVA_MF, n_c_GSVA_CC, n_c_GSVA_RE)
n_c_GSVA_TOTAL.g <- n_c_GSVA_TOTAL[grepl(paste(glyco_terms, collapse = "|"), rownames(n_c_GSVA_TOTAL), ignore.case = TRUE), ]

n_uc_GSVA_BP <- pull_gsva(gsva = gsva_BP$res, decidetests = gsva_BP$decidetests,
                         coeff = "neoplasia_UC", samples = c(neoplasia_samples, qUC_samples))
n_uc_GSVA_MF <- pull_gsva(gsva = gsva_MF$res, decidetests = gsva_MF$decidetests,
                         coeff = "neoplasia_UC", samples = c(neoplasia_samples, qUC_samples))
n_uc_GSVA_CC <- pull_gsva(gsva = gsva_CC$res, decidetests = gsva_CC$decidetests,
                         coeff = "neoplasia_UC", samples = c(neoplasia_samples, qUC_samples))
n_uc_GSVA_RE <- pull_gsva(gsva = gsva_RE$res, decidetests = gsva_RE$decidetests,
                         coeff = "neoplasia_UC", samples = c(neoplasia_samples, qUC_samples))
n_uc_GSVA_TOTAL <- rbind(n_uc_GSVA_BP, n_uc_GSVA_MF, n_uc_GSVA_CC, n_uc_GSVA_RE)
n_uc_GSVA_TOTAL.g <- n_uc_GSVA_TOTAL[grepl(paste(glyco_terms, collapse = "|"), rownames(n_uc_GSVA_TOTAL), ignore.case = TRUE), ]
save(n_uc_GSVA_TOTAL, file = "~/projects/tesis/output/integration/n_uc_GSVA_TOTAL.RData")

# GSVA heatmap ----
heatmap_colors <- paletteer_c("pals::coolwarm", 100, direction = 1)
phenotype_colors <- paletteer_d("lisa::JosefAlbers", 2)
uc_c_colors <- setNames(phenotype_colors, c("C", "qUC"))
n_c_colors <- setNames(phenotype_colors, c("neoplasia", "C"))
n_uc_colors <- setNames(phenotype_colors, c("neoplasia", "qUC"))

colnames(study_design) <- c("sampleID", "phenotype")
rownames(study_design) <- NULL
make_heatmap(data = uc_c_GSVA_TOTAL, phencols = uc_c_colors,
             path = "output/GSE37283/gsva_plots/separate/", filename = "uc_c_gsva.pdf")
make_heatmap(data = uc_c_GSVA_TOTAL.g[,2:ncol(uc_c_GSVA_TOTAL.g)], phencols = uc_c_colors,
             path = "output/GSE37283/gsva_plots/separate/", filename = "uc_c_glyco_gsva.pdf")
make_heatmap(data = n_c_GSVA_TOTAL, phencols = n_c_colors,
             path = "output/GSE37283/gsva_plots/separate/", filename = "n_c_gsva.pdf")
make_heatmap(data = n_c_GSVA_TOTAL.g[,2:ncol(n_c_GSVA_TOTAL.g)], phencols = n_c_colors,
             path = "output/GSE37283/gsva_plots/separate/", filename = "n_c_glyco_gsva.pdf")
make_heatmap(data = n_uc_GSVA_TOTAL, phencols = n_uc_colors,
             path = "output/GSE37283/gsva_plots/separate/", filename = "n_uc_gsva.pdf")
make_heatmap(data = n_uc_GSVA_TOTAL.g[,2:ncol(n_uc_GSVA_TOTAL.g)], phencols = n_uc_colors,
             path = "output/GSE37283/gsva_plots/separate/", filename = "n_uc_glyco_gsva.pdf")

pdfs <- list.files(path = "output/GSE37283/gsva_plots/separate/", full.names = TRUE)
qpdf::pdf_combine(input = pdfs[grepl("uc_c_", pdfs)],
                  output = "output/GSE37283/gsva_plots/uc_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("n_c_", pdfs)],
                  output = "output/GSE37283/gsva_plots/n_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("n_uc_", pdfs)],
                  output = "output/GSE37283/gsva_plots/n_uc_combined.pdf")

# Save Excel sheets ----
#old
openxlsx::write.xlsx(list("all_genes" = uc_c_toptable,
                          "all_significant" = uc_c_toptable.f,
                          "all_glyco" = uc_c_toptable.g,
                          "significant_glyco" = uc_c_toptable.fg,
                          "full_gsea" = uc_c_gsea,
                          "glyco_gsea" = uc_c_gsea.g),
                     file = "output/GSE37283/tables/uc_c_tables.xlsx", rowNames = TRUE)

openxlsx::write.xlsx(list("all_genes" = n_c_toptable,
                          "all_significant" = n_c_toptable.f,
                          "all_glyco" = n_c_toptable.g,
                          "significant_glyco" = n_c_toptable.fg,
                          "full_gsea" = n_c_gsea,
                          "glyco_gsea" = n_c_gsea.g),
                     file = "output/GSE37283/tables/n_c_tables.xlsx", rowNames = TRUE)

openxlsx::write.xlsx(list("all_genes" = n_uc_toptable,
                          "all_significant" = n_uc_toptable.f,
                          "all_glyco" = n_uc_toptable.g,
                          "significant_glyco" = n_uc_toptable.fg,
                          "full_gsea" = n_uc_gsea,
                          "glyco_gsea" = n_uc_gsea.g),
                     file = "output/GSE37283/tables/n_uc_tables.xlsx", rowNames = TRUE)

# new
openxlsx::write.xlsx(list("all_genes" = uc_c_toptable[,!colnames(uc_c_toptable) %in% "rank"],
                          "all_glyco" = uc_c_toptable.g[,!colnames(uc_c_toptable.g) %in% "rank"],
                          "full_gsea" = uc_c_gsea,
                          "glyco_gsea" = uc_c_gsea.g),
                     file = "~/projects/tesis/final_outputs/uc_control/GSE37283.xlsx", rowNames = TRUE)

openxlsx::write.xlsx(list("all_genes" = n_c_toptable[,!colnames(n_c_toptable) %in% "rank"],
                          "all_glyco" = n_c_toptable.g[,!colnames(n_c_toptable.g) %in% "rank"],
                          "full_gsea" = n_c_gsea,
                          "glyco_gsea" = n_c_gsea.g),
                     file = "~/projects/tesis/final_outputs/uc_activa_vs_quiescente/neo_vs_control.xlsx", rowNames = TRUE)

openxlsx::write.xlsx(list("all_genes" = n_uc_toptable[,!colnames(n_uc_toptable) %in% "rank"],
                          "all_glyco" = n_uc_toptable.g[,!colnames(n_uc_toptable.g) %in% "rank"],
                          "full_gsea" = n_uc_gsea,
                          "glyco_gsea" = n_uc_gsea.g),
                     file = "~/projects/tesis/final_outputs/uc_activa_vs_quiescente/neo_vs_quc.xlsx", rowNames = TRUE)
