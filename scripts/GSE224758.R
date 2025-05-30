# GSE224758 "Digby-Bell" cohort - bulk RNAseq
# Load packages -----
library(biomaRt)
library(Biobase)
library(clusterProfiler)
library(cowplot)
library(DT)
library(edgeR)
library(enrichplot)
library(ggrepel)
library(gplots)
library(gridExtra)
library(GSEABase)
library(GSVA)
library(gt)
library(limma)
library(matrixStats)
library(msigdbr)
library(plotly)
library(shiny)
library(tidyverse)
library(openxlsx)
library(qpdf)
source("scripts/functions.R")

# Set up data ----
logfc <- 1
pval <- 0.05

load("data/rdata/glycogenes.RData")
load("data/rdata/glycogenes2.RData")

# Import raw counts to R, and capture sample names
expression <- read.csv("data/raw/GSE224758/GSE224758_read_counts.csv")
rownames(expression) <- expression$geneID
expression <- expression[,!colnames(expression) == "geneID"]
expression <- expression[rowMeans(expression) != 0,]
expression <- cpm(expression, log = TRUE)
expression.t <- expression |> as_tibble(rownames = "geneID")
expression.tp <- expression.t |>
  pivot_longer(cols = colnames(expression.t[, -1]),
               names_to = "samples",
               values_to = "expression")
sample_names <- colnames(expression)

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

plot.r <- grid.arrange(eplot.r, density_plot.r, ncol = 2,
                       top = "Raw expression data (unfiltered, non-normalized)")

# Import metadata file to create study_design table
study_design <-
  t(read_tsv("data/raw/GSE224758/GSE224758_series_matrix.txt", skip = 26)) %>%
  as_tibble() %>%
  dplyr::select(1, 10) %>%
  dplyr::slice(-1) %>%
  dplyr::rename(accession = V1, phenotype = V10) %>%
  mutate(phenotype = if_else(phenotype == "condition: healthy", "H", "UC")) %>%
  add_column(sample = sample_names)
phenotypes <- factor(study_design$phenotype)

# Map gene IDs to gene names
listEnsembl()
datasets <- listDatasets(useEnsembl(biomart = "genes", mirror = "www"))
ensembl_connection <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl_connection)
filters <- listFilters(ensembl_connection)

ID_to_genes <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = counts$geneID,
  mart = ensembl_connection)

expression <- merge(expression, ID_to_genes, by.x = "row.names", by.y = "ensembl_gene_id", all.x = FALSE)
expression <- expression[!expression$external_gene_name == "",]
expression <- expression %>% group_by(external_gene_name) %>%
  summarise(across(-c(Row.names), mean)) %>% ungroup()
expression <- data.frame(expression)
rownames(expression) <- expression$external_gene_name
expression <- expression[,colnames(expression) != "external_gene_name"]
table(is.na(rownames(expression)))
geneIDs <- rownames(expression)

# Filter expression data ----
# Divide sample names according to phenotypes
UC_samples <- study_design$sample[study_design$phenotype == "UC"]
H_samples <- study_design$sample[study_design$phenotype == "H"]

# Keep genes with cpm > 2 in all samples of each group
UC_filtered_genes <- rowSums(expression[, UC_samples] > 2) == length(UC_samples)
H_filtered_genes <- rowSums(expression[, H_samples] > 2) == length(H_samples)
UC_filtered_genes <- names(UC_filtered_genes[UC_filtered_genes])
H_filtered_genes <- names(H_filtered_genes[H_filtered_genes])
filtered_genes <- unique(c(UC_filtered_genes, H_filtered_genes))
expression <- expression[filtered_genes,]
geneIDs <- rownames(expression)

# Plot filtered data
expression.t <- expression |> as_tibble(rownames = "geneID")
expression.tp <- expression.t |>
  pivot_longer(cols = colnames(expression.t[, -1]),
               names_to = "samples",
               values_to = "expression")

density_plot.f <- ggplot(data = expression.tp, aes(x = expression, color = samples)) +
  geom_density(adjust = 1.5) +
  theme_classic()

eplot.f <- ggplot(expression.tp) +
  aes(x = samples, y = expression, fill = samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", geom = "point", shape = 95, size = 10, color = "black", show.legend = FALSE) +
  labs(y = "log2 CPM", x = "sample") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_bw()

plot.f <- grid.arrange(eplot.f, density_plot.f, ncol = 2,
                       top = "Filtered expression data (filtered, non-normalized)")

# TMM normalization: adjusts for library sizes ----
counts <- read.csv("data/raw/GSE224758/GSE224758_read_counts.csv")
rownames(counts) <- counts$geneID
counts <- counts[,-1]
counts <- counts[rowMeans(counts) != 0,]
counts <- merge(counts, ID_to_genes, by.x = "row.names", by.y = "ensembl_gene_id", all.x = FALSE)
counts <- counts[!counts$external_gene_name == "",]
counts <- counts %>% group_by(external_gene_name) %>% summarise(across(-c(Row.names), mean)) %>% ungroup()
counts <- data.frame(counts)
rownames(counts) <- counts$external_gene_name
counts <- counts[,colnames(counts) != "external_gene_name"]
counts.f <- counts[filtered_genes,]
DGElist.fn <- calcNormFactors(DGEList(counts.f), method = "TMM")

# Plot normalized data
expression <- cpm(DGElist.fn, log = TRUE)
expression.t <- expression |> as_tibble(rownames = "geneID")
expression.tp <- expression.t |>
  pivot_longer(cols = colnames(expression.t[, -1]),
               names_to = "samples",
               values_to = "expression")

density_plot.n <- ggplot(data = expression.tp, aes(x = expression, color = samples)) +
  geom_density(adjust = 1.5) +
  theme_classic()

eplot.n <- ggplot(expression.tp) +
  aes(x = samples, y = expression, fill = samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", geom = "point", shape = 95, size = 10, color = "black", show.legend = FALSE) +
  labs(y = "log2 CPM", x = "sample") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_bw()

plot.n <- grid.arrange(eplot.n, density_plot.n, ncol = 2,
                       top = "Normalized expression data (filtered, TMM-normalized)")

plot.c <- grid.arrange(plot.r, plot.f, plot.n)

ggsave(plot = plot.c, filename = "quality_combined.pdf", device = "pdf", dpi = 300,
       height = 16, width = 12, path = "output/GSE224758/quality_plots")

# PCA ----
pca_expression <- expression # expression[, colnames(expression) != "GSM1162246"]
pca_samples <- colnames(pca_expression)
pca_phenotypes <- study_design$phenotype[study_design$sample %in% pca_samples]
pca <- prcomp(t(pca_expression), scale. = T, center = T, retx = T)
eigenvalues <- pca$sdev^2
pc_percentages <- round((eigenvalues) / sum(eigenvalues) * 100, 1)
PCs <- paste0("PC", 1:length(pc_percentages))
pc_percentages <- data.frame(Component = PCs, Percentage = pc_percentages)
pca_data <- as_tibble(pca$x)

pca_plot <- ggplot(pca_data) +
  aes(x = PC1, y = PC2, color = pca_phenotypes) +
  geom_point(size = 3) +
  scale_color_manual(values = c("cornflowerblue", "coral1")) +
  xlab(paste0("PC1 (", pc_percentages[1, 2], "%", ")")) +
  ylab(paste0("PC2 (", pc_percentages[2, 2], "%", ")")) +
  # geom_label_repel(aes(label = sample_names), size = 3) +
  theme_bw()

print(pca_plot)
ggsave(path = "output/GSE224758/quality_plots/",
       filename = "pca.pdf",
       pca_plot, device = "pdf", width = 6, height = 5, dpi = 300)

# Linear model ----
lm_dmatrix <- model.matrix(~ 0 + phenotypes)
colnames(lm_dmatrix) <- levels(phenotypes)
rownames(lm_dmatrix) <- sample_names
master_lm <- lmFit(expression, lm_dmatrix)
lm_cmatrix <- makeContrasts(UC_control = UC - H, levels = lm_dmatrix)
lm <- contrasts.fit(master_lm, lm_cmatrix)
lm_stats <- eBayes(lm)

# TopTable and DEGs -----
# Make toptables
toptable <- make_toptable("UC_control", "UC - control")
#toptable.f <- toptable %>% dplyr::filter(adj.P.Val <= pval)
toptable.g <- toptable %>% dplyr::filter(row.names(toptable) %in% glycogenes2)
#toptable.fg <- toptable.f %>% dplyr::filter(row.names(toptable.f) %in% glycogenes)

# Make and save volcano plots----
volcano(toptable, p = pval, lfc = logfc, title = "DEGs in UC vs. control",
        path = "output/GSE224758/volcano_plots/separate/", filename = "uc_c_volcano.pdf")
volcano(toptable.g, p = pval, lfc = logfc, title = "GlycoDEGs in UC vs. control",
        path = "output/GSE224758/volcano_plots/separate/", filename = "uc_c_volcano_glyco.pdf")

# Save combined .pdf according to comparison in title
pdfs <- list.files(path = "output/GSE224758/volcano_plots/separate/", full.names = TRUE)
qpdf::pdf_combine(input = pdfs[grepl("uc_c_", pdfs)],
                  output = "output/GSE224758/volcano_plots/uc_c_combined.pdf")

# Perform GSEA (gene set enrichment analysis) ----
load("data/rdata/c5_hs.RData")
load("data/rdata/reactome_hs.RData")

# Construct an ordered named vector to perform GSEA
ranked_DEGs <- toptable$rank
names(ranked_DEGs) <- rownames(toptable)

# Run GSEA
load("data/rdata/glyco_terms.RData")
uc_c_gsea <- run_gsea(ranked_DEGs, "UC vs Control", c5_hs, reactome_hs)
uc_c_gsea.g <- uc_c_gsea %>% dplyr::filter(grepl(paste(glyco_terms, collapse = "|"), as.character(ID), ignore.case = TRUE))

# Make and save bubble plots
gsea_bubble(uc_c_gsea, comp = "uc_c", path = "output/GSE224758/gsea_plots/separate/", n = 30)
gsea_bubble(uc_c_gsea.g, comp = "uc_c_glyco", path = "output/GSE224758/gsea_plots/separate/", n = 30)

# Save combined .pdf according to comparison in title
pdfs <- list.files(path = "output/GSE224758/gsea_plots/separate/", full.names = TRUE)
qpdf::pdf_combine(input = pdfs[grepl("uc_c_", pdfs)],
                  output = "output/GSE224758/gsea_plots/uc_c_combined.pdf")

# GSVA: pathway enrichment per sample: patient sub-populations ----
# Format genesets for GSVA
GO_BP <- c5_hs |> dplyr::filter(grepl("GOBP", gs_name))
GO_BP <- split(x = GO_BP$gene_symbol, f = GO_BP$gs_name)

GO_MF <- c5_hs |> dplyr::filter(grepl("GOMF", gs_name))
GO_MF <- split(x = GO_MF$gene_symbol, f = GO_MF$gs_name)

GO_CC <- c5_hs |>  dplyr::filter(grepl("GOCC", gs_name))
GO_CC <- split(x = GO_CC$gene_symbol, f = GO_CC$gs_name)

GO_RE <- reactome_hs |>  dplyr::filter(grepl("REACTOME", gs_name))
GO_RE <- split(x = GO_RE$gene_symbol, f = GO_RE$gs_name)

# Run GSVA
gsva_BP <- run_gsva(exp = expression, go = GO_BP, dmatrix = lm_dmatrix, cmatrix = lm_cmatrix, lfc = 0.1, p = 0.00001)  
gsva_MF <- run_gsva(exp = expression, go = GO_MF, dmatrix = lm_dmatrix, cmatrix = lm_cmatrix, lfc = 0.1, p = 0.00001)  
gsva_CC <- run_gsva(exp = expression, go = GO_CC, dmatrix = lm_dmatrix, cmatrix = lm_cmatrix, lfc = 0.1, p = 0.00001)  
gsva_RE <- run_gsva(exp = expression, go = GO_RE, dmatrix = lm_dmatrix, cmatrix = lm_cmatrix, lfc = 0.1, p = 0.00001)  
summary(gsva_BP$decidetests)
summary(gsva_MF$decidetests)
summary(gsva_CC$decidetests)
summary(gsva_RE$decidetests)

# Pull out enriched GSVA results --> (gsva, decidetests, coeff, samples)
uc_c_GSVA_BP <- pull_gsva(gsva = gsva_BP$res, decidetests = gsva_BP$decidetests,
                          coeff = 1, samples = c(UC_samples, H_samples))
uc_c_GSVA_MF <- pull_gsva(gsva = gsva_MF$res, decidetests = gsva_MF$decidetests,
                          coeff = 1, samples = c(UC_samples, H_samples))
uc_c_GSVA_CC <- pull_gsva(gsva = gsva_CC$res, decidetests = gsva_CC$decidetests,
                          coeff = 1, samples = c(UC_samples, H_samples))
uc_c_GSVA_RE <- pull_gsva(gsva = gsva_RE$res, decidetests = gsva_RE$decidetests,
                          coeff = 1, samples = c(UC_samples, H_samples))
uc_c_GSVA_TOTAL <- rbind(uc_c_GSVA_BP, uc_c_GSVA_MF, uc_c_GSVA_CC, uc_c_GSVA_RE)
uc_c_GSVA_TOTAL.g <- uc_c_GSVA_TOTAL[grepl(paste(glyco_terms, collapse = "|"), rownames(uc_c_GSVA_TOTAL), ignore.case = TRUE), ]

# GSVA heatmap ----
# Create color palettes
heatmap_colors <- paletteer_c("pals::coolwarm", 100, direction = 1)
phenotype_colors <- paletteer_d("lisa::JosefAlbers", 2)
uc_c_colors <- setNames(phenotype_colors, c("control", "H"))

# Make heatmaps
#make_heatmap <- function(data, phencols, path, filename) {
make_heatmap(data = uc_c_GSVA_TOTAL, phencols = uc_c_colors,
             path = "output/GSE224758/gsva_plots/separate/", filename = "uc_c_gsva.pdf")
make_heatmap(data = uc_c_GSVA_TOTAL.g[,2:ncol(uc_c_GSVA_TOTAL.g)], phencols = uc_c_colors,
             path = "output/GSE224758/gsva_plots/separate/", filename = "uc_c_glyco_gsva.pdf")

# Save Excel sheets ----
#old
openxlsx::write.xlsx(list("all_genes" = toptable,
                          "all_significant" = toptable.f,
                          "all_glyco" = toptable.g,
                          "significant_glyco" = toptable.fg,
                          "full_gsea" = uc_c_gsea,
                          "glyco_gsea" = uc_c_gsea.g),
                     file = "output/GSE224758/tables/uc_c_tables.xlsx", rowNames = TRUE)
#new
openxlsx::write.xlsx(list("all_genes" = toptable[,!colnames(toptable) %in% "rank"],
                          "all_glyco" = toptable.g[,!colnames(toptable.g) %in% "rank"],
                          "gsea" = uc_c_gsea,
                          "glyco_gsea" = uc_c_gsea.g),
                     file = "~/projects/tesis/final_outputs/uc_control/GSE224758.xlsx", rowNames = TRUE)
