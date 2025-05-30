# GSE3629 "Watanabe" cohort - bulk microarray
# Load packages -----
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
library(hgu133plus2.db)
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
load("data/rdata/glycogenes2.RData")

# Set up data ----
load("data/rdata/glycogenes.RData")
load("data/rdata/glyco_terms.RData")
logfc <- 2
pval <- 0.05
read <- read_tsv("data/raw/GSE3629/GSE3629_series_matrix.txt", skip = 30, col_names = TRUE)
read <- as.data.frame(read)
rownames(read) <- read[,1]
read[,1] <- NULL
metadata <- read[1:30,]
expression <- read[32:(nrow(read)-1),]
colnames(expression) <- expression[1,]
expression <- expression[-1,]
study_design <- metadata[c("!Sample_description", "!Sample_geo_accession"),]
study_design <- t(study_design)
rownames(study_design) <- NULL
colnames(study_design) <- c("phenotype", "sampleID")
study_design <- as_tibble(study_design)
unique(study_design$phenotype)
sporadic_ca_samples <- study_design$sampleID[study_design$phenotype == "Sporadic Ca"]
uc_associated_ca_samples <- study_design$sampleID[study_design$phenotype == "UC-associated Ca"]
uc_ca_samples <- study_design$sampleID[study_design$phenotype == "UC-Ca"]
uc_non_ca_samples <- study_design$sampleID[study_design$phenotype == "UC-NonCa"]
study_design <- study_design[!study_design$sampleID %in% uc_associated_ca_samples,]
phenotypes <- factor(study_design$phenotype)
levels(phenotypes)
sample_names <- study_design$sampleID

expression <- mutate_all(expression, function(x) as.numeric(as.character(x)))

expression$ID_REF <- mapIds(hgu133plus2.db,
                            keys = row.names(expression),
                            column = "SYMBOL",
                            keytype = "PROBEID",
                            multiVals = "first")
expression <- expression[!is.na(expression$ID_REF),]
anyDuplicated(expression$ID_REF)
expression <- expression %>%
  group_by(ID_REF) %>%
  summarise(across(everything(), mean))
expression <- as.data.frame(expression)
rownames(expression) <- expression$ID_REF
expression$ID_REF <- NULL
expression <- expression[,!colnames(expression) %in% uc_associated_ca_samples]

# filter without eliminating 0
# Keep genes with cpm > 2 in all samples of each group
sporadic_ca_genes <- rowSums(expression[, sporadic_ca_samples] > 1) == length(sporadic_ca_samples)
uc_ca_genes <- rowSums(expression[, uc_ca_samples] > 1) == length(uc_ca_samples)
uc_non_ca_genes <- rowSums(expression[, uc_non_ca_samples] > 1) == length(uc_non_ca_samples)

sporadic_ca_genes <- names(sporadic_ca_genes[sporadic_ca_genes])
uc_ca_genes <- names(uc_ca_genes[uc_ca_genes])
uc_non_ca_genes <- names(uc_non_ca_genes[uc_non_ca_genes])

filtered_genes <- unique(c(sporadic_ca_genes, uc_ca_genes, uc_non_ca_genes))
expression <- expression[filtered_genes,]
geneIDs <- rownames(expression)
# Replace 0 values with a small value in order to log2 transform without infinite
expression[expression == 0] <- 0.00001
log2(expression[is.infinite(rowSums(log2(expression))), ])
expression <- log2(expression)
expression <- na.omit(expression)
expression <- expression[!is.infinite(rowSums(expression)),]

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

ggsave(plot = density_plot.r, filename = "density_raw.pdf", device = "pdf", dpi = 300,
       height = 7, width = 14, path = "output/GSE3629/quality_plots")

ggsave(plot = eplot.r, filename = "quality_raw.pdf", device = "pdf", dpi = 300,
       height = 10, width = 30, path = "output/GSE3629/quality_plots")

# Quantile normalization
expression <- normalizeBetweenArrays(expression)
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

ggsave(plot = density_plot.n, filename = "density_normalized.pdf", device = "pdf", dpi = 300,
       height = 7, width = 14, path = "output/GSE3629/quality_plots")

ggsave(plot = eplot.n, filename = "quality_normalized.pdf", device = "pdf", dpi = 300,
       height = 10, width = 49, path = "output/GSE3629/quality_plots")

# PCA test ----
pca_expression <- expression
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
  scale_color_manual(values = c("cornflowerblue", "orchid", "aquamarine")) +
  xlab(paste0("PC1 (", pc_percentages[1, 2], "%", ")")) +
  ylab(paste0("PC2 (", pc_percentages[2, 2], "%", ")")) +
  #geom_label_repel(aes(label = sample_names), size = 3) +
  theme_bw()

print(pca_plot)
ggsave(path = "output/GSE3629/quality_plots/",
       filename = "pca_after.pdf",
       pca_plot, device = "pdf", width = 6, height = 5, dpi = 300)

# Linear model ----
phenotypes <- factor(study_design$phenotype)

lm_dmatrix <- model.matrix(~ 0 + phenotypes)
colnames(lm_dmatrix) <- levels(phenotypes)
rownames(lm_dmatrix) <- sample_names
master_lm <- lmFit(expression, lm_dmatrix)
lm_cmatrix <- makeContrasts(CAC_UC = CAC - UC,
                            CAC_sporadic = CAC - sporadic_ca,
                            sporadic_UC = sporadic_ca - UC,
                            levels = lm_dmatrix)
lm <- contrasts.fit(master_lm, lm_cmatrix)
lm_stats <- eBayes(lm)

# TopTable to view DEGs and filter glycogenes -----
cac_uc_toptable <- make_toptable("CAC_UC", "CAC - UC")
cac_sporadic_toptable <- make_toptable("CAC_sporadic", "CAC - sporadic Ca")
sporadic_uc_toptable <- make_toptable("sporadic_UC", "sporadic Ca - UC")
# 
# cac_uc_toptable.f <- cac_uc_toptable %>% dplyr::filter(adj.P.Val <= pval)
# cac_sporadic_toptable.f <- cac_sporadic_tobtable %>% dplyr::filter(adj.P.Val <= pval)
# sporadic_uc_toptable.f <- sporadic_uc_toptable %>% dplyr::filter(adj.P.Val <= pval)

cac_uc_toptable.g <- cac_uc_toptable %>% dplyr::filter(row.names(cac_uc_toptable) %in% glycogenes2)
cac_sporadic_toptable.g <- cac_sporadic_tobtable %>% dplyr::filter(row.names(cac_sporadic_tobtable) %in% glycogenes2)
sporadic_uc_toptable.g <- sporadic_uc_toptable %>% dplyr::filter(row.names(sporadic_uc_toptable) %in% glycogenes2)

# cac_uc_toptable.fg <- cac_uc_toptable.f %>% dplyr::filter(row.names(cac_uc_toptable.f) %in% glycogenes2)
# cac_sporadic_toptable.fg <- cac_sporadic_tobtable.f %>% dplyr::filter(row.names(cac_sporadic_tobtable.f) %in% glycogenes2)
# sporadic_uc_toptable.fg <- sporadic_uc_toptable.f %>% dplyr::filter(row.names(sporadic_uc_toptable.f) %in% glycogenes2)

# Volcano Plot: logFC/p-value ----
logfc <- 2
pval <- 0.05
data <- as_tibble(cac_uc_toptable.g, rownames = "geneID")
degs <- (data[(abs(data$logFC) > logfc) & (data$adj.P.Val < pval),])$geneID
v <- ggplot(data) +
  geom_point(size = 3, color = ifelse(!data$geneID %in% degs, "grey", ifelse(data$logFC > 0, "red3", "blue3"))) +
  aes(y = -log10(adj.P.Val), x = logFC, text = paste("Symbol:", geneID), customdata = geneID) +
  geom_hline(yintercept = -log10(pval), linetype = "longdash", colour = "grey", linewidth = 1) +
  geom_vline(xintercept = logfc, linetype = "longdash", colour = "#BE684D", linewidth = 1) +
  geom_vline(xintercept = -logfc, linetype = "longdash", colour = "#2C467A", linewidth = 1) +
  geom_text_repel(data = subset(data, geneID %in% degs), aes(label = geneID), vjust = -1, size = 5) +
  theme_bw(base_size = 15)

plot(v)
volcano(cac_uc_toptable, p = pval, lfc = logfc, title = "DEGs in CAC vs UC",
        path = "output/GSE3629/volcano_plots/separate/", filename = "cac_uc_volcano.pdf")
volcano(cac_uc_toptable.g, p = pval, lfc = logfc, title = "GlycoDEGs in CAC vs UC",
        path = "output/GSE3629/volcano_plots/separate/", filename = "cac_uc_volcano_glyco.pdf")
volcano(cac_sporadic_toptable, p = pval, lfc = logfc, title = "DEGs in CAC vs. sporadic Ca",
        path = "output/GSE3629/volcano_plots/separate/", filename = "cac_sporadic_volcano.pdf")
volcano(cac_sporadic_toptable.g, p = pval, lfc = logfc, title = "GlycoDEGs in CAC vs. sporadic Ca",
        path = "output/GSE3629/volcano_plots/separate/", filename = "cac_sporadic_volcano_glyco.pdf")
volcano(sporadic_uc_toptable, p = pval, lfc = logfc, title = "DEGs in sporadic Ca vs. UC",
        path = "output/GSE3629/volcano_plots/separate/", filename = "sporadic_uc_volcano.pdf")
volcano(sporadic_uc_toptable.g, p = pval, lfc = logfc, title = "GlycoDEGs in sporadic Ca vs. UC",
        path = "output/GSE3629/volcano_plots/separate/", filename = "sporadic_uc_volcano_glyco.pdf")

pdfs <- list.files(path = "output/GSE3629/volcano_plots/separate/", full.names = TRUE)
qpdf::pdf_combine(input = pdfs[grepl("cac_uc_", pdfs)],
                  output = "output/GSE3629/volcano_plots/cac_uc_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("cac_sporadic_", pdfs)],
                  output = "output/GSE3629/volcano_plots/cac_sporadic_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("sporadic_uc_", pdfs)],
                  output = "output/GSE3629/volcano_plots/sporadic_uc_combined.pdf")

# Perform GSEA (gene set enrichment analysis) ----
load("data/rdata/c5_hs.RData")
load("data/rdata/reactome_hs.RData")
load("data/rdata/glyco_terms.RData")

ranked_cac_uc_DEGs <- cac_uc_toptable$rank
ranked_cac_sporadic_DEGs <- cac_sporadic_toptable$rank
ranked_sporadic_uc_DEGs <- sporadic_uc_toptable$rank
names(ranked_cac_uc_DEGs) <- rownames(cac_uc_toptable)
names(ranked_cac_sporadic_DEGs) <- rownames(cac_sporadic_toptable)
names(ranked_sporadic_uc_DEGs) <- rownames(sporadic_uc_toptable)

cac_uc_gsea <- run_gsea(ranked_cac_uc_DEGs, "CAC vs. UC", c5_hs, reactome_hs)
cac_uc_gsea.g <- cac_uc_gsea %>% dplyr::filter(grepl(paste(glyco_terms, collapse = "|"), as.character(ID), ignore.case = TRUE))
cac_sporadic_gsea <- run_gsea(ranked_cac_sporadic_DEGs, "CAC vs. sporadic Cancer", c5_hs, reactome_hs)
cac_sporadic_gsea.g <- cac_sporadic_gsea %>% dplyr::filter(grepl(paste(glyco_terms, collapse = "|"), as.character(ID), ignore.case = TRUE))
sporadic_uc_gsea <- run_gsea(ranked_sporadic_uc_DEGs, "sporadic vs. UC", c5_hs, reactome_hs)
sporadic_uc_gsea.g <- sporadic_uc_gsea %>% dplyr::filter(grepl(paste(glyco_terms, collapse = "|"), as.character(ID), ignore.case = TRUE))

gsea_bubble(cac_uc_gsea, comp = "cac_uc", path = "output/GSE3629/gsea_plots/separate/", n = 30)
gsea_bubble(cac_uc_gsea.g, comp = "cac_uc_glyco", path = "output/GSE3629/gsea_plots/separate/", n = 30)
gsea_bubble(cac_sporadic_gsea, comp = "cac_sporadic", path = "output/GSE3629/gsea_plots/separate/", n = 30)
gsea_bubble(cac_sporadic_gsea.g, comp = "cac_sporadic_glyco", path = "output/GSE3629/gsea_plots/separate/", n = 30)
gsea_bubble(sporadic_uc_gsea, comp = "sporadic_uc", path = "output/GSE3629/gsea_plots/separate/", n = 30)
gsea_bubble(sporadic_uc_gsea.g, comp = "sporadic_uc_glyco", path = "output/GSE3629/gsea_plots/separate/", n = 30)

pdfs <- list.files(path = "output/GSE3629/gsea_plots/separate/", full.names = TRUE)
qpdf::pdf_combine(input = pdfs[grepl("cac_uc_", pdfs)],
                  output = "output/GSE3629/gsea_plots/cac_uc_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("cac_sporadic_", pdfs)],
                  output = "output/GSE3629/gsea_plots/cac_sporadic_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("sporadic_uc_", pdfs)],
                  output = "output/GSE3629/gsea_plots/sporadic_uc_combined.pdf")

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
cac_uc_GSVA_BP <- pull_gsva(gsva = gsva_BP$res, decidetests = gsva_BP$decidetests,
                          coeff = "CAC_UC", samples = c(uc_ca_samples, uc_non_ca_samples))
cac_uc_GSVA_MF <- pull_gsva(gsva = gsva_MF$res, decidetests = gsva_MF$decidetests,
                          coeff = "CAC_UC", samples = c(uc_ca_samples, uc_non_ca_samples))
cac_uc_GSVA_CC <- pull_gsva(gsva = gsva_CC$res, decidetests = gsva_CC$decidetests,
                          coeff = "CAC_UC", samples = c(uc_ca_samples, uc_non_ca_samples))
cac_uc_GSVA_RE <- pull_gsva(gsva = gsva_RE$res, decidetests = gsva_RE$decidetests,
                          coeff = "CAC_UC", samples = c(uc_ca_samples, uc_non_ca_samples))
cac_uc_GSVA_TOTAL <- rbind(cac_uc_GSVA_BP, cac_uc_GSVA_MF, cac_uc_GSVA_CC, cac_uc_GSVA_RE)
cac_uc_GSVA_TOTAL.g <- cac_uc_GSVA_TOTAL[grepl(paste(glyco_terms, collapse = "|"), rownames(cac_uc_GSVA_TOTAL), ignore.case = TRUE), ]
cacrc_pan_GSVA_TOTAL <- cac_uc_GSVA_TOTAL

save(cacrc_pan_GSVA_TOTAL, file = "~/projects/tesis/output/integration/cacrc_pan_GSVA_TOTAL.RData")

cac_sporadic_GSVA_BP <- pull_gsva(gsva = gsva_BP$res, decidetests = gsva_BP$decidetests,
                         coeff = "CAC_sporadic", samples = c(uc_ca_samples, sporadic_ca_samples))
cac_sporadic_GSVA_MF <- pull_gsva(gsva = gsva_MF$res, decidetests = gsva_MF$decidetests,
                         coeff = "CAC_sporadic", samples = c(uc_ca_samples, sporadic_ca_samples))
cac_sporadic_GSVA_CC <- pull_gsva(gsva = gsva_CC$res, decidetests = gsva_CC$decidetests,
                         coeff = "CAC_sporadic", samples = c(uc_ca_samples, sporadic_ca_samples))
cac_sporadic_GSVA_RE <- pull_gsva(gsva = gsva_RE$res, decidetests = gsva_RE$decidetests,
                         coeff = "CAC_sporadic", samples = c(uc_ca_samples, sporadic_ca_samples))
cac_sporadic_GSVA_TOTAL <- rbind(cac_sporadic_GSVA_BP, cac_sporadic_GSVA_MF, cac_sporadic_GSVA_CC, cac_sporadic_GSVA_RE)
cac_sporadic_GSVA_TOTAL.g <- cac_sporadic_GSVA_TOTAL[grepl(paste(glyco_terms, collapse = "|"), rownames(cac_sporadic_GSVA_TOTAL), ignore.case = TRUE), ]

# GSVA heatmap ----
heatmap_colors <- paletteer_c("pals::coolwarm", 100, direction = 1)
phenotype_colors <- paletteer_d("lisa::JosefAlbers", 2)
cac_uc_colors <- setNames(phenotype_colors, c("CAC", "UC"))
cac_sporadic_colors <- setNames(phenotype_colors, c("CAC", "sporadic_ca"))

make_heatmap(data = cac_uc_GSVA_TOTAL, phencols = cac_uc_colors,
             path = "output/GSE3629/gsva_plots/separate/", filename = "cac_uc_gsva.pdf")
make_heatmap(data = cac_uc_GSVA_TOTAL.g[,2:ncol(cac_uc_GSVA_TOTAL.g)], phencols = cac_uc_colors,
             path = "output/GSE3629/gsva_plots/separate/", filename = "cac_uc_glyco_gsva.pdf")
make_heatmap(data = cac_sporadic_GSVA_TOTAL, phencols = cac_sporadic_colors,
             path = "output/GSE3629/gsva_plots/separate/", filename = "cac_sporadic_gsva.pdf")
make_heatmap(data = cac_sporadic_GSVA_TOTAL.g[,2:ncol(cac_sporadic_GSVA_TOTAL.g)], phencols = cac_sporadic_colors,
             path = "output/GSE3629/gsva_plots/separate/", filename = "cac_sporadic_glyco_gsva.pdf")

pdfs <- list.files(path = "output/GSE3629/gsva_plots/separate/", full.names = TRUE)
qpdf::pdf_combine(input = pdfs[grepl("cac_uc_", pdfs)],
                  output = "output/GSE3629/gsva_plots/cac_uc_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("cac_sporadic_", pdfs)],
                  output = "output/GSE3629/gsva_plots/cac_sporadic_combined.pdf")

# Save Excel sheets ----
openxlsx::write.xlsx(list("all_genes" = cac_uc_toptable,
                          "all_significant" = cac_uc_toptable.f,
                          "all_glyco" = cac_uc_toptable.g,
                          "significant_glyco" = cac_uc_toptable.fg,
                          "full_gsea" = cac_uc_gsea,
                          "glyco_gsea" = cac_uc_gsea.g),
                     file = "output/GSE3629/tables/cac_uc_tables.xlsx", rowNames = TRUE)

openxlsx::write.xlsx(list("all_genes" = cac_sporadic_toptable,
                          "all_significant" = cac_sporadic_toptable.f,
                          "all_glyco" = cac_sporadic_toptable.g,
                          "significant_glyco" = cac_sporadic_toptable.fg,
                          "full_gsea" = cac_sporadic_gsea,
                          "glyco_gsea" = cac_sporadic_gsea.g),
                     file = "output/GSE3629/tables/cac_sporadic_tables.xlsx", rowNames = TRUE)

openxlsx::write.xlsx(list("all_genes" = sporadic_uc_toptable,
                          "all_significant" = sporadic_uc_toptable.f,
                          "all_glyco" = sporadic_uc_toptable.g,
                          "significant_glyco" = sporadic_uc_toptable.fg,
                          "full_gsea" = sporadic_uc_gsea,
                          "glyco_gsea" = sporadic_uc_gsea.g),
                     file = "output/GSE3629/tables/sporadic_uc_tables.xlsx", rowNames = TRUE)

save(cac_uc_toptable, file = "~/projects/tesis/data/pancolitis/cacrc_uc_toptable.RData")
