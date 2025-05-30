# GSE31106 "Tang" mouse cohort - bulk microarray
# Pooleé displasia
# Choose macros ----
pval <- 0.05
lfc <- 2

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
library(mouse4302.db)
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
load("data/rdata/glycogenes_mouse.RData")
load("data/rdata/glycogenes_mouse2.RData")
load("data/rdata/HtoM.RData")
read <- read_tsv("data/raw/GSE31106/GSE31106_series_matrix.txt",
                     skip = 28, col_names = FALSE) %>% as_tibble()

metadata <- read[1:32, ]
metadata <- t(metadata)
colnames(metadata) <- metadata[1, ]
metadata <- metadata[-1, ]

# Read expression data
expression <- read[35:nrow(read),]
expression <- expression[-nrow(expression),]
colnames(expression) <- expression[1, ]
expression <- expression[-1,]
expression <- expression %>% column_to_rownames(var = "ID_REF")
geneIDs <- mapIds(mouse4302.db,
                  keys = rownames(expression),
                  column = "SYMBOL",
                  keytype = "PROBEID",
                  multiVals = "first")
geneIDs <- geneIDs[!is.na(geneIDs)]
anyDuplicated(geneIDs) # check repeated rownames 
expression$geneID <- geneIDs[rownames(expression)]
expression <- expression[!is.na(expression$geneID),]

# Collapse rows into the same genename, taking average of expression
expression <- as.matrix(expression)
rownames(expression) <- expression[,"geneID"]
expression <- avereps(expression) # average values per rowname
anyDuplicated(rownames(expression)) # check repeated rownames
expression <- expression[,-19]
geneIDs <- rownames(expression) # Capture gene names
expression <- apply(expression, 2, as.numeric)
rownames(expression) <- geneIDs
expression <- log2(expression)

# Make study_design table
study_design <- as_tibble(metadata[,c(2,13)])
colnames(study_design) <- c("sampleID", "phenotype")
study_design <- study_design %>%
  mutate(phenotype = case_when(grepl("normal", phenotype) ~ "control",
                               grepl("inflamed", phenotype) ~ "inflamed",
                               grepl("low", phenotype) ~ "dysplasia",
                               grepl("high", phenotype) ~ "dysplasia",
                               grepl("adenocarcinoma", phenotype) ~ "adenocarcinoma"))

# Capture sample IDs and phenotypes
sample_names <- study_design$sampleID
phenotypes <- factor(study_design$phenotype)
adenocarcinoma_samples <- study_design$sampleID[study_design$phenotype == "adenocarcinoma"]
dysplasia_samples <- study_design$sampleID[study_design$phenotype == "dysplasia"]
control_samples <- study_design$sampleID[study_design$phenotype == "control"]
inflammation_samples <- study_design$sampleID[study_design$phenotype == "inflamed"]
phenotypes <- factor(study_design$phenotype)

# Plot data quality
expression.t <- expression |> as_tibble(rownames = "geneID")
expression.tp <- expression.t |>
  pivot_longer(cols = colnames(expression.t[, -1]),
               names_to = "samples",
               values_to = "expression")

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
                       top = "Expression data analyzed with MAS 5.0 using Affymetrix default settings and global scaling (target intensity = 500).")
ggsave(plot = plot.r, filename = "quality.pdf", device = "pdf", dpi = 300,
       height = 7, width = 12, path = "output/GSE31106/quality_plots")

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
  scale_color_manual(values = c("cornflowerblue", "coral1", "orchid", "aquamarine")) +
  xlab(paste0("PC1 (", pc_percentages[1, 2], "%", ")")) +
  ylab(paste0("PC2 (", pc_percentages[2, 2], "%", ")")) +
  #geom_label_repel(aes(label = sample_names), size = 3) +
  theme_bw()

print(pca_plot)
ggsave(path = "output/GSE31106/quality_plots/",
       filename = "pca.pdf",
       pca_plot, device = "pdf", width = 6, height = 5, dpi = 300)

# Linear model ----
lm_dmatrix <- model.matrix(~ 0 + phenotypes)
colnames(lm_dmatrix) <- levels(phenotypes)
rownames(lm_dmatrix) <- sample_names
master_lm <- lmFit(expression, lm_dmatrix)
unique(phenotypes)
lm_cmatrix <- makeContrasts(inflamed_control = inflamed - control,
                            dysplasia_control = dysplasia - control,
                            dysplasia_inflamed = dysplasia - inflamed,
                            adenocarcinoma_control = adenocarcinoma - control,
                            adenocarcinoma_inflamed = adenocarcinoma - inflamed,
                            adenocarcinoma_dysplasia = adenocarcinoma - dysplasia,
                            levels = lm_dmatrix)
lm <- contrasts.fit(master_lm, lm_cmatrix)
lm_stats <- eBayes(lm)

# TopTable to view DEGs and filter glycogenes -----
i_c_toptable <- make_toptable("inflamed_control", "inflamed - control")
d_c_toptable <- make_toptable("dysplasia_control", "dysplasia - control")
d_i_toptable <- make_toptable("dysplasia_inflamed", "dysplasia - inflamed")
a_c_toptable <- make_toptable("adenocarcinoma_control", "adenocarcinoma - control")
a_i_toptable <- make_toptable("adenocarcinoma_inflamed", "adenocarcinoma - inflamed")
a_d_toptable <- make_toptable("adenocarcinoma_dysplasia", "adenocarcinoma - dysplasia")

i_c_toptable.f <- i_c_toptable %>% dplyr::filter(adj.P.Val <= pval)
d_c_toptable.f <- d_c_toptable %>% dplyr::filter(adj.P.Val <= pval)
d_i_toptable.f <- d_i_toptable %>% dplyr::filter(adj.P.Val <= pval)
a_c_toptable.f <- a_c_toptable %>% dplyr::filter(adj.P.Val <= pval)
a_i_toptable.f <- a_i_toptable %>% dplyr::filter(adj.P.Val <= pval)
a_d_toptable.f <- a_d_toptable %>% dplyr::filter(adj.P.Val <= pval)

i_c_toptable.g <- i_c_toptable %>% dplyr::filter(row.names(i_c_toptable) %in% glycogenes_mouse)
d_c_toptable.g <- d_c_toptable %>% dplyr::filter(row.names(d_c_toptable) %in% glycogenes_mouse)
d_i_toptable.g <- d_i_toptable %>% dplyr::filter(row.names(d_i_toptable) %in% glycogenes_mouse)
a_c_toptable.g <- a_c_toptable %>% dplyr::filter(row.names(a_c_toptable) %in% glycogenes_mouse)
a_i_toptable.g <- a_i_toptable %>% dplyr::filter(row.names(a_i_toptable) %in% glycogenes_mouse)
a_d_toptable.g <- a_d_toptable %>% dplyr::filter(row.names(a_d_toptable) %in% glycogenes_mouse)

i_c_toptable.fg <- i_c_toptable.f %>% dplyr::filter(row.names(i_c_toptable.f) %in% glycogenes_mouse)
d_c_toptable.fg <- d_c_toptable.f %>% dplyr::filter(row.names(d_c_toptable.f) %in% glycogenes_mouse)
d_i_toptable.fg <- d_i_toptable.f %>% dplyr::filter(row.names(d_i_toptable.f) %in% glycogenes_mouse)
a_c_toptable.fg <- a_c_toptable.f %>% dplyr::filter(row.names(a_c_toptable.f) %in% glycogenes_mouse)
a_i_toptable.fg <- a_i_toptable.f %>% dplyr::filter(row.names(a_i_toptable.f) %in% glycogenes_mouse)
a_d_toptable.fg <- a_d_toptable.f %>% dplyr::filter(row.names(a_d_toptable.f) %in% glycogenes_mouse)

# Volcano Plot: logFC/p-value ----
volcano(i_c_toptable, p = pval, lfc = logfc, title = "DEGs in inflamed vs. control",
        path = "output/GSE31106/volcano_plots/separate/", filename = "i_c_volcano.pdf")
volcano(i_c_toptable.g, p = pval, lfc = logfc, title = "GlycoDEGs in inflamed vs. control",
        path = "output/GSE31106/volcano_plots/separate/", filename = "i_c_volcano_glyco.pdf")
volcano(d_c_toptable, p = pval, lfc = logfc, title = "DEGs in dysplasia vs. control",
        path = "output/GSE31106/volcano_plots/separate/", filename = "d_c_volcano.pdf")
volcano(d_c_toptable.g, p = pval, lfc = logfc, title = "GlycoDEGs in dysplasia vs. control",
        path = "output/GSE31106/volcano_plots/separate/", filename = "d_c_volcano_glyco.pdf")
volcano(d_i_toptable, p = pval, lfc = logfc, title = "DEGs in dysplasia vs. inflamed",
        path = "output/GSE31106/volcano_plots/separate/", filename = "d_i_volcano.pdf")
volcano(d_i_toptable.g, p = pval, lfc = logfc, title = "GlycoDEGs in dysplasia vs. inflamed",
        path = "output/GSE31106/volcano_plots/separate/", filename = "d_i_volcano_glyco.pdf")
volcano(a_c_toptable, p = pval, lfc = logfc, title = "DEGs in adenocarcinoma vs. control",
        path = "output/GSE31106/volcano_plots/separate/", filename = "a_c_volcano.pdf")
volcano(a_c_toptable.g, p = pval, lfc = logfc, title = "GlycoDEGs in adenocarcinoma vs. control",
        path = "output/GSE31106/volcano_plots/separate/", filename = "a_c_volcano_glyco.pdf")
volcano(a_i_toptable, p = pval, lfc = logfc, title = "DEGs in adenocarcinoma vs. inflamed",
        path = "output/GSE31106/volcano_plots/separate/", filename = "a_i_volcano.pdf")
volcano(a_i_toptable.g, p = pval, lfc = logfc, title = "GlycoDEGs in adenocarcinoma vs. inflamed",
        path = "output/GSE31106/volcano_plots/separate/", filename = "a_i_volcano_glyco.pdf")
volcano(a_d_toptable, p = pval, lfc = logfc, title = "DEGs in adenocarcinoma vs. dysplasia",
        path = "output/GSE31106/volcano_plots/separate/", filename = "a_d_volcano.pdf")
volcano(a_d_toptable.g, p = pval, lfc = logfc, title = "GlycoDEGs in adenocarcinoma vs. dysplasia",
        path = "output/GSE31106/volcano_plots/separate/", filename = "a_d_volcano_glyco.pdf")


# not significant
volcano(a_i_toptable, p = 10, lfc = 0, title = "DEGs in adenocarcinoma vs. inflamed",
        path = "output/GSE31106/volcano_plots/separate/", filename = "a_i_volcano_notsig.pdf")
# save
pdfs <- list.files(path = "output/GSE31106/volcano_plots/separate/", full.names = TRUE)
qpdf::pdf_combine(input = pdfs[grepl("i_c_", pdfs)],
                  output = "output/GSE31106/volcano_plots/i_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("d_c_", pdfs)],
                  output = "output/GSE31106/volcano_plots/d_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("d_i_", pdfs)],
                  output = "output/GSE31106/volcano_plots/d_i_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("a_c_", pdfs)],
                  output = "output/GSE31106/volcano_plots/a_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("a_i_", pdfs)],
                  output = "output/GSE31106/volcano_plots/a_i_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("a_d_", pdfs)],
                  output = "output/GSE31106/volcano_plots/a_d_combined.pdf")

# Access up-to-date collections/signatures with human gene IDs----
msigdbr_species()
msig_mm_collections <- msigdbr(species = "Mus musculus") |>
  dplyr::distinct(gs_cat, gs_subcat) |>
  dplyr::arrange(gs_cat, gs_subcat)

c5_mm <- msigdbr(species = "Mus musculus", category = "C5") |>
  dplyr::select(gs_name, gene_symbol)
save(c5_mm, file = "data/rdata/c5_mm.RData")

reactome_mm <- msigdbr(species = "Mus musculus", category = "C2") |>
  dplyr::select(gs_name, gene_symbol) |>
  filter(grepl("REACTOME", gs_name))
save(reactome_mm, file = "data/rdata/reactome_mm.RData")

# Perform GSEA (gene set enrichment analysis) ----
load("data/rdata/c5_mm.RData")
load("data/rdata/reactome_mm.RData")
load("data/rdata/glyco_terms.RData")

ranked_i_c_DEGs <- i_c_toptable$rank
names(ranked_i_c_DEGs) <- rownames(i_c_toptable)
ranked_d_c_DEGs <- d_c_toptable$rank
names(ranked_d_c_DEGs) <- rownames(d_c_toptable)
ranked_d_i_DEGs <- d_i_toptable$rank
names(ranked_d_i_DEGs) <- rownames(d_i_toptable)
ranked_a_c_DEGs <- a_c_toptable$rank
names(ranked_a_c_DEGs) <- rownames(a_c_toptable)
ranked_a_i_DEGs <- a_i_toptable$rank
names(ranked_a_i_DEGs) <- rownames(a_i_toptable)
ranked_a_d_DEGs <- a_d_toptable$rank
names(ranked_a_d_DEGs) <- rownames(a_d_toptable)

i_c_gsea <- run_gsea(ranked_i_c_DEGs, "Inflamed vs Control", c5_mm, reactome_mm)
i_c_gsea.g <- i_c_gsea %>% dplyr::filter(grepl(paste(glyco_terms, collapse = "|"), as.character(ID), ignore.case = TRUE))
d_c_gsea <- run_gsea(ranked_d_c_DEGs, "Dysplasia vs Control", c5_mm, reactome_mm)
d_c_gsea.g <- d_c_gsea %>% dplyr::filter(grepl(paste(glyco_terms, collapse = "|"), as.character(ID), ignore.case = TRUE))
d_i_gsea <- run_gsea(ranked_d_i_DEGs, "Dysplasia vs Inflamed", c5_mm, reactome_mm)
d_i_gsea.g <- d_i_gsea %>% dplyr::filter(grepl(paste(glyco_terms, collapse = "|"), as.character(ID), ignore.case = TRUE))
a_c_gsea <- run_gsea(ranked_a_c_DEGs, "Adenocarcinoma vs Control", c5_mm, reactome_mm)
a_c_gsea.g <- a_c_gsea %>% dplyr::filter(grepl(paste(glyco_terms, collapse = "|"), as.character(ID), ignore.case = TRUE))
a_i_gsea <- run_gsea(ranked_a_i_DEGs, "Adenocarcinoma vs Inflamed", c5_mm, reactome_mm)
a_i_gsea.g <- a_i_gsea %>% dplyr::filter(grepl(paste(glyco_terms, collapse = "|"), as.character(ID), ignore.case = TRUE))
a_d_gsea <- run_gsea(ranked_a_d_DEGs, "Adenocarcinoma vs Dysplasia", c5_mm, reactome_mm)
a_d_gsea.g <- a_d_gsea %>% dplyr::filter(grepl(paste(glyco_terms, collapse = "|"), as.character(ID), ignore.case = TRUE))

gsea_bubble(i_c_gsea, comp = "i_c", path = "output/GSE31106/gsea_plots/separate/", n = 30)
gsea_bubble(i_c_gsea.g, comp = "i_c_glyco", path = "output/GSE31106/gsea_plots/separate/", n = 30)
gsea_bubble(d_c_gsea, comp = "d_c", path = "output/GSE31106/gsea_plots/separate/", n = 30)
gsea_bubble(d_c_gsea.g, comp = "d_c_glyco", path = "output/GSE31106/gsea_plots/separate/", n = 30)
gsea_bubble(d_i_gsea, comp = "d_i", path = "output/GSE31106/gsea_plots/separate/", n = 30)
gsea_bubble(d_i_gsea.g, comp = "d_i_glyco", path = "output/GSE31106/gsea_plots/separate/", n = 30)
gsea_bubble(a_c_gsea, comp = "a_c", path = "output/GSE31106/gsea_plots/separate/", n = 30)
gsea_bubble(a_c_gsea.g, comp = "a_c_glyco", path = "output/GSE31106/gsea_plots/separate/", n = 30)
gsea_bubble(a_i_gsea, comp = "a_i", path = "output/GSE31106/gsea_plots/separate/", n = 30)
gsea_bubble(a_i_gsea.g, comp = "a_i_glyco", path = "output/GSE31106/gsea_plots/separate/", n = 30)
gsea_bubble(a_d_gsea, comp = "a_d", path = "output/GSE31106/gsea_plots/separate/", n = 30)
gsea_bubble(a_d_gsea.g, comp = "a_d_glyco", path = "output/GSE31106/gsea_plots/separate/", n = 30)

pdfs <- list.files(path = "output/GSE31106/gsea_plots/separate/", full.names = TRUE)
qpdf::pdf_combine(input = pdfs[grepl("i_c_", pdfs)],
                  output = "output/GSE31106/gsea_plots/i_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("d_c_", pdfs)],
                  output = "output/GSE31106/gsea_plots/d_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("d_i_", pdfs)],
                  output = "output/GSE31106/gsea_plots/d_i_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("a_c_", pdfs)],
                  output = "output/GSE31106/gsea_plots/a_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("a_i_", pdfs)],
                  output = "output/GSE31106/gsea_plots/a_i_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("a_d_", pdfs)],
                  output = "output/GSE31106/gsea_plots/a_d_combined.pdf")

# GSVA: pathway enrichment per sample: patient sub-populations ----
GO_BP <- c5_mm |> dplyr::filter(grepl("GOBP", gs_name))
GO_BP <- split(x = GO_BP$gene_symbol, f = GO_BP$gs_name)

GO_MF <- c5_mm |> dplyr::filter(grepl("GOMF", gs_name))
GO_MF <- split(x = GO_MF$gene_symbol, f = GO_MF$gs_name)

GO_CC <- c5_mm |>  dplyr::filter(grepl("GOCC", gs_name))
GO_CC <- split(x = GO_CC$gene_symbol, f = GO_CC$gs_name)

GO_RE <- reactome_mm |>  dplyr::filter(grepl("REACTOME", gs_name))
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
# i_c_
i_c_GSVA_BP <- pull_gsva(gsva = gsva_BP$res, decidetests = gsva_BP$decidetests,
                         coeff = "inflamed_control", samples = c(inflammation_samples, control_samples))
i_c_GSVA_MF <- pull_gsva(gsva = gsva_MF$res, decidetests = gsva_MF$decidetests,
                         coeff = "inflamed_control", samples = c(inflammation_samples, control_samples))
i_c_GSVA_CC <- pull_gsva(gsva = gsva_CC$res, decidetests = gsva_CC$decidetests,
                         coeff = "inflamed_control", samples = c(inflammation_samples, control_samples))
i_c_GSVA_RE <- pull_gsva(gsva = gsva_RE$res, decidetests = gsva_RE$decidetests,
                         coeff = "inflamed_control", samples = c(inflammation_samples, control_samples))
i_c_GSVA_TOTAL <- rbind(i_c_GSVA_BP, i_c_GSVA_MF, i_c_GSVA_CC, i_c_GSVA_RE)
i_c_GSVA_TOTAL.g <- i_c_GSVA_TOTAL[grepl(paste(glyco_terms, collapse = "|"), rownames(i_c_GSVA_TOTAL), ignore.case = TRUE), ]

# d_c_
d_c_GSVA_BP <- pull_gsva(gsva = gsva_BP$res, decidetests = gsva_BP$decidetests,
                         coeff = "dysplasia_control", samples = c(dysplasia_samples, control_samples))
d_c_GSVA_MF <- pull_gsva(gsva = gsva_MF$res, decidetests = gsva_MF$decidetests,
                         coeff = "dysplasia_control", samples = c(dysplasia_samples, control_samples))
d_c_GSVA_CC <- pull_gsva(gsva = gsva_CC$res, decidetests = gsva_CC$decidetests,
                         coeff = "dysplasia_control", samples = c(dysplasia_samples, control_samples))
d_c_GSVA_RE <- pull_gsva(gsva = gsva_RE$res, decidetests = gsva_RE$decidetests,
                         coeff = "dysplasia_control", samples = c(dysplasia_samples, control_samples))
d_c_GSVA_TOTAL <- rbind(d_c_GSVA_BP, d_c_GSVA_MF, d_c_GSVA_CC, d_c_GSVA_RE)
d_c_GSVA_TOTAL.g <- d_c_GSVA_TOTAL[grepl(paste(glyco_terms, collapse = "|"), rownames(d_c_GSVA_TOTAL), ignore.case = TRUE), ]

# d_i_
d_i_GSVA_BP <- pull_gsva(gsva = gsva_BP$res, decidetests = gsva_BP$decidetests,
                         coeff = "dysplasia_inflamed", samples = c(dysplasia_samples, inflammation_samples))
d_i_GSVA_MF <- pull_gsva(gsva = gsva_MF$res, decidetests = gsva_MF$decidetests,
                         coeff = "dysplasia_inflamed", samples = c(dysplasia_samples, inflammation_samples))
d_i_GSVA_CC <- pull_gsva(gsva = gsva_CC$res, decidetests = gsva_CC$decidetests,
                         coeff = "dysplasia_inflamed", samples = c(dysplasia_samples, inflammation_samples))
d_i_GSVA_RE <- pull_gsva(gsva = gsva_RE$res, decidetests = gsva_RE$decidetests,
                         coeff = "dysplasia_inflamed", samples = c(dysplasia_samples, inflammation_samples))
d_i_GSVA_TOTAL <- rbind(d_i_GSVA_BP, d_i_GSVA_MF, d_i_GSVA_CC, d_i_GSVA_RE)
d_i_GSVA_TOTAL.g <- d_i_GSVA_TOTAL[grepl(paste(glyco_terms, collapse = "|"), rownames(d_i_GSVA_TOTAL), ignore.case = TRUE), ]

# a_c_
a_c_GSVA_BP <- pull_gsva(gsva = gsva_BP$res, decidetests = gsva_BP$decidetests,
                         coeff = "adenocarcinoma_control", samples = c(adenocarcinoma_samples, control_samples))
a_c_GSVA_MF <- pull_gsva(gsva = gsva_MF$res, decidetests = gsva_MF$decidetests,
                         coeff = "adenocarcinoma_control", samples = c(adenocarcinoma_samples, control_samples))
a_c_GSVA_CC <- pull_gsva(gsva = gsva_CC$res, decidetests = gsva_CC$decidetests,
                         coeff = "adenocarcinoma_control", samples = c(adenocarcinoma_samples, control_samples))
a_c_GSVA_RE <- pull_gsva(gsva = gsva_RE$res, decidetests = gsva_RE$decidetests,
                         coeff = "adenocarcinoma_control", samples = c(adenocarcinoma_samples, control_samples))
a_c_GSVA_TOTAL <- rbind(a_c_GSVA_BP, a_c_GSVA_MF, a_c_GSVA_CC, a_c_GSVA_RE)
a_c_GSVA_TOTAL.g <- a_c_GSVA_TOTAL[grepl(paste(glyco_terms, collapse = "|"), rownames(a_c_GSVA_TOTAL), ignore.case = TRUE), ]

# a_i_
a_i_GSVA_BP <- pull_gsva(gsva = gsva_BP$res, decidetests = gsva_BP$decidetests,
                         coeff = "adenocarcinoma_inflamed", samples = c(adenocarcinoma_samples, inflammation_samples))
a_i_GSVA_MF <- pull_gsva(gsva = gsva_MF$res, decidetests = gsva_MF$decidetests,
                         coeff = "adenocarcinoma_inflamed", samples = c(adenocarcinoma_samples, inflammation_samples))
a_i_GSVA_CC <- pull_gsva(gsva = gsva_CC$res, decidetests = gsva_CC$decidetests,
                         coeff = "adenocarcinoma_inflamed", samples = c(adenocarcinoma_samples, inflammation_samples))
a_i_GSVA_RE <- pull_gsva(gsva = gsva_RE$res, decidetests = gsva_RE$decidetests,
                         coeff = "adenocarcinoma_inflamed", samples = c(adenocarcinoma_samples, inflammation_samples))
a_i_GSVA_TOTAL <- rbind(a_i_GSVA_BP, a_i_GSVA_MF, a_i_GSVA_CC, a_i_GSVA_RE)
a_i_GSVA_TOTAL.g <- a_i_GSVA_TOTAL[grepl(paste(glyco_terms, collapse = "|"), rownames(a_i_GSVA_TOTAL), ignore.case = TRUE), ]

# a_d_
a_d_GSVA_BP <- pull_gsva(gsva = gsva_BP$res, decidetests = gsva_BP$decidetests,
                         coeff = "adenocarcinoma_dysplasia", samples = c(adenocarcinoma_samples, dysplasia_samples))
a_d_GSVA_MF <- pull_gsva(gsva = gsva_MF$res, decidetests = gsva_MF$decidetests,
                         coeff = "adenocarcinoma_dysplasia", samples = c(adenocarcinoma_samples, dysplasia_samples))
a_d_GSVA_CC <- pull_gsva(gsva = gsva_CC$res, decidetests = gsva_CC$decidetests,
                         coeff = "adenocarcinoma_dysplasia", samples = c(adenocarcinoma_samples, dysplasia_samples))
a_d_GSVA_RE <- pull_gsva(gsva = gsva_RE$res, decidetests = gsva_RE$decidetests,
                         coeff = "adenocarcinoma_dysplasia", samples = c(adenocarcinoma_samples, dysplasia_samples))
a_d_GSVA_TOTAL <- rbind(a_d_GSVA_BP, a_d_GSVA_MF, a_d_GSVA_CC, a_d_GSVA_RE)
a_d_GSVA_TOTAL.g <- a_d_GSVA_TOTAL[grepl(paste(glyco_terms, collapse = "|"), rownames(a_d_GSVA_TOTAL), ignore.case = TRUE), ]

# GSVA heatmap ----
heatmap_colors <- paletteer_c("pals::coolwarm", 100, direction = 1)
phenotype_colors <- paletteer_d("lisa::JosefAlbers", 2)
i_c_colors <- setNames(phenotype_colors, c("inflamed", "control"))
d_c_colors <- setNames(phenotype_colors, c("dysplasia", "control"))
d_i_colors <- setNames(phenotype_colors, c("dysplasia", "inflamed"))
a_c_colors <- setNames(phenotype_colors, c("adenocarcinoma", "control"))
a_i_colors <- setNames(phenotype_colors, c("adenocarcinoma", "inflamed"))
a_d_colors <- setNames(phenotype_colors, c("adenocarcinoma", "dysplasia"))

### print for thesis draft
heatmap_colors <- paletteer_c("pals::coolwarm", 100, direction = 1)
a_i_GSVA_BP.g <- a_i_GSVA_BP[grepl(paste(glyco_terms, collapse = "|"), rownames(a_i_GSVA_BP), ignore.case = TRUE), ]
hm <- pheatmap(a_i_GSVA_BP.g,
         show_colnames = T,
         show_rownames = T,
         fontsize_col = 10,
         scale = "row",
         treeheight_row = 0,
         fontsize_row = 10,
         color = heatmap_colors,
         annotation_colors = list(phenotype = a_i_colors),
         annotation_col = column_to_rownames(study_design, "sampleID"),
         breaks = seq(-2, 2,length.out = 101),
         cutree_rows = 2)

row_tree <- hm$tree_row
clusters <- cutree(row_tree, k = 2)
df <- data.frame(cluster = clusters, pathway = row_tree$labels) %>% arrange(cluster)
df$pathway <- gsub("^GOBP_", "", df$pathway)
write_csv(df, file = "a_i_gsva.csv")

#i wanna keep the most relevant values
test <- rowMeans(abs(a_i_GSVA_BP))
names(test) <- gsub("^GOBP_", "", names(test))
test <- sort(test, decreasing = TRUE)
top_terms <- test[1:50]
View(df[df$pathway %in% names(top_terms),])
###


make_heatmap(data = i_c_GSVA_TOTAL, phencols = i_c_colors,
             path = "output/GSE31106/gsva_plots/separate/", filename = "i_c_gsva.pdf")
make_heatmap(data = i_c_GSVA_TOTAL.g[,2:ncol(i_c_GSVA_TOTAL.g)], phencols = i_c_colors,
             path = "output/GSE31106/gsva_plots/separate/", filename = "i_c_glyco_gsva.pdf")

make_heatmap(data = d_c_GSVA_TOTAL, phencols = d_c_colors,
             path = "output/GSE31106/gsva_plots/separate/", filename = "d_c_gsva.pdf")
make_heatmap(data = d_c_GSVA_TOTAL.g[,2:ncol(d_c_GSVA_TOTAL.g)], phencols = d_c_colors,
             path = "output/GSE31106/gsva_plots/separate/", filename = "d_c_glyco_gsva.pdf")

make_heatmap(data = d_i_GSVA_TOTAL, phencols = d_i_colors,
             path = "output/GSE31106/gsva_plots/separate/", filename = "d_i_gsva.pdf")
make_heatmap(data = d_i_GSVA_TOTAL.g[,2:ncol(d_i_GSVA_TOTAL.g)], phencols = d_i_colors,
             path = "output/GSE31106/gsva_plots/separate/", filename = "d_i_glyco_gsva.pdf")

make_heatmap(data = a_c_GSVA_TOTAL, phencols = a_c_colors,
             path = "output/GSE31106/gsva_plots/separate/", filename = "a_c_gsva.pdf")
make_heatmap(data = a_c_GSVA_TOTAL.g[,2:ncol(a_c_GSVA_TOTAL.g)], phencols = a_c_colors,
             path = "output/GSE31106/gsva_plots/separate/", filename = "a_c_glyco_gsva.pdf")

make_heatmap(data = a_i_GSVA_TOTAL, phencols = a_i_colors,
             path = "output/GSE31106/gsva_plots/separate/", filename = "a_i_gsva.pdf")
make_heatmap(data = a_i_GSVA_TOTAL.g[,2:ncol(a_i_GSVA_TOTAL.g)], phencols = a_i_colors,
             path = "output/GSE31106/gsva_plots/separate/", filename = "a_i_glyco_gsva.pdf")

make_heatmap(data = a_d_GSVA_TOTAL, phencols = a_d_colors,
             path = "output/GSE31106/gsva_plots/separate/", filename = "a_d_gsva.pdf")
make_heatmap(data = a_d_GSVA_TOTAL.g[,2:ncol(a_d_GSVA_TOTAL.g)], phencols = a_d_colors,
             path = "output/GSE31106/gsva_plots/separate/", filename = "a_d_glyco_gsva.pdf")

pdfs <- list.files(path = "output/GSE31106/gsva_plots/separate/", full.names = TRUE)
qpdf::pdf_combine(input = pdfs[grepl("i_c_", pdfs)],
                  output = "output/GSE31106/gsva_plots/i_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("d_c_", pdfs)],
                  output = "output/GSE31106/gsva_plots/d_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("d_i_", pdfs)],
                  output = "output/GSE31106/gsva_plots/d_i_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("a_c_", pdfs)],
                  output = "output/GSE31106/gsva_plots/a_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("a_i_", pdfs)],
                  output = "output/GSE31106/gsva_plots/a_i_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("a_d_", pdfs)],
                  output = "output/GSE31106/gsva_plots/a_d_combined.pdf")

# Save Excel sheets ----
openxlsx::write.xlsx(list("all_genes" = i_c_toptable,
                          "all_significant" = i_c_toptable.f,
                          "all_glyco" = i_c_toptable.g,
                          "significant_glyco" = i_c_toptable.fg,
                          "full_gsea" = i_c_gsea,
                          "glyco_gsea" = i_c_gsea.g),
                     file = "output/GSE31106/tables/i_c_tables.xlsx", rowNames = TRUE)

openxlsx::write.xlsx(list("all_genes" = d_c_toptable,
                          "all_significant" = d_c_toptable.f,
                          "all_glyco" = d_c_toptable.g,
                          "significant_glyco" = d_c_toptable.fg,
                          "full_gsea" = d_c_gsea,
                          "glyco_gsea" = d_c_gsea.g),
                     file = "output/GSE31106/tables/d_c_tables.xlsx", rowNames = TRUE)

openxlsx::write.xlsx(list("all_genes" = d_i_toptable,
                          "all_significant" = d_i_toptable.f,
                          "all_glyco" = d_i_toptable.g,
                          "significant_glyco" = d_i_toptable.fg,
                          "full_gsea" = d_i_gsea,
                          "glyco_gsea" = d_i_gsea.g),
                     file = "output/GSE31106/tables/d_i_tables.xlsx", rowNames = TRUE)

openxlsx::write.xlsx(list("all_genes" = a_c_toptable,
                          "all_significant" = a_c_toptable.f,
                          "all_glyco" = a_c_toptable.g,
                          "significant_glyco" = a_c_toptable.fg,
                          "full_gsea" = a_c_gsea,
                          "glyco_gsea" = a_c_gsea.g),
                     file = "output/GSE31106/tables/a_c_tables.xlsx", rowNames = TRUE)

openxlsx::write.xlsx(list("all_genes" = a_i_toptable,
                          "all_significant" = a_i_toptable.f,
                          "all_glyco" = a_i_toptable.g,
                          "significant_glyco" = a_i_toptable.fg,
                          "full_gsea" = a_i_gsea,
                          "glyco_gsea" = a_i_gsea.g),
                     file = "output/GSE31106/tables/a_i_tables.xlsx", rowNames = TRUE)

openxlsx::write.xlsx(list("all_genes" = a_d_toptable,
                          "all_significant" = a_d_toptable.f,
                          "all_glyco" = a_d_toptable.g,
                          "significant_glyco" = a_d_toptable.fg,
                          "full_gsea" = a_d_gsea,
                          "glyco_gsea" = a_d_gsea.g),
                     file = "output/GSE31106/tables/a_d_tables.xlsx", rowNames = TRUE)
