# GSE47908 "Bjerrum" cohort - bulk microarray
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
library(ggvenn)
library(gplots)
library(gridExtra)
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

# Set up data ----
# Choose pvalue and logFC
logfc <- 2
pval <- 0.05

# Import glycogenes
load("data/rdata/glycogenes.RData")
load("data/rdata/glycogenes2.RData")
load("~/projects/tesis/output/integration/glycogenes2.RData")

# Import metadata
metadata <- read_tsv("data/raw/GSE47908/GSE47908_series_matrix.txt",
                     skip = 32, col_names = FALSE) %>% t %>% as_tibble()
colnames(metadata) <- metadata[1, ]
metadata <- metadata[-1, ]

# Check the data processing
processing <- unlist(metadata[1, 19])

# Read expression data (.CEL files) and perform RMA normalization
setwd("~/projects/tesis/data/raw/GSE47908")
expression <- exprs(rma(ReadAffy())) # rma returns expression in log2 scale
setwd("~/projects/tesis")

# Alternatively, read expression file
expression <- read.csv("data/raw/GSE47908/GSE47908_series_matrix.txt", header = T, skip = 33, sep = "")
expression <- expression[34:nrow(expression) - 1,]
rownames(expression) <- expression[,1]
expression <- expression[,-1]
probeIDs <- rownames(expression)
expression <- sapply(expression, as.numeric)
rownames(expression) <- probeIDs
rownames(expression) <- mapIds(hgu133plus2.db,
                               keys = rownames(expression),
                               column = "SYMBOL",
                               keytype = "PROBEID",
                               multiVals = "first")

#colnames(expression) <- sub("_.*", "", colnames(expression))
expression <- expression[!is.na(rownames(expression)),]
anyDuplicated(rownames(expression)) # check repeated rownames 
expression <- avereps(expression) # average values per rowname

# After PCA, realize I need to remove a sample
# expression <- expression[, colnames(expression) != "GSM1162246"]

# Make study_design table
study_design <- select(metadata, 33, 11)
colnames(study_design) <- c("sampleID", "phenotype")
study_design <- study_design %>%
  mutate(phenotype = case_when(grepl("control", phenotype) ~ "control",
                               grepl("dysplasia", phenotype) ~ "dysplasia",
                               grepl("left", phenotype) ~ "L_UC",
                               grepl("pancolitis", phenotype) ~ "pancolitis")) # |> filter(!sampleID == "GSM1162246")
#study_design <- study_design[!study_design$sampleID %in% c("GSM1162246", "GSM1162279"),]
table(study_design$phenotype)
save(study_design, file = "data/study_designs/bjerrum_study_design.RData")

# Remove sample outliers after PCA
to_remove <- c("GSM1162279", "GSM1162246")
study_design <- study_design[!study_design$sampleID %in% to_remove,]

# Capture sample IDs and phenotypes
sample_names <- study_design$sampleID
phenotypes <- factor(study_design$phenotype)
expression <- expression[, colnames(expression) %in% sample_names]
colnames(expression)
# Remove pancolitis from expression data expression <- expression[, sample_names]

# Plot data quality ----
expression.t <- expression |> as_tibble(rownames = "geneID")
expression.tp <- expression.t |>
  pivot_longer(cols = colnames(expression.t[, -1]),
               names_to = "samples",
               values_to = "expression")

# Plot data density to check if need to remove lowly expressed genes
density_plot <- ggplot(data=expression.tp, aes(x=expression, color=samples)) +
  geom_density(adjust=1.5) +
  theme_classic()

# Plot normalized data
eplot <- ggplot(expression.tp) +
  aes(x = samples, y = expression, fill = samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median",
               geom = "point",
               shape = 95,
               size = 10,
               color = "black",
               show.legend = FALSE) +
  labs(y = "log2 fluorescence values",
       x = "sample",
       title = "Expression data") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_bw()

ggsave("output/GSE47908/quality_plots/density_plot.pdf", density_plot, device = pdf, width = 8, height = 5, dpi = 300)
ggsave("output/GSE47908/quality_plots/data_quality.pdf", eplot, device = pdf, width = 8, height = 5, dpi = 300)

# PCA test ----
pca_expression <- expression[, ! colnames(expression) %in% c("GSM1162246", "GSM1162279")]
#raw_expression <- expression
#pca_expression <- raw_expression
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
  geom_point(size = 5) +
  labs(color = "Grupo biológico") +
  scale_color_manual(values = c("cornflowerblue", "coral1", "orchid", "darkgreen"),
                     labels = c("Control", "Displasia", "L_UC", "Pancolitis")) +
  xlab(paste0("PC1 (", pc_percentages[1, 2], "%", ")")) +
  ylab(paste0("PC2 (", pc_percentages[2, 2], "%", ")")) +
  #geom_label_repel(aes(label = sample_names), size = 3) +
  theme_bw(base_size = 18)
print(pca_plot)

pca_plot2 <- ggplot(pca_data) +
  aes(x = PC1, y = PC2, color = pca_phenotypes) +
  geom_point(size = 5) +
  labs(color = "Grupo biológico") +
  scale_color_manual(values = c("cornflowerblue", "coral1", "orchid", "darkgreen"),
                     labels = c("Control", "Displasia", "L_UC", "Pancolitis")) +
  xlab(paste0("PC1 (", pc_percentages[1, 2], "%", ")")) +
  ylab(paste0("PC2 (", pc_percentages[2, 2], "%", ")")) +
  #geom_label_repel(aes(label = sample_names, size = 3)) +
  theme_bw(base_size = 18)

print(pca_plot2)

ggsave(path = "output/GSE47908/quality_plots/",
       filename = "pca_before_sample_removal.pdf",
       pca_plot, device = "pdf", width = 6, height = 5, dpi = 300)

pca_samples <- study_design[study_design$phenotype != "control",]$sampleID
pca_expression <- expression[,colnames(expression) %in% pca_samples]
pca_phenotypes <- study_design$phenotype[study_design$sampleID %in% pca_samples]
pca <- prcomp(t(pca_expression), scale. = T, center = T, retx = T)
eigenvalues <- pca$sdev^2
pc_percentages <- round((eigenvalues) / sum(eigenvalues) * 100, 1)
PCs <- paste0("PC", 1:length(pc_percentages))
pc_percentages <- data.frame(Component = PCs, Percentage = pc_percentages)
pca_data <- as_tibble(pca$x)

ggplot(pca_data) +
  aes(x = PC1, y = PC2, color = pca_phenotypes) +
  geom_point(size = 5) +
  labs(color = "Grupo biológico") +
  scale_color_manual(values = c("cornflowerblue", "coral1", "orchid"),
                     labels = c("Displasia", "L-UC", "Pancolitis")) +
  xlab(paste0("PC1 (", pc_percentages[1, 2], "%", ")")) +
  ylab(paste0("PC2 (", pc_percentages[2, 2], "%", ")")) +
  theme_bw(base_size = 18)
print(pca_plot)

# Find out which genes influence the PC1 more
rotation <- pca$rotation[,colnames(pca$rotation) == "PC1"]
rotation.sorted <- sort(rotation, decreasing = TRUE) 
top_contributors <- names(sort(abs(rotation), decreasing = TRUE)[1:10])
rotation_top <- rotation.sorted[top_contributors]
down_loadings <- names(rotation_top[rotation_top < 0])
up_loadings <- names(rotation_top[rotation_top > 0])

# Biplot
pca_small <- pca
pca_small$rotation <- pca$rotation[top_contributors,]
biplot(pca_small, scale = 0, col = c("white", "deeppink3"))
top_contributors

ora_down <- gprofiler2::gost(down_loadings, organism = "hsapiens", sources = c("GO"), significant = TRUE)
ora_up <- gprofiler2::gost(up_loadings, organism = "hsapiens", sources = c("GO"), significant = TRUE)
ora_down.res <- ora_down$result[ora_down$result$term_size < 500,]
ora_up.res <- ora_up$result[ora_up$result$term_size < 500,]
ora_down_filtered <- list(result = ora_down.res, meta = ora_down$meta)
ora_up_filtered <- list(result = ora_up.res, meta = ora_up$meta)
class(ora_down_filtered) <- "gostres"
class(ora_up_filtered) <- "gostres"
gprofiler2::gostplot(ora_down_filtered, interactive = T, capped = T)
gprofiler2::gostplot(ora_up_filtered, interactive = T, capped = T)

rotation.df <- data.frame(valor = rotation.sorted, posicion = c(1:21358))

ggplot(rotation.df, aes(x = valor, y = posicion)) +
       geom_point(size = 1, colour = "deeppink3") +
  theme_bw(base_size = 15)

# GSEA of the rotation genes
pc1_gsea_reactome <- GSEA(rotation.sorted, TERM2GENE = reactome_hs, verbose = TRUE, eps = 0, minGSSize = 5, maxGSSize = 500)
pc1_gsea_c5 <- GSEA(rotation.sorted, TERM2GENE = c5_hs, verbose = TRUE, eps = 0, minGSSize = 5, maxGSSize = 500)
pc1_gobp_gsea.res <- pc1_gsea_c5@result %>% dplyr::filter(grepl("GOBP", ID))
pc1_gobp_gsea.res.f <- pc1_gobp_gsea.res[pc1_gobp_gsea.res$p.adjust < 0.00000001,]
nrow(pc1_gobp_gsea.res.f)
selection <- c("GOBP_REGULATION_OF_CELL_CYCLE_PHASE_TRANSITION",
               "GOBP_MITOTIC_NUCLEAR_DIVISION",
               "GOBP_PROTEIN_LOCALIZATION_TO_NUCLEUS",
               "GOBP_CELL_KILLING",
               "GOBP_CHROMOSOME_SEGREGATION",
               "GOBP_REGULATION_OF_MITOTIC_CELL_CYCLE",
               "GOBP_DNA_REPLICATION",
               "GOBP_NEGATIVE_REGULATION_OF_CELL_CYCLE",
               "GOBP_MONOCARBOXYLIC_ACID_CATABOLIC_PROCESS",
               "GOBP_FATTY_ACID_BETA_OXIDATION",
               "GOBP_OXIDATIVE_PHOSPHORYLATION",
               "GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY",
               "GOBP_LEUKOCYTE_MIGRATION",
               "GOBP_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS",
               "GOBP_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN",
               "GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_CELL_CELL_ADHESION",
               "GOBP_ACTIVATION_OF_IMMUNE_RESPONSE",
               "GOBP_IMMUNE_RESPONSE_REGULATING_SIGNALING_PATHWAY",
               "GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY",
               "GOBP_LEUKOCYTE_CELL_CELL_ADHESION",
               "GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION",
               "GOBP_LEUKOCYTE_MEDIATED_IMMUNITY",
               "GOBP_T_CELL_ACTIVATION",
               "GOBP_REGULATION_OF_LYMPHOCYTE_ACTIVATION",
               "GOBP_ADAPTIVE_IMMUNE_RESPONSE")
gsea_data <- pc1_gobp_gsea.res.f[rownames(pc1_gobp_gsea.res.f) %in% selection,]

pc1_gsea.res <- rbind(pc1_gsea_reactome@result, pc1_gsea_c5@result)
gsea_data <- pc1_gsea.res %>%
  #dplyr::filter(grepl("GOBP", ID)) %>%
  dplyr::mutate(NES_category = ifelse(NES > 0, "Positive", "Negative"))
gsea_data$comparison <- "PC1 GSEA"

ggplot(gsea_data, aes(x = comparison, y = reorder(ID, -log10(p.adjust)))) +
      geom_point(aes(size = setSize, color = NES_category, alpha = -log10(p.adjust))) +
      scale_color_manual(values = c("Positive" = "red", "Negative" = "blue")) +
      labs(title = paste("GOBP", gsea_data$comparison), 
           x = element_blank(),
           y = "Pathway",
           size = "Gene Set Size",
           color = "NES Direction",
           alpha = "-log10(p.adjust)") +
      theme_bw(base_size = 15) +
      theme(plot.title = element_text(hjust = 1, size = 14, face = "bold"))

openxlsx::write.xlsx(pc1_gsea.res, file = "~/projects/tesis/final_outputs/displasia_vs_uc/GSEA_PC1.xlsx", rowNames = TRUE)

# Linear model ----
lm_dmatrix <- model.matrix(~ 0 + phenotypes)
colnames(lm_dmatrix) <- levels(phenotypes)
rownames(lm_dmatrix) <- sample_names

# Master lineal model
master_lm <- lmFit(expression, lm_dmatrix)

# Contrast matrix (pairwise comparisons: UC/C, D/UC)
lm_cmatrix <- makeContrasts(UC_control = L_UC - control,
                            dysplasia_UC = dysplasia - L_UC,
                            dysplasia_control = dysplasia - control,
                            pancolitis_uc = pancolitis - L_UC,
                            dysplasia_pancolitis = dysplasia - pancolitis,
                            levels = lm_dmatrix)

# Linear model of contrasts and bayesian statistics
lm <- contrasts.fit(master_lm, lm_cmatrix)
lm_stats <- eBayes(lm)

# TopTable and DEGs -----
# Make toptables
uc_c_toptable <- make_toptable("UC_control", "left-side UC - control")
d_c_toptable <- make_toptable("dysplasia_control", "dysplasia - control")
d_uc_toptable <- make_toptable("dysplasia_UC", "dysplasia - left-side UC")
pan_uc_toptable <- make_toptable("pancolitis_uc", "pancolitis - left-side UC")
dys_pan_toptable <- make_toptable("dysplasia_pancolitis", "dysplasia - pancolitis")

# Filter significant DEGs
# uc_c_toptable.f <- uc_c_toptable %>% dplyr::filter(adj.P.Val <= pval)
# d_c_toptable.f <- d_c_toptable %>% dplyr::filter(adj.P.Val <= pval)
# d_uc_toptable.f <- d_uc_toptable %>% dplyr::filter(adj.P.Val <= pval)
# pan_uc_toptable.f <- pan_uc_toptable %>% dplyr::filter(adj.P.Val <= pval)

# Filter glycogenes
uc_c_toptable.g <- uc_c_toptable %>% dplyr::filter(row.names(uc_c_toptable) %in% glycogenes2)
d_c_toptable.g <- d_c_toptable %>% dplyr::filter(row.names(d_c_toptable) %in% glycogenes2)
d_uc_toptable.g <- d_uc_toptable %>% dplyr::filter(row.names(d_uc_toptable) %in% glycogenes2)
pan_uc_toptable.g <- pan_uc_toptable %>% dplyr::filter(row.names(pan_uc_toptable) %in% glycogenes2)
dys_pan_toptable.g <- dys_pan_toptable %>% dplyr::filter(row.names(dys_pan_toptable) %in% glycogenes2)

# Filter significant glycogenes
# uc_c_toptable.fg <- uc_c_toptable.f %>% dplyr::filter(row.names(uc_c_toptable.f) %in% glycogenes)
# d_c_toptable.fg <- d_c_toptable.f %>% dplyr::filter(row.names(d_c_toptable.f) %in% glycogenes)
# d_uc_toptable.fg <- d_uc_toptable.f %>% dplyr::filter(row.names(d_uc_toptable.f) %in% glycogenes)
# pan_uc_toptable.fg <- pan_uc_toptable.f %>% dplyr::filter(row.names(pan_uc_toptable.f) %in% glycogenes)

# Make and save volcano plots----
logfc <- 3
pval
data <- as_tibble(dys_pan_toptable.g, rownames = "geneID")
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
ggplotly(v)

volcano(uc_c_toptable, p = pval, lfc = logfc, title = "DEGs in left-sided UC vs. control",
        path = "output/GSE47908/volcano_plots/separate/", filename = "uc_c_volcano.pdf")
volcano(d_c_toptable, p = pval, lfc = logfc, title = "DEGs in dysplasia vs. control",
        path = "output/GSE47908/volcano_plots/separate/", filename = "d_c_volcano.pdf")
volcano(d_uc_toptable, p = pval, lfc = logfc, title = "DEGs in dysplasia vs. left-sided UC",
        path = "output/GSE47908/volcano_plots/separate/", filename = "d_uc_volcano.pdf")
volcano(pan_uc_toptable, p = pval, lfc = logfc, title = "DEGs in pancolitis vs. left-sided UC",
        path = "output/GSE47908/volcano_plots/separate/", filename = "pan_uc_volcano.pdf")
volcano(uc_c_toptable.g, p = pval, lfc = logfc, title = "GlycoDEGs in left-sided UC vs. control",
        path = "output/GSE47908/volcano_plots/separate/", filename = "uc_c_volcano_glyco.pdf")
volcano(d_c_toptable.g, p = pval, lfc = logfc, title = "GlycoDEGs in dysplasia vs. control",
        path = "output/GSE47908/volcano_plots/separate/", filename = "d_c_volcano_glyco.pdf")
volcano(d_uc_toptable.g, p = pval, lfc = logfc, title = "GlycoDEGs in dysplasia vs. left-sided UC",
        path = "output/GSE47908/volcano_plots/separate/", filename = "d_uc_volcano_glyco.pdf")
volcano(pan_uc_toptable.g, p = pval, lfc = logfc, title = "GlycoDEGs in pancolitis vs. left-sided UC",
        path = "output/GSE47908/volcano_plots/separate/", filename = "pan_uc_volcano_glyco.pdf")

# Save combined .pdf according to comparison in title
pdfs <- list.files(path = "output/GSE47908/volcano_plots/separate/", full.names = TRUE)
qpdf::pdf_combine(input = pdfs[grepl("uc_c_", pdfs)],
                  output = "output/GSE47908/volcano_plots/uc_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("d_c_", pdfs)],
                  output = "output/GSE47908/volcano_plots/d_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("d_uc_", pdfs)],
                  output = "output/GSE47908/volcano_plots/d_uc_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("pan_uc_", pdfs)],
                  output = "output/GSE47908/volcano_plots/pan_uc_combined.pdf")

# ORA pancolitis uc ----
upreg <- rownames(pan_uc_toptable.f[pan_uc_toptable.f$direction == "up",])
downreg <- rownames(pan_uc_toptable.f[pan_uc_toptable.f$direction == "down",])
ora_up <- gprofiler2::gost(upreg, organism = "hsapiens", sources = c("GO"), significant = TRUE)
ora_down <- gprofiler2::gost(downreg, organism = "hsapiens", sources = c("GO"), significant = TRUE)
ora_up.res <- ora_up$result[ora_up$result$term_size < 500,]
ora_down.res <- ora_down$result[ora_down$result$term_size < 500,]
gprofiler2::gostplot(ora_up, interactive = T, capped = T)
gprofiler2::gostplot(ora_down, interactive = T, capped = T)

# Plot filtered results! Make custom gostres object
ora_up2 <- list(result = ora_up.res, meta = ora_up$meta)
class(ora_up2) <- "gostres"
gprofiler2::gostplot(ora_up2, interactive = T, capped = T)

ora_down2 <- list(result = ora_down.res, meta = ora_down$meta)
class(ora_down2) <- "gostres"
gprofiler2::gostplot(ora_down2, interactive = T, capped = T)

View(ora_down.res)

# Perform GSEA (gene set enrichment analysis) ----
# Import msigdb collections
load("data/rdata/c5_hs.RData")
load("data/rdata/reactome_hs.RData")

# Construct an ordered named vector to perform GSEA
ranked_uc_c_DEGs <- uc_c_toptable$rank
ranked_d_c_DEGs <- d_c_toptable$rank
ranked_d_uc_DEGs <- d_uc_toptable$rank
ranked_d_pan_DEGS <- dys_pan_toptable$rank
names(ranked_uc_c_DEGs) <- rownames(uc_c_toptable)
names(ranked_d_c_DEGs) <- rownames(d_c_toptable)
names(ranked_d_uc_DEGs) <- rownames(d_uc_toptable)
names(ranked_d_pan_DEGS) <- rownames(dys_pan_toptable)

# Run GSEA
load("data/rdata/glyco_terms.RData")

uc_c_gsea <- run_gsea(ranked_uc_c_DEGs, "Left-side UC vs Control", c5_hs, reactome_hs)
uc_c_gsea.g <- uc_c_gsea %>% dplyr::filter(grepl(paste(glyco_terms, collapse = "|"), as.character(ID), ignore.case = TRUE))

d_c_gsea <- run_gsea(ranked_d_c_DEGs, "Dysplasia vs Control", c5_hs, reactome_hs)
d_c_gsea.g <- d_c_gsea %>% dplyr::filter(grepl(paste(glyco_terms, collapse = "|"), as.character(ID), ignore.case = TRUE))

d_uc_gsea <- run_gsea(ranked_d_uc_DEGs, "Dysplasia vs Left-side UC", c5_hs, reactome_hs)
d_uc_gsea.g <- d_uc_gsea %>% dplyr::filter(grepl(paste(glyco_terms, collapse = "|"), as.character(ID), ignore.case = TRUE))

d_pan_gsea <- run_gsea(ranked_d_pan_DEGS, "Dysplasia vs Pancolitis", c5_hs, reactome_hs)
d_pan_gsea.g <- d_pan_gsea %>% dplyr::filter(grepl(paste(glyco_terms, collapse = "|"), as.character(ID), ignore.case = TRUE))

# Make and save bubble plots
gsea_bubble(uc_c_gsea, comp = "uc_c", path = "output/GSE47908/gsea_plots/separate/", n = 30)
gsea_bubble(uc_c_gsea.g, comp = "uc_c_glyco", path = "output/GSE47908/gsea_plots/separate/", n = 30)
gsea_bubble(d_c_gsea, comp = "d_c", path = "output/GSE47908/gsea_plots/separate/", n = 30)
gsea_bubble(d_c_gsea.g, comp = "d_c_glyco", path = "output/GSE47908/gsea_plots/separate/", n = 30)
gsea_bubble(d_uc_gsea, comp = "d_uc", path = "output/GSE47908/gsea_plots/separate/", n = 30)
gsea_bubble(d_uc_gsea.g, comp = "d_uc_glyco", path = "output/GSE47908/gsea_plots/separate/", n = 30)

# Save combined .pdf according to comparison in title
pdfs <- list.files(path = "output/GSE47908/gsea_plots/separate/", full.names = TRUE)
qpdf::pdf_combine(input = pdfs[grepl("uc_c_", pdfs)],
                  output = "output/GSE47908/gsea_plots/uc_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("d_c_", pdfs)],
                  output = "output/GSE47908/gsea_plots/d_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("d_uc_", pdfs)],
                  output = "output/GSE47908/gsea_plots/d_uc_combined.pdf")


d_uc_gsea.g$source <- ifelse(grepl("^GOBP_", rownames(d_uc_gsea.g)), "GOBP",
                            ifelse(grepl("^GOCC_", rownames(d_uc_gsea.g)), "GOCC",
                                   ifelse(grepl("^REACTOME_", rownames(d_uc_gsea.g)), "REACTOME", 
                                          ifelse(grepl("^GOMF_", rownames(d_uc_gsea.g)), "GOMF", "HP"))))
my_gsea <- d_uc_gsea.g %>% mutate(NES_category = ifelse(NES > 0, "Positivo", "Negativo"))
my_gsea <- my_gsea[!grepl("^HP_", rownames(my_gsea)),]
my_gsea <- my_gsea[1:10,]
my_gsea$comparison <- ""
ggplot(my_gsea, aes(x = comparison, y = reorder(ID, -log10(p.adjust)))) +
  geom_point(aes(size = setSize, color = NES_category, alpha = -log10(p.adjust))) +
  scale_color_manual(values = c("Positivo" = "red", "Negativo" = "blue")) +
  labs(title = "GSEA: displasia vs. CU-I (glicobiológico)", 
       x = element_blank(),
       y = "Vía",
       size = "Tamaño",
       color = "Sentido NES",
       alpha = "-log10(p-ajustado)") +
  facet_wrap(~ my_gsea$source, nrow = 1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 1, size = 14, face = "bold"))

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
gsva_BP <- run_gsva(exp = expression, go = GO_BP, dmatrix = lm_dmatrix, cmatrix = lm_cmatrix, lfc = 0.5)  
gsva_MF <- run_gsva(exp = expression, go = GO_MF, dmatrix = lm_dmatrix, cmatrix = lm_cmatrix, lfc = 0.5)  
gsva_CC <- run_gsva(exp = expression, go = GO_CC, dmatrix = lm_dmatrix, cmatrix = lm_cmatrix, lfc = 0.5)  
gsva_RE <- run_gsva(exp = expression, go = GO_RE, dmatrix = lm_dmatrix, cmatrix = lm_cmatrix, lfc = 0.5)  
summary(gsva_BP$decidetests)
summary(gsva_MF$decidetests)
summary(gsva_CC$decidetests)
summary(gsva_RE$decidetests)

# Extract samples per phenotype
c_samples <- study_design$sampleID[study_design$phenotype == "control"]
uc_samples <- study_design$sampleID[study_design$phenotype == "L_UC"]
d_samples <- study_design$sampleID[study_design$phenotype == "dysplasia"]
pan_samples <- study_design$sampleID[study_design$phenotype == "pancolitis"]

# Pull out enriched GSVA results --> (gsva, decidetests, coeff, samples)
uc_c_GSVA_BP <- pull_gsva(gsva = gsva_BP$res, decidetests = gsva_BP$decidetests,
                          coeff = "UC_control", samples = c(uc_samples, c_samples))
uc_c_GSVA_MF <- pull_gsva(gsva = gsva_MF$res, decidetests = gsva_MF$decidetests,
                          coeff = "UC_control", samples = c(uc_samples, c_samples))
uc_c_GSVA_CC <- pull_gsva(gsva = gsva_CC$res, decidetests = gsva_CC$decidetests,
                          coeff = "UC_control", samples = c(uc_samples, c_samples))
uc_c_GSVA_RE <- pull_gsva(gsva = gsva_RE$res, decidetests = gsva_RE$decidetests,
                          coeff = "UC_control", samples = c(uc_samples, c_samples))
uc_c_GSVA_TOTAL <- rbind(uc_c_GSVA_BP, uc_c_GSVA_MF, uc_c_GSVA_CC, uc_c_GSVA_RE)
uc_c_GSVA_TOTAL.g <- uc_c_GSVA_TOTAL[grepl(paste(glyco_terms, collapse = "|"), rownames(uc_c_GSVA_TOTAL), ignore.case = TRUE), ]

d_c_GSVA_BP <- pull_gsva(gsva = gsva_BP$res, decidetests = gsva_BP$decidetests,
                          coeff = "dysplasia_control", samples = c(d_samples, c_samples))
d_c_GSVA_MF <- pull_gsva(gsva = gsva_MF$res, decidetests = gsva_MF$decidetests,
                          coeff = "dysplasia_control", samples = c(d_samples, c_samples))
d_c_GSVA_CC <- pull_gsva(gsva = gsva_CC$res, decidetests = gsva_CC$decidetests,
                          coeff = "dysplasia_control", samples = c(d_samples, c_samples))
d_c_GSVA_RE <- pull_gsva(gsva = gsva_RE$res, decidetests = gsva_RE$decidetests,
                          coeff = "dysplasia_control", samples = c(d_samples, c_samples))
d_c_GSVA_TOTAL <- rbind(d_c_GSVA_BP, d_c_GSVA_MF, d_c_GSVA_CC, d_c_GSVA_RE)
d_c_GSVA_TOTAL.g <- d_c_GSVA_TOTAL[grepl(paste(glyco_terms, collapse = "|"), rownames(d_c_GSVA_TOTAL), ignore.case = TRUE), ]

d_uc_GSVA_BP <- pull_gsva(gsva = gsva_BP$res, decidetests = gsva_BP$decidetests,
                         coeff = "dysplasia_UC", samples = c(d_samples, uc_samples))
d_uc_GSVA_MF <- pull_gsva(gsva = gsva_MF$res, decidetests = gsva_MF$decidetests,
                         coeff = "dysplasia_UC", samples = c(d_samples, uc_samples))
d_uc_GSVA_CC <- pull_gsva(gsva = gsva_CC$res, decidetests = gsva_CC$decidetests,
                         coeff = "dysplasia_UC", samples = c(d_samples, uc_samples))
d_uc_GSVA_RE <- pull_gsva(gsva = gsva_RE$res, decidetests = gsva_RE$decidetests,
                         coeff = "dysplasia_UC", samples = c(d_samples, uc_samples))
d_uc_GSVA_TOTAL <- rbind(d_uc_GSVA_BP, d_uc_GSVA_MF, d_uc_GSVA_CC, d_uc_GSVA_RE)
d_uc_GSVA_TOTAL.g <- d_uc_GSVA_TOTAL[grepl(paste(glyco_terms, collapse = "|"), rownames(d_uc_GSVA_TOTAL), ignore.case = TRUE), ]
save(d_uc_GSVA_TOTAL, file = "~/projects/tesis/output/integration/d_uc_GSVA_TOTAL.RData")

d_pan_GSVA_BP <- pull_gsva(gsva = gsva_BP$res, decidetests = gsva_BP$decidetests,
                          coeff = "dysplasia_pancolitis", samples = c(d_samples, pan_samples))
d_pan_GSVA_MF <- pull_gsva(gsva = gsva_MF$res, decidetests = gsva_MF$decidetests,
                          coeff = "dysplasia_pancolitis", samples = c(d_samples, pan_samples))
d_pan_GSVA_CC <- pull_gsva(gsva = gsva_CC$res, decidetests = gsva_CC$decidetests,
                          coeff = "dysplasia_pancolitis", samples = c(d_samples, pan_samples))
d_pan_GSVA_RE <- pull_gsva(gsva = gsva_RE$res, decidetests = gsva_RE$decidetests,
                          coeff = "dysplasia_pancolitis", samples = c(d_samples, pan_samples))
d_pan_GSVA_TOTAL <- rbind(d_pan_GSVA_BP, d_pan_GSVA_MF, d_pan_GSVA_CC, d_pan_GSVA_RE)
save(d_pan_GSVA_TOTAL, file = "~/projects/tesis/output/integration/d_pan_GSVA_TOTAL.RData")

pancolitis_uc_GSVA_BP <- pull_gsva(gsva = gsva_BP$res, decidetests = gsva_BP$decidetests,
                                   coeff = "pancolitis_uc", samples = c(pan_samples, uc_samples))

class(d_uc_GSVA_TOTAL)

# GSVA heatmap ----
# Create color palettes
heatmap_colors <- paletteer_c("pals::coolwarm", 100, direction = 1)
phenotype_colors <- paletteer_d("lisa::JosefAlbers", 2)
uc_c_colors <- setNames(phenotype_colors, c("control", "L_UC"))
d_c_colors <- setNames(phenotype_colors, c("dysplasia", "control"))
d_uc_colors <- setNames(phenotype_colors, c("dysplasia", "L_UC"))
pan_uc_colors <- setNames(phenotype_colors, c("pancolitis", "L_UC"))

# Make heatmaps
#make_heatmap <- function(data, phencols, path, filename) {
make_heatmap(data = uc_c_GSVA_TOTAL, phencols = uc_c_colors,
             path = "output/GSE47908/gsva_plots/separate/", filename = "uc_c_gsva.pdf")
make_heatmap(data = d_c_GSVA_TOTAL, phencols = d_c_colors,
             path = "output/GSE47908/gsva_plots/separate/", filename = "d_c_gsva.pdf")
make_heatmap(data = d_uc_GSVA_TOTAL, phencols = d_uc_colors,
             path = "output/GSE47908/gsva_plots/separate/", filename = "d_uc_gsva.pdf")
make_heatmap(data = uc_c_GSVA_TOTAL.g[,2:ncol(uc_c_GSVA_TOTAL.g)], phencols = uc_c_colors,
             path = "output/GSE47908/gsva_plots/separate/", filename = "uc_c_glyco_gsva.pdf")
make_heatmap(data = d_c_GSVA_TOTAL.g[,2:ncol(d_c_GSVA_TOTAL.g)], phencols = d_c_colors,
             path = "output/GSE47908/gsva_plots/separate/", filename = "d_c_glyco_gsva.pdf")
make_heatmap(data = d_uc_GSVA_TOTAL.g, phencols = d_uc_colors,
             path = "output/GSE47908/gsva_plots/separate/", filename = "d_uc_glyco_gsva.pdf")

selected_terms <- c("GOBP_POSITIVE_REGULATION_OF_T_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL",
                    "GOBP_T_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL",
                    "GOBP_ISOTYPE_SWITCHING_TO_IGA_ISOTYPES",
                    "GOBP_REGULATION_OF_MITOTIC_SISTER_CHROMATID_SEGREGATION",
                    "GOBP_ESTABLISHMENT_OF_SISTER_CHROMATID_COHESION",
                    "GOBP_CELLULAR_RESPONSE_TO_INSULIN_LIKE_GROWTH_FACTOR_STIMULUS",
                    "GOBP_MITOTIC_DNA_REPLICATION_CHECKPOINT_SIGNALING",
                    "GOBP_MITOTIC_SPINDLE_MIDZONE_ASSEMBLY",
                    "GOBP_MINUS_END_DIRECTED_ORGANELLE_TRANSPORT_ALONG_MICROTUBULE",
                    "GOBP_PROTEIN_LOCALIZATION_TO_NUCLEAR_BODY",
                    "GOBP_PROTEIN_LOCALIZATION_TO_NUCLEOPLASM",
                    "GOBP_POSITIVE_REGULATION_OF_ESTABLISHMENT_OF_PROTEIN_LOCALIZATION_TO_TELOMERE",
                    "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_HYPOXIA",
                    "GOBP_HEPATOCYTE_APOPTOTIC_PROCESS",
                    "GOBP_NEGATIVE_REGULATION_OF_EXECUTION_PHASE_OF_APOPTOSIS",
                    "GOBP_NEGATIVE_REGULATION_OF_HEPATOCYTE_APOPTOTIC_PROCESS",
                    "GOBP_NEGATIVE_REGULATION_OF_OXIDATIVE_STRESS_INDUCED_NEURON_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
                    "GOBP_NEURON_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_OXIDATIVE_STRESS",
                    "GOBP_POSITIVE_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_ABSENCE_OF_LIGAND",
                    "GOBP_REGULATION_OF_RETROGRADE_TRANSPORT_ENDOSOME_TO_GOLGI",
                    "GOBP_VESICLE_FUSION_WITH_GOLGI_APPARATUS",
                    "GOBP_BMP_SIGNALING_PATHWAY_INVOLVED_IN_HEART_DEVELOPMENT",
                    "GOBP_CELL_MIGRATION_INVOLVED_IN_KIDNEY_DEVELOPMENT",
                    "GOBP_NEURAL_PLATE_MORPHOGENESIS",
                    "GOBP_NEURAL_PLATE_DEVELOPMENT",
                    "GOBP_MITOTIC_DNA_REPLICATION_CHECKPOINT_SIGNALING",
                    "GOBP_NEGATIVE_REGULATION_OF_DNA_DAMAGE_CHECKPOINT",
                    "GOBP_POSITIVE_REGULATION_OF_CELL_CYCLE_CHECKPOINT",
                    "GOBP_POSITIVE_REGULATION_OF_G0_TO_G1_TRANSITION")

my_gsva <- abs(pancolitis_uc_GSVA_BP)
rowsums <- rowSums(my_gsva)
rowsums <- sort(rowsums, decreasing = T)
selected_terms <- names(rowsums[1:25])

# make heatmaps
gsva_plot <- pheatmap(
  #pancolitis_uc_GSVA_BP,
  pancolitis_uc_GSVA_BP[rownames(pancolitis_uc_GSVA_BP) %in% selected_terms,],
         show_colnames = T,
         show_rownames = T,
         #fontsize_col = 8,
         scale = "row",
         treeheight_row = 0,
         #treeheight_col = 18,
         #fontsize_row = 6,
         border_color = "white",
         cutree_rows = 2,
         cutree_cols = 2,
         fontsize = 11,
         color = heatmap_colors,
         annotation_colors = list(phenotype = pan_uc_colors),
         annotation_col = column_to_rownames(study_design[!study_design$phenotype %in% c("control", "dysplasia"),], "sampleID"))

# Save combined .pdf according to comparison in title
pdfs <- list.files(path = "output/GSE47908/gsva_plots/separate/", full.names = TRUE)
qpdf::pdf_combine(input = pdfs[grepl("uc_c_", pdfs)],
                  output = "output/GSE47908/gsva_plots/uc_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("d_c_", pdfs)],
                  output = "output/GSE47908/gsva_plots/d_c_combined.pdf")
qpdf::pdf_combine(input = pdfs[grepl("d_uc_", pdfs)],
                  output = "output/GSE47908/gsva_plots/d_uc_combined.pdf")

# ORA displasia/uc----

lista <- unique(c("SLC26A2", "HMGCS2", "ADH1C", "GUCA2B", "CHP2", "MS4A12", "TRPM6",
           "DUOX2", "PI3", "LCN2", "TIMP1", "MMP9", "CXCL9","DMBT1","S100A9",
           "OLFM4","CCL2","SLC6A14","REG1A","DUOX2","REG1B","DUOXA2","CXCL1","MMP3",
           "TNIP3", "ADH1C","UGT2A3", "ABCG2", "CLDN8", "PCK1","AQP8"))
lista_glyco <- lista[lista %in% glycogenes2]

ora_up_genes <- rownames(d_uc_toptable.f[d_uc_toptable.f$direction == "up",])
ora_down_genes <- rownames(d_uc_toptable.f[d_uc_toptable.f$direction == "down",])

ora_up <- gprofiler2::gost(ora_up_genes, organism = "hsapiens", sources = c("GO:BP"), significant = TRUE)
ora_up.res <- ora_up$result[ora_up$result$term_size < 500,]
gprofiler2::gostplot(ora_up, interactive = T, capped = T)

ora_down <- gprofiler2::gost(ora_down_genes, organism = "hsapiens", sources = c("GO:BP"), significant = TRUE)
ora_down.res <- ora_down$result[ora_down$result$term_size < 500,]
gprofiler2::gostplot(ora_down, interactive = T, capped = T)

# Further dys vs pan or l-uc ----
d_uc_sig <- d_uc_toptable %>% dplyr::filter(adj.P.Val <= pval)
table(d_uc_sig$direction)
d_uc_sig.up <- rownames(d_uc_sig[d_uc_sig$direction == "up",])
d_uc_sig.down <- rownames(d_uc_sig[d_uc_sig$direction == "down",])
d_uc_sig.up.g <- d_uc_sig.up[d_uc_sig.up %in% glycogenes2]
d_uc_sig.down.g <- d_uc_sig.down[d_uc_sig.down %in% glycogenes2]

d_pan_sig <- dys_pan_toptable %>% dplyr::filter(adj.P.Val <= pval)
table(d_pan_sig$direction)
d_pan_sig.up <- rownames(d_pan_sig[d_pan_sig$direction == "up",])
d_pan_sig.down <- rownames(d_pan_sig[d_pan_sig$direction == "down",])
d_pan_sig.up.g <- d_pan_sig.up[d_pan_sig.up %in% glycogenes2]
d_pan_sig.down.g <- d_pan_sig.down[d_pan_sig.down %in% glycogenes2]

gene_lists <- list(LUC_DOWN = d_uc_sig.down, LUC_UP = d_uc_sig.up, PAN_UP = d_pan_sig.up, PAN_DOWN = d_pan_sig.down)
ggvenn(gene_lists, fill_color = c("blue", "red", "red3", "blue3"), text_size = 7, show_percentage = F)

gene_lists <- list(LUC_DOWN = d_uc_sig.down.g, LUC_UP = d_uc_sig.up.g, PAN_UP = d_pan_sig.up.g, PAN_DOWN = d_pan_sig.down.g)
ggvenn(gene_lists, fill_color = c("blue", "red", "red3", "blue3"), text_size = 7, show_percentage = F)

both_down <- intersect(d_uc_sig.down, d_pan_sig.down)
both_up <- intersect(d_uc_sig.up, d_pan_sig.up)
intersect_three <- d_uc_sig.down[d_uc_sig.up %in% d_pan_sig.down]

ora_down <- gprofiler2::gost(both_down, organism = "hsapiens", sources = c("GO"), significant = TRUE)
ora_up <- gprofiler2::gost(both_up, organism = "hsapiens", sources = c("GO"), significant = TRUE)
ora_down.res <- ora_down$result[ora_down$result$term_size < 500,]
ora_up.res <- ora_up$result[ora_up$result$term_size < 500,]
ora_down_filtered <- list(result = ora_down.res, meta = ora_down$meta)
ora_up_filtered <- list(result = ora_up.res, meta = ora_up$meta)
class(ora_down_filtered) <- "gostres"
class(ora_up_filtered) <- "gostres"
gprofiler2::gostplot(ora_down_filtered, interactive = T, capped = T)
gprofiler2::gostplot(ora_up_filtered, interactive = T, capped = T)

# Save Excel sheets ----
#old
openxlsx::write.xlsx(list("all_genes" = uc_c_toptable,
                          "all_significant" = uc_c_toptable.f,
                          "all_glyco" = uc_c_toptable.g,
                          "significant_glyco" = uc_c_toptable.fg,
                          "full_gsea" = uc_c_gsea,
                          "glyco_gsea" = uc_c_gsea.g),
                     file = "output/GSE47908/tables/uc_c_tables.xlsx", rowNames = TRUE)

openxlsx::write.xlsx(list("all_genes" = d_c_toptable,
                          "all_significant" = d_c_toptable.f,
                          "all_glyco" = d_c_toptable.g,
                          "significant_glyco" = d_c_toptable.fg,
                          "full_gsea" = d_c_gsea,
                          "glyco_gsea" = d_c_gsea.g),
                     file = "output/GSE47908/tables/d_c_tables.xlsx", rowNames = TRUE)

openxlsx::write.xlsx(list("all_genes" = d_uc_toptable,
                          "all_significant" = d_uc_toptable.f,
                          "all_glyco" = d_uc_toptable.g,
                          "significant_glyco" = d_uc_toptable.fg,
                          "full_gsea" = d_uc_gsea,
                          "glyco_gsea" = d_uc_gsea.g),
                     file = "output/GSE47908/tables/d_uc_tables.xlsx", rowNames = TRUE)

openxlsx::write.xlsx(list("all_genes" = pan_uc_toptable,
                          "all_significant" = pan_uc_toptable.f,
                          "all_glyco" = pan_uc_toptable.g,
                          "significant_glyco" = d_uc_toptable.fg),
                     file = "pan_uc.xlsx", rowNames = TRUE)

openxlsx::write.xlsx(list("all_genes" = dys_pan_toptable,
                          "full_gsea" = d_pan_gsea),
                     file = "dys_pan_uc.xlsx", rowNames = TRUE)

# Save pancolitis/lUC R object
save(pan_uc_toptable, file = "~/projects/tesis/data/pancolitis/pan_uc_toptable.RData")
save(dys_pan_toptable, file = "~/projects/tesis/data/pancolitis/dys_pan_toptable.RData")

#new
openxlsx::write.xlsx(list("all_genes" = uc_c_toptable[,!colnames(uc_c_toptable) %in% "rank"],
                          "all_glyco" = uc_c_toptable.g[,!colnames(uc_c_toptable.g) %in% "rank"],
                          "full_gsea" = uc_c_gsea,
                          "glyco_gsea" = uc_c_gsea.g),
                     file = "~/projects/tesis/final_outputs/uc_control/GSE47908.xlsx", rowNames = TRUE)

openxlsx::write.xlsx(list("all_genes" = d_uc_toptable[,!colnames(d_uc_toptable) %in% "rank"],
                          "all_glyco" = d_uc_toptable.g[,!colnames(d_uc_toptable.g) %in% "rank"],
                          "full_gsea" = d_uc_gsea,
                          "glyco_gsea" = d_uc_gsea.g),
                     file = "~/projects/tesis/final_outputs/displasia_vs_uc/dis_vs_luc.xlsx", rowNames = TRUE)

openxlsx::write.xlsx(list("all_genes" = dys_pan_toptable[,!colnames(dys_pan_toptable) %in% "rank"],
                          "all_glyco" = dys_pan_toptable.g[,!colnames(dys_pan_toptable.g) %in% "rank"],
                          "full_gsea" = d_pan_gsea,
                          "glyco_gsea" = d_pan_gsea.g),
                     file = "~/projects/tesis/final_outputs/displasia_vs_uc/dis_vs_pan.xlsx", rowNames = TRUE)
