# Define functions to use

# make_toptable
make_toptable <- function(coefficient, comparison) {
  toptable <- topTable(lm_stats, adjust = "BH", coef = coefficient, number = 90000) |>
    dplyr::select(logFC, P.Value, adj.P.Val) |>
    dplyr::mutate(rank = ifelse(logFC > 0, -log10(adj.P.Val), log10(adj.P.Val))) |>
    dplyr::mutate(direction = ifelse(logFC > 0, "up", "down")) |>
    dplyr::arrange(desc(rank))
  toptable$comparison <- comparison
  return(toptable)
}

# volcano
volcano <- function(table, p = 0.05, lfc = 1.5, title = "Volcano plot", path, filename, show = FALSE) {
  
  table <- as_tibble(table, rownames = "geneID")
  degs <- table %>% dplyr::filter(adj.P.Val <= p & abs(logFC) >= lfc)
  degs <- degs$geneID
  
  significant <- table %>% dplyr::filter(adj.P.Val <= pval)
  significant_up <- sum(significant$direction == "up")
  significant_down <- sum(significant$direction == "down")
  
  plot <- ggplot(table) +
    aes(y = -log10(adj.P.Val), x = logFC, text = paste("Symbol:", geneID), customdata = geneID) +
    geom_hline(yintercept = -log10(p), linetype = "longdash", colour = "grey", linewidth = 1) +
    geom_vline(xintercept = lfc, linetype = "longdash", colour = "#BE684D", linewidth = 1) +
    geom_vline(xintercept = -lfc, linetype = "longdash", colour = "#2C467A", linewidth = 1) +
    geom_point(size = 2, color = ifelse(!table$geneID %in% degs, "grey", ifelse(table$logFC > 0, "red3", "blue3"))) +
    geom_text_repel(data = subset(table, geneID %in% degs), aes(label = geneID), vjust = -1) +
    labs(title = title, subtitle = paste0("p-value = ", p, ", logFC = ", lfc),
         caption = paste0("Regardless of logFC, significant DEGs: ", significant_up,
                          " (upregulated), ", significant_down, " (downregulated).")) +
    theme_bw()
  
  if (show == TRUE) {
    print(plot)
  }
  ggsave(filename = filename, plot = plot, device = "pdf", path = path, width = 15, height = 10, dpi = 300)
}

# run_gsea
run_gsea <- function(ranked_vector, comparison, c5, reactome, min = 5, max = 500) {
  
  # Run GSEA
  c5_result <- (GSEA(ranked_vector, TERM2GENE = c5, verbose = TRUE, eps = 0,
                     minGSSize = min, maxGSSize = max))@result
  
  reactome_result <- (GSEA(ranked_vector, TERM2GENE = reactome, verbose = TRUE, eps = 0,
                           minGSSize = min, maxGSSize = max))@result
  
  result <- rbind(c5_result, reactome_result)
  result$comparison <- comparison
  result <- result %>% dplyr::arrange(p.adjust, desc(abs(NES)))
  result$geneset_size <- paste(min, "-", max)
  
  return(result)
}

# gsea_bubble
gsea_bubble <- function(input, comp, n = 35, path, show = FALSE) {
  for (i in c("GOBP", "GOCC", "GOMF", "HP", "REACTOME")) {

    data <- input %>%
      dplyr::filter(grepl(i, ID)) %>%
      dplyr::mutate(NES_category = ifelse(NES > 0, "Positive", "Negative"))
    
    if (n > nrow(data)) {
      n <- nrow(data)
    }
    
    data <- data[1:n,]
    filename <- paste0(comp, "_GSEA_", i, ".pdf")
    height <- nrow(data) * 0.3
    width <- max(nchar(data$ID)) * 0.1 + 2
    if (height < 7) {
      height <- 7
    }

    plot <- ggplot(data, aes(x = comparison, y = reorder(ID, -log10(p.adjust)))) +
      geom_point(aes(size = setSize, color = NES_category, alpha = -log10(p.adjust))) +
      scale_color_manual(values = c("Positive" = "red", "Negative" = "blue")) +
      labs(title = paste(i, data$comparison), 
           subtitle = paste("Geneset size range:", data$geneset_size),
           x = element_blank(),
           y = "Pathway",
           size = "Gene Set Size",
           color = "NES Direction",
           alpha = "-log10(p.adjust)") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 1, size = 14, face = "bold"))

    if (show == TRUE) {
      print(plot)
    }
    ggsave(filename = filename, plot = plot, device = "pdf",
           height = height, width = width, path = path, dpi = 300)
  }
}

# gsea_bubble2
gsea_bubble2 <- function(input, n = 50, set = "GOBP") {
  data <- uc_c_gsea[1:n,] %>%
    dplyr::filter(grepl(set, ID)) %>%
    dplyr::mutate(NES_category = ifelse(NES > 0, "Positive", "Negative"))
  if (n > nrow(data)) {
    n <- nrow(data) 
  }
  plot <- ggplot(data, aes(x = NES, y = -log10(p.adjust), color = NES)) +
    geom_point(aes(size = setSize, color = NES_category, alpha = -log10(p.adjust))) +
    scale_color_manual(values = c("Positive" = "red", "Negative" = "blue")) +
    geom_text_repel(aes(label = ID), size = 2,colour = "black") +
    guides(color = "none", alpha = "none") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 1, size = 14, face = "bold"))
  print(plot)
}

# run_gsva
run_gsva <- function(exp, go, min = 5, max = 500, dmatrix, cmatrix, p = 0.05, lfc = 1.5) {
  gsva <- GSVA::gsva(param = gsvaParam(exp, go, maxDiff = TRUE, minSize = min, maxSize = max))
  
  fit <- lmFit(gsva, dmatrix)
  fit2 <- contrasts.fit(contrasts = cmatrix, fit = fit)
  stats <- eBayes(fit2)
  
  # Identify the differentially enriched gene sets
  gsva_decidetests <- decideTests(stats, method = "global", adjust.method = "BH", p.value = p, lfc = lfc)
  
  return(list(res = gsva, stats = stats, decidetests = gsva_decidetests))
}

# pull enriched gsva
pull_gsva <- function (gsva, decidetests, coeff, samples) {
  enriched <- gsva[decidetests[, coeff] != 0,]
  enriched <- enriched[,samples]
  return(enriched)
}

# make heatmaps
make_heatmap <- function(data, phencols, path, filename) {
  if (nrow(data) > 550) {
    data <- data[1:550,]
  }
  heatmap <- pheatmap(data,
                      show_colnames = T,
                      fontsize_col = 6,
                      scale = "row",
                      treeheight_row = 0,
                      fontsize_row = 6,
                      color = heatmap_colors,
                      annotation_colors = list(phenotype = phencols),
                      annotation_col = column_to_rownames(study_design, "sampleID"))
  height <- nrow(data) * 0.3
  if (height > 49) {
    height <- 49
  }
  width <- max(nchar(rownames(data))) * 0.1 + 2
  ggsave(filename = filename, plot = heatmap, device = "pdf", height = height,
         width = width, path = path, dpi = 300)
}

# breaks for heatmaps
make_heatmap_b <- function(data, show_rows = T, manual_break = 0, krows = 1, p = F, pmat) {
  
  cols <- colorRampPalette(c("blue", "white", "red"))(100)
  
  if(manual_break == 0){
    max_abs <- max(abs(min(data)), abs(max(data)))
    breaks <- seq(-max_abs, max_abs, length.out = 101)
  } else {
    breaks <- seq(-manual_break, manual_break, length.out = 101)
  }
  
  if(!p){
    heatmap <- pheatmap(data, color = cols, treeheight_col = 0, cutree_rows = krows,
                      border_color = "black",
                      treeheight_row = 0, breaks = breaks, show_rownames = show_rows)
  } else {
    heatmap <- pheatmap(data, color = cols, treeheight_col = 0, cutree_rows = krows,
                        border_color = "black",display_numbers = pmat,
                        treeheight_row = 0, breaks = breaks, show_rownames = show_rows)
  }
  return(heatmap)
}
