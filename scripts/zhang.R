# UC / C - Zhang et. al: 10.3389/fimmu.2023.1095098
library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggrepel)
source("scripts/functions.R")

metadata <- t(read.xlsx(xlsxFile = "data/processed/zhang/uc_c_zhang.xlsx", sheet = "summary_information", rowNames = T))
data <- read.xlsx(xlsxFile = "data/processed/zhang/uc_c_zhang.xlsx", sheet = "deg_table", rowNames = T)[,c(1,5)]

upregulated <- rownames(data)[data$logFC > 0]
downregulated <- rownames(data)[data$logFC < 0]

v <- ggplot(data) +
  aes(y = -log10(adj.P.Val), x = logFC, text = "row.names") +
  geom_hline(yintercept = -log10(0.05), linetype = "longdash", colour = "grey", linewidth = 1) +
  geom_vline(xintercept = 2, linetype = "longdash", colour = "#BE684D", linewidth = 1) +
  geom_vline(xintercept = -2, linetype = "longdash", colour = "#2C467A", linewidth = 1) +
  geom_point(size = 2, color = ifelse(data$logFC > 0, "red", "blue")) +
  geom_text_repel(data = data, aes(label = rownames(data)), vjust = -1) +
  labs(title = "UC/C DEGs", subtitle = paste0("p-value = ", 0.05, ", logFC = ", 2),
       caption = paste0("Regardless of logFC, significant DEGs: ", length(upregulated),
                        " (upregulated), ", length(downregulated), " (downregulated).")) +
  theme_bw()

ggsave(filename = "volcano.pdf", plot = v, device = "pdf", path = "output/zhang/", width = 16, height = 10, dpi = 300)

load("data/rdata/glycogenes.RData")
glyco_data <- data[rownames(data) %in% glycogenes,]

g <- ggplot(glyco_data) +
  aes(y = -log10(adj.P.Val), x = logFC, text = "row.names") +
  geom_hline(yintercept = -log10(0.05), linetype = "longdash", colour = "grey", linewidth = 1) +
  geom_vline(xintercept = 2, linetype = "longdash", colour = "#BE684D", linewidth = 1) +
  geom_vline(xintercept = -2, linetype = "longdash", colour = "#2C467A", linewidth = 1) +
  geom_point(size = 2, color = ifelse(glyco_data$logFC > 0, "red", "blue")) +
  geom_text_repel(data = glyco_data, aes(label = rownames(glyco_data)), vjust = -1) +
  labs(title = "UC/C glycoDEGs", subtitle = paste0("p-value = ", 0.05, ", logFC = ", 2),
       caption = paste0("Regardless of logFC, significant DEGs: ", length(upregulated),
                        " (upregulated), ", length(downregulated), " (downregulated).")) +
  theme_bw()

ggsave(filename = "glyco_volcano.pdf", plot = g, device = "pdf", path = "output/zhang/", width = 16, height = 10, dpi = 300)

openxlsx::write.xlsx(list("all" = data , "glycogenes" = glyco_data),
                     file = "output/zhang/uc_c_tables.xlsx", rowNames = TRUE)
