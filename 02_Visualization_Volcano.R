# Load libraries
library(ggplot2)
library(dplyr)

# Load data 
# res_df <- read.csv("DESeq2_full_results.csv") 

# 1. Define significance categories
res_df <- res_df %>%
  mutate(
    significance = case_when(
      log2FC > 1.5 & padj < 0.05 ~ "Up-regulated",
      log2FC < -1.5 & padj < 0.05 ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

# 2. Create labels for legend
sig_counts <- res_df %>% count(significance)
sig_labels <- paste0(sig_counts$significance, " (", format(sig_counts$n, big.mark=","), ")")
names(sig_labels) <- sig_counts$significance

# 3. Generate Volcano Plot
p1 <- ggplot(res_df, aes(x = log2FC, y = -log10(padj), color = significance)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(
    values = c("Down-regulated" = "#18A1CD", "Not significant" = "black", "Up-regulated" = "#E91E63"),
    labels = sig_labels
  ) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  theme_classic(base_size = 14) +
  labs(
    x = "log2(Fold Change)",
    y = "-log10(Adj. p-value)",
    color = "Significance"
  ) +
  annotate("text", x = min(res_df$log2FC), y = max(-log10(res_df$padj), na.rm=TRUE),
           label = "(a)", size = 6, fontface = "bold")

print(p1)
# ggsave("Volcano_Plot.png", p1, width = 8, height = 6)

message("Volcano plot generated.")