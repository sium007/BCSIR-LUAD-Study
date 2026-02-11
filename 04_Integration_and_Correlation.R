# Load required package
library(ggstatsplot)

# 1. Read input data
# Format: sample_id | methylation | expression
luad_df <- read.csv("LUAD_methylation_expression.csv")

# 2. Spearman correlation test (Base R)
spearman_test <- cor.test(
  luad_df$methylation,
  luad_df$expression,
  method = "spearman",
  conf.level = 0.95
)

# Print statistics
print(paste("Spearman's rho:", spearman_test$estimate))
print(paste("P-value:", spearman_test$p.value))

# 3. Visualization using ggstatsplot
corr_plot <- ggstatsplot::ggscatterstats(
  data = luad_df,
  x = methylation,
  y = expression,
  type = "nonparametric",        # Spearman rank correlation
  conf.level = 0.95,
  marginal = TRUE,
  xlab = "DMC methylation level",
  ylab = "DElncRNA expression",
  title = "Correlation between DMC methylation and DElncRNA expression",
  subtitle = paste0("Spearman’s rank correlation (LUAD, n = ", nrow(luad_df), ")")
)

print(corr_plot)
# ggsave("Correlation_Plot.png", corr_plot, width = 8, height = 8)

message("Correlation analysis completed.")