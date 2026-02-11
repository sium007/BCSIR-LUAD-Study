# Load required libraries
library(survival)
library(survminer)
library(ggplot2)

# 1. Load data
clinical_data <- read.csv('clinical_data.csv')
rnaseq_data <- read.csv('rnaseq_data.csv')

# Merge clinical and RNAseq data
merged_data <- merge(clinical_data, rnaseq_data, by = 'sample_id')

# 2. Define gene of interest
gene_of_interest <- 'FAM83A-AS2' # Change to your gene, e.g., 'AC012213.1'

# 3. Stratify patients (High vs Low Expression)
# Using median cutoff
expression_cutoff <- median(merged_data[[gene_of_interest]])
merged_data$group <- ifelse(merged_data[[gene_of_interest]] > expression_cutoff, 
                            "High Expression", "Low Expression")

# 4. Fit Kaplan-Meier curves
fit <- survfit(Surv(time_to_event, event_of_interest) ~ group, data = merged_data)

# 5. Plot using ggsurvplot (Better visualization than base ggplot)
km_plot <- ggsurvplot(
  fit,
  data = merged_data,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.col = "strata",
  linetype = "strata",
  surv.median.line = "hv",
  ggtheme = theme_classic(),
  palette = c("#E91E63", "#18A1CD"),
  legend.labs = c("High Expression", "Low Expression"),
  xlab = "Time (days)",
  ylab = "Survival Probability"
)

print(km_plot)
# ggsave("KM_Plot.png", print(km_plot), width = 8, height = 6)

message("Survival analysis completed.")