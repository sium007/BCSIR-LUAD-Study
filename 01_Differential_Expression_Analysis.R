# Load required library
library(DESeq2)

# 1. Load input data
counts_data <- read.csv("counts_data.csv", row.names = 1)
colData <- read.csv("colData.csv", row.names = 1)

# Ensure sample names match exactly
stopifnot(all(colnames(counts_data) == rownames(colData)))

# 2. Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = colData,
  design = ~ Phenotype
)

# 3. Filter low-count genes

dds <- dds[rowSums(counts(dds)) >= 10, ]

# 4. Set reference level
dds$Phenotype <- relevel(dds$Phenotype, ref = "normal")

# 5. Differential expression analysis
dds <- DESeq(dds)

# 6. Extract DE results (BH-adjusted p-values)
res <- results(
  dds,
  contrast = c("Phenotype", "Tumor", "normal"),
  alpha = 0.05
)

# Save full DESeq2 results
write.csv(as.data.frame(res), "DESeq2_full_results.csv")

# 7. Variance Stabilizing Transformation (VST) for downstream analysis
vsd <- vst(dds, blind = FALSE)

# Extract VST-normalized expression matrix
vsd_mat <- assay(vsd)

# Save VST matrix
write.csv(vsd_mat, "DESeq2_VST_expression_matrix.csv")

message("Differential expression analysis completed successfully.")