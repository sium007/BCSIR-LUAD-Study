# NOTE:
# The ChAMP workflow was used here for IDAT import, preprocessing, normalization,
# and differential methylation analysis.
# The package `IlluminaHumanMethylation450kanno.ilmn12.hg19` is loaded for
# compatibility with the ChAMP preprocessing framework.
# However, the final CpG coordinate mapping and promoter assignment reported
# in the manuscript were performed separately using the hg38-based annotation
# resource `HM450_hg38_manifest_gencode_v36`.

# Load required library
if (!requireNamespace("Hmisc", quietly = TRUE)) install.packages("Hmisc")
library(Hmisc)

# 1. Read input data
data <- read.csv("correlation_network.csv")

# 2. Calculate Spearman correlation matrix
cor_matrix <- rcorr(as.matrix(data), type = "spearman")

# Extract correlation coefficients and p-values
correlation <- cor_matrix$r
p_values <- cor_matrix$P

# 3. Flatten matrix to pairwise list
results <- data.frame(Variable1 = character(), Variable2 = character(),
                      Correlation = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

for (i in 1:(ncol(correlation) - 1)) {
  for (j in (i + 1):ncol(correlation)) {
    results <- rbind(results, data.frame(
      Variable1 = colnames(correlation)[i],
      Variable2 = colnames(correlation)[j],
      Correlation = correlation[i, j],
      p_value = p_values[i, j]
    ))
  }
}

# 4. Save results
write.csv(results, "correlation_results.csv", row.names = FALSE)


message("Network correlation analysis completed.")
