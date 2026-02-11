# Load required libraries
suppressPackageStartupMessages({
  library(ChAMP)
  library(minfi)
  library(limma)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})

# 1. Load IDAT files
# Directory should contain paired Red/Green IDAT files and a sample sheet
myLoad <- champ.load(
  directory = "IDAT_files/",
  arraytype = "450K"
)

# 2. Probe filtering (QC)
# - Removes low-quality probes (detection p-value > 0.01)
# - Removes probes with low bead count (< 3 in > 5% samples)
# - Filters non-CpG, SNP-associated, multi-hit, and sex chromosome probes
myFilt <- champ.filter(
  beta       = myLoad$beta,
  M          = myLoad$M,
  pd         = myLoad$pd,
  detP       = myLoad$detP,
  beadcount  = myLoad$beadcount,
  detPcut        = 0.01,
  beadCutoff     = 0.05,
  filterNoCG     = TRUE,
  filterSNPs     = TRUE,
  filterMultiHit = TRUE,
  filterXY       = TRUE
)

# 3. Normalization (BMIQ)
# Corrects for Type I/II probe bias
myNorm <- champ.norm(
  beta = myFilt$beta,
  method = "BMIQ",
  arraytype = "450K"
)

# 4. Identify Differentially Methylated Probes (DMPs)
# Using limma-based linear modeling
myDMP <- champ.DMP(
  beta  = myNorm,
  pheno = myFilt$pd$Sample_Group,
  method = "limma",
  adjPVal = 0.05
)

# Extract tumor vs normal results
DMP_Tumor_vs_Normal <- myDMP$Tumor_to_Normal

# 5. Identify Differentially Methylated Regions (DMRs) - Optional
# Using Bumphunter algorithm
myDMR <- champ.DMR(
  beta = myNorm,
  pheno = myFilt$pd$Sample_Group,
  method = "Bumphunter",
  arraytype = "450K",
  minProbes = 7,
  maxGap = 300,
  cutoff = 0.454,
  B = 1000
)

DMR_Tumor_vs_Normal <- myDMR$Tumor_to_Normal

# 6. Save results
write.csv(DMP_Tumor_vs_Normal, "DMP_Tumor_vs_Normal.csv", row.names = TRUE)
write.csv(DMR_Tumor_vs_Normal, "DMR_Bumphunter.csv", row.names = FALSE)

message("Differential methylation analysis completed successfully.")