# Load required libraries
library(glmnet)
library(survival)
library(survivalROC)
library(survminer)

# 1. Read and clean data
clinical_new <- read.csv("Survival.csv")
clinical_new <- clinical_new[complete.cases(clinical_new) & clinical_new$survival > 0, ]

# 2. Univariate Cox Regression
covariates <- c("FAM83A.AS2", "AC012213.1")
univ_models <- lapply(covariates, function(v) {
  coxph(as.formula(paste("Surv(survival, status) ~", v)), data = clinical_new)
})

# Extract results
univ_results <- lapply(univ_models, function(model) {
  s <- summary(model)
  beta <- signif(s$coef[, "coef"], 2)
  HR <- signif(s$coef[, "exp(coef)"], 2)
  CI.low <- signif(s$conf.int[, "lower .95"], 2)
  CI.high <- signif(s$conf.int[, "upper .95"], 2)
  pval <- signif(s$coef[, "Pr(>|z|)"], 2)
  c(beta, paste0(HR, " (", CI.low, "-", CI.high, ")"), pval)
})

univ_table <- t(as.data.frame(univ_results))
colnames(univ_table) <- c("Beta", "HR (95% CI)", "P value")
rownames(univ_table) <- covariates
print(univ_table)

# 3. Prepare data for LASSO
X <- as.matrix(clinical_new[, !(colnames(clinical_new) %in% c("survival", "status"))])
surv_obj <- Surv(time = clinical_new$survival, event = clinical_new$status)

# 4. CV-LASSO Cox Regression
cvfit <- cv.glmnet(X, surv_obj, family = "cox", alpha = 1, nfolds = 10, maxit = 1000)

# Plot CV curve
plot(cvfit)
title("LASSO Cross-Validation", line = 2.5)

# 5. Risk Score Calculation
risk_score <- as.vector(predict(cvfit, newx = X, s = "lambda.min", type = "link"))

# 6. Time-dependent ROC Analysis
# Determine prediction time (e.g., median survival)
predict_time <- median(clinical_new$survival) 
cat("Prediction Time (Median Survival):", predict_time, "\n")

roc_obj <- survivalROC(
  Stime = clinical_new$survival,
  status = clinical_new$status,
  marker = risk_score,
  predict.time = predict_time,
  method = "NNE",
  lambda = 0.1
)

# Plot ROC Curve
plot(roc_obj$FP, roc_obj$TP, type = "b", col = "green3", lwd = 2, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "False Positive Rate", ylab = "True Positive Rate", main = "Time-dependent ROC")
abline(0, 1, col = "red", lwd = 1, lty = 2)
text(0.6, 0.2, paste0("AUC = ", round(roc_obj$AUC, 2)), cex = 1.2, col = "blue")

message("Prognostic modeling completed.")