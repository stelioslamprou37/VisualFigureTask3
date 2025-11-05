
# Set data directory

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
})

cat("============================================================================\n")
cat("Forest Plot Analysis - LOCAL VERSION (Windows)\n")
cat("Paper: Discovery Oncology 2025\n")
cat("Dataset: GSE84437 (433 samples)\n")
cat("============================================================================\n\n")
data_file <- "C:/Users/steli/Desktop/Figure Reproduction Project - Copy/input_data/GSE84437_series_matrix.txt"

if(!file.exists(data_file)) {
  stop("ERROR: File not found at ", data_file)
}

gse84437 <- getGEO(filename = data_file, GSEMatrix = TRUE)

expr_data <- exprs(gse84437)
pdata <- pData(gse84437)
fdata <- fData(gse84437)

cat("  ✓ Loaded:", nrow(expr_data), "probes x", ncol(expr_data), "samples\n\n")

# ============================================================================
# Step 2: Extract survival data (CORRECTED COLUMNS)
# ============================================================================

cat("Step 2: Extracting survival data...\n")

# CORRECT COLUMNS:
survival_status <- as.numeric(pdata[["death:ch1"]])
survival_time <- as.numeric(pdata[["duration overall survival:ch1"]])

cat("  Samples:", length(survival_time), "\n")
cat("  Events (death=1):", sum(survival_status == 1, na.rm = TRUE), "\n")
cat("  Censored (death=0):", sum(survival_status == 0, na.rm = TRUE), "\n\n")

# ============================================================================
# Step 2: Extract survival data
# ============================================================================

cat("Step 2: Extracting survival data...\n")

survival_status <- as.numeric(pdata[["death:ch1"]])
survival_time <- as.numeric(pdata[["duration overall survival:ch1"]])

cat("  Samples:", length(survival_time), "\n")
cat("  Events:", sum(survival_status == 1, na.rm = TRUE), "\n\n")

# ============================================================================
# Step 3: Preprocessing
# ============================================================================

cat("Step 3: Preprocessing...\n")

# Filter low-variance probes
probe_vars <- apply(expr_data, 1, var, na.rm = TRUE)
expr_filtered <- expr_data[probe_vars > 0, ]

cat("  Probes after filtering:", nrow(expr_filtered), "\n")

if(max(expr_filtered, na.rm = TRUE) > 100) {
  expr_filtered <- log2(expr_filtered + 1)
}

expr_normalized <- normalize.quantiles(as.matrix(expr_filtered))
rownames(expr_normalized) <- rownames(expr_filtered)
colnames(expr_normalized) <- colnames(expr_filtered)

cat("  ✓ Complete\n\n")

# ============================================================================
# Step 4: Handle missing values
# ============================================================================

cat("Step 4: Handling missing values...\n")

n_missing <- sum(is.na(expr_normalized))
if(n_missing > 0) {
  for(i in 1:nrow(expr_normalized)) {
    row_mean <- mean(expr_normalized[i, ], na.rm = TRUE)
    expr_normalized[i, is.na(expr_normalized[i, ])] <- row_mean
  }
  cat("  Imputed", n_missing, "values\n\n")
} else {
  cat("  No missing values\n\n")
}

# ============================================================================
# Step 5: Cox regression
# ============================================================================

cat("Step 5: Cox regression...\n")

valid_samples <- !is.na(survival_time) & !is.na(survival_status) & survival_time > 0

survival_time_clean <- survival_time[valid_samples]
survival_status_clean <- survival_status[valid_samples]
expr_clean <- expr_normalized[, valid_samples]

cat("  Samples:", ncol(expr_clean), "\n")
cat("  Probes:", nrow(expr_clean), "\n")

# Top probes by variance
probe_vars_clean <- apply(expr_clean, 1, var, na.rm = TRUE)
top_probe_idx <- order(probe_vars_clean, decreasing = TRUE)[1:min(500, nrow(expr_clean))]

cox_results <- data.frame()

for(i in top_probe_idx) {
  probe_expr <- expr_clean[i, ]
  
  cox_data <- data.frame(
    time = survival_time_clean,
    status = survival_status_clean,
    expr = probe_expr
  )
  
  fit <- tryCatch(
    coxph(Surv(time, status) ~ expr, data = cox_data),
    error = function(e) NULL
  )
  
  if(!is.null(fit)) {
    coef_val <- coef(fit)[1]
    hr <- exp(coef_val)
    se_coef <- sqrt(vcov(fit)[1, 1])
    ci_lower <- exp(coef_val - 1.96 * se_coef)
    ci_upper <- exp(coef_val + 1.96 * se_coef)
    pval <- summary(fit)$coefficients[1, 5]
    
    cox_results <- rbind(cox_results, data.frame(
      probe = rownames(expr_clean)[i],
      HR = hr,
      CI_lower = ci_lower,
      CI_upper = ci_upper,
      p_value = pval
    ))
  }
}

cox_results <- cox_results[order(cox_results$p_value), ]

cat("  ✓ Complete - ", sum(cox_results$p_value < 0.05), "significant\n\n")

# ============================================================================
# Step 6: Forest plot
# ============================================================================

cat("Step 6: Creating forest plot...\n")

top_probes <- cox_results[cox_results$p_value < 0.05, ][1:10, ]
top_probes <- top_probes[!is.na(top_probes$HR), ]

if(nrow(top_probes) > 0) {
  
  probe_names <- substr(top_probes$probe, 1, 12)
  
  forest_data <- data.frame(
    Probe = probe_names,
    HR = top_probes$HR,
    CI_lower = top_probes$CI_lower,
    CI_upper = top_probes$CI_upper,
    P_value = top_probes$p_value,
    Significant = ifelse(top_probes$p_value < 0.05, "Yes", "No")
  )
  
  png("output_forest_plot.png", width = 10, height = 6, units = "in", res = 300)
  
  forest_data$Probe <- factor(forest_data$Probe, levels = rev(forest_data$Probe))
  
  p <- ggplot(forest_data, aes(x = HR, y = Probe, color = Significant)) +
    geom_point(size = 4) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, size = 1) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", size = 1) +
    scale_x_log10(breaks = c(0.5, 1, 2, 4, 8)) +
    scale_color_manual(values = c("Yes" = "#E41A1C", "No" = "#377EB8")) +
    labs(title = "Forest Plot: Survival-Associated Probes",
         x = "Hazard Ratio (95% CI)", y = "Probe ID", color = "p < 0.05") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  
  print(p)
  dev.off()
  
  cat("  ✓ Saved output_forest_plot.png\n\n")
  
  write.csv(cox_results, "output_cox_regression_results.csv", row.names = FALSE)
  write.csv(forest_data, "output_forest_plot_data.csv", row.names = FALSE)
  
  cat("  ✓ Saved CSV files\n\n")
}

cat("============================================================================\n")
cat("COMPLETE!\n")
cat("============================================================================\n")

# ============================================================================
# Summary
# ============================================================================

cat("============================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("============================================================================\n\n")

cat("Results:\n")
cat("  Samples analyzed:", ncol(expr_clean), "\n")
cat("  Total deaths:", sum(survival_status_clean), "\n")
cat("  Genes analyzed:", nrow(expr_clean), "\n")
cat("  Significant genes (p<0.05):", sum(cox_results$p_value < 0.05), "\n\n")

cat("Output files:\n")
cat("  1. output_forest_plot.png\n")
cat("  2. output_cox_regression_results.csv\n")
cat("  3. output_forest_plot_data.csv\n\n")

cat("============================================================================\n")
