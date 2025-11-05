#!/usr/bin/env Rscript

# ============================================================================
# Forest Plot Analysis: Cox Regression for Gastric Cancer
# Paper: Discovery Oncology 2025 - DOI: 10.1007/s12672-025-01907-7
# Dataset: GSE84437 (483 gastric cancer samples)
# Analysis: Cox proportional hazards regression with forest plot
# ============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
})

cat("============================================================================\n")
cat("Forest Plot Analysis: Gastric Cancer Cox Regression\n")
cat("Paper: Discovery Oncology 2025\n")
cat("Dataset: GSE84437 (483 samples)\n")
cat("============================================================================\n\n")

# ============================================================================
# Step 1: Load GSE84437 expression data
# ============================================================================

cat("Step 1: Loading GSE84437 from GEO...\n")

data_dir <- "./input_data"
if(!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

max_retries <- 3
retry_count <- 0
gse84437 <- NULL

while(retry_count < max_retries && is.null(gse84437)) {
  retry_count <- retry_count + 1

  cat("  Attempt", retry_count, "of", max_retries, "...\n")

  tryCatch({
    options(timeout = 300)
    gse84437 <- getGEO("GSE84437", GSEMatrix = TRUE, destdir = data_dir,
                       AnnotGPL = FALSE, getGPL = FALSE)

    if(length(gse84437) > 0) {
      gse84437 <- gse84437[[1]]
      cat("  ✓ Successfully downloaded\n\n")
    }

  }, error = function(e) {
    cat("  Error:", conditionMessage(e), "\n")
    if(retry_count < max_retries) {
      cat("  Retrying in 5 seconds...\n")
      Sys.sleep(5)
    }
  })
}

if(is.null(gse84437)) {
  cat("\n❌ ERROR: Could not download GSE84437\n")
  quit(status = 1)
}

expr_data <- exprs(gse84437)
pdata <- pData(gse84437)

cat("Data loaded:"), nrow(expr_data), "probes x", ncol(expr_data), "samples\n\n")

# ============================================================================
# Step 2: Extract survival and clinical data
# ============================================================================

cat("Step 2: Extracting survival data...\n")

survival_status <- as.numeric(pdata[["death:ch1"]])
survival_time <- as.numeric(pdata[["duration overall survival:ch1"]])

cat("  Samples:", length(survival_time), "\n")
cat("  Events (death=1):", sum(survival_status == 1, na.rm = TRUE), "\n")
cat("  Censored (death=0):", sum(survival_status == 0, na.rm = TRUE), "\n")
cat("  Median survival:", round(median(survival_time, na.rm = TRUE), 1), "months\n\n")

# ============================================================================
# Step 3: Preprocessing - Filter and normalize
# ============================================================================

cat("Step 3: Preprocessing...\n")

# Filter low-variance probes
probe_vars <- apply(expr_data, 1, var, na.rm = TRUE)
expr_filtered <- expr_data[probe_vars > 0, ]

cat("  Probes after filtering:", nrow(expr_filtered), "\n")

# Log transformation (if needed)
if(max(expr_filtered, na.rm = TRUE) > 100) {
  expr_filtered <- log2(expr_filtered + 1)
}

# Quantile normalization
expr_normalized <- normalize.quantiles(as.matrix(expr_filtered))
rownames(expr_normalized) <- rownames(expr_filtered)
colnames(expr_normalized) <- colnames(expr_filtered)

cat("  ✓ Normalization complete\n\n")

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
# Step 5: Univariate Cox regression
# ============================================================================

cat("Step 5: Cox regression analysis...\n")

# Filter for complete cases
valid_samples <- !is.na(survival_time) & !is.na(survival_status) & survival_time > 0

survival_time_clean <- survival_time[valid_samples]
survival_status_clean <- survival_status[valid_samples]
expr_clean <- expr_normalized[, valid_samples]

cat("  Samples for analysis:", ncol(expr_clean), "\n")
cat("  Probes analyzed:", nrow(expr_clean), "\n")

# Select top probes by variance
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
      p_value = pval,
      coef = coef_val
    ))
  }
}

cox_results <- cox_results[order(cox_results$p_value), ]

cat("  ✓ Cox regression complete\n")
cat("  Significant probes (p<0.05):", sum(cox_results$p_value < 0.05), "\n\n")

# ============================================================================
# Step 6: Select top probes for forest plot
# ============================================================================

cat("Step 6: Selecting top probes for forest plot...\n")

top_probes <- cox_results[cox_results$p_value < 0.05, ][1:10, ]
top_probes <- top_probes[!is.na(top_probes$HR), ]

cat("  Selected", nrow(top_probes), "top probes\n\n")

# ============================================================================
# Step 7: Generate forest plot
# ============================================================================

if(nrow(top_probes) > 0) {

  cat("Step 7: Generating forest plot...\n")

  # Shorten probe names for readability
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
    geom_point(size = 4, position = position_dodge(width = 0.5)) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), 
                    height = 0.2, position = position_dodge(width = 0.5), size = 1) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", size = 1) +
    scale_x_log10(breaks = c(0.5, 1, 2, 4, 8)) +
    scale_color_manual(values = c("Yes" = "#E41A1C", "No" = "#377EB8")) +
    labs(
      title = "Forest Plot: Survival-Associated Probes",
      subtitle = "Cox Regression Analysis (GSE84437 - 483 Samples)",
      x = "Hazard Ratio (95% CI) [log scale]",
      y = "Illumina Probe ID",
      color = "p < 0.05"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      panel.grid.major.x = element_line(color = "gray90")
    )

  print(p)
  dev.off()

  cat("  ✓ Saved output_forest_plot.png\n\n")

  # ============================================================================
  # Step 8: Export results
  # ============================================================================

  cat("Step 8: Exporting results...\n")

  write.csv(cox_results, "output_cox_regression_results.csv", row.names = FALSE)
  cat("  ✓ Saved output_cox_regression_results.csv\n")

  write.csv(forest_data, "output_forest_plot_data.csv", row.names = FALSE)
  cat("  ✓ Saved output_forest_plot_data.csv\n\n")

}

# ============================================================================
# Summary
# ============================================================================

cat("============================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("============================================================================\n\n")

cat("Data Summary:\n")
cat("  Total samples analyzed:", ncol(expr_clean), "\n")
cat("  Total events (deaths):", sum(survival_status_clean), "\n")
cat("  Total probes analyzed:", nrow(expr_clean), "\n")
cat("  Significant probes (p<0.05):", sum(cox_results$p_value < 0.05), "\n\n")

if(nrow(cox_results) > 5) {
  cat("Top 5 Prognostic Probes:\n")
  print(head(cox_results, 5))
  cat("\n")
}

cat("Output Files:\n")
cat("  1. output_forest_plot.png - Publication-ready forest plot (300 DPI)\n")
cat("  2. output_cox_regression_results.csv - All Cox regression results\n")
cat("  3. output_forest_plot_data.csv - Forest plot data points\n\n")

cat("Paper Reference:\n")
cat("  Journal: Discovery Oncology (Springer Nature)\n")
cat("  Published: February 2025\n")
cat("  DOI: 10.1007/s12672-025-01907-7\n")
cat("  Dataset: GSE84437 (483 gastric cancer samples)\n\n")

cat("============================================================================\n")
