suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(survival)
  library(ggplot2)
  library(dplyr)
})

cat("============================================================================\n")
cat("Panel A: Multivariate Cox Regression - Forest Plot\n")
cat("============================================================================\n\n")

# Load data
data_file <- "C:/Users/steli/Desktop/Figure Reproduction Project - Copy/input_data/GSE84437_series_matrix.txt"

gse84437 <- getGEO(filename = data_file, GSEMatrix = TRUE)
expr_data <- exprs(gse84437)
pdata <- pData(gse84437)
fdata <- fData(gse84437)

cat(" ✓ Loaded:", nrow(expr_data), "probes x", ncol(expr_data), "samples\n\n")

# ============================================================================
# Apply log2 transformation (CRITICAL - paper uses this!)
# ============================================================================

cat("Step 1: Applying log2 transformation to expression data...\n")

expr_data_log <- log2(expr_data + 1)  # +1 to avoid log(0)

cat("Step 2: Finding LMOD1 probe...\n")

lmod1_probe_id <- which(grepl("LMOD1", fdata$Symbol, ignore.case = TRUE))[1]
cat("Using LMOD1 probe:", rownames(fdata)[lmod1_probe_id], "\n\n")

# ============================================================================
# Extract clinical data
# ============================================================================

cat("Step 3: Extracting survival and clinical data...\n")

survival_status <- as.numeric(pdata[["death:ch1"]])
survival_time <- as.numeric(pdata[["duration overall survival:ch1"]])

survival_idx <- which(!is.na(survival_time) & !is.na(survival_status))
cat("Samples with survival data:", length(survival_idx), "\n")

age <- as.numeric(pdata[survival_idx, "age:ch1"])
gender_raw <- as.character(pdata[survival_idx, "Sex:ch1"])
t_stage <- as.numeric(gsub("[^0-9]", "", pdata[survival_idx, "ptstage:ch1"]))
n_stage <- as.numeric(gsub("[^0-9]", "", pdata[survival_idx, "pnstage:ch1"]))
lmod1_expr <- as.numeric(expr_data_log[lmod1_probe_id, survival_idx])  # Use LOG-TRANSFORMED data

cat("Unique gender values:", paste(unique(gender_raw), collapse = ", "), "\n\n")

# ============================================================================
# Prepare Cox data
# ============================================================================

cat("Step 4: Building Cox data frame...\n")

gender_numeric <- as.numeric(factor(gender_raw, levels = rev(unique(gender_raw)))) - 1

cox_data <- data.frame(
  time = survival_time[survival_idx],
  status = survival_status[survival_idx],
  LMOD1 = lmod1_expr,
  Age = age,
  Gender = gender_numeric,
  T_stage = t_stage,
  N_stage = n_stage,
  stringsAsFactors = FALSE
)

valid_rows <- complete.cases(cox_data)
cox_data <- cox_data[valid_rows, ]

cat("Valid samples:", nrow(cox_data), "\n")
cat("Events:", sum(cox_data$status == 1), "\n")
cat("Censored:", sum(cox_data$status == 0), "\n\n")

# Standardize continuous variables
cox_data$LMOD1 <- scale(cox_data$LMOD1)[,1]
cox_data$Age <- scale(cox_data$Age)[,1]

# ============================================================================
# Fit Cox model
# ============================================================================

cat("Step 5: Fitting multivariate Cox model...\n\n")

fit <- coxph(Surv(time, status) ~ LMOD1 + Age + Gender + T_stage + N_stage, data = cox_data)
summary_fit <- summary(fit)

print(summary_fit)
cat("\n")

# Extract results
coef_data <- summary_fit$coefficients
conf_data <- summary_fit$conf.int

forest_data <- data.frame(
  Variable = c("LMOD1", "Age", "Gender", "T stage", "N stage"),
  HR = coef_data[, "exp(coef)"],
  CI_lower = conf_data[, "lower .95"],
  CI_upper = conf_data[, "upper .95"],
  p_value = coef_data[, "Pr(>|z|)"]
)

forest_data$Significant <- ifelse(is.na(forest_data$p_value), "No", 
                                  ifelse(forest_data$p_value < 0.05, "Yes", "No"))

cat("Forest Plot Data:\n")
print(forest_data)
cat("\n")

# ============================================================================
# Create forest plot
# ============================================================================

cat("Step 6: Creating forest plot...\n")

forest_plot_data <- forest_data[!is.na(forest_data$HR), ]
forest_plot_data$Variable <- factor(forest_plot_data$Variable, 
                                    levels = rev(forest_plot_data$Variable))

p <- ggplot(forest_plot_data, aes(x = HR, y = Variable, color = Significant)) +
  geom_point(size = 5) +
  geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper), width = 0.25, linewidth = 1.2, 
                orientation = "y") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 1.2) +
  scale_color_manual(values = c("Yes" = "#E41A1C", "No" = "#4DAF4A")) +
  scale_x_continuous(limits = c(0.5, 2.5), breaks = seq(0.5, 2.5, 0.5)) +
  labs(
    x = "Hazard Ratio (95% CI)",
    y = "",
    color = "p < 0.05"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    legend.position = "right"
  )

pdf("output_forest_plot_panel_A.pdf", width = 10, height = 6)
print(p)
dev.off()

png("output_forest_plot_panel_A.png", width = 10, height = 6, units = "in", res = 300)
print(p)
dev.off()

cat(" ✓ Saved: output_forest_plot_panel_A.pdf\n")
cat(" ✓ Saved: output_forest_plot_panel_A.png\n\n")

write.csv(forest_data, "output_cox_results.csv", row.names = FALSE)

