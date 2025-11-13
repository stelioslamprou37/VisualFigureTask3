#!/usr/bin/env Rscript

# ============================================================================
# Forest Plot - Multivariate Cox Regression Results
# Reads: extracted_solution.csv (Cox regression results)
# Outputs: extracted_solution.jpg (forest plot)
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# Read Cox regression results from current directory
cox_data <- read.csv("extracted_solution.csv", stringsAsFactors = FALSE)

cat("Loaded Cox regression results:\n")
print(cox_data)
cat("\n")

# Prepare data for forest plot
cox_data <- cox_data %>%
  mutate(
    Variable = factor(Variable, levels = rev(Variable)),
    color = ifelse(Significant == "Yes", "Yes", "No")
  )

# Create forest plot
p <- ggplot(cox_data, aes(y = Variable, x = HR, xmin = CI_lower, xmax = CI_upper)) +
  # Add vertical reference line at HR = 1
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +

  # Add error bars (confidence intervals)
  geom_errorbarh(aes(color = color), height = 0.3, linewidth = 1.2) +

  # Add point estimates
  geom_point(aes(color = color), size = 4, shape = 19) +

  # Color scheme: red for significant, green for non-significant
  scale_color_manual(
    values = c("Yes" = "#E41A1C", "No" = "#4DAF4A"),
    labels = c("Yes", "No"),
    name = "p < 0.05"
  ) +

  # Axis labels and title
  labs(
    x = "Hazard Ratio (95% CI)",
    y = "",
    title = ""
  ) +

  # Theme customization
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    plot.margin = margin(20, 20, 20, 20)
  ) +

  # Set x-axis limits and breaks
  scale_x_continuous(
    limits = c(0.5, 2.5),
    breaks = seq(0.5, 2.5, 0.5)
  )

# Save plot
ggsave(
  "extracted_solution.jpg",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300,
  device = "jpeg"
)

cat("✓ Generated extracted_solution.jpg\n")
