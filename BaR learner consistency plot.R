library(ggplot2)
library(dplyr)
library(tidyr)

# Load dataset
base_path <- "~/Desktop/Research/Yao/HTE inference/code/Panning/July 2025"
plotDirectory = "~/Desktop/Research/Yao/HTE inference/code/Panning"

record_bar <- readRDS(file.path(base_path, "BaR_learner_consistency.rds"))
sample_sizes <- record_bar$n.seq

# Function to extract R^2 summaries
get_r2_summary <- function(record, method_label) {
  df <- as.data.frame(record$R2)
  colnames(df) <- paste0("n", sample_sizes)
  df_long <- df %>%
    pivot_longer(everything(), names_to = "sample_size", values_to = "R2") %>%
    mutate(
      sample_size = as.integer(sub("n", "", sample_size)),
      Method = method_label
    )
  df_long %>%
    group_by(sample_size, Method) %>%
    summarise(
      mean = mean(R2),
      sd = sd(R2),
      upper = mean + sd,
      lower = mean - sd,
      .groups = "drop"
    )
}

# Get R^2 summary for BaR-learner only
r2_bar <- get_r2_summary(record_bar, "BaR-learner")

# Define color and linetype maps
color_map <- c("BaR-learner" = "#8C1515")
linetype_map <- c("BaR-learner" = "solid")

# Plot R^2
g <- ggplot(r2_bar, aes(x = sample_size, y = mean, color = Method, linetype = Method)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 100, size = 0.7) +
  scale_x_continuous(breaks = seq(1, 5) * 300) +
  scale_color_manual(values = color_map) +
  scale_linetype_manual(values = linetype_map) +
  labs(
    x = "Sample size",
    y = expression(R^2),
    color = "Method",
    linetype = "Method"
  ) +
  scale_y_continuous(breaks = seq(2, 10) * 0.1, labels = seq(2, 10) * 0.1) +  ylim(0.3, 1.05) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none",
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# pdf(file = paste(plotDirectory, "/", "Bar_learner_consistency_R2", ".pdf", sep = ""), width = 3, height = 3)
print(g)
# dev.off()

# Function to extract beta error summaries
get_betaError_summary <- function(record, method_label) {
  df <- as.data.frame(record$betaError)
  colnames(df) <- paste0("n", sample_sizes)
  df_long <- df %>%
    pivot_longer(everything(), names_to = "sample_size", values_to = "beta_error") %>%
    mutate(
      sample_size = as.integer(sub("n", "", sample_size)),
      Method = method_label
    )
  df_long %>%
    group_by(sample_size, Method) %>%
    summarise(
      mean = mean(beta_error),
      sd = sd(beta_error),
      upper = mean + sd,
      lower = mean - sd,
      .groups = "drop"
    )
}

# Get beta error summary for BaR-learner only
betaError_bar <- get_betaError_summary(record_bar, "BaR-learner")

# Plot beta error
g <- ggplot(betaError_bar, aes(x = sample_size, y = mean, color = Method, linetype = Method)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 100, size = 0.7) +
  scale_color_manual(values = color_map) +
  scale_linetype_manual(values = linetype_map) +
  labs(
    x = "Sample size",
    y = expression(frac("||" * hat(beta) - beta * "||"[2]^2, "||" * beta * "||"[2]^2)),
    color = "Method",
    linetype = "Method"
  ) +
  scale_x_continuous(breaks = seq(1, 5) * 300) +
  ylim(-0.04, 0.5) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none",
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# pdf(file = paste(plotDirectory, "/", "Bar_learner_consistency_betaError", ".pdf", sep = ""), width = 3.5, height = 3)
print(g)
# dev.off()
