library(ggplot2)
library(dplyr)
library(reshape2)

# Settings and directory
settings <- c("default",
              "larger sample size",
              "larger noise",
              "even larger noise",
              "larger sample size larger noise",
              "null",
              "no marginalization", 
              "mu xgboost",
              "fewer for nuisance",
              "fewer for inference",
              "null more repeats",
              "no marginalization larger sample size",
              "no marginalization larger noise")
directory = "~/Desktop/Research/Yao/HTE inference/code/Panning" 

record_default <- readRDS(file.path(directory, "July 2025", paste(settings[1], ".rds", sep = "")))
inference_art_default <- as.data.frame(record_default$inference.proportion$ART)
Group.number <- ncol(inference_art_default)
colnames(inference_art_default) <- paste0("G", 1:Group.number)
inference_art_default$Method <- "ART (n = 500)"

# Load ART from "larger sample size"
record_largern <- readRDS(file.path(directory, "April 2025", paste(settings[2], ".rds", sep = "")))
inference_art_largern <- as.data.frame(record_largern$inference.proportion.before.throw.away$ART)
colnames(inference_art_largern) <- paste0("G", 1:Group.number)
inference_art_largern$Method <- "ART (n = 1000)"

# Combine ART datasets only
combined_data <- rbind(inference_art_default, inference_art_largern)
combined_long <- reshape2::melt(combined_data, id.vars = "Method")

combined_long$Method <- factor(
  combined_long$Method,
  levels = c("ART (n = 500)", "ART (n = 1000)")
)
combined_long$variable <- factor(combined_long$variable, levels = paste0("G", 1:Group.number))

# Compute mean and SD
summary_df <- combined_long %>%
  group_by(Method, variable) %>%
  summarise(
    mean = mean(value),
    sd = sd(value),
    .groups = "drop"
  )

# Plot bar chart with ±1 SD error bars
# pdf(paste0(directory, "/inference_proportion_", setting, "_no_throw_away.pdf"), width = 3, height = 3)
# pdf(paste0(directory, "/inference_proportion_", setting, ".pdf"), width = 3, height = 3)
ggplot(summary_df, aes(x = variable, y = mean, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.5, color = "black") +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    position = position_dodge(width = 0.6),
    width = 0.3
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.6) +
  coord_cartesian(ylim = c(0.35, 0.98)) +
  scale_fill_manual(
    values = c("ART (n = 500)"  = "lightsalmon1", "ART (n = 1000)" = "red3"),
    labels = c("n = 500", "n = 1000")
  ) +
  labs(
    x = NULL,
    y = "Inference fold proportion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom", 
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )
# dev.off()

