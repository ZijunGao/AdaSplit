# plot for comparing AdaSplit and sample splitting with unit selection (SSSSRT) (units with positive predicted treatment effects).
path = "~/Desktop/Research/Yao/HTE inference/code/Panning/April 2025"
plotDirectory = file.path("~/Desktop/Research/Yao/HTE inference/code/Panning/0630")

library(ggplot2)
library(reshape2)
library(dplyr)
directory = "~/Desktop/Research/Yao/HTE inference/code/Panning" 

setting = "smaller sample size"
record = readRDS(file.path(directory, "April 2025", paste(setting, ".rds", sep = "")))

inference_ate_ssssrt <- as.data.frame(record$pval$SSSSRT)
inference_ate_art <- as.data.frame(record$pval$ART)
Group.number = dim(inference_ate_ssssrt)[2]

colnames(inference_ate_ssssrt) <- paste("G", 1:Group.number, sep = "")
colnames(inference_ate_art) <- paste("G", 1:Group.number, sep = "")

inference_ate_art$Method <- "ART"
inference_ate_ssssrt$Method <- "SSSSRT"

combined_inference_ate <- rbind(inference_ate_ssssrt, inference_ate_art)
combined_inference_ate_long <- reshape2::melt(combined_inference_ate, id.vars = "Method")
combined_inference_ate_long$Method <- factor(
  combined_inference_ate_long$Method,
  levels = c("RT", "SSRT", "SSSSRT", "ART") 
)
combined_inference_ate_long$variable <- factor(combined_inference_ate_long$variable, levels = paste("G", 1:Group.number, sep = ""))

pdf(paste0(directory, "/P_value ", setting, ".pdf"), width = 3, height = 3)
ggplot(combined_inference_ate_long, aes(x = variable, y = value, fill = Method)) +
  geom_boxplot(
    position = position_dodge(width = 0.7),  # 0.5
    notch = FALSE, 
    outlier.size = 1,
    width = 0.6
  ) +
  # geom_hline(yintercept = 0.0, linetype = "dashed", color = "black", linewidth = 0.6) +  
  ylim(-0, 1.0) +
  scale_fill_manual(
    values = c("RT" = "#A3C1AD", "SSRT" = "gold", "SSSSRT" = "coral", "ART" = "#8C1515"),
    labels = c("RT" = "RT", "SSRT" = "RT (RandomSplit)", "SSSSRT" = expression("RT (RandomSplit with " *hat(tau) > 0 * ")"), "ART" = "RT (AdaSplit)")
  ) +
  labs(
    x = NULL,
    y = "P-values"
    # title = paste0("Comparison of RT, ART, SSRT (Maximum Inference Fold: ", (1 - proportion) * 100, "%)")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom", 
    legend.direction = "vertical", # vertical legend
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

dev.off()


