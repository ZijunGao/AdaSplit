library(ggplot2)
library(reshape2)
library(dplyr)

directory = "~/Desktop/Research/Yao/HTE inference/code/Panning" 

library(ggplot2)
library(reshape2)
library(dplyr)
directory = "~/Desktop/Research/Yao/HTE inference/code/Panning" 

settings = c("default",
             "larger sample size",
             "larger noise",
             "even larger noise",
             "larger sample size larger noise",
             "null",
             "no marginalization", 
             "mu xgboost",
             "fewer for nuisance",
             "fewer for inference",
             "null more repeats")

count = 2 # "default": 1; "larger sample size": 2
setting = settings[count]
record = readRDS(file.path(directory, "April 2025", paste(setting, ".rds", sep = "")))

inference_ate_ssrt <- as.data.frame(record$inference.ATE$SSRT)
inference_ate_art <- as.data.frame(record$inference.ATE$ART)
Group.number = dim(inference_ate_ssrt)[2]

colnames(inference_ate_ssrt) <- paste("G", 1:Group.number, sep = "")
colnames(inference_ate_art) <- paste("G", 1:Group.number, sep = "")

inference_ate_art$Method <- "ART"
inference_ate_ssrt$Method <- "SSRT"

combined_inference_ate <- rbind(inference_ate_ssrt, inference_ate_art)
combined_inference_ate_long <- reshape2::melt(combined_inference_ate, id.vars = "Method")
combined_inference_ate_long$Method <- factor(
  combined_inference_ate_long$Method,
  levels = c("RT", "SSRT", "ART") 
)
combined_inference_ate_long$variable <- factor(combined_inference_ate_long$variable, levels = paste("G", 1:Group.number, sep = ""))

pdf(paste0(directory, "/inference_ate_", setting, ".pdf"), width = 3, height = 3)
ggplot(combined_inference_ate_long, aes(x = variable, y = value, fill = Method)) +
  geom_boxplot(
    position = position_dodge(width = 0.5), 
    notch = FALSE, 
    outlier.size = 1,
    width = 0.6
  ) +
  geom_hline(yintercept = 0.0, linetype = "dashed", color = "black", linewidth = 0.6) +  
  ylim(-0.2, 1.4) +
  scale_fill_manual(
    values = c("RT" = "#A3C1AD", "SSRT" = "gold", "ART" = "#8C1515"),
    labels = c("RT" = "RT", "SSRT" = "RRT", "ART" = "ART")
  ) +
  labs(
    x = NULL,
    y = "Inference fold ATE"
    # title = paste0("Comparison of RT, ART, SSRT (Maximum Inference Fold: ", (1 - proportion) * 100, "%)")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom", 
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

dev.off()
