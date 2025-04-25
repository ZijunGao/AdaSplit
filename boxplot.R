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
             "null more repeats",
             "no marginalization larger sample size",
             "no marginalization larger noise")

  count = count + 1
  setting = settings[count]
  record = readRDS(file.path(directory, "April 2025", paste(setting, ".rds", sep = "")))
  p_values_df_rt <- as.data.frame(record$pval$RT)
  p_values_df_ssrt <- as.data.frame(record$pval$SSRT)
  p_values_df_art <- as.data.frame(record$pval$ART)
  Group.number = dim(p_values_df_rt)[2]
  
  colnames(p_values_df_rt) <- paste("G", 1:Group.number, sep = "")
  colnames(p_values_df_ssrt) <- paste("G", 1:Group.number, sep = "")
  colnames(p_values_df_art) <- paste("G", 1:Group.number, sep = "")
  
  p_values_df_rt$Method <- "RT"
  p_values_df_art$Method <- "ART"
  p_values_df_ssrt$Method <- "SSRT"
  
  combined_p_values_df <- rbind(p_values_df_rt, p_values_df_ssrt, p_values_df_art)
  combined_p_values_long <- reshape2::melt(combined_p_values_df, id.vars = "Method")
  combined_p_values_long$Method <- factor(
    combined_p_values_long$Method,
    levels = c("RT", "SSRT", "ART") 
  )
  combined_p_values_long$variable <- factor(combined_p_values_long$variable, levels = paste("G", 1:Group.number, sep = ""))
  
  pdf(paste0(directory, "/Comparison_setting_", setting, ".pdf"), width = 3, height = 3)
  ggplot(combined_p_values_long, aes(x = variable, y = value, fill = Method)) +
    geom_boxplot(
      position = position_dodge(width = 0.75), 
      notch = FALSE, 
      outlier.size = 1,
      width = 0.6
    ) +
    scale_fill_manual(
      values = c("RT" = "#A3C1AD", "SSRT" = "gold", "ART" = "#8C1515"),
      labels = c("RT" = "RT", "SSRT" = "RRT", "ART" = "ART")
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
      legend.title = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  dev.off()
  count


