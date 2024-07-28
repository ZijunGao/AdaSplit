# plot for sample size
path = "~/Desktop/Research/Yao/HTE inference/code/Panning/0630"
plotDirectory = file.path("~/Desktop/Research/Yao/HTE inference/code/Panning/0630")
setting = c("unbalanced")

library(ggplot2)

record = readRDS(file.path(path, paste(setting, "rds", sep = ".")))
alpha.seq = c(seq(1,25)/1000, seq(1, 20) / 20)
q = 0.2 # FDR level
index = 11
method = c("ORT", "ART (overfitting)", "ART")
color_values <- c("AIPW normalized" = "#d62728",
                  "AIPW" = "#ff9896",
                  "denoise normalized" = "#1f77b4", 
                  "denoise" = "#aec7e8")
# values = c("AIPW normalized" = "dark red", "AIPW" = "coral", "denoise normalized" = "blue", "denoise" = "light blue"), breaks = c("AIPW normalized", "AIPW", "denoise normalized", "denoise")
for(i in 1 : length(method)){  
  plot.data = data.frame(lapply(record$pValue[[i]], function(x){sapply(alpha.seq, function(y)(mean(x[, index] <= y)))}))
  plot.data$alpha = alpha.seq
  
  pdf(file = paste(plotDirectory, "/", setting, method[i], ".pdf", sep = ""), width = 3.8, height = 3.5)
  line.width = 0.75
  point.size = 0
  g = ggplot(plot.data, aes(x = alpha)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  # Add diagonal reference line
    geom_point(aes(y = AIPW.normalized, color = "AIPW normalized"), size = point.size) +
    geom_line(aes(y = AIPW.normalized, color = "AIPW normalized"), size = line.width) +  # AIPW.normalized curve in red
    geom_point(aes(y = AIPW, color = "AIPW"), size = point.size) +
    geom_line(aes(y = AIPW, color = "AIPW"), size = line.width) +  # AIPW curve in coral
    geom_point(aes(y = denoise.normalized, color = "denoise normalized"), size = point.size) +
    geom_line(aes(y = denoise.normalized, color = "denoise normalized"), size = line.width) +  # denoise.normalized curve in blue
    geom_point(aes(y = denoise, color = "denoise"), size = point.size) +
    geom_line(aes(y = denoise, color = "denoise"), size = line.width) +  # denoise curve in light blue
    scale_color_manual(values = color_values, breaks = c("AIPW normalized", "AIPW", "denoise normalized", "denoise")) +
    labs(
      x = "alpha",
      y = "CDF",
      color = "Method",
      title = method[i]) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +  # Custom y-axis breaks 
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.title.x = element_text(size = 14),  # Increase x-axis title size
      axis.title.y = element_text(size = 14),  # Increase y-axis title size
      axis.text.x = element_text(size = 13),   # Increase x-axis text size
      axis.text.y = element_text(size = 13),   # Increase y-axis text size
      legend.title = element_text(size = 14),  # Increase legend title size
      legend.text = element_text(size = 13)    # Increase legend text size
    )
  print(g)
  dev.off()
}


plot.data = data.frame(lapply(record$pValue$ART, function(x){sapply(alpha.seq, function(y)(mean(x[, index] <= y)))}))
plot.data$alpha = alpha.seq

pdf(file = paste(plotDirectory, "/", setting, "ART (overfitting)", ".pdf", sep = ""), width = 3.8, height = 3.5)
line.width = 0.75
point.size = 0
g = ggplot(plot.data, aes(x = alpha)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  # Add diagonal reference line
  geom_point(aes(y = AIPW.normalized, color = "AIPW normalized"), size = point.size) +
  geom_line(aes(y = AIPW.normalized, color = "AIPW normalized"), size = line.width) +  # AIPW.normalized curve in red
  geom_point(aes(y = AIPW, color = "AIPW"), size = point.size) +
  geom_line(aes(y = AIPW, color = "AIPW"), size = line.width) +  # AIPW curve in coral
  geom_point(aes(y = denoise.normalized, color = "denoise normalized"), size = point.size) +
  geom_line(aes(y = denoise.normalized, color = "denoise normalized"), size = line.width) +  # denoise.normalized curve in blue
  geom_point(aes(y = denoise, color = "denoise"), size = point.size) +
  geom_line(aes(y = denoise, color = "denoise"), size = line.width) +  # denoise curve in light blue
  scale_color_manual(values = color_values, breaks = c("AIPW normalized", "AIPW", "denoise normalized", "denoise")) +
  labs(
    x = "alpha",
    y = "CDF",
    color = "Method",
    title = "ART (overfitting)") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +  # Custom y-axis breaks 
  theme_bw() +
  theme(
    legend.background = element_rect(fill = "white", color = "black"), # Black box with white background
    legend.position = c(0.65, 0.3),  # none
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.x = element_text(size = 14),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),  # Increase y-axis title size
    axis.text.x = element_text(size = 13),   # Increase x-axis text size
    axis.text.y = element_text(size = 13),   # Increase y-axis text size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 13)    # Increase legend text size
  )
print(g)
dev.off()




