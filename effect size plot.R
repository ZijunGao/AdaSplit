# plot for effect size
path = "~/Desktop/Research/Yao/HTE inference/code/Panning/0630"
plotDirectory = file.path("~/Desktop/Research/Yao/HTE inference/code/Panning/0630")

setting = c("EffectSize")
record = readRDS(file.path(path, paste(setting, "rds", sep = ".")))
delta.seq = c(c(0.5, 1, 1.5, 2, 2.5))
q = 0.2 # FDR level

library(ggplot2)

plot.data <- data.frame(
  delta = delta.seq,
  RT_FDR = apply(record$FDP$RT, 2, mean),
  SSRT_FDR = apply(record$FDP$SSRT, 2, mean),
  ART_FDR = apply(record$FDP$ART, 2, mean),
  RT_power = apply(record$power$RT, 2, mean),
  SSRT_power = apply(record$power$SSRT, 2, mean),
  ART_power = apply(record$power$ART, 2, mean)
)

line.width = 2; point.size = 3
pdf(file = paste(plotDirectory, "/", setting, "FDR", ".pdf", sep = ""), width = 3.8, height = 3.5)
ggplot(plot.data, aes(x = delta)) +
  geom_hline(yintercept = q, linetype = "dashed", color = "black") +  # Add horizontal reference line
  geom_point(aes(y = RT_FDR, color = "RT"), size = point.size) +
  geom_line(aes(y = RT_FDR, color = "RT"), size = line.width) +  # RT curve in black
  geom_point(aes(y = SSRT_FDR, color = "RT(SS)"), size = point.size) +
  geom_line(aes(y = SSRT_FDR, color = "RT(SS)"), size = line.width) +  # SSRT curve in dark blue
  geom_point(aes(y = ART_FDR, color = "ART"), size = point.size) +  
  geom_line(aes(y = ART_FDR, color = "ART"), size = line.width ) +  # ART curve in dark red
  scale_color_manual( values = c("RT" = "#9467bd", "RT(SS)" = "#1f77b4", "ART" = "#d35400"), breaks = c("RT", "RT(SS)", "ART")) +
  labs(# title = "FDR",
    x = "Effect size",
    y = "FDR",
    color = "Method") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +  # Custom y-axis breaks 
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 14),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),  # Increase y-axis title size
    axis.text.x = element_text(size = 13),   # Increase x-axis text size
    axis.text.y = element_text(size = 13),   # Increase y-axis text size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 13)    # Increase legend text size
  )
dev.off()

pdf(file = paste(plotDirectory, "/", setting, "Power", ".pdf", sep = ""), width = 3.8, height = 3.5)
line.width = 2; point.size = 3
ggplot(plot.data, aes(x = delta)) +
  geom_point(aes(y = RT_power, color = "RT"), size = point.size) +
  geom_line(aes(y = RT_power, color = "RT"), size = line.width) +  # RT curve in black
  geom_point(aes(y = SSRT_power, color = "RT(SS)"), size = point.size) +
  geom_line(aes(y = SSRT_power, color = "RT(SS)"), size = line.width) +  # SSRT curve in dark blue
  geom_point(aes(y = ART_power, color = "ART"), size = point.size) +  
  geom_line(aes(y = ART_power, color = "ART"), size = line.width ) +  # ART curve in dark red
  scale_color_manual( values = c("RT" = "#9467bd", "RT(SS)" = "#1f77b4", "ART" = "#d35400"), breaks = c("RT", "RT(SS)", "ART")) +
  labs(# title = "Power",
    x = "Effect size",
    y = "Power",
    color = "Method") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +  # Custom y-axis breaks
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 14),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),  # Increase y-axis title size
    axis.text.x = element_text(size = 13),   # Increase x-axis text size
    axis.text.y = element_text(size = 13),   # Increase y-axis text size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 13)    # Increase legend text size
  )
dev.off()


