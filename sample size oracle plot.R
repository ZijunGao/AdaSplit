# plot for sample size
path = "~/Desktop/Research/Yao/HTE inference/code/Panning/0630"
plotDirectory = file.path("~/Desktop/Research/Yao/HTE inference/code/Panning/0630")


setting = c("SampleSize")
record = readRDS(file.path(path, paste(setting, "rds", sep = ".")))
n.seq = seq(400, 1200, by = 200)
q = 0.2 # FDR level

plot.data <- data.frame(
  n = n.seq,
  ORT_FDR = apply(record$FDP$ORT, 2, mean),
  DDRT_FDR = apply(record$FDP$DDRT, 2, mean),
  ART_FDR = apply(record$FDP$ART, 2, mean),
  ORT_power = apply(record$power$ORT, 2, mean),
  DDRT_power = apply(record$power$DDRT, 2, mean),
  ART_power = apply(record$power$ART, 2, mean)
)

library(ggplot2)
pdf(file = paste(plotDirectory, "/", setting, "FDR", "_Oracle.pdf", sep = ""), width = 3.8, height = 3.5)
line.width = 2; point.size = 3
ggplot(plot.data, aes(x = n)) +
  geom_hline(yintercept = q, linetype = "dashed", color = "black") +  # Add horizontal reference line
  geom_point(aes(y = ORT_FDR, color = "RT(Oracle)"), size = point.size) +
  geom_line(aes(y = ORT_FDR, color = "RT(Oracle)"), size = line.width) +  # RT curve in black
  geom_point(aes(y = DDRT_FDR, color = "RT(DD)"), size = point.size) +
  geom_line(aes(y = DDRT_FDR, color = "RT(DD)"), size = line.width) +  # SSRT curve in dark blue
  geom_point(aes(y = ART_FDR, color = "ART"), size = point.size) +  
  geom_line(aes(y = ART_FDR, color = "ART"), size = line.width ) +  # ART curve in dark red
  scale_color_manual(values = c("RT(Oracle)" = "#7f7f7f", "RT(DD)" = "#2ca02c", "ART" = "#d35400"), breaks = c("RT(Oracle)", "RT(DD)", "ART")) +
  labs(# title = "FDR",
    x = "Sample size",
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

line.width = 2; point.size = 3
pdf(file = paste(plotDirectory, "/", setting, "Power", "_Oracle.pdf", sep = ""), width = 3.8, height = 3.5)
ggplot(plot.data, aes(x = n)) +
  geom_point(aes(y = ORT_power, color = "RT(Oracle)"), size = point.size) +
  geom_line(aes(y = ORT_power, color = "RT(Oracle)"), size = line.width) +  # RT curve in black
  geom_point(aes(y = DDRT_power, color = "RT(DD)"), size = point.size) +
  geom_line(aes(y = DDRT_power, color = "RT(DD)"), size = line.width) +  # SSRT curve in dark blue
  geom_point(aes(y = ART_power, color = "ART"), size = point.size) +  
  geom_line(aes(y = ART_power, color = "ART"), size = line.width ) +  # ART curve in dark red
  scale_color_manual( values = c("RT(Oracle)" = "#7f7f7f" , "RT(DD)" = "#2ca02c", "ART" = "#d35400"), breaks = c("RT(Oracle)", "RT(DD)", "ART")) +
  labs(# title = "Power",
    x = "Sample size",
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



