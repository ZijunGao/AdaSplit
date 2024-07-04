# plot for sample size
path = "~/Desktop/Research/Yao/HTE inference/code/Panning/0630"
plotDirectory = file.path("~/Desktop/Research/Yao/HTE inference/code/Panning/0630")

setting = c("SampleSize")
record = readRDS(file.path(path, paste(setting, "rds", sep = ".")))
n.seq = seq(300, 1000, by = 100)
q = 0.2 # FDR level
plot.data <- data.frame(
  n = n.seq,
  RT_FDR = apply(record$FDP$RT, 2, mean),
  SSRT_FDR = apply(record$FDP$SSRT, 2, mean),
  ART_FDR = apply(record$FDP$ART, 2, mean),
  RT_power = apply(record$power$RT, 2, mean),
  SSRT_power = apply(record$power$SSRT, 2, mean),
  ART_power = apply(record$power$ART, 2, mean)
)

library(ggplot2)
pdf(file = paste(plotDirectory, "/", setting, "FDR", ".pdf", sep = ""), width = 3.5, height = 3.5)
line.width = 1; point.size = 2
ggplot(plot.data, aes(x = n)) +
  geom_point(aes(y = RT_FDR, color = "RT"), size = point.size) +
  geom_line(aes(y = RT_FDR, color = "RT"), size = line.width) +  # RT curve in black
  geom_point(aes(y = SSRT_FDR, color = "SSRT"), size = point.size) +
  geom_line(aes(y = SSRT_FDR, color = "SSRT"), size = line.width) +  # SSRT curve in dark blue
  geom_point(aes(y = ART_FDR, color = "ART"), size = point.size) +  
  geom_line(aes(y = ART_FDR, color = "ART"), size = line.width ) +  # ART curve in dark red
  geom_hline(yintercept = q, linetype = "dashed", color = "red") +  # Add horizontal reference line
  scale_color_manual( values = c("RT" = "black", "SSRT" = "dark blue", "ART" = "dark red"), breaks = c("RT", "SSRT", "ART")) +
  labs(# title = "FDR",
    x = "Sample size",
    y = "FDR",
    color = "Method") +
  ylim(0, 1) +
  theme_bw() 
dev.off()

line.width = 1; point.size = 2
pdf(file = paste(plotDirectory, "/", setting, "Power", ".pdf", sep = ""), width = 3.5, height = 3.5)
ggplot(plot.data, aes(x = n)) +
  geom_point(aes(y = RT_power, color = "RT"), size = point.size) +
  geom_line(aes(y = RT_power, color = "RT"), size = line.width) +  # RT curve in black
  geom_point(aes(y = SSRT_power, color = "SSRT"), size = point.size) +
  geom_line(aes(y = SSRT_power, color = "SSRT"), size = line.width) +  # SSRT curve in dark blue
  geom_point(aes(y = ART_power, color = "ART"), size = point.size) +  
  geom_line(aes(y = ART_power, color = "ART"), size = line.width ) +  # ART curve in dark red
  scale_color_manual( values = c("RT" = "black", "SSRT" = "dark blue", "ART" = "dark red"), breaks = c("RT", "SSRT", "ART")) +
  labs(# title = "Power",
    x = "Sample size",
    y = "Power",
    color = "Method") +
  ylim(0, 1) +
  theme_bw() 
dev.off()