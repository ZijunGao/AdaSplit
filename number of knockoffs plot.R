# plot for number of knockoffs
path = "~/Desktop/Research/Yao/HTE inference/code/Panning/0630"
plotDirectory = file.path("~/Desktop/Research/Yao/HTE inference/code/Panning/0630")

library(ggplot2)


setting = c("NumberKnockoff")
B.seq = c(seq(5,50, by = 5))

record = readRDS(file.path(path, paste(setting, "rds", sep = ".")))

plot.data = data.frame(
  B = B.seq,
  FDR = apply(record$FDP$ART, 2, mean),
  power =  apply(record$power$ART, 2, mean)
)

line.width = 2; point.size = 3
pdf(file = paste(plotDirectory, "/", setting, ".pdf", sep = ""), width =9, height = 3)

ggplot(plot.data, aes(x = B, y = power)) +
  geom_point(aes(y = power, color = "ART"), size = point.size) +
  geom_line(aes(y = power, color = "ART"), size = line.width) +  # Plot the curve
  geom_line(aes(y = rep(mean(record$power$ORT),10), color = "RT(Oracle)"), size = line.width) + 
  geom_point(aes(y = rep(mean(record$power$ORT),10), color = "RT(Oracle)"), size = point.size) +
  # Add horizontal reference line
  scale_color_manual(values = c("RT(Oracle)" = "#7f7f7f", "ART" = "#d35400"), breaks = c("RT(Oracle)", "ART")) +
  labs(# title = "Power ",
       x = "Number of assignments",
       y = "Power",
       color = "Method") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1))+
  scale_x_continuous(breaks = seq(5, 50, by = 5), limits = c(5, 50))+
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),  # Increase y-axis title size
    axis.text.x = element_text(size = 13.3),   # Increase x-axis text size
    axis.text.y = element_text(size = 13),   # Increase y-axis text size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 13)    # Increase legend text size
  )

dev.off()