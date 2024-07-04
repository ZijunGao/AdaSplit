# plot for number of knockoffs
path = "~/Desktop/Research/Yao/HTE inference/code/Panning/0630"
plotDirectory = file.path("~/Desktop/Research/Yao/HTE inference/code/Panning/0630")

library(ggplot2)


setting = c("NumberKnockoff")
B.seq = c(seq(1, 4), seq(5, 40, by = 5))

record = readRDS(file.path(path, paste(setting, "rds", sep = ".")))

plot.data = data.frame(
  B = B.seq,
  FDR = apply(record$FDP$ART, 2, mean),
  power =  apply(record$power$ART, 2, mean)
)

line.width = 1; point.size = 2
pdf(file = paste(plotDirectory, "/", setting, ".pdf", sep = ""), width = 3.5, height = 3.5)

ggplot(plot.data, aes(x = B, y = power)) +
  geom_point(aes(y = power, color = "ART"), size = point.size) +
  geom_line(aes(y = power, color = "ART"), size = line.width) +  # Plot the curve
  geom_hline(aes(yintercept = mean(record$power$ORT), color = "ORT"), size = line.width) +  # Add horizontal reference line
  scale_color_manual(values = c("ORT" = "orange", "ART" = "dark red")) +
  labs(# title = "Power ",
       x = "Number of knockoffs",
       y = "Power",
       color = "Method") +
  ylim(0, 1) + 
  theme_bw() 

dev.off()