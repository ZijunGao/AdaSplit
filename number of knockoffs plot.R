# plot for number of knockoffs
path = "~/Desktop/Research/Yao/HTE inference/code/Panning/0630"
plotDirectory = file.path("~/Desktop/Research/Yao/HTE inference/code/Panning/0630")
setting.seq = c("NumberKnockoffDenoise", "NumberKnockoffITE")

library(ggplot2)


setting.seq = c("NumberKnockoffDenoise", "NumberKnockoffITE")
B.seq = seq(1, 10)
for(setting in setting.seq){
  record = readRDS(file.path(path, paste(setting, "rds", sep = ".")))
  
  plot.data = data.frame(
    B = B.seq,
    FDR = apply(record$FDP$ART, 2, mean),
    power =  apply(record$power$ART, 2, mean)
  )
  pdf(file = paste(plotDirectory, "/", setting, ".pdf", sep = ""), width = 3.5, height = 3.5)
  
  ggplot(plot.data, aes(x = B, y = power)) +
    geom_line(color = "blue") +  # Plot the curve
    geom_hline(yintercept = mean(record$power$ORT), linetype = "dashed", color = "red") +  # Add horizontal reference line
    labs(title = "Power ",
         x = "Number of knockoffs",
         y = "Power") +
    theme_minimal()
  
  dev.off()
}