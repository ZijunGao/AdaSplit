# plot for SPRINT data
path = "~/Desktop/Research/Yao/HTE inference/code/Panning/July 2025"
plotDirectory = file.path("~/Desktop/Research/Yao/HTE inference/code/Panning")

library(ggplot2)
library(reshape2)
library(gridExtra)

myHeatmap = function(data, title = ""){
  p = ggplot(data, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("firebrick3", "lightcoral", "lightsalmon",  "beige", "lightblue", "deepskyblue4"),
    trans = "log",
    # breaks = c(0.05, 0.2, 0.5, 1),
    # labels = c("0.05", "0.2", "0.5", "1"),
    limits = c(0.021, 1)
  ) +
  geom_text(aes(label = round(value, 2))) +
  labs(title = title, x = "", y = "") +
  theme_minimal()
  return(p)
}

setting = c("SPRINT")
record = readRDS(file.path(path, paste(setting, "rds", sep = ".")))

# RT
plot.data.merged = rbind(c(array(record$pValue$RT, dim = c(2,2,2))[,2,]), c(array(record$pValue$RT, dim = c(2,2,2))[,1,]))
colnames(plot.data.merged) = c("Low Risk; BMI < 27", "High Risk; BMI < 27", "Low Risk; BMI > 27", "High Risk;BMI > 27")
rownames(plot.data.merged) = c("Senior", "Junior")
# plot.data.merged = plot.data.merged[c(2,1),]
plot.data.merged = t(plot.data.merged)[,c(2,1)]
heatmap.RT.merged = myHeatmap(data = melt(plot.data.merged), title = "RT")
# heatmap.RT.merged

# RRT
index = 1
plot.data.merged = rbind(c(array(record$pValue$SSRT[index,], dim = c(2,2,2))[,2,]), c(array(record$pValue$SSRT[index,], dim = c(2,2,2))[,1,]))
colnames(plot.data.merged) = c("Low Risk; BMI < 27", "High Risk; BMI < 27", "Low Risk; BMI > 27", "High Risk;BMI > 27")
rownames(plot.data.merged) = c("Senior", "Junior")
plot.data.merged = t(plot.data.merged)[,c(2,1)]
heatmap.RRT.merged = myHeatmap(data = melt(plot.data.merged), title = "RT (RandomSplit)")
# heatmap.RRT.merged


# ART
plot.data.merged = rbind(c(array(record$pValue$ART, dim = c(2,2,2))[,2,]), c(array(record$pValue$ART, dim = c(2,2,2))[,1,]))
colnames(plot.data.merged) = c("Low Risk; BMI < 27", "High Risk; BMI < 27", "Low Risk; BMI > 27", "High Risk;BMI > 27")
rownames(plot.data.merged) = c("Senior", "Junior")
plot.data.merged = t(plot.data.merged)[,c(2,1)]
heatmap.ART.merged = myHeatmap(data = melt(plot.data.merged), title = "RT (AdaSplit)")
# heatmap.ART.merged 


# DDRT
plot.data.merged = rbind(c(array(record$pValue$DDRT, dim = c(2,2,2))[,2,]), c(array(record$pValue$DDRT, dim = c(2,2,2))[,1,]))
colnames(plot.data.merged) = c("Low Risk; BMI < 27", "High Risk; BMI < 27", "Low Risk; BMI > 27", "High Risk;BMI > 27")
rownames(plot.data.merged) = c("Senior", "Junior")
plot.data.merged = t(plot.data.merged)[,c(2,1)]
heatmap.DDRT.merged = myHeatmap(data = melt(plot.data.merged), title = "RT (Double-dipping)")
# heatmap.DDRT.merged


# # Arrange 1
# grid.arrange(heatmap.RT.merged,
#              heatmap.RRT.merged,
#              heatmap.ART.merged,
#              heatmap.DDRT.merged, ncol = 2, nrow = 2)
# # Arrange 2
# grid.arrange(heatmap.RT.merged,
#              heatmap.RRT.merged,
#              heatmap.ART.merged, ncol = 3, nrow = 1)


method = "RT" # "RT", "RRT", "ART", "DDRT"
# pdf(file = paste(plotDirectory, "/", "SPRINT","_", method, "_Heatmap.pdf", sep = ""), width = 4, height = 4)
heatmap.RT.merged
# dev.off()


method = "RRT" # "RT", "RRT", "ART", "DDRT"
# pdf(file = paste(plotDirectory, "/", "SPRINT","_", method, "_Heatmap.pdf", sep = ""), width = 4, height = 4)
heatmap.RRT.merged
# dev.off()

method = "ART" # "RT", "RRT", "ART", "DDRT"
# pdf(file = paste(plotDirectory, "/", "SPRINT","_", method, "_Heatmap.pdf", sep = ""), width = 4, height = 4)
heatmap.ART.merged
# dev.off()

method = "DDRT" # "RT", "RRT", "ART", "DDRT"
heatmap.DDRT.merged

