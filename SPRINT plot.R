# plot for SPRINT data
path = "~/Desktop/Research/Yao/HTE inference/code/Panning/0630"
plotDirectory = file.path("~/Desktop/Research/Yao/HTE inference/code/Panning/0630")

library(ggplot2)
library(reshape2)
library(gridExtra)

myHeatmap = function(data, title = ""){
  p = ggplot(data, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradiennt(
    colors = c("red", "lightcoral", "grey", "lightblue"),
    trans = "log",
    breaks = c(1e-2, 1e-1, 1),
    labels = c("1e-2", "1e-1", "1"),
    limits = c(1e-3, 1)
  ) +
  geom_text(aes(label = round(value, 2))) +
  labs(title = title, x = "", y = "") +
  theme_minimal()
  return(p)
}

setting = c("SPRINT")
record = readRDS(file.path(path, paste(setting, "rds", sep = ".")))

# RT
plot.data = matrix(record$pValue$RT, nrow = 4, ncol = 4)
colnames(plot.data) = paste0('Age', c("<62", " 62-67", " 68-76", ">76"))
rownames(plot.data) = paste(rep(c("low BMI", "high BMI"), times = c(2, 2)), rep(c("Male", "Female"), 2), sep = "/")

# Convert the matrix to a data frame
data_melt = melt(plot.data)
heatmapRT = myHeatmap(data = data_melt, title = "RT")

# SSRT
plot.data = matrix(record$pValue$SSRT, nrow = 4, ncol = 4)
colnames(plot.data) = paste0('Age', c("<62", " 62-67", " 68-76", ">76"))
rownames(plot.data) = paste(rep(c("low BMI", "high BMI"), times = c(2, 2)), rep(c("Male", "Female"), 2), sep = "/")

# Convert the matrix to a data frame
data_melt = melt(plot.data)
heatmapSSRT = myHeatmap(data = data_melt, title = "SSRT")

# ART
plot.data = matrix(record$pValue$ART, nrow = 4, ncol = 4)
colnames(plot.data) = paste0('Age', c("<62", " 62-67", " 68-76", ">76"))
rownames(plot.data) = paste(rep(c("low BMI", "high BMI"), times = c(2, 2)), rep(c("Male", "Female"), 2), sep = "/")

# Convert the matrix to a data frame
data_melt = melt(plot.data)
heatmapART = myHeatmap(data = data_melt, title = "ART")

# DDRT
plot.data = matrix(record$pValue$DDRT, nrow = 4, ncol = 4)
colnames(plot.data) = paste0('Age', c("<62", " 62-67", " 68-76", ">76"))
rownames(plot.data) = paste(rep(c("low BMI", "high BMI"), times = c(2, 2)), rep(c("Male", "Female"), 2), sep = "/")

# Convert the matrix to a data frame
data_melt = melt(plot.data)
heatmapDDRT = myHeatmap(data = data_melt, title = "DDRT")


# Arrange the two heatmaps side by side
# pdf(file = paste(plotDirectory, "/", setting, "Heatmap.pdf", sep = ""), width = 7.5, height = 7)
grid.arrange(heatmapRT, heatmapSSRT, heatmapART, heatmapDDRT, ncol = 2, nrow = 2)
# dev.off()