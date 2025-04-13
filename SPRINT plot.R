# plot for SPRINT data
path = "~/Desktop/Research/Yao/HTE inference/code/Panning/April 2025"
plotDirectory = file.path("~/Desktop/Research/Yao/HTE inference/code/Panning/April 2025")

library(ggplot2)
library(reshape2)
library(gridExtra)

myHeatmap = function(data, title = ""){
  p = ggplot(data, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("red", "lightcoral", "lightgrey", "lightblue"),
    trans = "log",
    breaks = c(0.1, 0.4, 1),
    labels = c("0.1", "0.4", "1"),
    limits = c(0.005, 1)
  ) +
  geom_text(aes(label = round(value, 2))) +
  labs(title = title, x = "", y = "") +
  theme_minimal()
  return(p)
}

setting = c("SPRINT")
record = readRDS(file.path(path, paste(setting, "rds", sep = ".")))

# X$High.risk.group = (X$RISK10YRS > 18) # 18 
# X$Senior.group = (X$AGE > 70) # 70
# X$Obese.group = (X$BMI > 27) # 27
# Group = X$Obese.group * 4 + X$Senior.group * 2 + X$High.risk.group

# RT
plot.data.merged = rbind(c(array(record$pValue$RT, dim = c(2,2,2))[,,1]), c(array(record$pValue$RT, dim = c(2,2,2))[,,2]))
colnames(plot.data.merged) = c("Young; Low Risk", "Senior; Low Risk", "Young; High Risk", "Senior; High Risk")
rownames(plot.data.merged) = c("BMI < 27", "BMI > 27")
# plot.data.merged = plot.data.merged[c(2,1),]
plot.data.merged = t(plot.data.merged)
heatmap.RT.merged = myHeatmap(data = melt(plot.data.merged), title = "RT")

plot.data = array(record$pValue$RT, dim = c(2,2,2))[,,1]
colnames(plot.data) = c("Young", "Senior")
rownames(plot.data) = c("Low Risk", "High Risk")
heatmap.RT.normal = myHeatmap(data = melt(plot.data), title = "RT; Normal")

plot.data = array(record$pValue$RT, dim = c(2,2,2))[,,2]
colnames(plot.data) = c("Young", "Senior")
rownames(plot.data) = c("Low Risk", "High Risk")
heatmap.RT.overweight = myHeatmap(data = melt(plot.data), title = "RT; Overweight")

# SSRT
index = 1
plot.data.merged = rbind(c(array(record$pValue$SSRT[index, ], dim = c(2,2,2))[,,1]), c(array(record$pValue$SSRT[index, ], dim = c(2,2,2))[,,2]))
colnames(plot.data.merged) = c("Young; Low Risk", "Senior; Low Risk", "Young; High Risk", "Senior; High Risk")
rownames(plot.data.merged) = c("BMI < 27", "BMI > 27")
# plot.data.merged = plot.data.merged[c(2,1),]
plot.data.merged = t(plot.data.merged)
heatmap.SSRT.merged = myHeatmap(data = melt(plot.data.merged), title = "SSRT")

plot.data = array(record$pValue$SSRT[index, ], dim = c(2,2,2))[,,1]
colnames(plot.data) = c("Young", "Senior")
rownames(plot.data) = c("Low Risk", "High Risk")
heatmap.SSRT.normal = myHeatmap(data = melt(plot.data), title = "SSRT; Normal")

plot.data = array(record$pValue$SSRT[index, ], dim = c(2,2,2))[,,2]
colnames(plot.data) = c("Young", "Senior")
rownames(plot.data) = c("Low Risk", "High Risk")
heatmap.SSRT.overweight = myHeatmap(data = melt(plot.data), title = "SSRT; Overweight")

# ART
plot.data.merged = rbind(c(array(record$pValue$ART, dim = c(2,2,2))[,,1]), c(array(record$pValue$ART, dim = c(2,2,2))[,,2]))
colnames(plot.data.merged) = c("Young; Low Risk", "Senior; Low Risk", "Young; High Risk", "Senior; High Risk")
rownames(plot.data.merged) = c("BMI < 27", "BMI > 27")
# plot.data.merged = plot.data.merged[c(2,1),]
plot.data.merged = t(plot.data.merged)
heatmap.ART.merged = myHeatmap(data = melt(plot.data.merged), title = "ART")

plot.data = array(record$pValue$ART, dim = c(2,2,2))[,,1]
colnames(plot.data) = c("Young", "Senior")
rownames(plot.data) = c("Low Risk", "High Risk")
heatmap.ART.normal = myHeatmap(data = melt(plot.data), title = "ART; Normal")

plot.data = array(record$pValue$ART, dim = c(2,2,2))[,,2]
colnames(plot.data) = c("Young", "Senior")
rownames(plot.data) = c("Low Risk", "High Risk")
heatmap.ART.overweight = myHeatmap(data = melt(plot.data), title = "ART; Overweight")

# DDRT
plot.data.merged = rbind(c(array(record$pValue$DDRT, dim = c(2,2,2))[,,1]), c(array(record$pValue$DDRT, dim = c(2,2,2))[,,2]))
colnames(plot.data.merged) = c("Young; Low Risk", "Senior; Low Risk", "Young; High Risk", "Senior; High Risk")
rownames(plot.data.merged) = c("BMI < 27", "BMI > 27")
# plot.data.merged = plot.data.merged[c(2,1),]
plot.data.merged = t(plot.data.merged)
heatmap.DDRT.merged = myHeatmap(data = melt(plot.data.merged), title = "DDRT")

plot.data = array(record$pValue$DDRT, dim = c(2,2,2))[,,1]
colnames(plot.data) = c("Young", "Senior")
rownames(plot.data) = c("Low Risk", "High Risk")
heatmap.DDRT.normal = myHeatmap(data = melt(plot.data), title = "DDRT; Normal")

plot.data = array(record$pValue$DDRT, dim = c(2,2,2))[,,2]
colnames(plot.data) = c("Young", "Senior")
rownames(plot.data) = c("Low Risk", "High Risk")
heatmap.DDRT.overweight = myHeatmap(data = melt(plot.data), title = "DDRT; Overweight")

heatmap.ART.merged
# Arrange the two heatmaps side by side
# pdf(file = paste(plotDirectory, "/", setting, "Heatmap.pdf", sep = ""), width = 7.5, height = 7)
# Arrange 1
# grid.arrange(heatmap.RT.normal,
#              heatmap.SSRT.normal,
#              heatmap.ART.normal,
#              heatmap.DDRT.normal,
#              heatmap.RT.overweight,
#              heatmap.SSRT.overweight,
#              heatmap.ART.overweight,
#              heatmap.DDRT.overweight, ncol = 4, nrow = 2)
# Arrange 2
# grid.arrange(heatmap.ART.normal, 
#              heatmap.ART.overweight, ncol = 2, nrow = 1)
# dev.off()


# pdf(file = paste(plotDirectory, "/", setting, "Heatmap.pdf", sep = ""), width = 7.5, height = 7)
# Arrange 1
# grid.arrange(heatmap.RT.merged,
#              heatmap.SSRT.merged,
#              heatmap.ART.merged,
#              heatmap.DDRT.merged, ncol = 2, nrow = 2)
# Arrange 2
grid.arrange(heatmap.RT.merged,
             heatmap.SSRT.merged,
             heatmap.ART.merged, ncol = 3, nrow = 1)
# dev.off()
