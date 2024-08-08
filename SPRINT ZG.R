# Create a hypothetical two-stage enrichment trial based on the SPRINT trial.
# Subgroup: age groups: [0, 59], [60, 69], [70, 79], [80, 100].
rm(list = ls())
library(ggplot2)
library(reshape2)
library(patchwork)
library(reshape2)

source("~/Desktop/Research/Yao/HTE inference/code/Panning/helper Yao.R")
#helper Yao.R

# TODO: 
# pull codes
# have better prediction model

# train causal trees

# load data
baselineData = read.csv("/Users/zijungao/Desktop/Research/Trevor/validataionCausalEffect/realData/SPRINT/baseline.csv")
outcomeData = read.csv("/Users/zijungao/Desktop/Research/Trevor/validataionCausalEffect/realData/SPRINT/outcomes.csv")


baselineData$EVENT_PRIMARY = outcomeData$EVENT_PRIMARY
baselineData$T_PRIMARY = outcomeData$T_PRIMARY

# preprocessing
# remove NA
data = baselineData[which(apply(is.na(baselineData), 1, sum) == 0), ]


seed = 318 # 318
set.seed(seed) 

data$RACE4 <- factor( as.numeric(factor(data$RACE4)), levels = 1:4)
one_hot_feature <- model.matrix(~ RACE4 - 1, data)
colnames(one_hot_feature) <- paste("RACE", 1:4, sep="_")


#index = (data[,"EVENT_PRIMARY"]==1)


Y = ((data$T_PRIMARY + 00)*(1-data$EVENT_PRIMARY) + data$EVENT_PRIMARY * data$T_PRIMARY * 0) #data$T_PRIMARY*(1-data$EVENT_PRIMARY); + data$T_PRIMARY 
#(data$T_PRIMARY>800)*data$T_PRIMARY #*(1-data$EVENT_PRIMARY) # ((data$T_PRIMARY + 500) *(1-data$EVENT_PRIMARY) + data$EVENT_PRIMARY * data$T_PRIMARY)

#Y = log(Y+0.000001)
W = data$INTENSIVE
data =cbind(data[, setdiff(colnames(data), c("MASKID","INTENSIVE","T_PRIMARY","EVENT_PRIMARY"))],one_hot_feature)  #
#numeric_columns <- sapply(data, is.numeric)
#df_scaled <- data
#df_scaled[, numeric_columns] <- scale(data[, numeric_columns])
#data <- df_scaled 


X <- data

d = ncol(X)
#bind(data[, setdiff(colnames(data), c("MASKID","INTENSIVE","T_PRIMARY","EVENT_PRIMARY"))],one_hot_feature) #

#col = c("RISK10YRS", "EGFR","BMI","UMALCR","TRR","CHR","GLUR","NEWSITEID","DBP","SBP","SCREAT","HDL","AGE")   
#X = X[col]

X <- as.data.frame(lapply(X, function(col) {
  if (is.factor(col) || is.character(col)) {
    return(as.numeric(as.factor(col)))
  } else {
    return(as.numeric(col))
  }
}))



# potential outcomes
# data$EVENT_PRIMARY; data$T_PRIMARY


library(xgboost)
library(DiagrammeR)
n = length(Y)
p = 0.5
q = 0.25
test.stats.method = "AIPW"
nuisance.learner.method = "random forest" 
M = 200
B = 6



W.knockoff = matrix(rbinom(n * B, 1, p), ncol = B) 
W.aug = cbind(W, W.knockoff)
W.tilde = apply(W.aug, 1, mean) 
nuisance.hat = nuisance.learner(Y = Y, X = X, prop = p, G = G, W = W.tilde, method = nuisance.learner.method, train.index = seq(1, n), test.index = seq(1, n))



Group.level.number = c(4,4) # c(5,5)
K = length(Group.level.number)

if (nuisance.learner.method == "random forest"){
  importance_scores <- importance(nuisance.hat$tau)
  top_k_features <- sort(importance_scores, decreasing = TRUE)[1:K]
  feature_list = rownames(data.frame(top_k_features))
  
}else if (nuisance.learner.method == "gradient boosting"){
  
  importance_matrix <- xgb.importance(model = nuisance.hat$tau)
  feature_list <- importance_matrix[1:K, ]$Feature
  
}
feature_list = c("UMALCR", "BMI")

cat("Top", K, "Important Features:\n",feature_list)




num_features = length(feature_list)

S <- matrix(nrow = n, ncol = num_features)

# Specify the number of quantiles to split

j = 0
for (feature in feature_list) {
  j = j + 1
  print(feature)
  num_grids = Group.level.number[j] + 1
  levels = seq(0, 1, length.out = num_grids)[-c(1,num_grids)]
  print(levels)
  S[, j] <- as.numeric(cut(X[, feature], c(-Inf, quantile(X[, feature],probs = levels), Inf)))
}


Group.number = prod(Group.level.number)
Group = S[, num_features]
for (i in 1:(num_features-1)) {
  Group = Group + (S[, i] - 1)*prod(Group.level.number[(i+1):num_features])
}

length(unique(Group))

G = model.matrix(~ factor(Group))[, -1] 


ate.org = c()
for  (g in sort(unique(Group))){
  index = (Group == g) *(W.tilde!=p)
  ate.org_g = mean(Y[(Group == g) & (W == 1)]) - mean(Y[(Group == g) & (W == 0)])
  ate.org = c(ate.org, ate.org_g)
}
print(ate.org)

ate = c()
for  (g in sort(unique(Group))){
  index = (Group == g) *(W.tilde!=p)
  ate_g = mean(nuisance.hat$tau.hat[index==1])
  ate = c(ate, ate_g)
}
print(ate)

sorted_numbers <- sort( ate, decreasing = TRUE)
sorted_indices <- order( ate, decreasing = TRUE)


#########################
feature_list = c("UMALCR", "BMI") # c('10-year CVD risk score','Urine Albumin/Creatinine ratio')


record = list()
record$pValue = list()
record$pValue$RT = record$pValue$SSRT = record$pValue$DDRT = record$pValue$ART = rep(0, ncol = Group.number) 




# Calculate test statistics for the original data
test.stats.value = test.stats.group(Y = Y, W = W, X = X, Group = Group, G = G, stats = test.stats.method, prop = p, mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat, W.tilde=W.tilde)

# Generate reference test statistics using permutations
test.stats.ref = t(replicate(M, test.stats.group(Y = Y, W = W.aug[cbind(seq(1, n), replicate(n, sample(1 + B, 1)))], X = X, G = G, 
                                                 Group = Group, stats = test.stats.method, prop = p, mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, 
                                                 mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat, W.tilde=W.tilde))) 

# Calculate p-values for each group
pval = sapply(seq(1, length(test.stats.value)), function(x) {
  permutation.p.value(stats = test.stats.value[x], stats.ref = test.stats.ref[, x])
})

record$pValue$ART = pval
record$R$ART = which(pval <= BH.threshold.Storey(pval = pval[sorted_indices], q = q))


# Other methods
# RT (baseline): standard RT
record$pValue$RT = RT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = M, seed =seed)$pval
record$R$RT = which(record$pValue$RT <= BH.threshold.Storey(pval = record$pValue$RT[sorted_indices], q = q))

# SSRT: sample-splitting RT  
#record$pValue$SSRT = SS(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = M, test.stats.method = test.stats.method, nuisance.learner.method = nuisance.learner.method, seed =seed)$pval
#record$R$SSRT = which(record$pValue$SSRT <= BH.threshold.Storey(pval = record$pValue$SSRT, q = q))

# DDRT: double-dipping RT
#record$pValue$DDRT = DD(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = M, test.stats.method = test.stats.method, nuisance.learner.method = nuisance.learner.method, seed =seed)$pval
#record$R$DDRT = which(record$pValue$DDRT <= BH.threshold.Storey(pval = record$pValue$DDRT[sorted_indices], q = q))



print(record)

# Load necessary libraries
library(ggplot2)
library(reshape2)

# Reshape data into a 4x4 matrix
data_matrix <- matrix(record$pValue$ART, nrow=4, byrow=TRUE)

# Convert the matrix to a data frame
df <- melt(data_matrix)

# Create the heatmap
ggplot(data = df, aes(Var2, Var1, fill = value)) + 
  geom_tile(color = "white") + 
  scale_fill_viridis_c() + 
  labs(title = "Heatmap of Given Data", x = "Column", y = "Row") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



#saveRDS(record, "~/Desktop/Research/Yao/HTE inference/code/Panning/0630/SPRINT.rds")




# Sort the p-values in each row from largest to smallest
matrix1 <- matrix(record$pValue$RT, nrow = Group.level.number[1], byrow = TRUE)
matrix2 <- matrix(record$pValue$ART, nrow = Group.level.number[1], byrow = TRUE)




# Melt the matrices for ggplot2
df1 <- melt(matrix1)
df2 <- melt(matrix2)

# Define a light color palette
lighter_colors <- colorRampPalette(c("lightcoral", "white", "azure2", "lightblue"))(200)

# Calculate the combined range of the p-values
combined_range <- range(c(df1$value, df2$value))


sorted_p1 = sort(record$pValue$RT)
df <- data.frame(
  Index = 1:length(sorted_p1),
  P_Value = sorted_p1
)


p1 <- ggplot(df1, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = lighter_colors, limits = combined_range) +
  geom_text(aes(label = round(value, 2)), color = "black", size = 5) +
  labs(
    title = "Group-wise P-values",
    x = feature_list[2],
    y = feature_list[1]
  ) +
  theme_minimal() +
  theme(
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.position = "none"
    # plot.margin = margin(t = 10, r = 60, b = 10, l = 10)
  )
p1 
#dev.off()


# Assume pval is already defined and sorted
lambda <- 0.5

# Estimate pi0 using Storey's method
pi0 <- min(1, mean(record$pValue$RT > lambda) / (1 - lambda))

# Sort p-values

# Create a data frame

# Add a column to indicate the color (red for the first 20, grey for the rest)
df$Color <- ifelse(df$Index <= length(record$R$RT), "brown1", "grey")

# Save the plot as a PNG file
#png(file = paste(plotDirectory, "/","BH_threshold.png", sep = ""), width =650, height =650, res = 150)

# Plot the data using ggplot2
p2<-ggplot(df, aes(x = Index, y = P_Value)) +
  geom_abline(slope = q/prod(Group.level.number)/pi0, intercept = 0, linetype = "dashed", color = "black", linewidth = 1.1) +  # Add BH threshold line
  geom_point(aes(color = Color), size = 4) +  # Use color based on the new column
  scale_color_identity() +  # Use colors specified in the Color column without a legend
  labs(
    title = expression("BH + Storey (" ~ hat(pi)[0] ~ "= )"),
    x = "Index",
    y = "Sorted P-value"
  ) +
  ylim(0, 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 17, hjust = 0.5), # Center the title
    # plot.title = element_text(size = 17, hjust = 0.5, margin = margin(b =10)), # Center the title
    axis.title = element_text(size = 17, color = "black"),
    # Increase axis label size and set color to black
    axis.text = element_text(size = 16, color = "black"),   # Increase axis tick label size and set color to black
    axis.ticks = element_line(color = "black"),             # Set axis ticks color to black
    legend.position = "none"
    # plot.margin = margin(t = 10, r = 10, b = 10, l = 60)
  )
#$annotate("text", x = max(df$Index) * 0.5, y = 0.1,
#         label = expression("BH + Storey (" ~ hat(pi)[0] ~ "= 0.72)"),
#         color = "black", size = 5, hjust = 0)





# Combine the plots with the legend on the left
combined_plot <- (p1 + p2) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(
    legend.position = "left",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.key.width = unit(0.25, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.title = element_blank(),  # Remove the legend title ("value")
    legend.text = element_text(size = 12)  # Increase legend text size
  )

# Open a PNG device
# png(file = paste(plotDirectory, "/","p_rt.png", sep = ""), width = 1800, height =  700, res = 150)

# Print the combined plot
print(combined_plot)

# Close the PNG device
# dev.off()









sorted_p2 = sort(record$pValue$ART)
df <- data.frame(
  Index = 1:length(sorted_p2),
  P_Value = sorted_p2
)


p1 <- ggplot(df2, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = lighter_colors, limits = combined_range) +
  geom_text(aes(label = round(value, 2)), color = "black", size = 5) +
  labs(
    title = "Group-wise P-values",
    x = feature_list[2],
    y = feature_list[1]
  ) +
  theme_minimal() +
  theme(
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.position = "none"
    # plot.margin = margin(t = 10, r = 60, b = 10, l = 10)
  )
#dev.off()



# Estimate pi0 using Storey's method
pi0 <- min(1, mean(record$pValue$ART > lambda) / (1 - lambda))

# Sort p-values
sorted_p2 <- sort(record$pValue$ART)

# Create a data frame


# Add a column to indicate the color (red for the first 20, grey for the rest)
df$Color <- ifelse(df$Index <= length(record$R$ART), "brown1", "grey")

# Save the plot as a PNG file
#png(file = paste(plotDirectory, "/","BH_threshold.png", sep = ""), width =650, height =650, res = 150)

# Plot the data using ggplot2
p2<-ggplot(df, aes(x = Index, y = P_Value)) +
  geom_abline(slope = q/prod(Group.level.number)/pi0, intercept = 0, linetype = "dashed", color = "black", linewidth = 1.1) +  # Add BH threshold line
  geom_point(aes(color = Color), size = 4) +  # Use color based on the new column
  scale_color_identity() +  # Use colors specified in the Color column without a legend
  labs(
    title = expression("BH + Storey (" ~ hat(pi)[0] ~ "= )"), #
    x = "Index",
    y = "Sorted P-value"
  ) +
  ylim(0, 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 17, hjust = 0.5), # Center the title
    # plot.title = element_text(size = 17, hjust = 0.5, margin = margin(b =10)), # Center the title
    axis.title = element_text(size = 17, color = "black"),
    # Increase axis label size and set color to black
    axis.text = element_text(size = 16, color = "black"),   # Increase axis tick label size and set color to black
    axis.ticks = element_line(color = "black"),             # Set axis ticks color to black
    legend.position = "none"
    # plot.margin = margin(t = 10, r = 10, b = 10, l = 60)
  )




# Combine the plots with the legend on the left
combined_plot <- (p1 + p2) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(
    legend.position = "left",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.key.width = unit(0.25, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.title = element_blank(),  # Remove the legend title ("value")
    legend.text = element_text(size = 12)  # Increase legend text size
  )

# Open a PNG device
# png(file = paste(plotDirectory, "/","p_art.png", sep = ""), width = 1800, height = 700, res = 150)

# Print the combined plot
print(combined_plot)

# Close the PNG device
# dev.off()
