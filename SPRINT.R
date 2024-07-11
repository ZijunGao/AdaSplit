# Create a hypothetical two-stage enrichment trial based on the SPRINT trial.
# Subgroup: age groups: [0, 59], [60, 69], [70, 79], [80, 100].
rm(list = ls())
source("~/Desktop/Research/Yao/HTE inference/code/Panning/helper.R")

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

set.seed(318) # 123
Group.level.number = c(4, 2, 2)
Group.number = prod(Group.level.number)
S = cbind(as.numeric(cut(data$AGE, quantile(c(-Inf, Inf, data$AGE)))), 
          as.numeric(cut(data$BMI, quantile(c(-Inf, Inf, data$BMI), probs = c(0, 1/2, 1)))),
          1 + as.numeric(data$FEMALE)) 
Group = (S[, 1] - 1) * Group.level.number[2] * Group.level.number[3] + (S[, 2] - 1) * Group.level.number[3] + S[, 3]; 
G = model.matrix(~ factor(Group))[, -1] # one-hot encoding of group membership
X = data[, c("DBP", "SCREAT", "BMI", "GLUR")] # covariates
W = data$INTENSIVE # treatment assignment
# potential outcomes
Y = data$T_PRIMARY # data$EVENT_PRIMARY; data$T_PRIMARY

p = 0.5
q = 0.2
test.stats.method = "AIPW"
nuisance.learner.method = "gradient boosting"
M = 2000
B = 6

record = list()
record$pValue = list()
record$pValue$RT = record$pValue$SSRT = record$pValue$DDRT = record$pValue$ART = rep(0, ncol = Group.number) 

# Testing
# RT (baseline): standard RT
record$pValue$RT = RT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = M)$pval
record$R$RT = which(record$pValue$RT <= BH.threshold(pval = record$pValue$RT, q = q))

# SSRT: sample-splitting RT  
record$pValue$SSRT = SS(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = M, test.stats.method = test.stats.method, nuisance.learner.method = nuisance.learner.method)$pval
record$R$SSRT = which(record$pValue$SSRT <= BH.threshold(pval = record$pValue$SSRT, q = q))

# DDRT: double-dipping RT
record$pValue$DDRT = DD(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = M, test.stats.method = test.stats.method, nuisance.learner.method = nuisance.learner.method)$pval
record$R$DDRT = which(record$pValue$DDRT <= BH.threshold(pval = record$pValue$DDRT, q = q))

# ART: augmented RT
record$pValue$ART = ART(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = M, test.stats.method = test.stats.method, nuisance.learner.method = nuisance.learner.method, B = B)$pval
record$R$ART = which(record$pValue$ART <= BH.threshold(pval = record$pValue$ART, q = q))

record

# saveRDS(record, "~/Desktop/Research/Yao/HTE inference/code/Panning/0630/SPRINT.rds")
