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

set.seed(1235) # 123
Group.level.number = c(4,2, 2)
Group.number = prod(Group.level.number)

S = cbind(as.numeric(cut(data$RISK10YRS, c(-Inf,11.66, 28.05, 35, Inf))),
          as.numeric(cut(data$HDL, c(-Inf,63.50, Inf))),
          as.numeric(cut(data$EGFR, c(-Inf,95.20, Inf))))

#S = cbind(as.numeric(cut(data$AGE, quantile(c(-Inf, Inf, data$AGE)))), 
#                  as.numeric(cut(data$BMI, quantile(c(-Inf, Inf, data$BMI)))))
#                                                    #, probs = c(0, 1/2, 1)))), 1 + as.numeric(data$FEMALE)) 

#S = cbind(as.numeric(cut(data$AGE, quantile(c(-Inf, Inf, data$AGE)))), 
#          as.numeric(cut(data$BMI, quantile(c(-Inf, Inf, data$BMI), probs = c(0, 1/2, 1)))), 
#          1 + as.numeric(data$FEMALE))
          
#(S[, 1] - 1) * max(Group.level.number) + S[, 2]

Group = (S[, 1] - 1) * Group.level.number[2] * Group.level.number[3] + (S[, 2] - 1) * Group.level.number[3] + S[, 3]; 

G = model.matrix(~ factor(Group))[, -1] # one-hot encoding of group membership
#X = data[, c("DBP", "SCREAT", "BMI", "GLUR")] # covariates
data$RACE4 <- factor( as.numeric(factor(data$RACE4)), levels = 1:4)
one_hot_feature <- model.matrix(~ RACE4 - 1, data)
colnames(one_hot_feature) <- paste("RACE", 1:4, sep="_")


X <- cbind(data[, setdiff(colnames(data), c("MASKID","INTENSIVE","EVENT_PRIMARY","T_PRIMARY"))], one_hot_feature)
W = data$INTENSIVE # treatment assignment

X <- as.data.frame(lapply(X, function(col) {
  if (is.factor(col) || is.character(col)) {
    return(as.numeric(as.factor(col)))
  } else {
    return(as.numeric(col))
  }
}))


# potential outcomes
Y = data$T_PRIMARY # data$EVENT_PRIMARY; data$T_PRIMARY
#X = data[, c("DBP", "SCREAT", "BMI", "GLUR")]
#X = data[,c("NEWSITEID","RISK10YRS","UMALCR","EGFR","TRR","BMI", "AGE","CHR","DBP","SBP", "GLUR","HDL")]



library(xgboost)
library(DiagrammeR)
n = length(Y)
p = 0.5
q = 0.2
test.stats.method = "AIPW"
nuisance.learner.method = "gradient boosting"
M = 2000
B = 6


W.tilde = apply(cbind(W,  matrix(rbinom(n * B, 1, p), ncol = B)), 1, mean)
nuisance.hat = nuisance.learner(Y = Y, X = X, prop = p, G = G, W = W.tilde, method = nuisance.learner.method, train.index = seq(1, n), test.index = seq(1, n))
#importance_matrix <- xgb.importance(model = nuisance.hat$tau)
#xgb.plot.importance(importance_matrix)

#library(DiagrammeR)
#gr <- xgb.plot.tree(model=nuisance.hat$tau, trees=0:1, render=FALSE)
#export_graph(gr, '~/Desktop/Research/Yao/HTE inference/code/Panning/0630/tree.pdf', width=3000, height=4000)


tree_dump <- xgb.model.dt.tree(model = nuisance.hat$tau)
print(tree_dump)

#specific_tree <- subset(tree_dump, Tree == 0)
#print(specific_tree)








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

#ART: augmented RT
record$pValue$ART = ART(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = M, test.stats.method = test.stats.method, nuisance.learner.method = nuisance.learner.method, B = B)$pval
record$R$ART = which(record$pValue$ART <= BH.threshold(pval = record$pValue$ART, q = q))

record

#saveRDS(record, "~/Desktop/Research/Yao/HTE inference/code/Panning/0630/SPRINT.rds")
