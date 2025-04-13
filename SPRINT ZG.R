# Create a hypothetical two-stage enrichment trial based on the SPRINT trial.
# Subgroup: age groups: [0, 59], [60, 69], [70, 79], [80, 100].

# Not making much sense
# Random

# Not linear.
rm(list = ls())
library(ggplot2)
library(reshape2)
library(patchwork)
library(reshape2)
library(xgboost)
library(DiagrammeR)


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
data$RACE4 <- factor( as.numeric(factor(data$RACE4)), levels = 1:4)
one_hot_feature <- model.matrix(~ RACE4 - 1, data)
colnames(one_hot_feature) <- paste("RACE", 1:4, sep="_")

# X = cbind(data[, setdiff(colnames(data), c("NEWSITEID", "MASKID","INTENSIVE","T_PRIMARY","EVENT_PRIMARY"))],one_hot_feature) #
X = cbind(data[, c("SBP", "DBP", "N_AGENTS", "SMOKE_3CAT", "ASPIRIN", "EGFR", "SCREAT", "AGE", "FEMALE", "CHR", "GLUR", "HDL", "TRR", "UMALCR", "BMI", "STATIN")], one_hot_feature[, dim(one_hot_feature)[2]-1])
X <- as.data.frame(lapply(X, function(col) {
  if (is.factor(col) || is.character(col)) {
    return(as.numeric(as.factor(col)))
  } else {
    return(as.numeric(col))
  }
}))
X = as.matrix(X)
W = data$INTENSIVE
# Y = 1 - data$EVENT_PRIMARY
# Y = data$T_PRIMARY #data$T_PRIMARY*(1-data$EVENT_PRIMARY); + data$T_PRIMARY
# Y = ((data$T_PRIMARY + 00)*(1-data$EVENT_PRIMARY) + data$EVENT_PRIMARY * data$T_PRIMARY * 0) #data$T_PRIMARY*(1-data$EVENT_PRIMARY); + data$T_PRIMARY
# Y = (data$T_PRIMARY>800)*data$T_PRIMARY # *(1-data$EVENT_PRIMARY) # ((data$T_PRIMARY + 500) *(1-data$EVENT_PRIMARY) + data$EVENT_PRIMARY * data$T_PRIMARY)
# Y = pmin(1000,data$T_PRIMARY)
Y = pmax(800,data$T_PRIMARY)
Y = log(Y+1) # log(Y+0.000001)
d = ncol(X)
n = nrow(X)
set.seed(318); random.index = sample(n, 2000); X = X[random.index,]; Y = Y[random.index]; W = W[random.index]; n = length(random.index)
Group.number = 1 # total number of groups
G = rep(0, n)
Group = rep(0, n)
Ex = rep(0.5, n)

M = 1000 # number of permutations
proportion = 0.5 # proportion of randomness in the nuisance fold
test.stats.method ="AIPW" #Test statistics

set.seed(318) 
mu.hat = nuisance.mu(Y,X)

tau.hat.active = nuisance.tau.active(Y = Y, X = X, Ex = Ex, W = W, mu = mu.hat, Group = Group, proportion = proportion) #estimate tau using active learning
train.index.active = tau.hat.active$train.index #the data points fitted to the nuisance
test.index.active = setdiff(1:n,train.index.active) #the data points used for inference


# Randomization tests
test.stats.value.active = test.stats.group(Y = Y[test.index.active], W = W[test.index.active], Group = Group[test.index.active], stats = test.stats.method, Ex = Ex[test.index.active], 
                                           mu0.hat = tau.hat.active$mu0[test.index.active], mu1.hat = tau.hat.active$mu1[test.index.active], mu.hat = mu.hat[test.index.active], 
                                           tau.hat = tau.hat.active$tau[test.index.active])

test.stats.ref.active = t(replicate(M, test.stats.group(Y = Y[test.index.active], W =  W[test.index.active], Group = Group[test.index.active], stats = test.stats.method, Ex = Ex[test.index.active], 
                                                        mu0.hat = tau.hat.active$mu0[test.index.active], mu1.hat = tau.hat.active$mu1[test.index.active], mu.hat = mu.hat[test.index.active], 
                                                        tau.hat = tau.hat.active$tau[test.index.active], test.index=test.index.active[test.index.active])))

if(dim(test.stats.ref.active)[1] == 1){test.stats.ref.active = t(test.stats.ref.active)}

# Calculate p-values for each group
p.value.active = (sapply(seq(1, length(test.stats.value.active)), function(x) { sum(test.stats.ref.active[, x]>=test.stats.value.active[x])/(M+1) + 1/(M+1)}))
#simulation_results_art[j,]


# Sample splitting
# set.seed(12) # 123: select 4000 units
train.index.ss = select.train(Group, proportion)  #sample(n, proportion*n) #sample(n, length(train.index.active))  
#sample(n, proportion*n)  #select.train(Group, proportion) 
test.index.ss = setdiff(1:n, train.index.ss)

# Fit the nuisance parameter
tau.hat.ss = nuisance.tau.ss(Y = Y, X = X, Ex = Ex, W = W, mu = mu.hat, train.index = train.index.ss)

# Randomization tests
test.stats.value = test.stats.group(Y = Y[test.index.ss], W = W[test.index.ss], Group = Group[test.index.ss], stats = test.stats.method, Ex = Ex[test.index.ss], 
                                    mu0.hat = tau.hat.ss$mu0[test.index.ss], mu1.hat = tau.hat.ss$mu1[test.index.ss], mu.hat = mu.hat[test.index.ss], tau.hat = tau.hat.ss$tau[test.index.ss])

test.stats.ref = t(replicate(M, test.stats.group(Y = Y[test.index.ss], W = W[test.index.ss], Group = Group[test.index.ss], stats = test.stats.method, Ex = Ex[test.index.ss], 
                                                 mu0.hat = tau.hat.ss$mu0[test.index.ss], mu1.hat = tau.hat.ss$mu1[test.index.ss], mu.hat = mu.hat[test.index.ss], 
                                                 tau.hat = tau.hat.ss$tau[test.index.ss], test.index = test.index.ss))) 
if(dim(test.stats.ref)[1] == 1){test.stats.ref = t(test.stats.ref)}

# Calculate p-values for each group
p.value.ss = (sapply(seq(1, length(test.stats.value)), function(x) { sum(test.stats.ref[, x]>=test.stats.value[x])/(M+1) + 1/(M+1)}))

# result
print("active splitting"); print(p.value.active)
print("sample splitting"); print(p.value.ss)
