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

# load data
baselineData = read.csv("/Users/zijungao/Desktop/Research/Trevor/validataionCausalEffect/realData/SPRINT/baseline.csv")
outcomeData = read.csv("/Users/zijungao/Desktop/Research/Trevor/validataionCausalEffect/realData/SPRINT/outcomes.csv")
baselineData$EVENT_PRIMARY = outcomeData$EVENT_PRIMARY
baselineData$T_PRIMARY = outcomeData$T_PRIMARY

# preprocessing
data = baselineData[which(apply(is.na(baselineData), 1, sum) == 0), ] # remove NA
data$RACE4 <- factor( as.numeric(factor(data$RACE4)), levels = 1:4)
one_hot_feature <- model.matrix(~ RACE4 - 1, data)
colnames(one_hot_feature) <- paste("RACE", 1:4, sep="_")
# X = cbind(data[, setdiff(colnames(data), c("NEWSITEID", "MASKID","INTENSIVE","T_PRIMARY","EVENT_PRIMARY"))],one_hot_feature) #
X = cbind(data[, c("RISK10YRS", "SBP", "DBP", "N_AGENTS", "SMOKE_3CAT", "ASPIRIN", "EGFR", "SCREAT", "AGE", "FEMALE", "CHR", "GLUR", "HDL", "TRR", "UMALCR", "BMI", "STATIN")], one_hot_feature[, dim(one_hot_feature)[2]-1]) # remove the last race group
X = as.data.frame(lapply(X, function(col) {
  if (is.factor(col) || is.character(col)) {
    return(as.numeric(as.factor(col)))
  } else {
    return(as.numeric(col))
  }
}))
W = data$INTENSIVE
# response
Y = 1 - data$EVENT_PRIMARY
# Y = data$T_PRIMARY #data$T_PRIMARY*(1-data$EVENT_PRIMARY); + data$T_PRIMARY
# Y = ((data$T_PRIMARY + 00)*(1-data$EVENT_PRIMARY) + data$EVENT_PRIMARY * data$T_PRIMARY * 0) #data$T_PRIMARY*(1-data$EVENT_PRIMARY); + data$T_PRIMARY
# Y = (data$T_PRIMARY>800)*data$T_PRIMARY # *(1-data$EVENT_PRIMARY) # ((data$T_PRIMARY + 500) *(1-data$EVENT_PRIMARY) + data$EVENT_PRIMARY * data$T_PRIMARY)
# Y = pmin(1000,data$T_PRIMARY)
# Y = pmax(800,data$T_PRIMARY)
# Y = log(Y+1) # log(Y+0.000001)
d = ncol(X)
n = nrow(X)

set.seed(318)
index = seq(1, n) # seq(1, n); sample(n, 4000)
X = X[index,]; Y = Y[index]; W = W[index]; n = length(index)
# Group.number = 1; G = rep(0, n); Group = rep(0, n) # a single group
Group.number = 8 # total number of groups
X$High.risk.group = (X$RISK10YRS > 18) # 18 
X$Senior.group = (X$AGE > 70) # 70
X$Obese.group = (X$BMI > 27) # 27
Group = X$Obese.group * 4 + X$Senior.group * 2 + X$High.risk.group
X$Senior.group = NULL; X$Obese.group = NULL; X$High.risk.group = NULL
table(Group)
X = as.matrix(X)
Ex = rep(0.5, n)

M = 2000 # number of permutations
proportion = 0.5 # proportion of randomness in the nuisance fold
test.stats.method ="AIPW" # test statistics

set.seed(318) # 318; 813; 123
# mu.hat = nuisance.mu(Y, X) # linear regression
mu.hat = nuisance.mu.xgboost(Y, X) # random forests

# Randomization test
set.seed(318) 
tau.hat.active = nuisance.tau.active(Y = Y, X = X, Ex = Ex, W = W, mu = mu.hat, Group = Group, proportion = proportion, robust = T) # estimate tau using active learning; robust = T, use power = -(TE^2)^2, and throw the point with the maximal power to the nuisance fold.
train.index.active = tau.hat.active$train.index # the data points fitted to the nuisance
test.index.active = setdiff(1:n,train.index.active) # the data points used for inference

# Randomization test
# Active Splitting
test.stats.value.active = test.stats.group(Y = Y[test.index.active], 
                                           W = W[test.index.active], 
                                           Group = Group[test.index.active], 
                                           stats = test.stats.method, Ex = Ex[test.index.active], 
                                           mu0.hat = tau.hat.active$mu0[test.index.active],
                                           mu1.hat = tau.hat.active$mu1[test.index.active], mu.hat = mu.hat[test.index.active], 
                                           tau.hat = tau.hat.active$tau[test.index.active])

test.stats.ref.active = t(replicate(M, test.stats.group(Y = Y[test.index.active],
                                                        W =  W[test.index.active], 
                                                        Group = Group[test.index.active], 
                                                        stats = test.stats.method, 
                                                        Ex = Ex[test.index.active], 
                                                        mu0.hat = tau.hat.active$mu0[test.index.active], 
                                                        mu1.hat = tau.hat.active$mu1[test.index.active], 
                                                        mu.hat = mu.hat[test.index.active], 
                                                        tau.hat = tau.hat.active$tau[test.index.active], 
                                                        test.index=test.index.active[test.index.active])))

if(dim(test.stats.ref.active)[1] == 1){test.stats.ref.active = t(test.stats.ref.active)}

# Calculate p-values for each group
p.value.active = sapply(seq(1, length(test.stats.value.active)), function(x) {sum(test.stats.ref.active[, x]>=test.stats.value.active[x])/(M+1) + 1/(M+1)})


set.seed(318)
m.pass = 1
p.value.ss = matrix(0, nrow = m.pass, ncol = Group.number)
for(i in 1 : m.pass){
  # Sample splitting
  train.index.ss = select.train(Group, proportion) # sample(n, proportion*n); sample(n, length(train.index.active))  
  test.index.ss = setdiff(1 : n, train.index.ss)
  
  # Fit the nuisance parameter
  tau.hat.ss = nuisance.tau.ss(Y = Y, X = X, Ex = Ex, W = W, mu = mu.hat, train.index = train.index.ss)
  
  # Randomization tests
  test.stats.value = test.stats.group(Y = Y[test.index.ss], 
                                      W = W[test.index.ss], 
                                       Group = Group[test.index.ss], 
                                      stats = test.stats.method, 
                                      Ex = Ex[test.index.ss], 
                                      mu0.hat = tau.hat.ss$mu0[test.index.ss], 
                                      mu1.hat = tau.hat.ss$mu1[test.index.ss], 
                                      mu.hat = mu.hat[test.index.ss], 
                                      tau.hat = tau.hat.ss$tau[test.index.ss])
  
  test.stats.ref = t(replicate(M, test.stats.group(Y = Y[test.index.ss], 
                                                   W = W[test.index.ss], 
                                                   Group = Group[test.index.ss], 
                                                   stats = test.stats.method, 
                                                   Ex = Ex[test.index.ss], 
                                                   mu0.hat = tau.hat.ss$mu0[test.index.ss], 
                                                   mu1.hat = tau.hat.ss$mu1[test.index.ss], 
                                                   mu.hat = mu.hat[test.index.ss], 
                                                   tau.hat = tau.hat.ss$tau[test.index.ss], test.index = test.index.ss))) 
  if(dim(test.stats.ref)[1] == 1){test.stats.ref = t(test.stats.ref)}
  
  # Calculate p-values for each group
  p.value.ss[i,] = sapply(seq(1, length(test.stats.value)), function(x) { sum(test.stats.ref[, x]>=test.stats.value[x])/(M+1) + 1/(M+1)})
  
  print(i)
}

# Double-dipping
set.seed(318)
train.index.dd = seq(1, n); test.index.dd = seq(1, n)

# Fit the nuisance parameter
tau.hat.dd = nuisance.tau.ss(Y = Y, X = X, Ex = Ex, W = W, mu = mu.hat, train.index = train.index.dd)

# Randomization tests
test.stats.value.dd = test.stats.group(Y = Y, 
                                    W = W, 
                                    Group = Group, 
                                    stats = test.stats.method, 
                                    Ex = Ex, 
                                    mu0.hat = tau.hat.dd$mu0, 
                                    mu1.hat = tau.hat.dd$mu1, 
                                    mu.hat = mu.hat, 
                                    tau.hat = tau.hat.dd$tau)

test.stats.ref.dd = t(replicate(M, test.stats.group(Y = Y, 
                                                 W = W, 
                                                 Group = Group, 
                                                 stats = test.stats.method, 
                                                 Ex = Ex, 
                                                 mu0.hat = tau.hat.dd$mu0, 
                                                 mu1.hat = tau.hat.dd$mu1, 
                                                 mu.hat = mu.hat, 
                                                 tau.hat = tau.hat.dd$tau, 
                                                 test.index = test.index.dd))) 
if(dim(test.stats.ref.dd)[1] == 1){test.stats.ref.dd = t(test.stats.ref.dd)}

# Calculate p-values for each group
p.value.dd = sapply(seq(1, length(test.stats.value.dd)), function(x) { sum(test.stats.ref.dd[, x]>=test.stats.value.dd[x])/(M+1) + 1/(M+1)})


# Baseline, not model-based
set.seed(318)
test.stats.value.baseline = test.stats.group(Y = Y, 
                                             W = W, 
                                             Group = Group, 
                                             stats = test.stats.method, 
                                             Ex = Ex, 
                                    
                                             mu0.hat = rep(0, n), 
                                             mu1.hat = rep(0, n), 
                                             mu.hat = rep(0, n), 
                                             tau.hat = rep(0, n))

test.stats.ref.baseline = t(replicate(M, test.stats.group(Y = Y, 
                                                          W = W, 
                                                          Group = Group, 
                                                          stats = test.stats.method, 
                                                          Ex = Ex, 
                                                 
                                                          mu0.hat = rep(0, n), 
                                                          mu1.hat = rep(0, n), 
                                                          mu.hat = rep(0, n), 
                                                 
                                                          tau.hat = rep(0, n), 
                                                          test.index = seq(1, n)))) 
if(dim(test.stats.ref.baseline)[1] == 1){test.stats.ref.baseline = t(test.stats.ref.baseline)}

# Calculate p-values for each group
p.value.baseline = (sapply(seq(1, length(test.stats.value.baseline)), function(x) { sum(test.stats.ref.baseline[, x]>=test.stats.value.baseline[x])/(M+1) + 1/(M+1)}))
  

# result
# X$High.risk.group = (X$RISK10YRS > 18) # 18 
# X$Senior.group = (X$AGE > 70) # 70
# X$Obese.group = (X$BMI > 27) # 27
# Group = X$Obese.group * 4 + X$Senior.group * 2 + X$High.risk.group
names(p.value.baseline) = names(p.value.dd) = names(p.value.active) = colnames(p.value.ss) = paste(rep(c("BMI < 27", "BMI > 27"),  times=c(4,4)), rep(c("junior", "senior"),  times=c(2,2)), c("low risk", "high risk"))
record = list(); record$pValue = list()
record$pValue$RT = p.value.baseline
record$pValue$SSRT = p.value.ss
record$pValue$DDRT = p.value.dd
record$pValue$ART = p.value.active

# FWER
q = 0.2 # 0.2
record$R$RT = sum(closing(p_val = record$pValue$RT, q = q, global.null.test = "Fisher"))
record$R$SSRT = sum(closing(p_val = record$pValue$SSRT[1,], q = q, global.null.test = "Fisher")) # use the first run
record$R$DDRT = sum(closing(p_val = record$pValue$DDRT, q = q, global.null.test = "Fisher"))
record$R$ART = sum(closing(p_val = record$pValue$ART, q = q, global.null.test = "Fisher"))

print(record)


# saveRDS(record, "~/Desktop/Research/Yao/HTE inference/code/Panning/July 2025/SPRINT.rds") # SPRINT.rds; SPRINT0421.rds

