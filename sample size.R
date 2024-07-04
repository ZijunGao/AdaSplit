# Use FRT to test subgroup treatment effects
# Sample size
rm(list = ls())
source("~/Desktop/Research/Yao/HTE inference/code/Panning/helper.R")

n.seq = seq(300, 1000, by = 100) # sample size
d = 5 # number of covariates
p = 0.5 # propensity score
Group.level.number = c(4, 4) # number of levels per group
Group.number = prod(Group.level.number) # total number of groups
beta0 = 1; beta = rep(1, d) * 1; theta = rep(1, 2)
delta = 1
sigma = 1 # error magnitude in generating Y(0)

test.stats.method = "AIPW + ITE" # "denoise", "ATE", "denoise + ATE", "AIPW", "ITE"; test statistic
B = 5 # number of knockoffs: seq(1, 10); c(1, 5, 10)
M = 400 # number of permutations
q = 0.2 # FDR level

setting = "SampleSize"

start.time = proc.time()
m = 400 # number of trials
record = list()
record$pValue = list()
record$pValue$ORT = record$pValue$RT = record$pValue$SSRT = record$pValue$DDRT = record$pValue$ART = array(0, dim = c(m, Group.number, length(n.seq))) # ORT: oracle RT; RT (baseline): standard RT; SSRT: sample-splitting RT; DDRT: double-dipping RT; ART: augmented RT; 
record$R = list(); record$R$ART = list(); for(j in seq(1, length(n.seq))){record$R$ART[[j]] = list()}; record$R$ORT = record$R$RT = record$R$SSRT = record$R$DDRT = record$R$ART
record$FDP = list();record$FDP$ORT = record$FDP$RT = record$FDP$SSRT = record$FDP$DDRT = record$FDP$ART = matrix(0, m, length(n.seq)) 
record$power = record$FDP

set.seed(318)
for(j in 1 : length(n.seq)){
  n = n.seq[j]
  for(i in 1:m){
    # generate data
    S = cbind(apply(rmultinom(n, 1, rep(1, Group.level.number[1])), 2, which.max),
              apply(rmultinom(n, 1, rep(1, Group.level.number[2])), 2, which.max)) # S used to define subgroups
    Group = (S[, 1] - 1) * Group.level.number[2] + S[, 2]; 
    G = model.matrix(~ factor(Group))[, -1] # one-hot encoding of group membership
    X = matrix(rnorm(n * d, 0, 1), nrow = n, ncol = d) # covariates
    W = rbinom(n, 1, p) # treatment assignment
    # potential outcomes
    # mu0 is linear in X and S
    Y0 = beta0 + X %*% beta + S %*% theta + rnorm(n, 0, sigma)
    # tau is linear in S and independent of X
    tau = delta * (S[, 1] >= Group.level.number[1]) * (S[, 2] >= (Group.level.number[2] - 1)) * ((setting != "HTE") + (setting == "HTE") * X[, 1]) #  (S[, 1] >= 3) * (S[, 2] >= 3)
    tau.group = sapply(seq(1, Group.number), function(x) {
      mean(tau[Group == x])
    }) # average treatment effect in each group.
    Y1 = Y0 + tau
    Y = Y1 * W + Y0 * (1 - W) # observed outcome
    # nuisance functions
    mu0 =  beta0 + X %*% beta + S %*% theta
    mu1 =  beta0 + X %*% beta + S %*% theta + tau
    mu = mu0 * (1 - p) + mu1 * p
    
    # inference
    # ORT (oracle): RT with true nuisance functions
    record$pValue$ORT[i,,j] = ORT(Y = Y, X = X, G = G, Group = Group, prop = p, test.stats.method = test.stats.method, mu0 = mu0, mu1 = mu1, mu = mu, tau = tau, M = M)$pval
    record$R$ORT[[j]][[i]] = which(record$pValue$ORT[i,,j] <= BH.threshold(pval = record$pValue$ORT[i,,j], q = q))
    record$FDP$ORT[i, j] = sum(tau.group[record$R$ORT[[j]][[i]]] == 0) / max(1, length(record$R$ORT[[j]][[i]]))
    record$power$ORT[i, j] = sum(tau.group[record$R$ORT[[j]][[i]]] != 0) / max(1, sum(tau.group != 0))
    
    # RT (baseline): standard RT
    record$pValue$RT[i,,j] = RT(Y = Y, X = X, G = G, Group = Group, prop = p, M = M)$pval
    record$R$RT[[j]][[i]] = which(record$pValue$RT[i,,j] <= BH.threshold(pval = record$pValue$RT[i,,j], q = q))
    record$FDP$RT[i, j] = sum(tau.group[record$R$RT[[j]][[i]]] == 0) / max(1, length(record$R$RT[[j]][[i]]))
    record$power$RT[i, j] = sum(tau.group[record$R$RT[[j]][[i]]] != 0) / max(1, sum(tau.group != 0))
    
    # SSRT: sample-splitting RT  
    record$pValue$SSRT[i,,j] = SS(Y = Y, X = X, G = G, Group = Group, prop = p, M = M, test.stats.method = test.stats.method)$pval
    record$R$SSRT[[j]][[i]] = which(record$pValue$SSRT[i,,j] <= BH.threshold(pval = record$pValue$SSRT[i,,j], q = q))
    record$FDP$SSRT[i, j] = sum(tau.group[record$R$SSRT[[j]][[i]]] == 0) / max(1, length(record$R$SSRT[[j]][[i]]))
    record$power$SSRT[i,j] = sum(tau.group[record$R$SSRT[[j]][[i]]] != 0) / max(1, sum(tau.group != 0))
    
    # DDRT: double-dipping RT
    record$pValue$DDRT[i,,j] = DD(Y = Y, X = X, G = G, Group = Group, prop = p, M = M, test.stats.method = test.stats.method)$pval
    record$R$DDRT[[j]][[i]] = which(record$pValue$DDRT[i,,j] <= BH.threshold(pval = record$pValue$DDRT[i,,j], q = q))
    record$FDP$DDRT[i, j] = sum(tau.group[record$R$DDRT[[j]][[i]]] == 0) / max(1, length(record$R$DDRT[[j]][[i]]))
    record$power$DDRT[i, j] = sum(tau.group[record$R$DDRT[[j]][[i]]] != 0) / max(1, sum(tau.group != 0))
    
    # ART: augmented RT
    record$pValue$ART[i,,j] = ART(Y = Y, X = X, G = G, Group = Group, prop = p, M = M, test.stats.method = test.stats.method, B = B)$pval
    record$R$ART[[j]][[i]] = which(record$pValue$ART[i,,j] <= BH.threshold(pval = record$pValue$ART[i,,j], q = q))
    record$FDP$ART[i, j] = sum(tau.group[record$R$ART[[j]][[i]]] == 0) / max(1, length(record$R$ART[[j]][[i]]))
    record$power$ART[i, j] = sum(tau.group[record$R$ART[[j]][[i]]] != 0) / max(1, sum(tau.group != 0))
  }
  print(j)
}

# results
result = list()
result$FDR = data.frame(lapply(record$FDP, function(x){apply(x, 2, mean)})); rownames(result$FDR) = n.seq
result$power = data.frame(lapply(record$power, function(x){apply(x, 2, mean)})); rownames(result$power) = n.seq
result


# saveRDS(record, file.path("~/Desktop/Research/Yao/HTE inference/code/Panning/0630", paste(setting, ".rds", sep= "")))

