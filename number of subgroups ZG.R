# Use FRT to test subgroup treatment effects
# Number of subgroups
rm(list = ls())
source("~/Desktop/Research/Yao/HTE inference/code/Panning/helper.R")

n = 800 # sample size
d = 5 # number of covariates
p = 0.4 # propensity score
Group.level.number = c(4, 4) # number of levels per group
Group.level.number2.seq = seq(4, 8)
beta0 = 1; beta = rep(1, d) * 1; theta = rep(1, 2)
delta = 1.5
sigma = 1 # error magnitude in generating Y(0)
nuisance.learner.method = "gradient boosting"
test.stats.method = "AIPW + ITE" # "denoise", "ATE", "denoise + ATE", "AIPW", "ITE"; test statistic
B = 6 # number of knockoffs: seq(1, 10); c(1, 5, 10)
M = 200 # number of permutations
q = 0.2 # FDR level

setting = "NumberSubgroup"

start.time = proc.time()
m = 100 # number of trials
record = list()
record$pValue = list()
record$pValue$ORT = record$pValue$RT = record$pValue$SSRT = record$pValue$DDRT = record$pValue$ART = array(0, dim = c(m, Group.level.number[1] * max(Group.level.number2.seq), length(Group.level.number2.seq))) # ORT: oracle RT; RT (baseline): standard RT; SSRT: sample-splitting RT; DDRT: double-dipping RT; ART: augmented RT; 
record$R = list(); record$R$ART = list(); for(j in seq(1, length(Group.level.number2.seq))){record$R$ART[[j]] = list()}; record$R$ORT = record$R$RT = record$R$SSRT = record$R$DDRT = record$R$ART
record$FDP = list();record$FDP$ORT = record$FDP$RT = record$FDP$SSRT = record$FDP$DDRT = record$FDP$ART = matrix(0, m, length(Group.level.number2.seq)) 
record$power = record$FDP

set.seed(318)
for(j in 1 : length(Group.level.number2.seq)){
  Group.level.number[2] = Group.level.number2.seq[j]
  Group.number = prod(Group.level.number) # total number of groups
  for(i in 1:m){
    # generate data
    X = matrix(rnorm(n * d, 0, 1), nrow = n, ncol = d) # covariates
    
    quantiles_1 <- quantile(c(-Inf,X[,d-1],Inf), probs = seq(0, 1, length.out = Group.level.number[1] + 1))
    quantiles_2 <- quantile(c(-Inf,X[,d-2],Inf), probs = seq(0, 1, length.out = Group.level.number[2] + 1))
    S = cbind(as.numeric(cut(X[,d-1], quantiles_1)), as.numeric(cut(X[,d-2], quantiles_2)))
    
    Group = (S[, 1] - 1) * Group.level.number[2] + S[, 2]; 
    G = model.matrix(~ factor(Group))[, -1] # one-hot encoding of group membership
    
    W = rbinom(n, 1, p) # treatment assignment
    
    # potential outcomes
    # mu0 is linear in X and S
    mu0 =  beta0 + X %*% beta + 2* X[, 1]^2 
    Y0 = mu0 + rnorm(n, 0, sigma)
    # tau is linear in S and independent of X #(setting == "HTE") * 
    tau = delta * (S[, 1] >= (Group.level.number[1]-1)) * (S[, 2] >= (Group.level.number[2] - 1)) * (1+ X[, 2]^2)
    tau.group = sapply(seq(1, Group.number), function(x) {
      mean(tau[Group == x])
    }) # average treatment effect in each group.
    Y1 = Y0 + tau
    Y = Y1 * W + Y0 * (1 - W) # observed outcome
    # nuisance functions
    mu1 = mu0 + tau
    mu = mu0 * (1 - p) + mu1 * p
    
    
    # inference
    # ORT (oracle): RT with true nuisance functions
    record$pValue$ORT[i,seq(1, Group.number),j] = ORT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = test.stats.method, mu0 = mu0, mu1 = mu1, mu = mu, tau = tau, M = M)$pval
    record$R$ORT[[j]][[i]] = which(record$pValue$ORT[i,seq(1, Group.number),j] <= BH.threshold(pval = record$pValue$ORT[i,seq(1, Group.number),j], q = q))
    record$FDP$ORT[i, j] = sum(tau.group[record$R$ORT[[j]][[i]]] == 0) / max(1, length(record$R$ORT[[j]][[i]]))
    record$power$ORT[i, j] = sum(tau.group[record$R$ORT[[j]][[i]]] != 0) / max(1, sum(tau.group != 0))
    
    # RT (baseline): standard RT
    record$pValue$RT[i,seq(1, Group.number),j] = RT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = M)$pval
    record$R$RT[[j]][[i]] = which(record$pValue$RT[i,seq(1, Group.number),j] <= BH.threshold(pval = record$pValue$RT[i,seq(1, Group.number),j], q = q))
    record$FDP$RT[i, j] = sum(tau.group[record$R$RT[[j]][[i]]] == 0) / max(1, length(record$R$RT[[j]][[i]]))
    record$power$RT[i, j] = sum(tau.group[record$R$RT[[j]][[i]]] != 0) / max(1, sum(tau.group != 0))
    
    # SSRT: sample-splitting RT  
    record$pValue$SSRT[i,seq(1, Group.number),j] = SS(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = M,nuisance.learner.method = nuisance.learner.method, test.stats.method = test.stats.method)$pval
    record$R$SSRT[[j]][[i]] = which(record$pValue$SSRT[i,seq(1, Group.number),j] <= BH.threshold(pval = record$pValue$SSRT[i,seq(1, Group.number),j], q = q))
    record$FDP$SSRT[i, j] = sum(tau.group[record$R$SSRT[[j]][[i]]] == 0) / max(1, length(record$R$SSRT[[j]][[i]]))
    record$power$SSRT[i,j] = sum(tau.group[record$R$SSRT[[j]][[i]]] != 0) / max(1, sum(tau.group != 0))
    
    # DDRT: double-dipping RT
    record$pValue$DDRT[i,seq(1, Group.number),j] = DD(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = M,nuisance.learner.method = nuisance.learner.method, test.stats.method = test.stats.method)$pval
    record$R$DDRT[[j]][[i]] = which(record$pValue$DDRT[i,seq(1, Group.number),j] <= BH.threshold(pval = record$pValue$DDRT[i,seq(1, Group.number),j], q = q))
    record$FDP$DDRT[i, j] = sum(tau.group[record$R$DDRT[[j]][[i]]] == 0) / max(1, length(record$R$DDRT[[j]][[i]]))
    record$power$DDRT[i, j] = sum(tau.group[record$R$DDRT[[j]][[i]]] != 0) / max(1, sum(tau.group != 0))
    
    # ART: augmented RT
    record$pValue$ART[i,seq(1, Group.number),j] = ART(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = M,nuisance.learner.method = nuisance.learner.method, test.stats.method = test.stats.method, B = B)$pval
    record$R$ART[[j]][[i]] = which(record$pValue$ART[i,seq(1, Group.number),j] <= BH.threshold(pval = record$pValue$ART[i,seq(1, Group.number),j], q = q))
    record$FDP$ART[i, j] = sum(tau.group[record$R$ART[[j]][[i]]] == 0) / max(1, length(record$R$ART[[j]][[i]]))
    record$power$ART[i, j] = sum(tau.group[record$R$ART[[j]][[i]]] != 0) / max(1, sum(tau.group != 0))
  }
  print(j)
}

# results
result = list()
result$FDR = data.frame(lapply(record$FDP, function(x){apply(x, 2, mean)})); rownames(result$FDR) = Group.level.number2.seq
result$power = data.frame(lapply(record$power, function(x){apply(x, 2, mean)})); rownames(result$power) = Group.level.number2.seq
result


# saveRDS(record, file.path("~/Desktop/Research/Yao/HTE inference/code/Panning/0630", paste(setting, ".rds", sep= "")))

