# Use FRT to test subgroup treatment effects
# TODO: find a case where denoise > AIPW.normalized; test out the test statistic selection method
rm(list = ls())
source("~/Desktop/Research/Yao/HTE inference/code/Panning/helper ZG.R")

n = 400 # sample size
d = 5 # number of covariates
p = 0.5 # propensity score
Group.level.number = c(3, 3) # number of levels per group
Group.number = prod(Group.level.number) # total number of groups
beta0 = 1; beta = rep(1, d) * 1; theta = rep(1, 2)
delta = 1 # 1
sigma = 1 # error magnitude in generating Y(0)

test.stats.method = "denoise" # "denoise", "ATE", "denoise + ATE", "AIPW", "ITE"; test statistic
B = 3 # number of knockoffs
M = 400 # number of permutations
q = 0.2 # FDR level

# setting
setting = "default" # "default"

start.time = proc.time()
m = 40 # number of trials
record = list()
record$pValue = list()
record$pValue$ORT = record$pValue$RT = record$pValue$SSRT = record$pValue$DDRT = record$pValue$ART = matrix(0, nrow = m, ncol = Group.number) # ORT: oracle RT; RT (baseline): standard RT; SSRT: sample-splitting RT; DDRT: double-dipping RT; ART: augmented RT; 
record$R = list(); record$R$ORT = record$R$RT = record$R$SSRT = record$R$DDRT = record$R$ART = list()
record$FDP = list(); record$FDP$ORT = record$FDP$RT = record$FDP$SSRT = record$FDP$DDRT = record$FDP$ART = rep(0, m)
record$power = record$FDP

set.seed(318)
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
  tau = delta * (S[, 1] >= Group.level.number[1]) * (S[, 2] >= Group.level.number[2]) * (1 + X[, 1]) # delta * (S[, 1] >= Group.level.number[1]) * (S[, 2] >= Group.level.number[2])
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
  record$pValue$ORT[i,] = ORT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = test.stats.method, mu0 = mu0, mu1 = mu1, mu = mu, tau = tau, M = M)$pval
  record$R$ORT[[i]] = which(record$pValue$ORT[i,] <= BH.threshold(pval = record$pValue$ORT[i,], q = q))
  record$FDP$ORT[i] = sum(tau.group[record$R$ORT[[i]]] == 0) / max(1, length(record$R$ORT[[i]]))
  record$power$ORT[i] = sum(tau.group[record$R$ORT[[i]]] != 0) / max(1, sum(tau.group != 0))
  
  # RT (baseline): standard RT
  record$pValue$RT[i,] = RT(Y = Y, X = X, W = W, G = G, Group = Group, prop = p, M = M)$pval
  record$R$RT[[i]] = which(record$pValue$RT[i,] <= BH.threshold(pval = record$pValue$RT[i,], q = q))
  record$FDP$RT[i] = sum(tau.group[record$R$RT[[i]]] == 0) / max(1, length(record$R$RT[[i]]))
  record$power$RT[i] = sum(tau.group[record$R$RT[[i]]] != 0) / max(1, sum(tau.group != 0))
  
  # SSRT: sample-splitting RT  
  record$pValue$SSRT[i,] = SS(Y = Y, X = X, W = W, G = G, Group = Group, prop = p, M = M, test.stats.method = test.stats.method)$pval
  record$R$SSRT[[i]] = which(record$pValue$SSRT[i,] <= BH.threshold(pval = record$pValue$SSRT[i,], q = q))
  record$FDP$SSRT[i] = sum(tau.group[record$R$SSRT[[i]]] == 0) / max(1, length(record$R$SSRT[[i]]))
  record$power$SSRT[i] = sum(tau.group[record$R$SSRT[[i]]] != 0) / max(1, sum(tau.group != 0))
  
  # DDRT: double-dipping RT
  record$pValue$DDRT[i,] = DD(Y = Y, X = X, W = W, G = G, Group = Group, prop = p, M = M, test.stats.method = test.stats.method)$pval
  record$R$DDRT[[i]] = which(record$pValue$DDRT[i,] <= BH.threshold(pval = record$pValue$DDRT[i,], q = q))
  record$FDP$DDRT[i] = sum(tau.group[record$R$DDRT[[i]]] == 0) / max(1, length(record$R$DDRT[[i]]))
  record$power$DDRT[i] = sum(tau.group[record$R$DDRT[[i]]] != 0) / max(1, sum(tau.group != 0))
  
  # ART: augmented RT
  record$pValue$ART[i,] = ART(Y = Y, X = X, W = W, G = G, Group = Group, prop = p, M = M, test.stats.method = test.stats.method, B = B)$pval
  record$R$ART[[i]] = which(record$pValue$ART[i,] <= BH.threshold(pval = record$pValue$ART[i,], q = q))
  record$FDP$ART[i] = sum(tau.group[record$R$ART[[i]]] == 0) / max(1, length(record$R$ART[[i]]))
  record$power$ART[i] = sum(tau.group[record$R$ART[[i]]] != 0) / max(1, sum(tau.group != 0))
  
  if(i %% 10 == 0){print(i)}
}

# results
record$setting = setting
result = rbind(lapply(record$FDP, mean),
               lapply(record$FDP, sd),
               lapply(record$power, mean),
               lapply(record$power, sd))
row.names(result) = c("FDR", "FDR sd", "power", "power sd")
print(result)
end.time = proc.time()
print(end.time[3] - start.time[3])


# saveRDS(record, file.path("~/Desktop/Research/Yao/HTE inference/code/Panning/0630", paste(setting, "rds", sep= ".")))


