# Use FRT to test subgroup treatment effects
# Number of knockoffs
rm(list = ls())
source("~/Desktop/Research/Yao/HTE inference/code/Panning/helper.R")

n = 800 # sample size
d = 5 # number of covariates
p = 0.5 # propensity score
Group.level.number = c(4, 5) # number of levels per group
Group.number = prod(Group.level.number) # total number of groups
beta0 = 1; beta = rep(1, d) * 1; theta = rep(1, 2)
delta = 1.5
sigma = 1 # error magnitude in generating Y(0)

nuisance.learner.method = "gradient boosting"
test.stats.method = "AIPW + ITE" # "denoise", "ATE", "denoise + ATE", "AIPW", "ITE"; test statistic
B.seq = c(seq(6,60, by = 6)) # number of knockoffs: seq(1, 10); c(1, 5, 10)
M = 200 # number of permutations
q = 0.2 # FDR level

setting = "NumberKnockoff"
if(setting == "NumberKnockoff"){} 

start.time = proc.time()
m = 100 # number of trials
record = list()
record$pValue = list()
record$pValue$ORT = matrix(0, nrow = m, ncol = Group.number) 
record$pValue$ART = array(0, dim = c(m, Group.number, length(B.seq))) # ORT: oracle RT; RT (baseline): standard RT; SSRT: sample-splitting RT; DDRT: double-dipping RT; ART: augmented RT; 
record$R = list(); record$R$ORT = list(); record$R$ART = list(); for(j in seq(1, length(B.seq))){record$R$ART[[j]] = list()}
record$FDP = list(); record$FDP$ORT = rep(0, m)  
record$FDP$ART = matrix(0, m, length(B.seq)) 
record$power = record$FDP

set.seed(318)
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
  Y0 = beta0 + X %*% beta + S %*% theta + rnorm(n, 0, sigma)
  # tau is linear in S and independent of X
  tau = delta * (S[, 1] >= Group.level.number[1]-1) * (S[, 2] >= (Group.level.number[2] - 1)) * ((setting != "HTE") + (setting == "HTE") * X[, 1]) #  (S[, 1] >= 3) * (S[, 2] >= 3)
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
  
  # ART: augmented RT
  for(j in seq(1, length(B.seq))){
    B = B.seq[j]
    record$pValue$ART[i,,j] = ART(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = M, nuisance.learner.method = nuisance.learner.method, test.stats.method = test.stats.method, B = B)$pval
    record$R$ART[[j]][[i]] = which(record$pValue$ART[i,,j] <= BH.threshold(pval = record$pValue$ART[i,,j], q = q))
    record$FDP$ART[i,j] = sum(tau.group[record$R$ART[[j]][[i]]] == 0) / max(1, length(record$R$ART[[j]][[i]]))
    record$power$ART[i,j] = sum(tau.group[record$R$ART[[j]][[i]]] != 0) / max(1, sum(tau.group != 0))
  }
  
  if(i %% 10 == 0){print(i)}
}

# results
result = rbind(c(mean(record$FDP$ORT), apply(record$FDP$ART, 2, mean)),
               c(sd(record$FDP$ORT), apply(record$FDP$ART, 2, sd)),
               c(mean(record$power$ORT), apply(record$power$ART, 2, mean)),        c(sd(record$power$ORT), apply(record$power$ART, 2, sd)))
rownames(result) = c("FDR", "FDR sd", "power", "power sd")
colnames(result) = c("ORT", paste("B = ", B.seq, sep = ""))
print(result)
end.time = proc.time()
print(end.time[3] - start.time[3])

# saveRDS(record, file.path("~/Desktop/Research/Yao/HTE inference/code/Panning/0630", paste(setting, ".rds", sep= "")))

