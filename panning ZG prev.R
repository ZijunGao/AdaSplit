# Deprecated at June 30
# Use FRT to test subgroup treatment effects

# helper functions
BH.threshold = function(pval, q = 0.1){
  # preprocess
  m = length(pval)
  pval = sort(pval, decreasing = FALSE)
  
  # compute the threshold of the BH procedure
  FDR.hat = m * pval/seq(1, m)
  pval.index = which(FDR.hat <= q)
  if(length(pval.index) == 0){return(-1e6)}
  threshold = pval[max(pval.index)]
  return(threshold)
}

# observations
# The smaller the sample size, the more power gained by ART.
# n = c(200, 400, 800)
# The smaller magnitude of HTE, the more power gained by ART.
# delta = c(1/4, 1/3, 1/2)
# The smaller magnitude of baseline, the more power gained by ART.
# hyperparameters

n = 800 # sample size
d = 5 # number of covariates
p = 0.5 # propensity score
M = 100 # number of repetitions
q = 0.2 # FDR level
Group.number = 9
sigma = 0.1 

m = 10 # number of trials
record = list()
record$pValue = list()
record$pValue$RT = record$pValue$SSRT = record$pValue$ART = matrix(0, nrow = m, ncol = Group.number) # RT (baseline): standard RT; SSRT: sample-splitting RT; ART: augmented RT
record$R = list(); record$R$RT = record$R$SSRT = record$R$ART = list()
record$FDP = list(); record$FDP$RT = record$FDP$SSRT = record$FDP$ART = rep(0, m) 
record$power = record$FDP

set.seed(318)
for(i in 1:m){
  # generate data
  S = cbind(apply(rmultinom(n, 1, rep(1, 3)), 2, which.max),
            apply(rmultinom(n, 1, rep(1, 3)), 2, which.max)) # S used to define subgroups
  Group = (S[,1] - 1) * 3 + S[,2]; 
  G = model.matrix(~ factor(Group))[,-1] # one-hot encoding of Group
  Group.number = length(unique(Group))
  X = matrix(rnorm(n * d, 0, 1), nrow = n, ncol = d) # covariates
  W = rbinom(n, 1, p) # treatment assignment
  # potential outcomes
  # mu0 is linear in X and S
  # tau is linear in S and independent of X
  beta0 = 1; beta = rep(1, d) * 1; theta = rep(1, 2); delta = 1/2 # 1/3
  Y0 = beta0 + X %*% beta + S %*% theta + rnorm(n, 0, sigma)
  tau = delta * (S[, 1] > 2) * (S[, 2] > 2)
  tau.group = sapply(seq(1, Group.number), function(x) {
    mean(tau[Group == x])
  }) # average treatment effect in each group; previously (S[,1] - S[, 2])
  Y1 = Y0 + tau
  Y = Y1 * W + Y0 * (1 - W) # observed outcome
  
  # inference
  # RT (baseline): standard RT
  T.RT = rep(0, Group.number)
  for (k in 1: Group.number){
    T.RT[k] = mean(Y[(W == 1) & (Group == k)]) - mean(Y[(W == 0) & (Group == k)])
  }
  T.RT.ref = matrix(0, M, Group.number)
  for(j in 1 : M){
    W.ref = rbinom(n, 1, p) 
    for (k in 1 : Group.number){
      T.RT.ref[j, k] = mean(Y[(W.ref == 1) & (Group == k)]) - mean(Y[(W.ref == 0) & (Group == k)])
    }
  }
  record$pValue$RT[i,] = sapply(seq(1, Group.number), function(x){return(permutation.p.value(abs(T.RT[x]), abs(T.RT.ref[,x])))}) 
  record$R$RT[[i]] = which(record$pValue$RT[i,] <= BH.threshold(pval = record$pValue$RT[i,], q = q))
  record$FDP$RT[i] = sum(tau.group[record$R$RT[[i]]] == 0) / max(1, length(record$R$RT[[i]]))
  record$power$RT[i] = sum(tau.group[record$R$RT[[i]]] != 0) / sum(tau.group != 0)
  
  # SSRT: sample-splitting RT  
  # no further sample splitting
  data.train = data.frame(Y, X, G, W - 0.5, (W-0.5) * X, (W-0.5) * G)
  data = data.frame(Y, X, G, 0 - 0.5, (0 - 0.5) * X, (0 - 0.5) * G)
  colnames(data) = colnames(data.train)
  nuisance.index = sample(n, n/2)
  nuisance.model = lm(Y ~ ., data = data.train[nuisance.index,])
  mu0.hat = predict(nuisance.model, newdata = data[-nuisance.index, ])
  
  # nuisance function estimator fixed
  T.SSRT = rep(0, Group.number)
  for (k in 1: Group.number){
    T.SSRT[k] = mean((Y[-nuisance.index] - mu0.hat)[(W[-nuisance.index] == 1) & (Group[-nuisance.index] == k)]) -  mean((Y[-nuisance.index] - mu0.hat)[(W[-nuisance.index] == 0) & (Group[-nuisance.index] == k)]) # T.SSRT[k] may be NA
  }
  T.SSRT.ref = matrix(0, M, Group.number)
  for(j in 1 : M){
    W.ref = rbinom(n, 1, p) 
    for (k in 1: Group.number){
      data.ref = data.frame(Y, X, G, 0 - 0.5, (0 - 0.5) * X, (0 - 0.5) * G)
      colnames(data.ref) = colnames(data)
      mu0.hat.ref = predict(nuisance.model, newdata = data.ref[-nuisance.index,])
      T.SSRT.ref[j, k] =  mean((Y[-nuisance.index] - mu0.hat.ref)[(W.ref[-nuisance.index] == 1) & (Group[-nuisance.index] == k)]) - mean((Y[-nuisance.index] - mu0.hat.ref)[(W.ref[-nuisance.index] == 0) & (Group[-nuisance.index] == k)])
    }
  }
  record$pValue$SSRT[i,] = sapply(seq(1, Group.number), function(x){return(permutation.p.value(abs(T.SSRT[x]), abs(T.SSRT.ref[,x])))}) 
  record$R$SSRT[[i]] = which(record$pValue$SSRT[i,] <= BH.threshold(pval = record$pValue$SSRT[i,], q = q))
  record$FDP$SSRT[i] = sum(tau.group[record$R$SSRT[[i]]] == 0) / max(1, length(record$R$SSRT[[i]]))
  record$power$SSRT[i] = sum(tau.group[record$R$SSRT[[i]]] != 0) / sum(tau.group != 0)
  
  # ART: augmented RT
  # Only use one knockoff copy.
  W.knockoff = rbinom(n, 1, p)
  W.tilde = (W.knockoff + W) / 2
  data.tilde = data.frame(Y, X, G, W.tilde - 0.5, (W.tilde - 0.5) * X, (W.tilde - 0.5) * G)
  colnames(data.tilde) = colnames(data)
  nuisance.model = lm(Y ~ ., data = data.tilde)
  mu0.hat = predict(nuisance.model, data)
  T.ART = rep(0, Group.number)
  for(k in 1: Group.number){
    T.ART[k] = mean((Y - mu0.hat)[(W == 1) & (Group == k)]) -  mean((Y - mu0.hat)[(W == 0) & (Group == k)])
  }
  T.ART.ref = matrix(0, M, Group.number)
  for(j in 1 : M){
    swap.index = which(rbinom(n, 1, 0.5) == 1)
    W.ref = W; W.ref[swap.index] = W.knockoff[swap.index]
    for (k in 1: Group.number){
      data.ref = data.frame(Y, X, G, 0 - 0.5, (0 - 0.5) * X, (0 - 0.5) * G)
      colnames(data.ref) = colnames(data)
      mu0.hat.ref = predict(nuisance.model, newdata = data.ref)
      T.ART.ref[j, k] =  mean((Y - mu0.hat.ref)[(W.ref == 1) & (Group == k)]) -  mean((Y - mu0.hat.ref)[(W.ref == 0) & (Group == k)])
    }
  }
  record$pValue$ART[i,] = sapply(seq(1, Group.number), function(x){return(permutation.p.value(abs(T.ART[x]), abs(T.ART.ref[,x])))}) 
  record$R$ART[[i]] = which(record$pValue$ART[i,] <= BH.threshold(pval = record$pValue$ART[i,], q = q))
  record$FDP$ART[i] = sum(tau.group[record$R$ART[[i]]] == 0) / max(1, length(record$R$ART[[i]]))
  record$power$ART[i] = sum(tau.group[record$R$ART[[i]]] != 0) / sum(tau.group != 0)
  
  if(i %% 10 == 0){print(i)}
}

# results
result = rbind(lapply(record$FDP, mean),
               lapply(record$FDP, sd),
               lapply(record$power, mean),
               lapply(record$power, sd))
row.names(result) = c("FDR", "FDR sd", "power", "power sd")
print(result)
