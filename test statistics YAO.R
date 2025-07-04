# Use FRT to test subgroup treatment effects
# Sample size
rm(list = ls())
source("~/Desktop/Research/Yao/HTE inference/code/Panning/helper.R")

n.seq = c(800) #seq(400, 1200, by = 200) # sample size
d = 5 # number of covariates
p = 0.2 # propensity score
Group.level.number = c(4, 5) # number of levels per group
Group.number = prod(Group.level.number) # total number of groups
beta0 = 1; beta = rep(1, d) * 1; theta = rep(1, 2)
delta = 1.5
sigma = 1  # error magnitude in generating Y(0)

nuisance.learner.method = "gradient boosting"
test.stats.method = "AIPW + ITE" # "denoise", "ATE", "denoise + ATE", "AIPW", "ITE"; test statistic
B = 6 # number of knockoffs: seq(1, 10); c(1, 5, 10)
M = 200 # number of permutations
q = 0.2 # FDR level

setting = "Statistics"

start.time = proc.time()
m = 200 # number of trials
record = list()
record$pValue = list()
record$pValue$ART1 = record$pValue$ART2 = record$pValue$ART3 = record$pValue$ART4 = record$pValue$ART5 = array(0, dim = c(m, Group.number, length(n.seq))) # ORT: oracle RT; RT (baseline): standard RT; SSRT: sample-splitting RT; DDRT: double-dipping RT; ART: augmented RT; 
record$R = list(); 
record$R$ART = list(); 
for(j in seq(1, length(n.seq))){record$R$ART[[j]] = list()}; 
record$R$ART1 = record$R$ART2 = record$R$ART3 = record$R$ART4 = record$R$ART5 = record$R$ART
record$FDP = list();
record$FDP$ART1 = record$FDP$ART2 = record$FDP$ART3 = record$FDP$ART4 = record$FDP$ART5 = matrix(0, m, length(n.seq)) 
record$power = record$FDP


set.seed(318)


for(j in 1 : length(n.seq)){
  n = n.seq[j]
  
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
    # tau is linear in S and independent of X
    #tau = delta * (S[, 1] >= (Group.level.number[1]-1)) * (S[, 2] >= (Group.level.number[2] - 1)) * (1+ grepl("HTE", setting)* X %*% beta2) #  (S[, 1] >= 3) * (S[, 2] >= 3)
    tau = delta * (S[, 1] >= (Group.level.number[1]-1)) * (S[, 2] >= (Group.level.number[2] - 1)) * (1+ X[, 2]^2 )
    tau.group = sapply(seq(1, Group.number), function(x) {
      mean(tau[Group == x])
    }) # average treatment effect in each group.
    Y1 = Y0 + tau
    Y = Y1 * W + Y0 * (1 - W) # observed outcome
    # nuisance functions
    mu1 = mu0 + tau
    mu = mu0 * (1 - p) + mu1 * p
    
    
    record$pValue$ART2[i,,j] = ART(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = M, 
                                   nuisance.learner.method = nuisance.learner.method, test.stats.method = "denoise", B = B)$pval
    record$R$ART2[[j]][[i]] = which(record$pValue$ART2[i,,j] <= BH.threshold(pval = record$pValue$ART2[i,,j], q = q))
    record$FDP$ART2[i, j] = sum(tau.group[record$R$ART2[[j]][[i]]] == 0) / max(1, length(record$R$ART2[[j]][[i]]))
    record$power$ART2[i, j] = sum(tau.group[record$R$ART2[[j]][[i]]] != 0) / max(1, sum(tau.group != 0))
    
    record$pValue$ART3[i,,j] = ART(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = M, 
                                   nuisance.learner.method = nuisance.learner.method, test.stats.method = "AIPW", B = B)$pval
    record$R$ART3[[j]][[i]] = which(record$pValue$ART3[i,,j] <= BH.threshold(pval = record$pValue$ART3[i,,j], q = q))
    record$FDP$ART3[i, j] = sum(tau.group[record$R$ART3[[j]][[i]]] == 0) / max(1, length(record$R$ART3[[j]][[i]]))
    record$power$ART3[i, j] = sum(tau.group[record$R$ART3[[j]][[i]]] != 0) / max(1, sum(tau.group != 0))
    
    
  }
  print(j)
}

# results
result = list()
result$FDR = data.frame(lapply(record$FDP, function(x){apply(x, 2, mean)})); rownames(result$FDR) = n.seq
result$power = data.frame(lapply(record$power, function(x){apply(x, 2, mean)})); rownames(result$power) = n.seq
result


# saveRDS(record, file.path("~/Desktop/Research/Yao/HTE inference/code/Panning/0630", paste(setting, ".rds", sep= "")))

