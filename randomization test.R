# Use FRT to test subgroup treatment effects
rm(list = ls())
source("~/Desktop/Research/Yao/HTE inference/code/Panning/helper.R")

n = 50 # sample size; 400
d = 5 # number of covariates
p = 0.2 # propensity score
Group.level.number = c(4, 4) # number of levels per group
Group.number = prod(Group.level.number) # total number of groups
beta0 = 1; beta = rep(1, d); theta = rep(1, 2)
delta = 0.5
sigma = 1 # error magnitude in generating Y(0)

M = 400 # number of permutations
q = 0.2 # FDR level

# setting
setting = "HTE" # "constant treatment effect", "HTE"

start.time = proc.time() 
m = 1000 # number of trials
record = list()
record$pValue = list()
record$pValue$denoise = record$pValue$denoise.normalized = record$pValue$AIPW = record$pValue$AIPW.normalized = matrix(0, nrow = m, ncol = 1) 

set.seed(318)
for(i in 1:m){
  # generate data
  S = cbind(apply(rmultinom(n, 1, rep(1, Group.level.number[1])), 2, which.max),
            apply(rmultinom(n, 1, rep(1, Group.level.number[2])), 2, which.max)) # S used to define subgroups
  Group = (S[, 1] - 1) * Group.level.number[2] + S[, 2]; 
  G = model.matrix(~ factor(Group))[, -1] # one-hot encoding of group membership
  Group = Group * 0 + 1
  X = matrix(rnorm(n * d, 0, 1), nrow = n, ncol = d) # covariates
  W = rbinom(n, 1, p) # treatment assignment
  # potential outcomes
  # mu0 is linear in X and S
  Y0 = beta0 + X %*% beta + S %*% theta + rnorm(n, 0, sigma)
  # tau is linear in S and independent of X
  if(setting == "constant treatment effect"){
    tau = delta * rep(1, n)
  }else if(setting == "HTE"){
    tau = delta * (X[,1])^2
  }
  Y1 = Y0 + tau
  Y = Y1 * W + Y0 * (1 - W) # observed outcome
  # nuisance functions
  mu0 =  beta0 + X %*% beta + S %*% theta
  mu1 =  beta0 + X %*% beta + S %*% theta + tau
  mu = mu0 * (1 - p) + mu1 * p
  
  # inference
  # AIPW
  record$pValue$AIPW[i,] = ORT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "AIPW", mu0 = mu0, mu1 = mu1, mu = mu, tau = tau, M = M)$pval
  
  # AIPW normalized
  record$pValue$AIPW.normalized[i,] = ORT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "AIPW normalized", mu0 = mu0, mu1 = mu1, mu = mu, tau = tau, M = M)$pval
  
  # denoise
  record$pValue$denoise[i,] = ORT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "denoise", mu0 = mu0, mu1 = mu1, mu = mu, tau = tau, M = M)$pval
  
  # denoise normalized
  record$pValue$denoise.normalized[i,] = ORT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "denoise normalized", mu0 = mu0, mu1 = mu1, mu = mu, tau = tau, M = M)$pval

  if(i %% 10 == 0){print(i)}
}

# results
record$setting = setting
par(mfrow = c(1,1))
# lapply(record$pValue, hist)
cdf.curve = data.frame(lapply(record$pValue, function(x){sapply(seq(1, 20)/20, function(y)(mean(x <= y)))}))
matplot(cdf.curve, main = setting, x = seq(1, 20)/20, xlab = "alpha", type = "l", lty = 1); legend("bottomright", colnames(cdf.curve), col = seq(1, 4), lty = 1)

# saveRDS(record, file.path("~/Desktop/Research/Yao/HTE inference/code/Panning/0630", paste(setting, "rds", sep= ".")))


