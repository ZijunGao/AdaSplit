# Use FRT to test subgroup treatment effects
rm(list = ls())
source("~/Desktop/Research/Yao/HTE inference/code/Panning/helper ZG.R")
# n too large: pnorm(-sqrt(n)) super close to zero
# n too small: asymptotics has not kicked in

# TODO: add variation quantification of the comparison

# HTE2: n not too large, var(tau) too large, asymptotics have not kicked in. Denoise appears to be more powerful than AIPW.normalized because denoise's test statistic is quite noisy and the absolute value tends to be relatively large. In fact, if we change the two tailed test to one tailed test, then the CDF of denoise's p-value stays below that of AIPW.normalized except in a neighborhood of zero.

# Estiamted distribution of the test statistics.
  # Reasonably well when the aymptotics have kicked in. If n is too small, or ATE small but the treatment is highly variant, the estimator does not work very well. 
  # Works better for B small, since the error is proportional to \sqrt{B}.
  # Useful if the asymptotic power of the test with some test statistic is challenging to derive and the test statistic is a function of the moments of g(Z_i, Y_i, X_i), for some known function g.

n = 100 # sample size is small because it is just for one group
d = 5 # number of covariates
p = 0.5 # propensity score; 0.5, 0.2
Group.level.number = c(1) # number of levels per group
Group.number = prod(Group.level.number) # total number of groups
beta0 = 1; beta =rep(1,d); theta = rep(1, 2)
delta = 0.5 # 1
sigma = 1  # error magnitude in generating Y(0); 1

nuisance.learner.method = "gradient boosting" # "gradient boosting", "gradient boosting early
B = 6
M = 200 # number of permutations; 400
q = 0.10 # FDR level

# setting
setting = "unbalanced" # "balanced", "unbalanced"
if(setting == "unbalanced"){p = 0.2}

start.time = proc.time() 
m = 200 # number of trials
record = list()
record$pValue = list()
record$pValue$ORT$denoise = record$pValue$ORT$denoise.normalized = record$pValue$ORT$AIPW = record$pValue$ORT$AIPW.normalized = matrix(0, nrow = m, ncol = Group.number)
record$pValue$ART.early.stopping = record$pValue$ART  = record$pValue$ORT

# set.seed(318)
for(i in 1:m){
  
  # generate data
  X = matrix(runif(n * d, -1, 1), nrow = n, ncol = d) # covariates
  
  quantiles_1 <- quantile(c(-Inf,X[,d],Inf), probs = seq(0, 1, length.out = Group.level.number[1] + 1))
  #quantiles_2 <- quantile(c(-Inf,X[,d-2],Inf), probs = seq(0, 1, length.out = Group.level.number[2] + 1))
  S = as.numeric(cut(X[,d-1], quantiles_1)) #cbind(as.numeric(cut(X[,d-1], quantiles_1)), as.numeric(cut(X[,d-2], quantiles_2)))
  
  Group = S-1 #(S[, 1] - 1) * Group.level.number[2] + S[, 2]; 
  G = Group # model.matrix(~ factor(Group))[, -1] # one-hot encoding of group membership
  
  W = rbinom(n, 1, p) # treatment assignment
  
  # potential outcomes
  mu0 = (X %*% beta)
  Y0 = mu0
  tau = 3  + 10* rnorm(n, 0, sigma)
  tau.group = sapply(seq(1, Group.number), function(x) {
    mean(tau[Group == x])
  }) # average treatment effect in each group.
  Y1 = mu0 + tau 
  Y = Y1 * W + Y0 * (1 - W) # observed outcome
  # nuisance functions
  mu1 = mu0 + tau
  mu = mu0 * (1 - p) + mu1 * p
  
  
  # inference
  tau.hat = tau
  mu.hat = mu # mu 
  mu1.hat = mu.hat + (1 - p) * tau.hat
  mu0.hat = mu.hat - p * tau.hat
  # ORT
  record$pValue$ORT$AIPW[i,] = ORT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "AIPW", mu0 = mu0.hat, mu1 = mu1.hat, mu = mu.hat, tau = tau.hat, M = M)$pval

  record$pValue$ORT$denoise[i,] = ORT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "denoise", mu0 = mu0.hat, mu1 = mu1.hat, mu = mu.hat, tau = tau.hat, M = M)$pval
  
  record$pValue$ORT$AIPW.normalized[i,] = ORT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "AIPW normalized", mu0 = mu0.hat, mu1 = mu1.hat, mu = mu.hat, tau = tau.hat, M = M)$pval
  
  record$pValue$ORT$denoise.normalized[i,] = ORT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "denoise normalized", mu0 = mu0.hat, mu1 = mu1.hat, mu = mu.hat, tau = tau.hat, M = M)$pval
  
  # ART
  # gradient boosting, overfitting
  # record$pValue$ART$AIPW[i,] = ART(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "AIPW", M = M, nuisance.learner.method = "gradient boosting", B = B)$pval
  # 
  # record$pValue$ART$denoise[i,] = ART(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "denoise", M = M, nuisance.learner.method = "gradient boosting", B = B)$pval
  # 
  # record$pValue$ART$AIPW.normalized[i,] = ART(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "AIPW normalized", M = M, nuisance.learner.method = "gradient boosting", B = B)$pval
  # 
  # record$pValue$ART$denoise.normalized[i,] = ART(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "denoise normalized", M = M, nuisance.learner.method = "gradient boosting", B = B)$pval

  # gradient boosting, early stopping
  # record$pValue$ART.early.stopping$AIPW[i,] = ART(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "AIPW", M = M, nuisance.learner.method = "gradient boosting early stopping", B = B)$pval

  #record$pValue$ART.early.stopping$denoise[i,] = ART(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "denoise", M = M, nuisance.learner.method = "gradient boosting early stopping", B = B)$pval

  #record$pValue$ART.early.stopping$AIPW.normalized[i,] = ART(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "AIPW normalized", M = M, nuisance.learner.method = "gradient boosting early stopping", B = B)$pval

  #record$pValue$ART.early.stopping$denoise.normalized[i,] = ART(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "denoise normalized", M = M, nuisance.learner.method = "gradient boosting early stopping", B = B)$pval

  if(i %% 10 == 0){print(i)}
}

# results
record$setting = setting
index = Group.number - 0
par(mfrow = c(1, 1))
# lapply(record$pValue, hist)
for(i in 1:1){
  cdf.curve = data.frame(lapply(record$pValue[[i]], function(x){sapply(seq(1, 100)/100, function(y)(mean(x[, index] <= y)))}))
  matplot(cdf.curve, main = names(record$pValue)[i], x = seq(1, 100)/100, xlab = "alpha", type = "l", lty = 1, col = seq(2, 5), ylim = c(0, 1)); legend("bottomright", colnames(cdf.curve), col = seq(2, 5), lty = 1); abline(a = 0, b = 1, lty = 3)
}

# distribution of randomization p-values
# hist(record$pValue$denoise, breaks = 100)
# abline(v = pnorm(-sqrt(n)), "denoise.p.value") # pnorm(-sqrt(n)) could be too small to plot 


# saveRDS(record, file.patht("~/Desktop/Research/Yao/HTE inference/code/Panning/0630", paste(setting, "rds", sep= ".")))

