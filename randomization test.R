# Use FRT to test subgroup treatment effects
rm(list = ls())
source("~/Desktop/Research/Yao/HTE inference/code/Panning/helper ZG.R")
# n too large: pnorm(-sqrt(n)) super close to zero
# n too small: asymptotics has not kicked in

# HTE2: n not too large, var(tau) too large, asymptotics have not kicked in. Denoise appears to be more powerful than AIPW.normalized because denoise's test statistic is quite noisy and the absolute value tends to be relatively large. In fact, if we change the two tailed test to one tailed test, then the CDF of denoise's p-value stays below that of AIPW.normalized except in a neighborhood of zero.

# Estiamted distribution of the test statistics.
  # Reasonably well when the aymptotics have kicked in. If n is too small, or ATE small but the treatment is highly variant, the estimator does not work very well. 
  # Works better for B small, since the error is proportional to \sqrt{B}.
  # Useful if the asymptotic power of the test with some test statistic is challenging to derive and the test statistic is a function of the moments of g(Z_i, Y_i, X_i), for some known function g.

n = 40 # sample size; 400
d = 5 # number of covariates
p = 0.2 # propensity score
Group.level.number = c(4, 4) # number of levels per group
Group.number = prod(Group.level.number) # total number of groups
beta0 = 1; beta = rep(1, d); theta = rep(1, 2)
delta = 1
sigma = 0.5 # error magnitude in generating Y(0); 1

B = 2
M = 1000 # number of permutations
q = 0.2 # FDR level

# setting
setting = "HTE2" # "default", "HTE", "HTE2", "HTE3", "misspecified"

start.time = proc.time() 
m = 40 # number of trials
record = list()
record$pValue = list()
record$pValue$denoise = record$pValue$denoise.normalized = record$pValue$AIPW = record$pValue$AIPW.normalized = matrix(0, nrow = m, ncol = 1)
# record$empirical.test.stats = record$empirical.test.stats = record$theoretical.test.stats = record$theoretical.ref.var = record$estimated.test.stats = record$estimated.ref.var = record$pValue
  
# set.seed(318)
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
  if(setting == "default"){
    tau = delta * rep(1, n)
  }else if(setting == "HTE"){
    tau = delta * (X[,1] + 1) # delta * (X[,1] + 1); delta * (X[,1])^2
  }else if(setting == "HTE2"){
    tau = delta * (10 * X[,1] + 0) # delta * (10 * X[,1] + 1)
  }else if(setting == "HTE3"){
    tau = delta * (X[,1])^2
  }
  Y1 = Y0 + tau
  Y = Y1 * W + Y0 * (1 - W) # observed outcome
  # nuisance functions
  mu0 =  beta0 + X %*% beta + S %*% theta
  mu1 =  beta0 + X %*% beta + S %*% theta + tau
  mu = mu0 * (1 - p) + mu1 * p
  
  # inference
  mu1.hat = mu1 # mu1 + rnorm(n, 0, 1)
  mu0.hat = mu0 # mu0 + rnorm(n, 0, 1)
  tau.hat = mu1.hat - mu0.hat; mu.hat = p * mu1.hat + (1 - p) * mu0.hat 
  # AIPW
  record$pValue$AIPW[i,] = ORT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "AIPW", mu0 = mu0.hat, mu1 = mu1.hat, mu = mu.hat, tau = tau.hat, M = M)$pval
  # record$empirical.test.stats[i] = recor
  
  # AIPW normalized
  record$pValue$AIPW.normalized[i,] = ORT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "AIPW normalized", mu0 = mu0.hat, mu1 = mu1.hat, mu = mu.hat, tau = tau.hat, M = M)$pval
  
  # denoise
  record$pValue$denoise[i,] = ORT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "denoise", mu0 = mu0.hat, mu1 = mu1.hat, mu = mu.hat, tau = tau.hat, M = M)$pval
  
  # denoise normalized
  record$pValue$denoise.normalized[i,] = ORT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "denoise normalized", mu0 = mu0.hat, mu1 = mu1.hat, mu = mu.hat, tau = tau.hat, M = M)$pval

  if(i %% 10 == 0){print(i)}
}

# results
record$setting = setting
par(mfrow = c(1,1))
# lapply(record$pValue, hist)
cdf.curve = data.frame(lapply(record$pValue, function(x){sapply(seq(1, 100)/100, function(y)(mean(x <= y)))}))
matplot(cdf.curve, main = setting, x = seq(1, 100)/100, xlab = "alpha", type = "l", lty = 1, col = seq(2, 5), ylim = c(0, 1)); legend("bottomright", colnames(cdf.curve), col = seq(2, 5), lty = 1); abline(a = 0, b = 1, lty = 3)

# distribution of randomization p-values
# hist(record$pValue$denoise, breaks = 100)
# abline(v = pnorm(-sqrt(n)), "denoise.p.value") # pnorm(-sqrt(n)) could be too small to plot 

# Comparison of the theoretical and the empirical value of the test statistic and the variance of the associated reference distribution
# Only work for homoscedastic Y_i(j); test statistic without taking the absolute value
# denoise
test = ORT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "denoise", mu0 = mu0, mu1 = mu1, mu = mu, tau = tau, M = M)
tau.hat = tau; mu.hat = mu; mu1.hat = mu1; mu0.hat = mu0
(result.denoise = data.frame(empirical.ref.mean = mean(test$test.stats.ref),
           theoretical.ref.mean = 0,
           estimated.ref.mean = 0,
           empirical.ref.var = var(test$test.stats.ref) * n, 
           theoretical.ref.var = mean(tau^2) + sigma^2/p/(1-p), 
           estimated.ref.var = p * (1 - p) * mean(((Y - mu.hat) / p + (Y - mu.hat) / (1 - p))^2), 
           empirical.test.stats= test$test.stats, 
           theoretical.test.stats = mean(tau),
           estimated.test.stats = mean(tau.hat)))
# AIPW normalized
test = ORT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, test.stats.method = "AIPW normalized", mu0 = mu0, mu1 = mu1, mu = mu, tau = tau, M = M)
W.knockoff = matrix(rbinom(n * B, 1, p), ncol = B) 
first.moment.hat = sum(apply(cbind(W, W.knockoff), 2, function(x){mean(x * ((Y - mu1.hat) / p + mu1.hat - mu0.hat) + (1-x) * (-(Y - mu0.hat) / (1-p) + mu1.hat - mu0.hat))})) # mean(mu1.hat - mu0.hat)
second.moment.hat = sum(apply(cbind(W, W.knockoff), 2, function(x){mean(x * ((Y - mu1.hat) / p + mu1.hat - mu0.hat)^2 + (1-x) * (-(Y - mu0.hat) / (1-p) + mu1.hat - mu0.hat)^2)})) - B * mean(p * ((Y - mu1.hat) / p + mu1.hat - mu0.hat)^2 + (1-p) * (-(Y - mu0.hat) / (1-p) + mu1.hat - mu0.hat)^2)
(result.AIPW.normalized = data.frame(empirical.ref.mean = mean(test$test.stats.ref),
           theoretical.ref.mean = 0,
           estimated.ref.mean = 0,
           empirical.ref.var = var(test$test.stats.ref) * n,
           theoretical.ref.var = 1,
           estimated.ref.var = 1,
           empirical.test.stats = test$test.stats,
           theoretical.test.stats = mean(tau) / sqrt(sigma^2/p/(1-p)),
           estimated.test.stats =  mean(tau) / sqrt(max(1e-6, second.moment.hat - first.moment.hat^2))))

(result.comparison = data.frame(empirical.ratio = ((result.denoise$empirical.test.stats - result.denoise$empirical.ref.mean)^2 / result.denoise$empirical.ref.var) / ((result.AIPW.normalized$empirical.test.stats - result.AIPW.normalized$empirical.ref.mean)^2 / result.AIPW.normalized$empirical.ref.var),
          theoretical.ratio = ((result.denoise$theoretical.test.stats - result.denoise$theoretical.ref.mean)^2 / result.denoise$theoretical.ref.var) / ((result.AIPW.normalized$theoretical.test.stats - result.AIPW.normalized$theoretical.ref.mean)^2 / result.AIPW.normalized$theoretical.ref.var),
          estimated.ratio = ((result.denoise$estimated.test.stats - result.denoise$estimated.ref.mean)^2 / result.denoise$estimated.ref.var) / ((result.AIPW.normalized$estimated.test.stats - result.AIPW.normalized$estimated.ref.mean)^2 / result.AIPW.normalized$estimated.ref.var)))



# saveRDS(record, file.patht("~/Desktop/Research/Yao/HTE inference/code/Panning/0630", paste(setting, "rds", sep= ".")))

