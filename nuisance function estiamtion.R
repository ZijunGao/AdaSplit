# evaluate the quality of nuisance function estimators

n = 1000 # sample size; 40
d = 5 # number of covariates
p.seq = seq(0.2, 0.8, by = 0.1) # propensity score
Group.level.number = c(4, 4) # number of levels per group
Group.number = prod(Group.level.number) # total number of groups
beta0 = 1; beta = rep(1, d); theta = rep(1, 2)
delta = 2
sigma = 1.5 # error magnitude in generating Y(0); 1

M = 400 # number of permutations
q = 0.2 # FDR level

m = 100
record = list(); record$cor = record$mu.R2 = record$tau.R2 = record$ratio = matrix(0, nrow = m, ncol = length(p.seq))
for(j in 1 : length(p.seq)){
  p = p.seq[j]
  for(i in 1 : m){
  S = cbind(apply(rmultinom(n, 1, rep(1, Group.level.number[1])), 2, which.max),
            apply(rmultinom(n, 1, rep(1, Group.level.number[2])), 2, which.max)) # S used to define subgroups
  Group = (S[, 1] - 1) * Group.level.number[2] + S[, 2]; 
  G = model.matrix(~ factor(Group))[, -1] # one-hot encoding of group membership
  Group = Group * 0 + 1
  X = matrix(rnorm(n * d, 0, 1), nrow = n, ncol = d) # covariates
  # X = matrix(runif(n * d, 0, 1), nrow = n, ncol = d)
  W = rbinom(n, 1, p) # treatment assignment
  # potential outcomes
  # tau is linear in S and independent of X
  tau = delta * (X[,1])^2
  mu0 =  beta0 + X %*% beta + S %*% theta
  mu1 =  beta0 + X %*% beta + S %*% theta + tau
  mu = mu0 * (1 - p) + mu1 * p
  # mu0 is linear in X and S
  Y0 = mu0 + rnorm(n, 0, sigma)
  Y1 = Y0 + tau
  Y = Y1 * W + Y0 * (1 - W) # observed outcome

  test = nuisance.learner(Y = Y, X = X, prop = p, G = G, W = W, method = "gradient boosting early stopping", train.index = seq(1, n / 2), test.index = seq(n/2 + 1, n)) # half training, half testing
  
  record$mu.R2[i, j] = 1 - mean((test$mu.hat - mu[seq(n/2 + 1, n)])^2) / mean((mu[seq(n/2 + 1, n)])^2)
  record$tau.R2[i, j] = 1 - mean((test$tau.hat - tau[seq(n/2 + 1, n)])^2) / mean((tau[seq(n/2 + 1, n)])^2)
  record$cor[i, j] = cor(test$mu.hat - mu[seq(n/2 + 1, n)], test$tau.hat - tau[seq(n/2 + 1, n)])
  record$ratio[i, j] = mean((test$mu0.hat * p - mu0[seq(n/2 + 1, n)] * p + test$mu1.hat * (1 - p) - mu1[seq(n/2 + 1, n)] * (1-p))^2) / mean((test$mu.hat - mu[seq(n/2 + 1, n)])^2)
  }
  print(j)
}

par(mfrow = c(2,2))
plot(p.seq, apply(record$tau.R2, 2, mean), type = "l", main = "tau R2")
plot(p.seq, apply(record$mu.R2, 2, mean), type = "l",main = "mu R2")
plot(p.seq, apply(record$cor, 2, mean), type = "l",main = "cor")
plot(p.seq, apply(record$ratio, 2, mean), type = "l",main = "error ratio")

    
    
    