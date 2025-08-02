# Verify the consistency of Bar-learner
# lambda is chosen by cross-validation
rm(list = ls())
source("~/Desktop/Research/Yao/HTE inference/code/Panning/helper Yao.R")
directory = "~/Desktop/Research/Yao/HTE inference/code/Panning/April 2025" 

n.seq = 300 * seq(1, 5) # 300 * seq(1, 5); sequence of sample sizes
num_features = 5
Group.number = 5 # total number of groups
delta = 1 # effect size; 1, 0
sigma = 1 # std of noise: 1
M = 1000 # number of permutations
proportion = 0.5 # proportion of randomness in the nuisance fold
test.stats.method ="AIPW" # test statistics
n_trial = 100 # number of runs; 2; 100
verbose = F # if true, print out intermediate results
q = 0.2 # FWER level
marginalize = T # if marginalize
mu.learner = "linear"

# Store the experimental results
record = list()
record$betaError.R = record$R2.R = record$betaError = record$R2 = matrix(0, nrow = n_trial, ncol = length(n.seq))
record$lambda.seq = list()
set.seed(318)

# Data for evaluating the nuisance estimation performance
n.val = 10^4
X.val = matrix(runif(n.val*num_features, 0, 1), nrow = n.val, ncol = num_features) 
quantiles_1 <- quantile(c(-Inf,X.val[,1],Inf), probs = seq(0, 1, length.out = Group.number[1] + 1))
S.val = as.numeric(cut(X.val[,1], quantiles_1))
Group.val = S.val - 1
G.val = model.matrix(~ factor(Group.val))[, -1]
Ex.val = rep(0.5,n.val)
X.val[,2] = (X.val[,2]>0.75); X.val[,3] = (X.val[,3]>0.25)
mu.val = X.val %*%rnorm(num_features, 1, 1)
tau.val = (X.val - 0.5) %*% rep(delta, num_features) + 0.5 * delta

mu0.val <- mu.val - Ex.val * tau.val
mu1.val <- mu0.val + tau.val

W.val <- sapply(Ex.val, function(p) rbinom(1, size = 1, prob = p))
Y.val = mu0.val + W.val*tau.val + rnorm(n.val, 0, sigma)

for(n.index in seq(1, length(n.seq))){
  record$lambda.seq[[n.index]] = list()
  n = n.seq[n.index]
  for (j in 1:n_trial){
    # Data generation
    X = matrix(runif(n*num_features, 0, 1), nrow = n, ncol = num_features) 
    quantiles_1 <- quantile(c(-Inf,X[,1],Inf), probs = seq(0, 1, length.out = Group.number[1] + 1))
    S = as.numeric(cut(X[,1], quantiles_1))
    Group = S - 1
    G = model.matrix(~ factor(Group))[, -1]
    Ex = rep(0.5,n)
    X[,2] = (X[,2]>0.75)
    mu = X %*%rnorm(num_features, 1, 1)
    X[,3] = (X[,3]>0.25); tau = (X - 0.5) %*% rep(delta, num_features) + 0.5 * delta; beta = rep(0, num_features + 1); beta[1] = -0.5 * delta * num_features + 0.5 * delta; beta[-1] = delta
    mu0 <- mu - Ex * tau
    mu1 <- mu0 + tau
    W <- sapply(Ex, function(p) rbinom(1, size = 1, prob = p))
    Y = mu0 + W*tau + rnorm(n, 0, sigma)
    
    
    if(mu.learner == "linear"){mu.hat = nuisance.mu(Y,X)
    }else if(mu.learner == "xgboost"){mu.hat = nuisance.mu.xgboost(Y,X)}
    
    # BaR-learner
    tau.hat.active = nuisance.tau.active(Y = Y, X = X, Ex = Ex, W = W, mu = mu.hat, Group = Group, proportion = proportion, lambda = 1, lambda.update.period = 100000, loss_threshold = -Inf, marginalize = T, initial.proportion = 0.2, mixture.period = 2, weighting = T) #estimate tau using active learning
    record$lambda.seq[[n.index]][[j]] = tau.hat.active$lambda.seq
    v = sum((tau.val-cbind(1, X.val) %*% tau.hat.active$beta)**2)/sum((tau.val - mean(tau.val))**2)
    cat("R^2 of AdaSplit + BaR-learner:", 1- v, "\n")
    record$R2[j, n.index] = 1-v
    record$betaError[j, n.index] = sum((tau.hat.active$beta - beta)^2)/ sum(beta^2)
    
    # R-learner
    # Sample splitting
    train.index.ss = select.train(Group, proportion)
    test.index.ss = setdiff(1:n, train.index.ss)
    # Fit the nuisance parameter
    tau.hat.ss = nuisance.tau.ss(Y = Y, X = X, Ex = Ex, W = W, mu = mu.hat, train.index = train.index.ss, marginalize = F)
    v = sum((tau.val-cbind(1, X.val) %*% tau.hat.ss$beta)**2)/sum((tau.val - mean(tau.val))**2)
    cat("R^2 of SS + R-learner:", 1 - v, "\n")
    record$R2.R[j, n.index] = 1-v
    record$betaError.R[j, n.index] = sum((tau.hat.ss$beta - beta)^2) / sum(beta^2)
    
    if(j %% 10 == 0){print(j)}
  }
  print(n)
}

# visualization

# results
record$n.seq = n.seq

data.frame(sample.size = n.seq,
           BaR.Ada.R2 = apply(record$R2, 2, mean),  
           R.SS.R2 = apply(record$R2.R, 2, mean),  
           BaR.Ada.BetaError = apply(record$betaError, 2, mean), 
           R.SS.BetaError = apply(record$betaError.R, 2, mean))


# saveRDS(record,"~/Desktop/Research/Yao/HTE inference/code/Panning/July 2025/BaR_learner_consistency_mixture.rds") # BaR_learner_consistency; BaR_learner_consistency_mixture
