# Verify the consistency of Bar-learner
# lambda is chosen by cross-validation

source("~/Desktop/Research/Yao/HTE inference/code/Panning/helper Yao.R")

set.seed(318)
setting = 2
# n = 300 # sample size
n.seq = c(200, 400, 800) # c(200, 400, 600, 800, 1000) # c(200, 400, 800, 1600) # sequence of sample sizes
Group.number = 5 # total number of groups
delta = 0.5 # effect size
sigma = 0.2 # std of noise
n_trial = 20 # number of runs
proportion = 0.2 # 0.5; proportion of randomness in the nuisance fold
test.stats.method ="AIPW" # test statistics
num_features = 5

# Store the experimental results
record = list()
record$betaError.R = record$R2.R = record$betaError = record$R2 = matrix(0, nrow = n_trial, ncol = length(n.seq))
record$lambda.seq = list()
for(n.index in seq(1, length(n.seq))){
  record$lambda.seq[[n.index]] = list()
  n = n.seq[n.index]
  for (j in 1:n_trial){
    
    X = matrix(runif(n*num_features, 0, 1), nrow = n, ncol = num_features) 
    quantiles_1 <- quantile(c(-Inf,X[,1],Inf), probs = seq(0, 1, length.out = Group.number[1] + 1))
    S = as.numeric(cut(X[,1], quantiles_1))
    
    Group = S - 1
    G = model.matrix(~ factor(Group))[, -1]
    Ex = rep(0.5,n)
    X[,2] = (X[,2]>0.75)
    mu = X %*%rnorm(num_features, 1, 1)
    tau = X[,setting]*delta
    beta = rep(0, num_features + 1); beta[setting + 1] = delta
    
    mu0 <- mu - Ex * tau
    mu1 <- mu0 + tau
    
    W <- sapply(Ex, function(p) rbinom(1, size = 1, prob = p))
    Y = mu0 + W*tau + rnorm(n, 0, sigma)
    
    mu.hat = nuisance.mu(Y,X)
    
    # BaR-learner
    tau.hat.active = nuisance.tau.active(Y = Y, X = X, Ex = Ex, W = W, mu = mu.hat, Group = Group, proportion = proportion, lambda.update.period = 50, loss_threshold = -Inf, marginalize = T) #estimate tau using active learning
    record$lambda.seq[[n.index]][[j]] = tau.hat.active$lambda.seq

    train.index.active = tau.hat.active$train.index #the data points fitted to the nuisance
    test.index.active = setdiff(1:n,train.index.active) #the data points used for inference
    
    v = sum((tau[test.index.active]-tau.hat.active$tau[test.index.active])**2)/sum((tau[test.index.active]-mean(tau[test.index.active]))**2)
    cat("R^2 of AdaSplit + BaR-learner:", 1- v, "\n")
    record$R2[j, n.index] = 1-v
    record$betaError[j, n.index] = sum((tau.hat.active$beta - beta)^2)/ sum(beta^2)
    
    # R-learner
    # Sample splitting
    train.index.ss = select.train(Group, proportion)
    test.index.ss = setdiff(1:n, train.index.ss)
    # Fit the nuisance parameter
    tau.hat.ss = nuisance.tau.ss(Y = Y, X = X, Ex = Ex, W = W, mu = mu.hat, train.index = train.index.ss, marginalize = F)
    v = sum((tau[test.index.ss]-tau.hat.ss$tau[test.index.ss])**2)/sum((tau[test.index.ss]-mean(tau[test.index.ss]))**2)
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
           BaR.R2 = apply(record$R2, 2, mean),  
           R.R2 = apply(record$R2.R, 2, mean),  
           BaR.BetaError = apply(record$betaError, 2, mean), 
           R.BetaError = apply(record$betaError.R, 2, mean))


# saveRDS(record,"~/Desktop/Research/Yao/HTE inference/code/Panning/April 2025/BaR_learner_consistency.rds")


