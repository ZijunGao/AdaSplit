# helper functions for panning
library("xgboost")
library(randomForest)
library(caret)
library(ranger)
library(pracma)
library(Matrix)


# BH
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


BH.threshold.Storey = function(pval, q = 0.1, lambda = 0.5){
  # preprocess
  m = length(pval)
  pval = sort(pval, decreasing = FALSE)
  
  # Estimate pi0 using Storey's method
  pi0 = min(1, mean(pval > lambda) / (1 - lambda))
  
  # compute the threshold of the BH procedure with Storey correction
  FDR.hat = pi0 * m * pval / seq(1, m)
  pval.index = which(FDR.hat <= q)
  if(length(pval.index) == 0){return(-1e6)}
  threshold = pval[max(pval.index)]
  return(threshold)
}



permutation.p.value = function(stats, stats.ref){
  if(is.na(stats)){return(1)}
  return((sum(stats <= stats.ref[!is.na(stats.ref)]) + 1) / (length(stats.ref[!is.na(stats.ref)]) + 1))
}




# nuisance learner
nuisance.learner = function(Y, X = NULL, prop = NULL, G = NULL, W = NULL, method = "linear", train.index = NULL, test.index = NULL, ...){

  
  if(method == "gradient boosting"){

    X_ = matrix(X[train.index])
    Y_ = Y[train.index]
    W_ = (W-prop)[train.index]
    params = list(objective = "reg:squarederror", eval_metric = "rmse",eta=0.01)
    nfold = 10
    nrounds =1000
    early_stopping_rounds = 25
    subsample =0.8
    
    nuisance.mu <- xgb.cv(
      data = xgb.DMatrix(data = X_, label = Y_),
      params = params, 
      nfold = nfold,
      nrounds = nrounds,
      subsample = subsample,
      early_stopping_rounds = early_stopping_rounds,
      verbose = FALSE,
      callbacks = list(cb.cv.predict(save_models = TRUE))
    )
    
  
    nuisance.tau <- xgb.cv(
      data = xgb.DMatrix(data = X_, label = (Y_ - nuisance.mu$pred)/W_, weight =W_^2),
      params = params, 
      nfold = nfold,
      nrounds = nrounds,
      subsample = subsample,
      early_stopping_rounds = early_stopping_rounds,
      verbose = FALSE,
      callbacks = list(cb.cv.predict(save_models = TRUE))
    )
    
    #tau.hat <- predict(nuisance.tau, xgb.DMatrix(data =  matrix(X[test.index])))
    
    
    
    models <- nuisance.mu$models
    mu.hat =  0 
    for (model in models) {
      new_data_dmatrix <- xgb.DMatrix(data =  matrix(X[test.index]))
      mu.hat = mu.hat + predict(model, new_data_dmatrix)/nfold
    }
    
    
    models <- nuisance.tau$models
    tau.hat =  0 
    for (model in models) {
      new_data_dmatrix <- xgb.DMatrix(data =  matrix(X[test.index]))
      tau.hat = tau.hat + predict(model, new_data_dmatrix)/nfold
    }
  
    #mu.hat <- predict(nuisance.mu, xgb.DMatrix(data =  matrix(X[test.index])))
  
    
    
    mu0.hat <- mu.hat - prop[test.index] * tau.hat
    mu1.hat <- mu0.hat + tau.hat
    
    
    return(result = list(mu0.hat = mu0.hat, mu1.hat = mu1.hat, mu.hat = mu.hat, tau.hat = tau.hat, tau=nuisance.tau)) #tau=nuisance.tau 
  }
  
  
  
}







test.stats = function(xwx,  tb, I, Y_full, mu_full, X_full, W_full, W_tilde_full, tau_full, Group_full, Y, W, X = NULL, G = NULL, stats = "denoise", prop = NULL, mu0.hat = NULL, mu1.hat = NULL, mu.hat = NULL, tau.hat = NULL, W.tilde = NULL) {
  
  
  
  prop  <- pmax(0.001, pmin(prop, 0.999))

  n = length(Y)
  if (is.null(prop)){
    n1 = max(1, sum(W))
    n0 = max(1, sum(1 - W))
  }else{
    n1 = n * prop
    n0 = n * (1 - prop)
  }


  if (stats == "denoise") {

    IF = W * (Y - mu.hat)/prop - (1 - W) * (Y - mu.hat)/(1-prop)
    value = mean(IF)
    
  }else if(stats == "AIPW"){

    IF = W * (Y - mu1.hat)/prop - (1 - W) * (Y - mu0.hat)/(1-prop) + tau.hat
    value = mean(IF)/sd(IF)
    
  }  
  
  # Return the computed statistic value
  return(value)
}



test.stats.group = function(Y, W, X = NULL, G = NULL, Group, stats = "denoise", prop = NULL, mu0.hat = NULL, mu1.hat = NULL, mu.hat = NULL, tau.hat = NULL, W.tilde = NULL) {
  # Get the unique levels of the grouping variable
  Group.level = sort(unique(Group))
  
  # Compute the test statistic value for each group level
  values = sapply(Group.level, function(x) {
    index = which(Group == x)
    test.stats(I = x, Y_full = Y, mu_full = mu.hat,  X_full = X, W_full = W, W_tilde_full = W.tilde, tau_full = tau.hat, Group_full = Group,
               Y = Y[index], W = W[index], X = X[index,], G = G[index,], 
               stats = stats, prop = prop[index], 
               mu0.hat = mu0.hat[index], 
               mu1.hat = mu1.hat[index], 
               mu.hat = mu.hat[index], 
               tau.hat = tau.hat[index] , 
               W.tilde = W.tilde[index])
  })
  
  # Return the computed values for each group level
  return(values)
}


test.stats.group_new = function(XWX = XWX, TB=TB, Y, W, X = NULL, G = NULL, Group, stats = "denoise", prop = NULL, mu0.hat = NULL, mu1.hat = NULL, mu.hat = NULL, tau.hat = NULL, W.tilde = NULL) {
  # Get the unique levels of the grouping variable
  Group.level = sort(unique(Group))
  
  # Compute the test statistic value for each group level
  values = sapply(Group.level, function(x) {
    index = which(Group == x)
    test.stats(xwx = XWX[[x]], tb = TB[[x]], I = x, Y_full = Y, mu_full = mu.hat,  X_full = X, W_full = W, W_tilde_full = W.tilde, tau_full = tau.hat, Group_full = Group,
               Y = Y[index], W = W[index], X = X[index,], G = G[index,], 
               stats = stats, prop = prop, 
               mu0.hat = mu0.hat[index], 
               mu1.hat = mu1.hat[index], 
               mu.hat = mu.hat[index], 
               tau.hat = tau.hat[index] , 
               W.tilde = W.tilde[index])
  })
  
  # Return the computed values for each group level
  return(values)
}





# oracle
  # Example:
  # set.seed(318)
  # ORT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, mu0 = mu0, mu1 = mu1, mu = mu, tau = tau, M = 100, test.stats.method = "denoise")$pval
ORT = function(Y, W, X, G, Group, prop = NULL, mu0, mu1, mu, tau, test.stats.method = "denoise",  treatment.assignment = "Bernoulli", M = 10, seed = 318) {
  # Convert X and G to matrix format
  X = as.matrix(X)
  G = as.matrix(G)
  
  # Get the number of observations
  n = length(Y)
  
  # Calculate test statistics for the original data
  test.stats.value = test.stats.group(Y = Y, W = W, X = X, Group = Group, G = G, stats = test.stats.method, prop = prop, mu0.hat = mu0, mu1.hat = mu1, mu.hat = mu, tau.hat = tau)
  
  # Generate reference test statistics using permutations
  test.stats.ref = t(replicate(M, test.stats.group(Y = Y, W = rbinom(n, 1, prop), X = X, G = G, Group = Group, stats = test.stats.method, prop = prop, mu0.hat = mu0, mu1.hat = mu1, mu.hat = mu, tau.hat = tau))) #   Each column represents a group
  
  # Calculate p-values for each group
  pval = sapply(seq(1, length(test.stats.value)), function(x) {
    permutation.p.value(stats = test.stats.value[x], stats.ref = test.stats.ref[, x])
  })
  
  # Return the p value, the test statistic, and the test statistic values used to compute the reference distribution
  return(list(pval = pval, test.stats = test.stats.value, test.stats.ref = test.stats.ref))
}


# standard RT
  # set.seed(318)
  # RT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = 100)$pval
RT = function(Y, W, X, G, Group, prop = NULL, treatment.assignment = "Bernoulli", M = 10, seed = 318){
  # Convert X and G to matrix format
  X = as.matrix(X)
  G = as.matrix(G)
  
  # Get the number of observations
  n = length(Y)

  # Calculate test statistics for the original data
  test.stats.value = test.stats.group(Y = Y, W = W, X = X, Group = Group, G = G, stats = "plain", prop = prop)
  
  # Generate reference test statistics using permutations
  test.stats.ref = t(replicate(M, test.stats.group(Y = Y, W = rbinom(n, 1, prop), X = X, G = G, Group = Group, stats = "plain", prop = prop))) #  Each column represents a group
  
  # Calculate p-values for each group
  pval = sapply(seq(1, length(test.stats.value)), function(x) {
    permutation.p.value(stats = test.stats.value[x], stats.ref = test.stats.ref[, x])
  })
  
  # Return the p value, the test statistic, and the test statistic values used to compute the reference distribution
  return(list(pval = pval, test.stats = test.stats.value, test.stats.ref = test.stats.ref))
}


# sample splitting
# TODO: other treatmnet assignment
  # Example:
  # set.seed(318)
  # SS(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = 100, test.stats.method = "denoise")$pval
SS = function(Y, W, X, G, Group, prop = NULL, nuisance.learner.method = "linear", test.stats.method = "denoise", treatment.assignment = "Bernoulli", M = 10, seed = 318){
  # Convert X and G to matrix format
  X = as.matrix(X)
  G = as.matrix(G)
  
  # Get the number of observations
  n = length(Y)
  
  nuisance.index = c()

  # Estimate the nuisance functions
  nuisance.index = sample(n, n * 1/2) # TODO 
  train.index = nuisance.index; test.index = setdiff(seq(1, n), nuisance.index)
  nuisance.hat = nuisance.learner(Y = Y, X = X, prop = prop, G = G, W = W, method = nuisance.learner.method, train.index = train.index, test.index = test.index)
  
  # Calculate test statistics for the original data
  test.stats.value = test.stats.group(Y = Y[test.index], W = W[test.index], X = X[test.index,], Group = Group[test.index], G = G[test.index,], stats = test.stats.method, prop = prop[test.index], mu0.hat = nuisance.hat$mu0.hat[test.index,], mu1.hat = nuisance.hat$mu1.hat[test.index,], mu.hat = nuisance.hat$mu.hat[test.index,], tau.hat = nuisance.hat$tau.hat[test.index,])
  
  # Generate reference test statistics using permutations
  test.stats.ref = t(replicate(M, test.stats.group(Y = Y[test.index], W =  sapply(Ex[test.index], function(p) rbinom(1, size = 1, prob = p)), X = X[test.index,], G = G[test.index,], Group = Group[test.index], stats = test.stats.method, prop = prop[test.index], mu0.hat = nuisance.hat$mu0.hat[test.index,], mu1.hat = nuisance.hat$mu1.hat[test.index,], mu.hat = nuisance.hat$mu.hat[test.index,], tau.hat = nuisance.hat$tau.hat[test.index,]))) #   Each column represents a group 

  # Calculate p-values for each group
  pval = sapply(seq(1, length(test.stats.value)), function(x) {
    permutation.p.value(stats = test.stats.value[x], stats.ref = test.stats.ref[, x])
  })
  
  # Return the p value, the test statistic, and the test statistic values used to compute the reference distribution
  return(list(pval = pval, test.stats = test.stats.value, test.stats.ref = test.stats.ref))
}










ART = function(Y, W, X, G, Group, prop = NULL, nuisance.learner.method = "linear", test.stats.method = "denoise", treatment.assignment = "Bernoulli", M = 10, B = 1, seed = 318){
  # Convert X and G to matrix format
  X = as.matrix(X)
  G = as.matrix(G)
  
  # Get the number of observations
  n = length(Y)
  
  # Estimate the nuisance functions
  W.knockoff = matrix(rbinom(n * B, 1, p), ncol = B) 
  
  W.aug = cbind(W, W.knockoff)
  W.tilde = apply(W.aug, 1, mean) # TODO more knockoffs
  nuisance.hat = nuisance.learner(Y = Y, X = X, prop = prop, G = G, W = W.tilde, method = nuisance.learner.method, train.index = seq(1, n), test.index = seq(1, n))
  

  
  # Calculate test statistics for the original data
  test.stats.value = test.stats.group_new(XWX = XWX, TB=TB, Y = Y, W = W, X = X, Group = Group, G = G, stats = test.stats.method, prop = prop, mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat, W.tilde = W.tilde)
  
  # Generate reference test statistics using permutations
  test.stats.ref = t(replicate(M, test.stats.group_new(XWX = XWX, TB=TB, Y = Y, W = W.aug[cbind(seq(1, n), replicate(n, sample(1 + B, 1)))], X = X, G = G, Group = Group, stats = test.stats.method, prop = prop, mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat, W.tilde= W.tilde))) #   Each column represents a group 

    # Calculate p-values for each group
    pval = sapply(seq(1, length(test.stats.value)), function(x) {
      permutation.p.value(stats = test.stats.value[x], stats.ref = test.stats.ref[, x])
    })
    
    return(list(pval = pval, test.stats = test.stats.value, test.stats.ref = test.stats.ref))
}




# double-dipping
# Example:
# set.seed(318)
# DD(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = 100, test.stats.method = "denoise")$pval
DD = function(Y, W, X, G, Group, prop = NULL, nuisance.learner.method = "linear", test.stats.method = "denoise", treatment.assignment = "Bernoulli", M = 10, seed = 318){
  # Convert X and G to matrix format
  X = as.matrix(X)
  G = as.matrix(G)
  
  # Get the number of observations
  n = length(Y)
  
  # Estimate the nuisance functions
  nuisance.hat = nuisance.learner(Y = Y, X = X, prop = prop, G = G, W = W, method = nuisance.learner.method, train.index = seq(1, n), test.index = seq(1, n))
  # overfitting
  # nuisance.hat = list(); nuisance.hat$mu0.hat = Y0 * (1 - W);  nuisance.hat$mu1.hat = Y1 * W;  nuisance.hat$tau.hat = nuisance.hat$mu1.hat - nuisance.hat$mu0.hat; nuisance.hat$mu.hat = (nuisance.hat$mu1.hat - nuisance.hat$mu0.hat)/2; 
  
  # Calculate test statistics for the original data
  test.stats.value = test.stats.group(Y = Y, W = W, X = X, Group = Group, G = G, stats = test.stats.method, prop = prop, mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat, W.tilde = W.tilde)
  
  # Generate reference test statistics using permutations
  test.stats.ref = t(replicate(M, test.stats.group(Y = Y, W = rbinom(n, 1, prop), X = X, G = G, Group = Group, stats = test.stats.method, prop = prop, mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat, W.tilde = W.tilde))) #   Each column represents a group 
  
  # Calculate p-values for each group
  pval = sapply(seq(1, length(test.stats.value)), function(x) {
    permutation.p.value(stats = test.stats.value[x], stats.ref = test.stats.ref[, x])
  })
  
  return(list(pval = pval, test.stats = test.stats.value, test.stats.ref = test.stats.ref))
}


