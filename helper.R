# helper functions for panning

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

# permutation p-value
# the larger the stats, the smaller the p-value
# permutation.p.value(2, c(1,NA,3))
permutation.p.value = function(stats, stats.ref){
  if(is.na(stats)){return(1)}
  return((sum(stats <= stats.ref[!is.na(stats.ref)]) + 1) / (length(stats.ref[!is.na(stats.ref)]) + 1))
}

# nuisance learner
nuisance.learner = function(Y, X = NULL, prop = NULL, G = NULL, W = NULL, method = "linear", train.index = NULL, test.index = NULL, ...){
  # TODO: other nuisance estimation method
  n = length(Y)
  data.train = data.frame(Y, X, G, W - prop, (W-prop) * X, (W-prop) * G)
  data.0 = data.frame(Y, X, G, 0 - prop, (0 - prop) * X, (0 - prop) * G); colnames(data.0) = colnames(data.train)
  data.1 = data.frame(Y, X, G, 1 - prop, (1 - prop) * X, (1 - prop) * G); colnames(data.1) = colnames(data.train)
  nuisance.model = lm(Y ~ ., data = data.train[train.index,])
  mu0.hat = predict(nuisance.model, newdata = data.0[test.index, ])
  mu1.hat = predict(nuisance.model, newdata = data.1[test.index, ])
  mu.hat = mu0.hat * (1 - prop) + mu1.hat * prop
  tau.hat = mu1.hat -  mu0.hat
  return(result = list(mu0.hat = mu0.hat, mu1.hat = mu1.hat, mu.hat = mu.hat, tau.hat = tau.hat)) 
}

nuisance.learner = function(Y, X = NULL, prop = NULL, G = NULL, W = NULL, method = "linear", train.index = NULL, test.index = NULL, ...) {
  # TODO: other nuisance estimation method
  
  # Get the number of observations
  n = length(Y)
  
  # Prepare training data
  data.train = data.frame(Y, X, G, W - prop, (W - prop) * X, (W - prop) * G)
  
  # Prepare test data for W = 0 and W = 1; the same structure as training data
  data.0 = data.frame(Y, X, G, 0 - prop, (0 - prop) * X, (0 - prop) * G)
  colnames(data.0) = colnames(data.train)  # Ensure column names match
  
  data.1 = data.frame(Y, X, G, 1 - prop, (1 - prop) * X, (1 - prop) * G)
  colnames(data.1) = colnames(data.train)  # Ensure column names match
  
  # Fit a linear model to the training data
  nuisance.model = lm(Y ~ ., data = data.train[train.index, ])
  
  # Predict outcomes for W = 0 and W = 1 using the fitted model
  mu0.hat = predict(nuisance.model, newdata = data.0[test.index, ])
  mu1.hat = predict(nuisance.model, newdata = data.1[test.index, ])
  
  # Compute the predicted marginal mean function
  mu.hat = mu0.hat * (1 - prop) + mu1.hat * prop
  
  # Compute the treatment effect estimate
  tau.hat = mu1.hat - mu0.hat
  
  # Return a list of estimated nuisance functions
  return(result = list(mu0.hat = mu0.hat, mu1.hat = mu1.hat, mu.hat = mu.hat, tau.hat = tau.hat))
}



# test statistics: reject large test statistic
# plain: abs(sum W'_i Y_i / sum W'_i - sum (1 - W'_i) Y_i / sum (1 - W'_i)
# denoise: sum W'_i(Y_i - mu.hat) / sum W'_i - sum (1 - W'_i) (Y_i - mu.hat) / sum (1 - W'_i)
# average treatment effect comparison (ATE): abs(sum W'_i Y_i / sum W'_i - sum (1 - W'_i) Y_i / sum (1 - W'_i) - mean(tau.hat))
# denoise + ATE: abs(sum W'_i (Y_i - mu.hat) / sum W'_i - sum (1 - W'_i) (Y_i - mu.hat) / sum (1 - W'_i) - mean(tau.hat))
# individual treatment effect comparison (ITE): sum_i abs(W'_i(Y_i - mu1.hat) + sum_i abs((1-W'_i)(Y_i - mu0.hat) 
# Example: 
#   mu0.hat = mu.hat = Y0; mu1.hat = Y1; tau.hat = mu1.hat - mu0.hat
#   test.stats(Y, W, X = X, G = G, "denoise", mu.hat)
#   test.stats(Y, W, X = X, G = G, "ATE", tau.hat)
#   test.stats(Y, W, X = X, G = G, "denoise + ATE", mu.hat, tau.hat)
#   test.stats(Y, sample(W), X = X, G = G, "ITE", mu0.hat, mu1.hat)
test.stats = function(Y, W, X = NULL, G = NULL, stats = "denoise", prop = NULL, mu0.hat = NULL, mu1.hat = NULL, mu.hat = NULL, tau.hat = NULL) {
  if (stats == "plain") {
    # Absolute value of the plain difference in means
    # value = abs(sum(W * Y) / max(1, sum(W)) - sum((1 - W) * Y) / max(1, sum(1 - W)))
    n = length(Y)
    value = abs(sum(W * Y) / (n * prop) - sum((1 - W) * Y) / (n * (1 - prop)))
  }else if (stats == "denoise") {
    # Absolute value of the difference in means with denoising using mu.hat
    value = abs(sum(W * (Y - mu.hat)) / max(1, sum(W)) - sum((1 - W) * (Y - mu.hat)) / max(1, sum(1 - W)))
    
  }else if (stats == "ATE") {
    # Absolute difference between the difference in means and the estimated ATE
    value = -abs(sum(W * Y) / max(1, sum(W)) - sum((1 - W) * Y) / max(1, sum(1 - W)) - mean(tau.hat))
  }else if (stats == "denoise + ATE") {
    # Negative absolute difference between the difference in means with denoising using mu.hat and the estimated ATE
    value = -abs(sum(W * (Y - mu.hat)) / max(1, sum(W)) - sum((1 - W) * (Y - mu.hat)) / max(1, sum(1 - W)) - mean(tau.hat))
  }else if(stats == "AIPW"){
    # value = abs(sum(W * (Y - mu1.hat)) / max(1, sum(W)) - sum((1 - W) * (Y - mu0.hat)) / max(1, sum(1 - W)) + mean(tau.hat))
    n = length(Y)
    value1  = abs(sum(W * (Y - mu1.hat)) / (n * prop) - sum((1 - W) * (Y - mu0.hat)) / (n * (1 - prop)) + mean(tau.hat)) / sqrt(1 / (n * prop) + 1 / (n * (1 - prop)))
  }else if (stats == "ITE") {
    # Average absolute difference between the outcome and the estimated nuisance function
    # value = -mean(abs(W * (Y - mu1.hat) + (1 - W) * (Y - mu0.hat)))
    value = mean(abs((W * (Y - mu1.hat)) - (1 - W) * (Y - mu0.hat) + tau.hat)) 
  }else if(stats == "AIPW + ITE"){
    # value1  = abs(sum(W * (Y - mu1.hat)) / max(1, sum(W)) - sum((1 - W) * (Y - mu0.hat)) / max(1, sum(1 - W)) + mean(tau.hat)) / sqrt(min(1, sum(W)) / max(1, sum(W)) +  min(1, sum(W)) / max(1, sum(1 - W)))
    n = length(Y)
    value1  = abs(sum(W * (Y - mu1.hat)) / (n * prop) - sum((1 - W) * (Y - mu0.hat)) / (n * (1 - prop)) + mean(tau.hat)) / sqrt(1 / (n * prop) + 1 / (n * (1 - prop)))
    # value2 = - mean(abs(W * (Y - mu1.hat) + (1 - W) * (Y - mu0.hat)))
    value2 = mean(abs((W * (Y - mu1.hat)) - (1 - W) * (Y - mu0.hat) + tau.hat)) 
    value = value1 + value2 
  }
  
  # Return the computed statistic value
  return(value)
}


# Example:
  # set.seed(318)
  # test.stats.group(Y, sample(W), X, G, Group, mu0.hat = Y0, mu1.hat = Y1, mu.hat = (Y0 + Y1)/2, tau.hat = Y1 - Y0,  stats = "ITE")
test.stats.group = function(Y, W, X = NULL, G = NULL, Group, stats = "denoise", prop = NULL, mu0.hat = NULL, mu1.hat = NULL, mu.hat = NULL, tau.hat = NULL) {
  # Get the unique levels of the grouping variable
  Group.level = sort(unique(Group))
  
  # Compute the test statistic value for each group level
  values = sapply(Group.level, function(x) {
    index = which(Group == x)
    test.stats(Y = Y[index], W = W[index], X = X[index,], G = G[index,], stats = stats, prop = prop, mu0.hat = mu0.hat[index], mu1.hat = mu1.hat[index], mu.hat = mu.hat[index], tau.hat = tau.hat[index])
  })
  
  # Return the computed values for each group level
  return(values)
}


# oracle
  # Example:
  # set.seed(318)
  # ORT(Y = Y, X = X, G = G, Group = Group, prop = p, mu0 = mu0, mu1 = mu1, mu = mu, tau = tau, M = 100, test.stats.method = "denoise")$pval
ORT = function(Y, X, G, Group, prop = NULL, mu0, mu1, mu, tau, test.stats.method = "denoise",  treatment.assignment = "Bernoulli", M = 10) {
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
  # RT(Y = Y, X = X, G = G, Group = Group, prop = p, M = 100)$pval
RT = function(Y, X, G, Group, prop = NULL, treatment.assignment = "Bernoulli", M = 10){
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

# double dipping
# Example:
# set.seed(318)
# DD(Y = Y, X = X, G = G, Group = Group, prop = p, M = 100, test.stats.method = "denoise")$pval
DD = function(Y, X, G, Group, prop = NULL, nuisance.learner.method = "linear", test.stats.method = "denoise", treatment.assignment = "Bernoulli", M = 10){
  # Convert X and G to matrix format
  X = as.matrix(X)
  G = as.matrix(G)
  
  # Get the number of observations
  n = length(Y)
  
  # Estimate the nuisance functions
  nuisance.hat = nuisance.learner(Y = Y, X = X, prop = prop, G = G, W = W, method = nuisance.learner.method, train.index = seq(1, n), test.index = seq(1, n)) 
  
  # Calculate test statistics for the original data
  test.stats.value = test.stats.group(Y = Y, W = W, X = X, Group = Group, G = G, stats = test.stats.method, prop = prop, mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat)

  # Generate reference test statistics using permutations
  test.stats.ref = t(replicate(M, test.stats.group(Y = Y, W = rbinom(n, 1, prop), X = X, G = G, Group = Group, stats = test.stats.method, prop = prop, mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat))) #  Each column represents a group 

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
  # SS(Y = Y, X = X, G = G, Group = Group, prop = p, M = 100, test.stats.method = "denoise")$pval
SS = function(Y, X, G, Group, prop = NULL, nuisance.learner.method = "linear", test.stats.method = "denoise", treatment.assignment = "Bernoulli", M = 10){
  # Convert X and G to matrix format
  X = as.matrix(X)
  G = as.matrix(G)
  
  # Get the number of observations
  n = length(Y)
  
  # Estimate the nuisance functions
  nuisance.index = sample(n, n * 0.5) # TODO 
  train.index = nuisance.index; test.index = setdiff(seq(1, n), nuisance.index)
  nuisance.hat = nuisance.learner(Y = Y, X = X, prop = prop, G = G, W = W, method = nuisance.learner.method, train.index = train.index, test.index = test.index)
  
  # Calculate test statistics for the original data
  test.stats.value = test.stats.group(Y = Y[test.index], W = W[test.index], X = X[test.index,], Group = Group[test.index], G = G[test.index,], stats = test.stats.method, prop = prop, mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat)
  
  # Generate reference test statistics using permutations
  test.stats.ref = t(replicate(M, test.stats.group(Y = Y[test.index], W = rbinom(length(test.index), 1, prop), X = X[test.index,], G = G[test.index,], Group = Group[test.index], stats = test.stats.method, prop = prop, mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat))) #   Each column represents a group 

  # Calculate p-values for each group
  pval = sapply(seq(1, length(test.stats.value)), function(x) {
    permutation.p.value(stats = test.stats.value[x], stats.ref = test.stats.ref[, x])
  })
  
  # Return the p value, the test statistic, and the test statistic values used to compute the reference distribution
  return(list(pval = pval, test.stats = test.stats.value, test.stats.ref = test.stats.ref))
}


# double-dipping
  # Example:
  # set.seed(318)
  # DD(Y = Y, X = X, G = G, Group = Group, prop = p, M = 100, test.stats.method = "denoise")$pval
DD = function(Y, X, G, Group, prop = NULL, nuisance.learner.method = "linear", test.stats.method = "denoise", treatment.assignment = "Bernoulli", M = 10){
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
  test.stats.value = test.stats.group(Y = Y, W = W, X = X, Group = Group, G = G, stats = test.stats.method, prop = prop, mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat)
  
  # Generate reference test statistics using permutations
  test.stats.ref = t(replicate(M, test.stats.group(Y = Y, W = rbinom(n, 1, prop), X = X, G = G, Group = Group, stats = test.stats.method, prop = prop, mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat))) #   Each column represents a group 

  # Calculate p-values for each group
  pval = sapply(seq(1, length(test.stats.value)), function(x) {
    permutation.p.value(stats = test.stats.value[x], stats.ref = test.stats.ref[, x])
  })
  
  return(list(pval = pval, test.stats = test.stats.value, test.stats.ref = test.stats.ref))
}

# panning
  # Example:
  # set.seed(318)
  # ART(Y = Y, X = X, G = G, Group = Group, prop = p, M = 100, test.stats.method = "ATE")$pval
ART = function(Y, X, G, Group, prop = NULL, nuisance.learner.method = "linear", test.stats.method = "denoise", treatment.assignment = "Bernoulli", M = 10, B = 1){
  # Convert X and G to matrix format
  X = as.matrix(X)
  G = as.matrix(G)
  
  # Get the number of observations
  n = length(Y)
  
  # Estimate the nuisance functions
  W.knockoff = matrix(rbinom(n * B, 1, p), ncol = B) # TODO
  W.aug = cbind(W, W.knockoff)
  W.tilde = apply(W.aug, 1, mean) # TODO more knockoffs
  nuisance.hat = nuisance.learner(Y = Y, X = X, prop = prop, G = G, W = W.tilde, method = nuisance.learner.method, train.index = seq(1, n), test.index = seq(1, n))
  
  # Calculate test statistics for the original data
  test.stats.value = test.stats.group(Y = Y, W = W, X = X, Group = Group, G = G, stats = test.stats.method, prop = prop, mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat)
  
  # Generate reference test statistics using permutations
  test.stats.ref = t(replicate(M, test.stats.group(Y = Y, W = W.aug[cbind(seq(1, n), replicate(n, sample(1 + B, 1)))], X = X, G = G, Group = Group, stats = test.stats.method, prop = prop, mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat))) #   Each column represents a group 

    # Calculate p-values for each group
    pval = sapply(seq(1, length(test.stats.value)), function(x) {
      permutation.p.value(stats = test.stats.value[x], stats.ref = test.stats.ref[, x])
    })
    
    return(list(pval = pval, test.stats = test.stats.value, test.stats.ref = test.stats.ref))
}


