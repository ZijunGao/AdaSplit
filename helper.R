# helper functions for panning

# Changes: Aug. 14
  # Replace mu.hat by mu in "gradient boosting early stopping"
  # Added sample splitting proportion rho (nuisance proportion) in function SS
  # replaced prop by rep(prop, n)
  # replaced prop by W.tilde for test statistic computation in ART
  # ? true tau and mu

library("xgboost")
library(caret)
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
  n = length(Y)
  if(length(prop) == 1){prop = rep(prop, n)}

  if(method == "linear"){
    data.train = data.frame(Y, X, G, W - prop, (W-prop) * X, (W-prop) * G)
    data.0 = data.frame(Y, X, G, 0 - prop, (0 - prop) * X, (0 - prop) * G); colnames(data.0) = colnames(data.train)
    data.1 = data.frame(Y, X, G, 1 - prop, (1 - prop) * X, (1 - prop) * G); colnames(data.1) = colnames(data.train)
    nuisance.model = lm(Y ~ ., data = data.train[train.index,])
    mu0.hat = predict(nuisance.model, newdata = data.0[test.index, ])
    mu1.hat = predict(nuisance.model, newdata = data.1[test.index, ])
    mu.hat = mu0.hat * (1 - prop[test.index]) + mu1.hat * prop[test.index]
    tau.hat = mu1.hat -  mu0.hat
    
    return(result = list(mu0.hat = mu0.hat, mu1.hat = mu1.hat, mu.hat = mu.hat, tau.hat = tau.hat)) 
    
  }else if(method == "gradient boosting"){
    
    data.full = data.frame(X, W - prop,Y)
    train.index_1 <- createDataPartition(data.full$Y[train.index], p = 0.9, list = FALSE) # 90% training data
    train.data <- data.full[train.index,]
    test.data <- data.full[test.index,]
    
    num_cols <- ncol(data.full)
    features = 1:(num_cols-2)
    # Define the data matrices
    dtrain = xgb.DMatrix(data = as.matrix(train.data[, features]), label = train.data$Y)
    # Define the parameters
    params = list(objective = "reg:squarederror", eval_metric = "rmse")


    # Fit the model with early stopping
    nrounds = 100 
    
    nuisance.mu <- xgb.train(
      #booster = "gblinear",
      params = params,
      data = dtrain,
      nrounds = nrounds,
      verbose = 0
    )
    
    
    # Make predictions on the training and validation data
    train.pred <- predict(nuisance.mu, newdata = dtrain)
    
    # Calculate residuals for training and validation data
    train.residual <- (train.data$Y - train.pred) / train.data[, num_cols-1] # train.data[, num_cols-1] could be zero

    dtest = xgb.DMatrix(data = as.matrix(test.data[,features]))
    

    # Define the data matrices
    dtrain <- xgb.DMatrix(data = as.matrix(train.data[, features]), label = train.residual)


    weights <- abs(train.data[, num_cols-1])**2
    setinfo(dtrain, "weight", weights)
    
    weighted_squared_error <- function(preds, dtrain) {
      labels <- getinfo(dtrain, "label")
      weights <- getinfo(dtrain, "weight")
      
      # Calculate the gradient and hessian
      grad <- weights * (preds - labels)
      hess <- weights
      
      return(list(grad = grad, hess = hess))
    }
    
    # Define the parameters for the second model
    params <- list(eval_metric = "rmse")
    
    # Fit the second model with early stopping
    nuisance.tau <- xgb.train(
      #booster = "gblinear",
      params = params,
      data = dtrain,
      nrounds = nrounds,
      verbose = 0,
      obj = weighted_squared_error
    )
    
    
    #dtest.0 = xgb.DMatrix(data = as.matrix(data.0[test.index, -1]))
    #dtest.1 = xgb.DMatrix(data = as.matrix(data.1[test.index, -1]))
    
    # Define the data matrix for the test data
    dtest <- xgb.DMatrix(data = as.matrix(test.data[, 1:(num_cols-2)]))
    
    # Make predictions for the test data using the first and second models
    mu.hat <- predict(nuisance.mu, newdata = dtest)
    tau.hat <- predict(nuisance.tau, newdata = dtest)
    
    # Calculate the final predictions for mu0.hat and mu1.hat
    mu0.hat <- mu.hat - prop[test.index] * tau.hat
    mu1.hat <- mu0.hat + tau.hat
    
    #mu0.hat = predict(nuisance.model, newdata = dtest.0)
    #mu1.hat = predict(nuisance.model, newdata = dtest.1)
    return(result = list(mu0.hat = mu0.hat, mu1.hat = mu1.hat, mu.hat = mu.hat, tau.hat = tau.hat,tau=nuisance.tau)) 
  }
  else if(method == "gradient boosting early stopping"){
    
    data.full = data.frame(X, W - prop,Y)
    train.index_1 <- createDataPartition(data.full$Y[train.index], p = 0.9, list = FALSE) # 90% training data
    train.data <- data.full[train.index,][train.index_1, ]
    val.data <- data.full[train.index,][-train.index_1, ]
    test.data <- data.full[test.index,]
    
    num_cols <- ncol(data.full)
    features = 1:(num_cols-2)
    # Define the data matrices
    dtrain = xgb.DMatrix(data = as.matrix(train.data[, features]), label = train.data$Y)
    dval = xgb.DMatrix(data = as.matrix(val.data[, features]), label = val.data$Y)
    # Define the parameters
    params = list(objective = "reg:squarederror", eval_metric = "rmse")
    
    #Watchlist to track performance on validation set
    watchlist = list(train = dtrain, eval = dval)
    
    # Fit the model with early stopping
    nrounds = 100  # Maximum number of boosting rounds
    early_stopping_rounds =10  # Stop early if there is no improvement
    
    
    nuisance.mu <- xgb.train(
      #booster = "gblinear",
      params = params,
      data = dtrain,
      nrounds = nrounds,
      #lambda = 1,
      watchlist = watchlist,
      early_stopping_rounds = early_stopping_rounds,
      verbose = 0
    )
    
    
    # Make predictions on the training and validation data
    # train.pred <- predict(nuisance.mu, newdata = dtrain)
    # val.pred <- predict(nuisance.mu, newdata = dval)
    # to be deleted start
    train.pred = mu[train.index][train.index_1]
    val.pred = mu[train.index][-train.index_1]
    # to be deleted end
    
    # Calculate residuals for training and validation data
    train.residual <- (train.data$Y - train.pred) / train.data[, num_cols-1] # train.data[, num_cols-1] could be zero
    val.residual <- (val.data$Y - val.pred) / val.data[, num_cols-1]
    
    
    dtest = xgb.DMatrix(data = as.matrix(test.data[,features]))
    
    
    # Define the data matrices
    dtrain <- xgb.DMatrix(data = as.matrix(train.data[, features]), label = train.residual)
    dval <- xgb.DMatrix(data = as.matrix(val.data[, features]), label = val.residual)
    
    watchlist = list(train = dtrain, eval = dval)
    weights <- abs(train.data[, num_cols-1])**2
    setinfo(dtrain, "weight", weights)
    
    weighted_squared_error <- function(preds, dtrain) {
      labels <- getinfo(dtrain, "label")
      weights <- getinfo(dtrain, "weight")
      
      # Calculate the gradient and hessian
      grad <- weights * (preds - labels)
      hess <- weights
      
      return(list(grad = grad, hess = hess))
    }
    
    # Define the parameters for the second model
    params <- list(eval_metric = "rmse")
    
    # Fit the second model with early stopping
    nuisance.tau <- xgb.train(
      #booster = "gblinear",
      params = params,
      data = dtrain,
      nrounds = nrounds,
      #lambda = 1,
      watchlist = watchlist,
      early_stopping_rounds = early_stopping_rounds,
      verbose = 0,
      obj = weighted_squared_error
    )
    
    
    #dtest.0 = xgb.DMatrix(data = as.matrix(data.0[test.index, -1]))
    #dtest.1 = xgb.DMatrix(data = as.matrix(data.1[test.index, -1]))
    
    # Define the data matrix for the test data
    dtest <- xgb.DMatrix(data = as.matrix(test.data[, 1:(num_cols-2)]))
    
    # Make predictions for the test data using the first and second models
    # mu.hat <- predict(nuisance.mu, newdata = dtest)
    mu.hat = mu[test.index]  # to be deleted
    tau.hat <- predict(nuisance.tau, newdata = dtest)
    
    # Calculate the final predictions for mu0.hat and mu1.hat
    mu0.hat <- mu.hat - prop[test.index] * tau.hat
    mu1.hat <- mu0.hat + tau.hat
    
    #mu0.hat = predict(nuisance.model, newdata = dtest.0)
    #mu1.hat = predict(nuisance.model, newdata = dtest.1)
    return(result = list(mu0.hat = mu0.hat, mu1.hat = mu1.hat, mu.hat = mu.hat, tau.hat = tau.hat,tau=nuisance.tau)) 
  }else if(method == "gradient boosting linear"){
    
    data.full = data.frame(X, W - prop,Y)
    train.index_1 <- createDataPartition(data.full$Y[train.index], p = 0.9, list = FALSE) # 90% training data
    train.data <- data.full[train.index,][train.index_1, ]
    val.data <- data.full[train.index,][-train.index_1, ]
    test.data <- data.full[test.index,]
    
    num_cols <- ncol(data.full)
    features = 1:(num_cols-2)
    # Define the data matrices
    dtrain = xgb.DMatrix(data = as.matrix(train.data[, features]), label = train.data$Y)
    dval = xgb.DMatrix(data = as.matrix(val.data[, features]), label = val.data$Y)
    # Define the parameters
    params = list(objective = "reg:squarederror", eval_metric = "rmse")
    
    #Watchlist to track performance on validation set
    watchlist = list(train = dtrain, eval = dval)
    
    # Fit the model with early stopping
    nrounds = 100  # Maximum number of boosting rounds

    
    nuisance.mu <- xgb.train(
      booster = "gblinear",
      params = params,
      data = dtrain,
      nrounds = nrounds,
      #lambda = 1,
      watchlist = watchlist,
      # early_stopping_rounds = early_stopping_rounds,
      verbose = 0
    )
    
    
    # Make predictions on the training and validation data
    # train.pred <- predict(nuisance.mu, newdata = dtrain)
    # val.pred <- predict(nuisance.mu, newdata = dval)
    # to be deleted start
    train.pred = mu[train.index][train.index_1]
    val.pred = mu[train.index][-train.index_1]
    # to be deleted end
    
    # Calculate residuals for training and validation data
    train.residual <- (train.data$Y - train.pred) / train.data[, num_cols-1] # train.data[, num_cols-1] could be zero
    val.residual <- (val.data$Y - val.pred) / val.data[, num_cols-1]
    
    
    dtest = xgb.DMatrix(data = as.matrix(test.data[,features]))
    
    
    # Define the data matrices
    dtrain <- xgb.DMatrix(data = as.matrix(train.data[, features]), label = train.residual)
    dval <- xgb.DMatrix(data = as.matrix(val.data[, features]), label = val.residual)
    
    watchlist = list(train = dtrain, eval = dval)
    weights <- abs(train.data[, num_cols-1])**2
    setinfo(dtrain, "weight", weights)
    
    weighted_squared_error <- function(preds, dtrain) {
      labels <- getinfo(dtrain, "label")
      weights <- getinfo(dtrain, "weight")
      
      # Calculate the gradient and hessian
      grad <- weights * (preds - labels)
      hess <- weights
      
      return(list(grad = grad, hess = hess))
    }
    
    # Define the parameters for the second model
    params <- list(eval_metric = "rmse")
    
    # Fit the second model with early stopping
    nuisance.tau <- xgb.train(
      booster = "gblinear",
      params = params,
      data = dtrain,
      nrounds = nrounds,
      #lambda = 1,
      watchlist = watchlist,
      # early_stopping_rounds = early_stopping_rounds,
      verbose = 0,
      obj = weighted_squared_error
    )
    
    
    #dtest.0 = xgb.DMatrix(data = as.matrix(data.0[test.index, -1]))
    #dtest.1 = xgb.DMatrix(data = as.matrix(data.1[test.index, -1]))
    
    # Define the data matrix for the test data
    dtest <- xgb.DMatrix(data = as.matrix(test.data[, 1:(num_cols-2)]))
    
    # Make predictions for the test data using the first and second models
    # mu.hat <- predict(nuisance.mu, newdata = dtest)
    mu.hat = mu[test.index]  # to be deleted
    tau.hat <- predict(nuisance.tau, newdata = dtest)
    
    # Calculate the final predictions for mu0.hat and mu1.hat
    mu0.hat <- mu.hat - prop[test.index] * tau.hat
    mu1.hat <- mu0.hat + tau.hat # 
    
    #mu0.hat = predict(nuisance.model, newdata = dtest.0)
    #mu1.hat = predict(nuisance.model, newdata = dtest.1)
    return(result = list(mu0.hat = mu0.hat, mu1.hat = mu1.hat, mu.hat = mu.hat, tau.hat = tau.hat,tau=nuisance.tau)) 
  }
  
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
  
  n = length(Y)
  if (is.null(prop)){
    n1 = max(1, sum(W))
    n0 = max(1, sum(1 - W))
  }else{
    n1 = n * prop
    n0 = n * (1 - prop)
  }
  

  if (stats == "plain") {
    # Absolute value of the plain difference in means
    # value = abs(sum(W * Y) / max(1, sum(W)) - sum((1 - W) * Y) / max(1, sum(1 - W)))
    #value = abs(sum(W * Y) / n1 - sum((1 - W) * Y) / n0)
    
    IF = W * (Y)/prop - (1 - W) * (Y)/(1-prop)
    value =  abs(mean(IF)/sd(IF))
    

  }else if (stats == "denoise" | stats == "denoise normalized") {
    # Absolute value of the difference in means with denoising using mu.hat
   #value = abs(sum(W * (Y - mu.hat)) / n1 - sum((1 - W) * (Y - mu.hat)) / n0)
    
    IF = W * (Y - mu.hat)/prop - (1 - W) * (Y - mu.hat)/(1-prop)
    value =  abs(mean(IF)/sd(IF))
    
  }else if (stats == "ATE") {
    # Absolute difference between the difference in means and the estimated ATE
    value = -abs(sum(W * Y) / n1 - sum((1 - W) * Y) / n0 - mean(tau.hat))
  }else if (stats == "denoise + ATE") {
    # Negative absolute difference between the difference in means with denoising using mu.hat and the estimated ATE
    value = -abs(sum(W * (Y - mu.hat)) / n1 - sum((1 - W) * (Y - mu.hat)) / n0 - mean(tau.hat))
  }else if(stats == "AIPW" | stats == "AIPW normalized"){
    IF = W * (Y - mu1.hat)/(prop + (prop == 0)) - (1 - W) * (Y - mu0.hat)/(1-prop + (prop == 1)) + tau.hat
    value =  abs(mean(IF)/sd(IF))
    #IF2 =  W * (Y - mu.hat)/prop - (1 - W) * (Y - mu.hat)/(1-prop)
    #value2 =  abs(mean(IF2)/sd(IF2))
    #value = max(value,value2)
  }else if (stats == "ITE") {
    # Average absolute difference between the outcome and the estimated nuisance function
    
    IF = W * (Y - mu1.hat)/prop - (1 - W) * (Y - mu0.hat)/(1-prop) + tau.hat
    value = mean(abs(IF))/sd(abs(IF))
    
  }else if(stats == "AIPW + ITE"){
    #value1  = abs(sum(W * (Y - mu1.hat)) / n1 - sum((1 - W) * (Y - mu0.hat)) / n0 + mean(tau.hat)) / sqrt(1 / n1 + 1 /n0)
    # value2 = - mean(abs(W * (Y - mu1.hat) + (1 - W) * (Y - mu0.hat)))
   # value2 = mean(abs((W * (Y - mu1.hat)) - (1 - W) * (Y - mu0.hat) + tau.hat)) 
    
    IF = W * (Y - mu1.hat)/prop - (1 - W) * (Y - mu0.hat)/(1-prop) 
    value = abs(mean(IF)/sd(IF)) + mean(abs(IF))/sd(abs(IF)) #value1 + value2 
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
  if(length(prop) == 1){prop = rep(prop, length(Y))}
  # Compute the test statistic value for each group level
  values = sapply(Group.level, function(x) {
    index = which(Group == x)
    test.stats(Y = Y[index], W = W[index], X = X[index,], G = G[index,], stats = stats, prop = prop[index], mu0.hat = mu0.hat[index], mu1.hat = mu1.hat[index], mu.hat = mu.hat[index], tau.hat = tau.hat[index])
  })
  
  # Return the computed values for each group level
  return(values)
}


# oracle
  # Example:
  # set.seed(318)
  # ORT(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, mu0 = mu0, mu1 = mu1, mu = mu, tau = tau, M = 100, test.stats.method = "denoise")$pval
ORT = function(Y, W, X, G, Group, prop = NULL, mu0, mu1, mu, tau, test.stats.method = "denoise",  treatment.assignment = "Bernoulli", M = 10) {
  # Convert X and G to matrix format
  X = as.matrix(X)
  G = as.matrix(G)
  
  # Get the number of observations
  n = length(Y)
  
  # create prop vector
  if(length(prop) == 1){prop = rep(prop, n)}
  
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
RT = function(Y, W, X, G, Group, prop = NULL, treatment.assignment = "Bernoulli", M = 10){
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
# DD(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = 100, test.stats.method = "denoise")$pval
DD = function(Y, W, X, G, Group, prop = NULL, nuisance.learner.method = "linear", test.stats.method = "denoise", treatment.assignment = "Bernoulli", M = 10){
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
  # SS(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = 100, test.stats.method = "denoise")$pval
SS = function(Y, W, X, G, Group, prop = NULL, nuisance.learner.method = "linear", test.stats.method = "denoise", treatment.assignment = "Bernoulli", M = 10, rho = 0.5){
  # Convert X and G to matrix format
  X = as.matrix(X)
  G = as.matrix(G)
  
  # Get the number of observations
  n = length(Y)
  
  # create prop vector
  if(length(prop) == 1){prop = rep(prop, n)}

  
  # Estimate the nuisance functions
  nuisance.index = sample(n, n * rho) # TODO 
  train.index = nuisance.index; test.index = setdiff(seq(1, n), nuisance.index)
  nuisance.hat = nuisance.learner(Y = Y, X = X, prop = prop, G = G, W = W, method = nuisance.learner.method, train.index = train.index, test.index = test.index) # here!!!
  
  # Calculate test statistics for the original data
  # test.stats.value = test.stats.group(Y = Y[test.index], W = W[test.index], X = X[test.index,], Group = Group[test.index], G = G[test.index,], stats = test.stats.method, prop = prop[test.index], mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat)
  test.stats.value = test.stats.group(Y = Y[test.index], W = W[test.index], X = X[test.index,], Group = Group[test.index], G = G[test.index,], stats = test.stats.method, prop = prop[test.index], mu0.hat = mu0[test.index], mu1.hat = mu1[test.index], mu.hat = mu[test.index], tau.hat = tau[test.index])
  
  # Generate reference test statistics using permutations
  # test.stats.ref = t(replicate(M, test.stats.group(Y = Y[test.index], W = rbinom(length(test.index), 1, prop), X = X[test.index,], G = G[test.index,], Group = Group[test.index], stats = test.stats.method, prop = prop[test.index], mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat))) #   Each column represents a group 
  test.stats.ref = t(replicate(M, test.stats.group(Y = Y[test.index], W = rbinom(length(test.index), 1, prop), X = X[test.index,], G = G[test.index,], Group = Group[test.index], stats = test.stats.method, prop = prop[test.index], mu0.hat = mu0[test.index], mu1.hat = mu1[test.index], mu.hat = mu[test.index], tau.hat = tau[test.index]))) #   Each column represents a group 
  
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
  # DD(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = 100, test.stats.method = "denoise")$pval
DD = function(Y, W, X, G, Group, prop = NULL, nuisance.learner.method = "linear", test.stats.method = "denoise", treatment.assignment = "Bernoulli", M = 10){
  # Convert X and G to matrix format
  X = as.matrix(X)
  G = as.matrix(G)
  
  # Get the number of observations
  n = length(Y)
  
  # create prop vector
  if(length(prop) == 1){prop = rep(prop, n)}
  
  
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
  # ART(Y = Y, W = W, X = X, G = G, Group = Group, prop = p, M = 100, test.stats.method = "ATE")$pval
ART = function(Y, W, X, G, Group, prop = NULL, nuisance.learner.method = "linear", test.stats.method = "denoise", treatment.assignment = "Bernoulli", M = 10, B = 1){
  # Convert X and G to matrix format
  X = as.matrix(X)
  G = as.matrix(G)
  
  # Get the number of observations
  n = length(Y)
  
  # create prop vector
  if(length(prop) == 1){prop = rep(prop, n)}
  
  
  # Estimate the nuisance functions
  if (is.null(prop)){
    W.knockoff = permute_within_groups(W, Group, B)
  }else{
    W.knockoff = matrix(rbinom(n * B, 1, prop), ncol = B) 
  }
  
  W.aug = cbind(W, W.knockoff)
  W.tilde = apply(W.aug, 1, mean) # TODO more knockoffs

  nuisance.hat = nuisance.learner(Y = Y, X = X, prop = prop, G = G, W = W.tilde, method = nuisance.learner.method, train.index = seq(1, n), test.index = seq(1, n))
  
  # Calculate test statistics for the original data
  # test.stats.value = test.stats.group(Y = Y, W = W, X = X, Group = Group, G = G, stats = test.stats.method, prop = W.tilde, mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat) # prop changed to W.tilde
  test.stats.value = test.stats.group(Y = Y, W = W, X = X, Group = Group, G = G, stats = test.stats.method, prop = W.tilde, mu0.hat =  mu0, mu1.hat = mu1, mu.hat = mu, tau.hat = tau) # prop changed to W.tilde
  
  # Generate reference test statistics using permutations
  # test.stats.ref = t(replicate(M, test.stats.group(Y = Y, W = W.aug[cbind(seq(1, n), replicate(n, sample(1 + B, 1)))], X = X, G = G, Group = Group, stats = test.stats.method, prop = W.tilde, mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat))) #   Each column represents a group; # prop changed to W.tilde
  test.stats.ref = t(replicate(M, test.stats.group(Y = Y, W = W.aug[cbind(seq(1, n), replicate(n, sample(1 + B, 1)))], X = X, G = G, Group = Group, stats = test.stats.method, prop = W.tilde, mu0.hat = mu0, mu1.hat = mu1, mu.hat = mu, tau.hat = tau))) #   Each column represents a group; # prop changed to W.tilde
    # Calculate p-values for each group
    pval = sapply(seq(1, length(test.stats.value)), function(x) {
      permutation.p.value(stats = test.stats.value[x], stats.ref = test.stats.ref[, x])
    })
    
    return(list(pval = pval, test.stats = test.stats.value, test.stats.ref = test.stats.ref))
}




permute_within_groups <- function(W, G, B) {
  # Ensure that W and G have the same length
  if (length(W) != length(G)) {
    stop("The length of W and G must be the same.")
  }
  
  # Initialize the matrix to store the original and permuted vectors
  n <- length(W)
  permutations_matrix <- matrix(NA, nrow = n, ncol = B + 1)
  
  # Store the original W as the first column
  permutations_matrix[, 1] <- W
  
  # Get the unique group labels
  unique_groups <- unique(G)
  
  for (b in 1:B) {
    # Create a copy of W to store the permuted values
    W_permuted <- W
    
    # Permute the elements of W within each group
    for (group in unique_groups) {
      # Get the indices of elements in the current group
      group_indices <- which(G == group)
      
      # Permute the elements within the current group
      W_permuted[group_indices] <- sample(W[group_indices])
    }
    
    # Store the permuted vector in the matrix
    permutations_matrix[, b + 1] <- W_permuted
  }
  
  return(permutations_matrix)
}


