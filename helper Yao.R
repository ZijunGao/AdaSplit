# helper functions for panning
library("xgboost")
library(randomForest)
library(caret)
library(ranger)
library(pracma)
library(Matrix)


nuisance.tau.ss = function(Y= NULL, X = NULL, Ex = NULL, W = NULL, mu = NULL, train.index = NULL, k=20, marginalize = T, lambda = 1, ...){
  
  # initial fit
  n = length(Y)
  IW = rep(1,n)
  X = as.matrix(X) # trun data.frame to matrx
  R = Y - mu
  A = knn.indices(X,Y,k)
  X_ = cbind(1, X[train.index,])
  R_ = (R/ (W-Ex))[train.index]
  W_ = diag(c(((W-Ex)**2)[train.index]))
  inv_XTWX = solve(t(X_) %*% W_ %*% X_ + 10e-10*diag(rep(1.0,dim(X_)[2])))
  XTWR = t(X_) %*% W_ %*% R_
  beta = inv_XTWX %*% XTWR
  tau = cbind(1, X) %*% beta
  
  Q = Posterior(train.index, tau, X, R, W, Ex,IW, weighting = FALSE)
  beta.imputed = Posterior_fit(train.index, X, R, W, Ex, Q, A, marginalize = marginalize, lambda = lambda)
  tau.imputed = cbind(1, X) %*% beta.imputed
  
  mu0.hat <- mu - Ex * tau.imputed
  mu1.hat <- mu0.hat + tau.imputed
  
  return(result = list(mu0 = mu0.hat, mu1 = mu1.hat, tau = tau.imputed, beta = beta.imputed))
}


Posterior_fit = function(train.index, X, R, W, Ex, Q, A, marginalize = T, lambda = 1){

  n = length(W)
  k = sum(A[1,]) 
  n.train = length(train.index)
  test.index = setdiff(1:n,train.index)
  
  
  XX <- cbind(1, rbind(X,X))
  W1 = W1_ = Ex**2 ; W2 = W2_ =  (1 - Ex)**2
  W1_[train.index] = W[train.index] * W1[train.index] 
  W1_[test.index] =  W1[test.index]*Q[test.index] * marginalize

  W2_[train.index] = (1-W[train.index]) * W2[train.index]
  W2_[test.index] = W2[test.index]*(1-Q[test.index]) * marginalize

  inv_XTWX = solve(t(XX) %*% diag(c(W1_,W2_)) %*% XX + 10e-10*diag(rep(1.0,dim(XX)[2])))
  beta = inv_XTWX %*% t(XX) %*% diag(c(W1_,W2_)) %*% c( R/(1 - Ex),  R/(0 - Ex))

  return(beta)
}


nuisance.tau.active = function(Y= NULL, X = NULL, Ex = NULL, W = NULL, mu = NULL, Group=NULL, proportion = NULL, groupwise=FALSE, weighting = TRUE, k= 10, p=0.1, robust = F, marginalize = T, lambda = 1, lambda.update.period = Inf, loss_threshold = NULL, initial.proportion = 0.05, ...){
  
  # Y is an n-dimensional outcome matrix
  # X is an n x d covariate matrix 
  # Ex is an n-dimensional treatment assignment probability vector
  # W is an observed treatment assignment vector
  # mu is an outcome regression model
  # Group is the group indicator
  # proportion is the max prop. of data points used for estimation
  # groupwise indicate if we will select points evenly across the groups
  # k is the number of nearest neighbors used in the knn model for correcting selection bias
  # p is the threshold for clipping the inverse probability weights
  # lambda balances the contributions to the loss from units with observed treatment assignments, denoted by loss(observed), and those with imputed treatment assignments, denoted by loss(imputed). Explicitly, loss = loss(observed) + lambda * loss(imputed). Default value is one.
  # lambda.update.period specifies the number of epochs between updates of lambda. The default is Inf, meaning that lambda is not updated.
  # loss_threshold is the cutoff below which the algorithm is allowed to stop. If the max_loss falls below this value, the iteration terminates.


  
  # Compute the sample size
  n = length(Y) 
  # Compute the k-nearest neighbors of every point; the j-th row of A indicates the neighbors of the j-th point
  A = knn.indices(X,Y,k)
  
  X = as.matrix(X) # turn data.frame to matrx

  # Compute the residuals
  R = Y - mu
  
  # Initialization
  n_0 = max(dim(X)[2]+5, round(initial.proportion*n)) # min(dim(X)[2]+5, round(initial.proportion*n))
  train.index = select.train(Group, n_0/n, random = F, X) # Select units in decreasing order of X_j^\top (X_{[n]}^{T}X_{[n]})^{-2}X_j
  test.index = setdiff(1:n,train.index) 
  n.train = length(train.index)
  n.test = length(test.index)
  prop.group.test = prop.group(Group,test.index)
  
  if (weighting){
    v = rep(0,n)
    v[test.index] = 1.0
    p.test = (A %*% v)/k
    p.train = 1 - p.test
    IW =  1 / pmax(p.train,p)
  }else{
    IW = rep(1,n)
  }

  # Compute OLS
  X_ = cbind(1, X[train.index,])
  R_ = (R/ (W-Ex))[train.index]
  W_ = diag(c(((W-Ex)**2 * IW)[train.index]))
  inv_XTWX = solve(t(X_) %*% W_ %*% X_ + 10e-10*diag(rep(1.0,dim(X_)[2])))
  XTWR = t(X_) %*% W_ %*% R_
  beta = inv_XTWX %*% XTWR
  tau = cbind(1, X) %*% beta
  

  
  
  #X_ = cbind(1, X)
  #inv_XTWX = solve(t(X_) %*% diag(c(((W-Ex)**2 ))) %*% X_ + 10e-10*diag(rep(1.0,dim(X_)[2])))
  #inv_XTWX_X = inv_XTWX %*% t(cbind(1, rbind(X,X)))
  
  # Compute Imputed OLS
  Q = Posterior(train.index, tau, X, R, W, Ex,IW, weighting)
  beta.imputed = Posterior_fit(train.index, X, R, W, Ex, Q, A, marginalize = marginalize, lambda = lambda)
  tau.imputed = cbind(1, X) %*% beta.imputed
  
  # Set the stopping rule
  loss <- 1;
  
  if(robust){
    window_size <- 10
    epoch.gap = 20 # 20
    if(is.null(loss_threshold)){loss_threshold <- 0.1}
    # lambda.update.period = 20
  }else{
    window_size <- 50
    epoch.gap = 1
    if(is.null(loss_threshold)){loss_threshold <- 0.01}
    # lambda.update.period = 20
  }
  
  
  loss_history <- rep(Inf, window_size)
  epoch <- 0
  stop = FALSE
  Group.idx = sort(unique(Group))
  lambda.seq = lambda

  
  while ( length(train.index) < ceiling(proportion * n) ) {
    # if(epoch %% 10 == 0){print(epoch)}
    epoch <- epoch + 1
    
    mu0.hat <- mu - Ex * tau.imputed
    mu1.hat <- mu0.hat + tau.imputed
    TE = (Y-mu1.hat)/Ex + (Y-mu0.hat)/(1-Ex)

    # Predict the CRT's power loss after moving the j-th point from inference to nuisance estimation
    power = numeric(n)
    

    power = abs(Ex - Q) *sign(tau.imputed)
    power = power[test.index]

    if(robust){power = (TE^2)[test.index]}
  
    
    # Find out the proportion of inference units in each group
    prop.group.test = prop.group(Group, Group.idx,test.index)
    
    # Find out which group has a large enough inference fold and compute the proportions
    B = prop.group.test > (1-proportion)
    Group.0 = Group.idx[B]

    if (groupwise){
      # If groupwise, randomly choose a group with enough proportion
      prop.group.test = prop.group.test[B]
      g.selected = Group.0[which(rmultinom(n = 1, size = 1, prob = prop.group.test/sum(prop.group.test))==1)]
      # Choose the best point in the group
      idx = Group[test.index] == g.selected
  
    }else{
      # Otherwise, select the best point from the groups with enough proportions
      idx = which(Group[test.index] %in% Group.0)
    }
    
    if(epoch %% Inf == 0){
      # new = sample(test.index[idx], 1)
      new = intersect(order(diag(X %*% solve(t(X) %*% X) %*% solve(t(X) %*% X) %*% t(X)), decreasing = T), test.index[idx])[1]

    }else{
      new = test.index[idx][which.min(power[idx])] # Add a new data point to the nuisance estimation fold
    }

    # Update the indices of the nuisance and inference folds
    train.index = c(train.index, new)
    test.index = setdiff(1:n,train.index)
    
    # Update the OLS model
    n.train = length(train.index); n.test = length(test.index)
    
    beta.copy = beta.imputed
    tau.copy = tau.imputed
    
    
    if (epoch %% epoch.gap ==0){
      

      
      if (weighting){
        
        X_ = cbind(1, X[train.index,])
        R_ = (R/ (W-Ex))[train.index]
        W_ = diag(c(((W-Ex)**2 * IW)[train.index]))
        inv_XTWX = solve(t(X_) %*% W_ %*% X_ + 10e-10*diag(rep(1.0,dim(X_)[2])))
        XTWR = t(X_) %*% W_ %*% R_
        beta = inv_XTWX %*% XTWR
        tau = cbind(1, X) %*% beta
        
        Q = Posterior(train.index, tau, X, R, W, Ex,IW, weighting)
        beta.imputed = Posterior_fit(train.index, X, R, W, Ex, Q, A, marginalize = marginalize, lambda = lambda)
        tau.imputed = cbind(1, X) %*% beta.imputed
        
      }
      
      Q = Posterior(train.index, tau.imputed, X, R, W, Ex, IW=rep(1,n),weighting=FALSE)
      
      beta.imputed = Posterior_fit(train.index, X, R, W, Ex, Q, A, marginalize = marginalize, lambda = lambda)
      tau.imputed = cbind(1, X) %*% beta.imputed
      
      # Compute the estimator change due to the new data point
      loss = sum((tau.imputed[test.index]-tau.copy[test.index])**2)/sum((tau.copy-mean(tau.copy[test.index]))**2)
      #loss = sum((beta.imputed - beta.copy)**2)/sum(beta.copy**2)
      
      loss_history <- c(loss_history[-1], loss)
      max_loss <- max(loss_history)
  
      if (epoch >= window_size) {
        if (max_loss < loss_threshold) { 
          stop =TRUE
          break
        }
      }
    
    
    }
    
    if (stop){break} 
    
    # Use 5-fold cross-validation to choose lambda (start)
    if(epoch %% lambda.update.period == 0){
      lambda.candidate = 5^(seq(-4, 0))
      shuffled = sample(train.index)
      folds = split(shuffled, cut(seq_along(shuffled), breaks = 5, labels = FALSE)) # Split the training data randomly into 5 equal folds
      cv.loss = matrix(0, 5, length(lambda.candidate))
      for(fold.index in 1 : 5){
        validation.index = folds[[fold.index]]
        for(lambda.index in seq(1, length(lambda.candidate))){
          new.train.index = sapply(setdiff(train.index, validation.index), function(x){x - sum(validation.index < x)})
          # print(train.index); print(validation.index); print(new.train.index)
          beta.temp = Posterior_fit(new.train.index, X[-validation.index,,drop = F], R[-validation.index], W[-validation.index], Ex[-validation.index], Q[-validation.index], A[-validation.index,,drop = F], marginalize = marginalize, lambda = lambda.candidate[lambda.index])
          cv.loss[fold.index, lambda.index] = sum((R[validation.index] - (W - Ex)[validation.index] * (cbind(1, X[validation.index,,drop = F]) %*% beta.temp))^2)
        }        
      }
      cv.loss.average = apply(cv.loss, 2, mean)
      cv.loss.se = apply(cv.loss, 2, sd) / sqrt(5) 
      # print(data.frame(loss = cv.loss.average, lse = cv.loss.average + cv.loss.se)) # To be deleted
      # print(cv.loss.average * 5 / length(train.index)) # To be deleted
      lambda = max(lambda.candidate[cv.loss.average < (cv.loss.average + cv.loss.se)[which.min(cv.loss.average)]])
      lambda.seq = c(lambda.seq, lambda)
    }
    # Use 5-fold cross-validation to choose lambda (end)
    
  }
  
  
  if(!robust){
    
    #compute the inference proportion of each group
    prop.group.test = prop.group(Group, Group.idx, test.index)
    #choose the groups with enough proportion
    prop.gap = prop.group.test - (1-proportion)
    Group.0 = Group.idx[prop.gap > 0]
    #collect the units in these groups
    idx = which(Group[test.index] %in% Group.0)
    

    test.index.move = NULL
    for (g in Group.0){
      size.g = sum(Group==g)
      idx.g = which(Group[test.index]==g)
      test.index.g = test.index[idx.g]
      num.negative =  sum(tau.imputed[test.index.g]<0)
      top_indices <- order(tau.imputed[test.index.g],decreasing=FALSE)[1:min(num.negative,floor(size.g*prop.gap[g+1]))]
 
      test.index.move = c(test.index.move,test.index.g[top_indices])
    }
    train.index.before.throw.away = train.index
    train.index = c(train.index, test.index.move)
    test.index = setdiff(1:n, train.index)
    
    
    

    
    if (weighting){
      
      X_ = cbind(1, X[train.index,])
      R_ = (R/ (W-Ex))[train.index]
      W_ = diag(c(((W-Ex)**2 * IW)[train.index]))
      inv_XTWX = solve(t(X_) %*% W_ %*% X_ + 10e-10*diag(rep(1.0,dim(X_)[2])))
      XTWR = t(X_) %*% W_ %*% R_
      beta = inv_XTWX %*% XTWR
      tau = cbind(1, X) %*% beta
      
      Q = Posterior(train.index, tau, X, R, W, Ex,IW,weighting)
      beta.imputed = Posterior_fit(train.index, X, R, W, Ex, Q, A, marginalize = marginalize, lambda = lambda)
      tau.imputed = cbind(1, X) %*% beta.imputed
      
    }
    
    Q = Posterior(train.index, tau.imputed, X, R, W, Ex, IW = rep(1,n),weighting=FALSE)
    beta.imputed = Posterior_fit(train.index, X, R, W, Ex, Q, A, marginalize = marginalize, lambda = lambda)
    tau.imputed = cbind(1, X) %*% beta.imputed
  }
  cat("Data used for the nuisance in ART (%):", 100*round(length(train.index)/n,6),"\n")
  
  return(result = list(mu0 = mu0.hat, mu1 = mu1.hat, tau = tau.imputed, train.index = train.index, beta = beta.imputed, lambda.seq = lambda.seq)) 
}











































test.stats = function(Y, W, stats = "AIPW", Ex = NULL, mu0.hat = NULL, mu1.hat = NULL, mu.hat = NULL, tau.hat = NULL) {
  
  
  
  n = length(Y)
  if (is.null(Ex)){
    n1 = max(1, sum(W))
    n0 = max(1, sum(1 - W))
  }else{
    n1 = n * Ex
    n0 = n * (1 - Ex)
  }
  

  if (stats == "denoise") {
    
    IF = W * (Y - mu.hat)/Ex - (1 - W) * (Y - mu.hat)/(1-Ex)
    value = mean(IF)/sd(IF)
    
  }else if(stats == "AIPW"){
    
    #IF = -(Y-mu.hat-(W-Ex)*tau.hat)**2
    
    IF = (W * (Y - mu1.hat)/Ex - (1 - W) * (Y - mu0.hat)/(1-Ex) + tau.hat) 
    
    value = mean(IF)#/sd(IF)
    
  }  
  
  # Return the computed statistic value
  return(value)
}





test.stats.group = function(Y, W, Group, stats = "AIPW", Ex = NULL, mu0.hat = NULL, mu1.hat = NULL, mu.hat = NULL, tau.hat = NULL, test.index = NULL) {
  # Get the unique levels of the grouping variable
  Group.level = sort(unique(Group))
  
  if (is.null(test.index)==FALSE){
    W = sapply(Ex, function(p) rbinom(1, size = 1, prob = p))
    #W[test.index] = W_[test.index]
  }
  
  # Compute the test statistic value for each group level
  values = sapply(Group.level, function(x) {
    index = which(Group == x)
    test.stats(Y = Y[index], W = W[index],
               stats = stats, Ex = Ex[index], 
               mu0.hat = mu0.hat[index], 
               mu1.hat = mu1.hat[index], 
               mu.hat = mu.hat[index], 
               tau.hat = tau.hat[index])
  })
  
  # Return the computed values for each group level
  return(values)
}


nuisance.mu = function(Y = NULL, X = NULL) {
  
  # Load required library
  library(glmnet)
  
  # Perform cross-validated ridge regression
  cv_model <- cv.glmnet(
    x = as.matrix(X), 
    y = Y, 
    alpha = 0,  # Ridge regression
    nfolds = 25
  )
  
  # Get out-of-sample predictions
  oos_pred <- predict(cv_model, as.matrix(X), s = "lambda.min")
  
  return(oos_pred)
}


nuisance.mu.xgboost = function(Y = NULL, X = NULL, nfold = 5, eta = 0.01, nrounds = 1000, early_stopping_rounds = 10) {
  
  # Load required library
  library(xgboost)
  library(caret)
  
  # Create folds for cross-validation
  folds <- createFolds(Y, k = nfold, list = TRUE, returnTrain = FALSE)
  
  # Initialize vector to store out-of-sample predictions
  oos_pred <- rep(NA, length(Y))
  
  # XGBoost parameters
  params <- list(
    objective = "reg:squarederror",
    eval_metric = "rmse",
    eta = eta
  )
  
  # Perform cross-validation manually
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(Y), test_idx)
    
    dtrain <- xgb.DMatrix(data = X[train_idx, ], label = Y[train_idx])
    dtest <- xgb.DMatrix(data = X[test_idx, ])
    
    model <- xgb.train(
      params = params,
      data = dtrain,
      nrounds = nrounds,
      watchlist = list(train = dtrain),
      early_stopping_rounds = early_stopping_rounds,
      verbose = 0
    )
    
    # Store out-of-sample predictions
    oos_pred[test_idx] <- predict(model, dtest)
  }
  
  return(oos_pred)
}


inv_update = function(inv_x, w, x, sign){
  
  # Sherman–Morrison–Woodbury update of inversion

  u <- sqrt(w) * c(1,x)
  inv_u = inv_x %*% u
  inv_update = inv_x - sign * (inv_u %*% t(inv_u)) / (1 + sign * t(u) %*% inv_u)[1]
  return(inv_update)
}

xr_update = function(xtwr, w, r, x, sign){
  
  # Projection update
  
  return(xtwr + sign* w * r *  c(1,x))
}




Posterior = function(train.index, tau, X, R, W, Ex, IW,weighting){
  
  # generate leave-one-out residual
  n = length(Ex) 
  res <- numeric(n) #rep(0,length(train.index))
  
  treated.index = which(W > 0.5) # treated units
  control.index = which(W < 0.5) # control units
  
  #for (i in 1:length(train.index)){
  #  res[i] = R[train.index[i]] - (W[train.index[i]]-Ex[train.index[i]])*tau[train.index[i]] 
  #}
  
  if (weighting){
    size = n
  }else{
    size = length((train.index))
  }
  
  res.y = R[train.index] -  (W[train.index]-Ex[train.index])*tau[train.index]
  res.mean = 0 #sum(res.y*IW[train.index])/n
  res.std =  sqrt(sum((res.y-res.mean)**2*IW[train.index])/size)
  

  z.1 = pmax(pmin((R - (1-Ex)*tau - res.mean)/res.std,4),-4)
  z.0 =  pmax(pmin((R - (0-Ex)*tau - res.mean)/res.std,4),-4)

  p.zxy.1 = Ex*dnorm(z.1)
  p.zxy.0 = (1-Ex)*dnorm(z.0)
  
  p.zxy.1 = (p.zxy.1)/(p.zxy.1 + p.zxy.0)
  
  return(p.zxy.1)
}  






knn.indices <- function(X,Y, k) {
  # X is an n x d matrix
  # k is the total number of "neighbors" to pick, including i itself
  # Returns an n x n matrix where row i has 1's for
  # the k closest points to i (including i)
  
  n <- nrow(X)
  # Compute pairwise distances (as.matrix to ensure it’s a dense matrix)
  D <- as.matrix(dist(X)) + as.matrix(dist(Y))
  
  # Initialize the result as an n x n matrix of 0s
  A <- matrix(0, nrow = n, ncol = n)
  
  for (i in seq_len(n)) {
    # Get all indices sorted by distance from i (including i)
    sorted_indices <- order(D[i, ])
    
    # Take the top k, which now includes 'i' itself
    knn <- sorted_indices[1:k]
    
    # Mark these k positions in row i as 1
    A[i, knn] <- 1
  }
  
  return(A)
}


prop.group <- function(Group, Group.idx,test.index){
  
  probs = NULL
  for (g in Group.idx){
    probs = c(probs, sum(Group[test.index] == g)/sum(Group == g))
  }
  return(probs)
}


select.train <- function(Group, prop, random = T, ...) {

  unique_groups <- unique(Group)
  
  # Initialize an empty vector for the indices
  training.indices <- NULL
  
  for (g in unique_groups) {
    # Indices for group g
    idx.g <- which(Group == g)
    
    # How many to pick in group g
    n.pick <- round(length(idx.g) * prop)
    
    if(random){
      # Randomly sample those indices
      chosen <- sample(idx.g, n.pick)
    }else{
      # Select in decreasing order of X_j^\top (X_{[n]}^{T}X_{[n]})^{-2}X_j
      chosen <- intersect(order(diag(X %*% solve(t(X) %*% X) %*% solve(t(X) %*% X) %*% t(X)), decreasing = T), idx.g)[seq(1, n.pick)]
    }
    
    # Combine into the overall result
    training.indices <- c(training.indices, chosen)
  }
  
  return(training.indices)
}





closing <- function(p_val, q = 0.05, global.null.test = c("Fisher", "Bonferroni", "Simes")) {
  global.null.test <- match.arg(global.null.test)
  
  K <- length(p_val)
  subsets <- unlist(lapply(1:K, function(k) combn(K, k, simplify = FALSE)), recursive = FALSE)
  
  # Compute global p-values for each subset using the specified combining method
  global_p <- sapply(subsets, function(S) {
    S_pvals <- p_val[S]
    k <- length(S)
    
    if (global.null.test == "Fisher") {
      stat <- -2 * sum(log(S_pvals))
      pchisq(stat, df = 2 * k, lower.tail = FALSE)
      
    } else if (global.null.test == "Bonferroni") {
      min_p <- min(S_pvals)
      min(1, k * min_p)
      
    } else if (global.null.test == "Simes") {
      sorted_p <- sort(S_pvals)
      simes_p <- min(sorted_p * k / seq_len(k))
      min(1, simes_p)
    }
  })
  
  # For each hypothesis j, find max global p over all subsets that include j
  adjusted_p <- rep(NA, K)
  for (j in 1:K) {
    relevant_indices <- which(sapply(subsets, function(S) j %in% S))
    adjusted_p[j] <- max(global_p[relevant_indices])
  }
  
  rejected <- adjusted_p <= q
  return(rejected)
}




