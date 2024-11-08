# helper functions for panning
library("xgboost")
library(randomForest)
library(caret)
library(ranger)
library(pracma)
library(Matrix)


# BH
BH.threshold = function(pval, q = 0.1){
  
  # BH procedure
  
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
  
  # BH procedure + Storey correction
  
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
  
  # Compute a permutation p-value given the observed stats and reference distribution
  
  if(is.na(stats)){return(1)}
  return((sum(stats <= stats.ref[!is.na(stats.ref)]) + 1) / (length(stats.ref[!is.na(stats.ref)]) + 1))
}



nuisance.mu = function(Y= NULL, X = NULL){
  
  # Estimate the expectation of Y given X using xgboost
  
  params = list(objective = "reg:squarederror", eval_metric = "rmse",eta=0.1)
  nfold = 10
  nrounds = 1000
  early_stopping_rounds = 10
  
  nuisance.mu <- xgb.cv(
    data = xgb.DMatrix(data = X, label = Y),
    params = params, 
    nfold = nfold,
    nrounds = nrounds,
    early_stopping_rounds = early_stopping_rounds,
    verbose = FALSE,
    callbacks = list(cb.cv.predict(save_models = TRUE))
  )
  
  nuisance.mu$pred
  
}


nuisance.tau.ss = function(Y= NULL, X = NULL, Ex = NULL, W = NULL, mu = NULL, train.index = NULL, ...){
  
  # Estimate CATE using sample splitting
  
  X_copied <- cbind(1, X[train.index,])
  inv_XTWX = solve(t(X_copied) %*% diag((W-Ex)[train.index]**2) %*% X_copied+ 10e-16*diag(rep(1.0,dim(X_copied)[2])))
  XTWR = t(X_copied) %*% diag((W-Ex)[train.index]**2) %*% ((Y[train.index] - mu[train.index]) / (W-Ex)[train.index])
  tau.hat <- cbind(1, X) %*% inv_XTWX %*% XTWR
  
  mu0.hat <- mu - Ex * tau.hat
  mu1.hat <- mu0.hat + tau.hat
  
  
  return(result = list(mu0 = mu0.hat, mu1 = mu1.hat, tau = tau.hat))
}




inv_update = function(inv_x, w, x, sign){
  
  # Sherman–Morrison–Woodbury update of inversion
  
  u <- sqrt(w) * c(1,x); inv_u = inv_x %*% u
  inv_update = inv_x - sign * (inv_u %*% t(inv_u)) / (1 + sign * t(u) %*% inv_u)[1]
  return(inv_update)
}

xr_update = function(xtwr, w, r, x, sign){
  
  # Projection update
  
  return(xtwr + sign* w * r *  c(1,x))
}





nuisance.tau.active = function(Y= NULL, X = NULL, Ex = NULL, W = NULL, mu = NULL, max_proportion = NULL,...){
  
  # Estimate CATE using active learning (using less data if possible)
  
  n = length(Y)
  
  R = Y - mu # residual of mu
  R1 = R / (1 - Ex) # residual of tau given z = 1
  R2 = R / (0 - Ex) # residual of tau given z = 0
  R3 = R/ (W-Ex) # true residual of tau given z
  
  W1 = (1 - Ex)**2 # weight when z = 1
  W2 = Ex**2 # weight when z = 0
  W3 = (W-Ex)**2 # true weight given z
  
  treated.index = which(W > 0.5) # treated units
  control.index = which(W < 0.5) # control units
  
  n_0 = dim(X)[2] + 1 # initial training set
  train.index = c(1:n_0) # training units
  test.index = setdiff(1:n,train.index) #testing units
  
  # initial fit
  X_ = cbind(1, X[train.index,])
  R_ = R3[train.index]
  W_ = diag(c(W3[train.index]))
  inv_XTWX = solve(t(X_) %*% W_ %*% X_ + 10e-16*diag(rep(1.0,dim(X_)[2])))
  XTWR = t(X_) %*% W_ %*% R_
  tau.hat <- cbind(1, X) %*% inv_XTWX %*% XTWR
  
  
  
  # generate leave-one-out residual
  tau.hat.val <- rep(0,n)
  res.val <-  rep(0,n)
  
  inv_XTWX_list = list()
  XTWR_list = list()
  
  for (i in train.index){
    
    # remove the i-th training data point
    inv_XTWX_list[[i]] = inv_update(inv_XTWX, W3[i], X[i,], -1) 
    XTWR_list[[i]] = xr_update(XTWR, W3[i], R3[i], X[i,], -1)
    
    tau.hat.val[i] = c(1, X[i,]) %*% inv_XTWX_list[[i]] %*% XTWR_list[[i]] # compute the leave-one-out prediction
    res.val[i] = R[i] - (W[i]-Ex[i])*tau.hat.val[i] # compute the leave-one-out residual
     
  }
  
  # estimate the prob. of z given x and y
  train.treated.index = intersect(train.index,treated.index) 
  mean.treated = mean(res.val[train.treated.index]) # treated mean
  std.treated = std(res.val[train.treated.index]) # treated std
  
  train.control.index = intersect(train.index,control.index)
  mean.control = mean(res.val[train.control.index]) # control mean
  std.control = std(res.val[train.control.index]) # treated std
  
  p.zxy.1 = Ex*dnorm((R - (1-Ex)*tau.hat - mean.treated)/std.treated)
  p.zxy.0 = (1-Ex)*dnorm((R - (0-Ex)*tau.hat - mean.control)/std.control)
  p.zxy = p.zxy.1 + p.zxy.0+0.000002
  p.zxy.1 = (p.zxy.1+0.000001)/p.zxy # prob. of z = 1 given x and y
  p.zxy.0 = (p.zxy.0+0.000001)/p.zxy # prob. of z = 0 given x and y
  

  
  active.sampling = function(){
    
    # # sample the point that minimize the entropy of p(z|x,y) in the remaining points
  
    entropy = c()
    
    num.test = length(test.index)
    test.sample_0 = test.index #sample(test.index,min(250,num.test))

    for (i in test.index){
      test.sample = setdiff(test.sample_0,i) # remove the i-th test point
      x = X[i,]
      w1 = W1[i]; w2 = W2[i]
      r1 = R1[i]; r2 = R2[i]
      
      # add the point x_i, y_i, z_i = 1
      inv_XTWX_1 = inv_update(inv_XTWX, w1, x, 1)
      XTWR_1 = xr_update(XTWR, w1, r1, x, 1)
      tau.hat.1 <- cbind(1, X[test.sample,]) %*% inv_XTWX_1 %*% XTWR_1
      
      # update the prob. of z given x and y
      p.zxy.i.1 = Ex[test.sample]*dnorm((R[test.sample] - (1-Ex[test.sample])*tau.hat.1 - mean.treated)/std.treated)
      p.zxy.i.0 = (1-Ex[test.sample])*dnorm((R[test.sample] - (0-Ex[test.sample])*tau.hat.1 - mean.control)/std.control)
      p.zxy.i = p.zxy.i.1 + p.zxy.i.0+0.000002; p.zxy.i.1 = (p.zxy.i.1+0.000001)/p.zxy.i; p.zxy.i.0 = (p.zxy.i.0+0.000001)/p.zxy.i
      
      # compute the entropy
      entropy.treated = - p.zxy.i.1*log(p.zxy.i.1) - p.zxy.i.0*log(p.zxy.i.0)
      
      # add x_i, y_i, z_i = 0
      inv_XTWX_2 = inv_update(inv_XTWX, w2, x, 1)
      XTWR_2 = xr_update(XTWR, w2, r2, x, 1)
      tau.hat.2 <- cbind(1, X[test.sample,]) %*% inv_XTWX_2 %*% XTWR_2
      
      # update the prob. of z given x and y
      p.zxy.i.1 = Ex[test.sample]*dnorm((R[test.sample] - (1-Ex[test.sample])*tau.hat.2 - mean.treated)/std.treated)
      p.zxy.i.0 = (1-Ex[test.sample])*dnorm((R[test.sample] - (0-Ex[test.sample])*tau.hat.2 - mean.control)/std.control)
      p.zxy.i = p.zxy.i.1 + p.zxy.i.0+0.000002; p.zxy.i.1 = (p.zxy.i.1+0.000001)/p.zxy.i; p.zxy.i.0 = (p.zxy.i.0+0.000001)/p.zxy.i
      
      # compute the entropy
      entropy.control = - p.zxy.i.1*log(p.zxy.i.1) - p.zxy.i.0*log(p.zxy.i.0)

      # sum over the predictive entropies      
      entropy_i = p.zxy.1[i]*mean(entropy.treated)+ p.zxy.0[i]*mean(entropy.control)
      
      entropy = c(entropy, entropy_i)
      
    }
    
    
    test.index[which.min(entropy)] # return the point that leaves minimum uncertainty
    
  }
  
  

  max_epochs <- round(max_proportion*n)-n_0  # maximum number of training epochs
  loss <- 1      # initial loss
  window_size <- 30      # window size for calculating the parameter change
  loss_threshold <- 0.005      # stopping threshold 
  loss_history <- rep(Inf, window_size)
  epoch <- 1
  
  
  while (epoch <= max_epochs) {
    
    
    idx_min = active.sampling() #actively sample a new data point
    x = X[idx_min,] #feature
    w3 = W3[idx_min] #weight
    r3 = R3[idx_min] #residual
    
    # update the leave-one-out fit with the new data point
    for (i in train.index){
      inv_XTWX_list[[i]] = inv_update(inv_XTWX_list[[i]], w3, x, 1)
      XTWR_list[[i]] = xr_update(XTWR_list[[i]], w3, r3, x, 1)
      tau.hat.val[i] = c(1, X[i,]) %*% inv_XTWX_list[[i]] %*% XTWR_list[[i]]
      res.val[i] = R[i] - (W[i]-Ex[i])*tau.hat.val[i]
    }
    
    # save the fit leaving out the new data point 
    inv_XTWX_list[[idx_min]] = inv_XTWX
    XTWR_list[[idx_min]] = XTWR
    beta_min = inv_XTWX_list[[idx_min]] %*% XTWR_list[[idx_min]]
    tau.hat.val[idx_min] = c(1, x) %*% beta_min
    res.val[idx_min] = R[idx_min] - (W[idx_min]-Ex[idx_min])*tau.hat.val[idx_min]
    
    # create the fit using the old and new data point
    inv_XTWX = inv_update(inv_XTWX, w3, x, 1)
    XTWR = xr_update(XTWR, w3, r3, x, 1)
    beta = inv_XTWX %*% XTWR
    tau.hat <- cbind(1, X) %*% beta
    
    # update the train and test index
    train.index = c(train.index,idx_min)
    test.index = setdiff(1:n,train.index)
    

    #update the distribution of z given x and y
    train.treated.index = intersect(train.index,treated.index)
    mean.treated = mean(res.val[train.treated.index])
    std.treated = std(res.val[train.treated.index])
  
    train.control.index = intersect(train.index,control.index)
    mean.control = mean(res.val[train.control.index])
    std.control = std(res.val[train.control.index])

    p.zxy.1 = Ex*dnorm((R - (1-Ex)*tau.hat - mean.treated)/std.treated)
    p.zxy.0 = (1-Ex)*dnorm((R - (0-Ex)*tau.hat - mean.control)/std.control)
    p.zxy = p.zxy.1 + p.zxy.0+0.000002; p.zxy.1 = (p.zxy.1+0.000001)/p.zxy; p.zxy.0 = (p.zxy.0+0.000001)/p.zxy
    
    
    # compute the parameter change due to the new data point
    loss = sum((beta - beta_min)**2)/sum((beta_min)**2)
    loss_history <- c(loss_history[-1], loss)
  
    if (epoch >= window_size) {

      ave_loss <- mean(loss_history)
      # Stop when the added points in the window do not lead to a large change on beta coefficients, and so tau and p.zxy.
      if (ave_loss < loss_threshold) {
        cat("Data used for the nuisance in ART (%):", 100*round(length(train.index)/n,6),"\n")
        cat("Change of parameter stopped in ART (%):", 100*round(ave_loss,6),"\n")
        break
      }
    }
    
    epoch <- epoch + 1
  }
  
  # we have sampled the most influential (i.e. representative) points in the sense that they reduce the uncertainty of p(z|x,y) of the others
  # compute the final model with the unknown treatment variables marginalized out (i.e. imputed).
  X_copied <- cbind(1, rbind(X,X))
  W1[train.index] = W[train.index] * W1[train.index]
  W1[test.index] = p.zxy.1[test.index] * W1[test.index]
  W2[train.index] = (1-W[train.index]) * W2[train.index]
  W2[test.index] = p.zxy.0[test.index] * W2[test.index]
  W_copied <- diag(c(W1,W2))
  R_copied <- c(R1,R2)

  inv_XTWX = solve(t(X_copied) %*% W_copied %*% X_copied + 10e-16*diag(rep(1.0,dim(X_copied)[2])))
  XTWR = t(X_copied) %*% W_copied %*% R_copied
  tau.hat <- cbind(1, X) %*% inv_XTWX %*% XTWR
  
  mu0.hat <- mu - Ex * tau.hat
  mu1.hat <- mu0.hat + tau.hat
  
  return(result = list(mu0 = mu0.hat, mu1 = mu1.hat, tau = tau.hat, train.index = train.index)) 
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
    
    IF = W * (Y - mu1.hat)/Ex - (1 - W) * (Y - mu0.hat)/(1-Ex) + tau.hat
    value = mean(IF)/sd(IF)
    
  }  
  
  # Return the computed statistic value
  return(value)
}



test.stats.group = function(Y, W, Group, stats = "AIPW", Ex = NULL, mu0.hat = NULL, mu1.hat = NULL, mu.hat = NULL, tau.hat = NULL, test.index = NULL) {
  # Get the unique levels of the grouping variable
  Group.level = sort(unique(Group))
  

  if (is.null(test.index)==FALSE){
    W[test.index] =  sapply(Ex[test.index], function(p) rbinom(1, size = 1, prob = p))
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

