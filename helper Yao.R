# helper functions for panning
library("xgboost")
library(randomForest)
library(caret)
library(ranger)
library(pracma)
library(Matrix)


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




nuisance.tau.ss = function(Y= NULL, X = NULL, Ex = NULL, W = NULL, mu = NULL, train.index = NULL, ...){
  
  # initial fit
  n = length(Y)
  R = Y - mu
  PW = rep(1,n)
  X_ = cbind(1, X[train.index,])
  R_ = (R/ (W-Ex))[train.index]
  W_ = diag(c(((W-Ex)**2)[train.index]))
  inv_XTWX = solve(t(X_) %*% W_ %*% X_ + 10e-4*diag(rep(1.0,dim(X_)[2])))
  XTWR = t(X_) %*% W_ %*% R_
  beta = inv_XTWX %*% XTWR
  tau = cbind(1, X) %*% beta
  
  #Q = Posterior(train.index, tau, X, R, W, Ex)
  #beta.imputed = Posterior_fit(train.index, inv_XTWX, XTWR, X, R, W, Ex, Q, PW)
  #tau.imputed = cbind(1, X) %*% beta.imputed
  
  mu0.hat <- mu - Ex * tau
  mu1.hat <- mu0.hat + tau
  
  return(result = list(mu0 = mu0.hat, mu1 = mu1.hat, tau = tau))
}


  
Posterior = function(train.index, tau, X, R, W, Ex){
  
  # generate leave-one-out residual
  n = length(Ex) 
  res <-  rep(0,n)

  treated.index = which(W > 0.5) # treated units
  control.index = which(W < 0.5) # control units
  
  for (i in train.index){
    res[i] = R[i] - (W[i]-Ex[i])*tau[i] 
  }
  
  train.treated.z = rep(0,n)
  train.treated.z[intersect(train.index,treated.index)] = 1 
  mean.treated = mean(train.treated.z[train.index]*res[train.index]/Ex[train.index]) 
  std.treated = mean(train.treated.z[train.index]*(res[train.index]-mean.treated)**2/Ex[train.index])

  train.control.z = rep(0,n)
  train.control.z[intersect(train.index,control.index)] = 1 
  mean.control = mean(train.control.z[train.index]*res[train.index]/(1-Ex[train.index])) 
  std.control = mean(train.control.z[train.index]*(res[train.index]-mean.control)**2/(1-Ex[train.index])) 
  
  p.zxy.1 = Ex*dnorm((R - (1-Ex)*tau - mean.treated)/std.treated)
  p.zxy.0 = (1-Ex)*dnorm((R - (0-Ex)*tau - mean.control)/std.control)
  p.zxy = p.zxy.1 + p.zxy.0
  p.zxy.1 = (p.zxy.1)/p.zxy # prob. of z = 1 given x and y
  p.zxy.0 = (p.zxy.0)/p.zxy # prob. of z = 0 given x and y
  
  return(list(p.zxy.0, p.zxy.1))
}  


Posterior_fit = function(train.index, inv_XTWX, XTWR, X, R, W, Ex, Q, PW){
  
  p.zxy.0 = Q[[1]]; p.zxy.1 = Q[[2]]
  test.index = setdiff(1:n,train.index)

  W1 = (1 - Ex)**2/PW
  W2 = Ex**2/PW
  
  XX <- cbind(1, rbind(X,X))
  W1_ = W1; W2_ = W2
  
  W1_[train.index] = W[train.index] * W1[train.index]
  W1_[test.index] = p.zxy.1[test.index] * W1[test.index]
  
  W2_[train.index] = (1-W[train.index]) * W2[train.index]
  W2_[test.index] = p.zxy.0[test.index] * W2[test.index]
  
  inv_XTWX_ = solve(t(XX) %*% diag(c(W1_,W2_)) %*% XX + 10e-4*diag(rep(1.0,dim(XX)[2])))
  beta = inv_XTWX_ %*% t(XX) %*% diag(c(W1_,W2_)) %*% c( R/(1 - Ex), R/(0 - Ex))
  return(beta)
}




nuisance.tau.active = function(Y= NULL, X = NULL, Ex = NULL, W = NULL, mu = NULL, Group=NULL, max_proportion = NULL,...){
  
  
  # Estimate CATE using active splitting (using less data if possible)
  n = length(Y) #sample size
  PW = rep(1,n)
  
  # Calculate diversity score
  Xt <- t(X)              
  inv_XtX <- solve(Xt %*% X)  
  s <- colSums(X)          
  Div <- (X %*% inv_XtX %*% s)**2
  
  # Calculate labels
  R = Y - mu # residual of mu
  
  prob = sqrt(Div)/sum(sqrt(Div))*round((dim(X)[2]*3 + 1))
  prob[prob>=1] = 1
  prob[prob< 0.05] =  0.05
  
  sample = rbinom(n = n, size = 1, prob)
  n_0 = sum(sample)
  train.index = which(sample>0.5) # training units
  test.index = setdiff(1:n,train.index) #testing units
  
  
  PW[train.index] = prob[train.index]
  PW[test.index] = 1-prob[test.index]

  X_ = cbind(1, X[train.index,])
  R_ = (R/ (W-Ex))[train.index]
  W_ = diag(c(((W-Ex)**2)[train.index])) /PW[train.index]
  inv_XTWX = solve(t(X_) %*% W_ %*% X_ + 10e-4*diag(rep(1.0,dim(X_)[2])))
  XTWR = t(X_) %*% W_ %*% R_
  beta = inv_XTWX %*% XTWR
  tau = cbind(1, X) %*% beta
  
  Q = Posterior(train.index, tau, X, R, W, Ex)
  #beta.imputed = Posterior_fit(train.index, inv_XTWX, XTWR, X, R, W, Ex, Q, PW)
  R.imputed = Q[[1]]*(R / (0 - Ex)) + Q[[2]]*(R / (1 - Ex))
  Rsqr.imputed = Q[[1]]*(R / (0 - Ex)-R.imputed)**2 + Q[[2]]*(R / (1 - Ex)-R.imputed)**2
  
  active.sampling = function(){
    
    DR.test = sqrt(Div[test.index]*Rsqr.imputed[test.index])
    prob = DR.test/sum(DR.test)
    prob[prob< 0.05] = 0.05
    
    n_sample = 0
    t = 0
    while (n_sample<1){
      sample = rbinom(n = length(DR.test), size = 1, prob = prob)
      n_sample = sum(sample)
      t = t+1
    }
    

    #prob[prob>0.975] = 0.975
    prob[sample>0.5] = prob[sample>0.5]*(1-prob[sample>0.5])**(t-1)
    prob[sample<0.5] = (1-prob[sample<0.5])**(t)
    selected = which(sample > 0.5) 
    return(list(selected,prob))
    
  }
  
  max_epochs <- round(max_proportion*n)-n_0  # maximum number of training epochs
  loss <- 1      # initial loss
  window_size <- 20      # window size for calculating the parameter change
  loss_threshold <- 0.001    # stopping threshold 
  loss_history <- rep(Inf, window_size)
  epoch <- 0
  
  weight = (W-Ex)**2
  residual = R/ (W-Ex)

  stop = FALSE
  while (length(train.index) <= 0.5*n) {

    out = active.sampling() 
    idx.active = out[[1]]
    PW[test.index] = PW[test.index] * out[[2]]
    
    for (i in idx.active){
      new = test.index[i]
      epoch <- epoch + 1
      x = X[new,]
      w = weight[new] /(PW[new])
      r = residual[new]
      inv_XTWX = inv_update(inv_XTWX, w, x, 1)
      XTWR = xr_update(XTWR, w, r, x, 1)

      beta.copy = beta
      beta = inv_XTWX %*% XTWR
      tau = cbind(1, X) %*% beta
      train.index = c(train.index, new)
      
      Q = Posterior(train.index, tau, X, R, W, Ex)
      R.imputed = Q[[1]]*(R / (0 - Ex)) + Q[[2]]*(R / (1 - Ex))
      Rsqr.imputed = Q[[1]]*(R / (0 - Ex)-R.imputed)**2 + Q[[2]]*(R / (1 - Ex)-R.imputed)**2
      
      # compute the parameter change due to the new data point
      loss = sum((beta - beta.copy)**2)/sum((beta.copy)**2)
      loss_history <- c(loss_history[-1], loss)
      
      
      if (epoch >= window_size) {
        ave_loss <- mean(loss_history)
        if (ave_loss < loss_threshold) { 
          stop =TRUE
          cat("Data used for the nuisance in ART (%):", 100*round(length(train.index)/n,6),"\n")
          cat("Change of parameter stopped in ART (%):", 100*round(ave_loss,6),"\n")
          break
        }
      }
    }
    test.index = setdiff(1:n,train.index)
    if (stop){break} 
    
  }
  
  # we have sampled the most influential (i.e. representative) points in the sense that they reduce the uncertainty of p(z|x,y) of the others
  # compute the final model with the unknown treatment variables marginalized out (i.e. imputed).
  beta.imputed = Posterior_fit(train.index, inv_XTWX, XTWR, X, R, W, Ex, Q, PW)
  tau.imputed = cbind(1, X) %*% beta.imputed
  mu0.hat <- mu - Ex * tau.imputed
  mu1.hat <- mu0.hat + tau.imputed
  
  return(result = list(mu0 = mu0.hat, mu1 = mu1.hat, tau = tau.imputed, train.index = train.index)) 
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
    W =  sapply(Ex, function(p) rbinom(1, size = 1, prob = p))
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

