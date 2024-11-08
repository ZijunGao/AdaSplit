source("~/Desktop/Research/Yao/HTE inference/code/Panning/helper Yao.R")

set.seed(318)

n = 500 #sample size
Group.number = 8 # total number of groups
delta = 1 # effect size
sigma = 3 # std of noise
M = 2000 # number of permutations
n_trial = 30 # number of runs
proportion = 0.5 # proportion of randomness in the nuisance fold
test.stats.method ="AIPW" #Test statistics

# Data generation
num_features = 20
X = matrix(runif(n*num_features, 0, 6), nrow = n, ncol = num_features) 
quantiles_1 <- quantile(c(-Inf,X[,1],Inf), probs = seq(0, 1, length.out = Group.number[1] + 1))
S = as.numeric(cut(X[,1], quantiles_1))

Group = S - 1
G = model.matrix(~ factor(Group))[, -1]
Ex = exp((X[,1]-3))/(1+exp((X[,1]-3)))
W <- sapply(Ex, function(p) rbinom(1, size = 1, prob = p))
mu0 =  1 + X %*% c(1:num_features)*0.1
Y0 = mu0 + rnorm(n, 0, sigma)
tau = 2 + delta * (X[,1])
Y1 = Y0 + tau 
Y = Y1 * W + Y0 * (1 - W) 
mu1 = mu0 + tau
mu = mu0 * (1 - Ex) + mu1 * Ex


# Store the experimental results
simulation_results_art <- matrix(0, nrow = n_trial, ncol = Group.number)
colnames(simulation_results_art) <- paste0("p", seq(1, Group.number))

simulation_results_ssrt <- matrix(0, nrow = n_trial, ncol = Group.number)
colnames(simulation_results_ssrt) <- paste0("p", seq(1, Group.number))




for (j in 1:n_trial){
  
  mu.hat = nuisance.mu(Y,X)

  #estimate mu using xgboost
  tau.hat.active = nuisance.tau.active(Y = Y, X = X, Ex = Ex, W = W, mu = mu.hat, max_proportion = 0.8) #estimate tau using active learning
  train.index.active = tau.hat.active$train.index #the data points fitted to the nuisance
  test.index.active = setdiff(1:n,train.index.active) #the data points used for inference
  
  #plot(X[test.index.active,1], tau.hat.active$tau[test.index.active], col = 'blue',pch = 17, main='ART-learner',ylim = c(-20,30))
  #points(X[test.index.active,1],tau[test.index.active], col = 'red', pch = 17)

  #print(mean(abs( tau.hat.active$tau[test.index.active]- tau[test.index.active])))
  
  # Randomization tests
  test.stats.value = test.stats.group(Y = Y, W = W, Group = Group, stats = test.stats.method, Ex = Ex, 
                                      mu0.hat = tau.hat.active$mu0, mu1.hat = tau.hat.active$mu1, mu.hat = mu.hat, tau.hat = tau.hat.active$tau)
  test.stats.value
  test.stats.ref = t(replicate(M, test.stats.group(Y = Y, W =  W, Group = Group, stats = test.stats.method, Ex = Ex, 
                                                   mu0.hat = tau.hat.active$mu0, mu1.hat = tau.hat.active$mu1, mu.hat = mu.hat, 
                                                   tau.hat = tau.hat.active$tau, test.index=test.index.active))) 
  
  # Calculate p-values for each group
  simulation_results_art[j,] = sapply(seq(1, length(test.stats.value)), function(x) {
    permutation.p.value(stats = test.stats.value[x], stats.ref = test.stats.ref[, x])
  })
  simulation_results_art[j,]
  
  # Sample splitting
  train.index = sample(n, n * proportion) 
  test.index = setdiff(1:n, train.index)
  
  # Fit the nuisance parameter
  tau.hat.ss = nuisance.tau.ss(Y = Y, X = X, Ex = Ex, W = W, mu = mu.hat, train.index = train.index)
  
  
  #plot(X[test.index,1],tau.hat.ss$tau[test.index], col = 'blue',pch = 17,,main='SSRT-learner',ylim = c(-20,30))
  #points(X[test.index,1],tau[test.index], col = 'red', pch = 17)
  
  #print(mean(abs(tau.hat.ss$tau[test.index]- tau[test.index])))
  
  
  # Randomization tests
  test.stats.value = test.stats.group(Y = Y, W = W, Group = Group, stats = test.stats.method, Ex = Ex, 
                                      mu0.hat = tau.hat.ss$mu0, mu1.hat = tau.hat.ss$mu1, mu.hat = mu.hat, tau.hat = tau.hat.ss$tau)
  
  test.stats.ref = t(replicate(M, test.stats.group(Y = Y, W =  W, Group = Group, stats = test.stats.method, Ex = Ex, 
                                                   mu0.hat = tau.hat.ss$mu0, mu1.hat = tau.hat.ss$mu1, mu.hat = mu.hat, 
                                                   tau.hat = tau.hat.ss$tau, test.index = test.index))) 
  
  # Calculate p-values for each group
  simulation_results_ssrt[j,] = sapply(seq(1, length(test.stats.value)), function(x) {
    permutation.p.value(stats = test.stats.value[x], stats.ref = test.stats.ref[, x])
  })
  
  round(simulation_results_ssrt[j,],3)
  round(simulation_results_art[j,] ,3)
  
  if (j>=2){
    cat(" ART:", format(round(colMeans(simulation_results_art[1:j,]),2),2), "\n")
    cat("SSRT:", format(round(colMeans(simulation_results_ssrt[1:j,]),2),2), "\n")
    cat("\n")
  }


}






# Generate the Boxplots
pdf(paste0("Comparison_plot_", proportion, ".pdf"), width = 13, height = 6)

p_values_df_art <- as.data.frame(simulation_results_art)
p_values_df_ssrt <- as.data.frame(simulation_results_ssrt)


colnames(p_values_df_art) <- paste("Group", 1:Group.number)
colnames(p_values_df_ssrt) <- paste("Group", 1:Group.number)


p_values_df_art$Method <- "ART"
p_values_df_ssrt$Method <- "SSRT"


combined_p_values_df <- rbind(p_values_df_art, p_values_df_ssrt)


combined_p_values_long <- reshape2::melt(combined_p_values_df, id.vars = "Method")


par(mar = c(5, 5, 4, 5))


group_spacing <- 2.5 
at_positions <- rep(1:Group.number, each = 2) * group_spacing + c(-0.5, 0.5)


boxplot(value ~ Method + variable, data = combined_p_values_long, 
        at = at_positions, col = c("lightblue", "grey"), notch = TRUE, xaxt = "n",
        ylab = "P-values", xlab = "", # Remove the x-axis label
        main = paste0("Comparison of ART and SSRT (Inference fold: ", proportion*100, "%)"),
        cex.axis = 1.5, 
        cex.lab = 1.5,   
        cex.main = 1.6)   


axis(1, at = 1:Group.number * group_spacing, labels = paste("Group", 1:Group.number),
     cex.axis = 1.5)  


legend("topright", legend = c("ART", "SSRT"), fill = c("lightblue", "grey"), 
       cex = 1.7, 
       pt.cex = 2,  
       bg = "white") 


dev.off()


