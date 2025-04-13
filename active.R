
#source("~/Desktop/Research/Yao/HTE inference/code/Panning/helper Yao.R")

source("/Users/yaozhang/Documents/GitHub/Panning/helper Yao.R")

set.seed(318)


setting = 2
n = 400 #sample size
Group.number = 5 # total number of groups
delta = 0.5 # effect size
sigma = 0.2 # std of noise
M = 1000 # number of permutations
n_trial = 20 # number of runs
proportion = 0.5 # proportion of randomness in the nuisance fold
test.stats.method ="AIPW" #Test statistics

# Data generation
num_features = 5

# Store the experimental results
simulation_results_art <- matrix(0, nrow = n_trial, ncol = Group.number)
colnames(simulation_results_art) <- paste0("p", seq(1, Group.number))

simulation_results_ssrt <- matrix(0, nrow = n_trial, ncol = Group.number)
colnames(simulation_results_ssrt) <- paste0("p", seq(1, Group.number))


aa = c()
bb = c()

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
  
  mu0 <- mu - Ex * tau
  mu1 <- mu0 + tau
  
  W <- sapply(Ex, function(p) rbinom(1, size = 1, prob = p))
  Y = mu0 + W*tau + rnorm(n, 0, sigma)


  mu.hat = nuisance.mu(Y,X)
  tau.hat.active = nuisance.tau.active(Y = Y, X = X, Ex = Ex, W = W, mu = mu.hat, Group = Group, proportion = proportion) #estimate tau using active learning
  train.index.active = tau.hat.active$train.index #the data points fitted to the nuisance
  test.index.active = setdiff(1:n,train.index.active) #the data points used for inference

  v = sum((tau[test.index.active]-tau.hat.active$tau[test.index.active])**2)/sum((tau[test.index.active]-mean(tau[test.index.active]))**2)
  cat("R^2 of ART:", 1- v, "\n")
  aa = c(aa, 1-v)

  # Randomization tests
  test.stats.value.active = test.stats.group(Y = Y[test.index.active], W = W[test.index.active], Group = Group[test.index.active], stats = test.stats.method, Ex = Ex[test.index.active], 
                                      mu0.hat = tau.hat.active$mu0[test.index.active], mu1.hat = tau.hat.active$mu1[test.index.active], mu.hat = mu.hat[test.index.active], 
                                      tau.hat = tau.hat.active$tau[test.index.active])

  test.stats.ref.active = t(replicate(M, test.stats.group(Y = Y[test.index.active], W =  W[test.index.active], Group = Group[test.index.active], stats = test.stats.method, Ex = Ex[test.index.active], 
                                                   mu0.hat = tau.hat.active$mu0[test.index.active], mu1.hat = tau.hat.active$mu1[test.index.active], mu.hat = mu.hat[test.index.active], 
                                                   tau.hat = tau.hat.active$tau[test.index.active], test.index=test.index.active[test.index.active])))
  

  
  # Calculate p-values for each group
  simulation_results_art[j,] = sapply(seq(1, length(test.stats.value.active)), function(x) { sum(test.stats.ref.active[, x]>=test.stats.value.active[x])/(M+1) + 1/(M+1)})
  #simulation_results_art[j,]


  # Sample splitting
  train.index.ss = select.train(Group, proportion)  #sample(n, proportion*n) #sample(n, length(train.index.active))  
    #sample(n, proportion*n)  #select.train(Group, proportion) 
  test.index.ss = setdiff(1:n, train.index.ss)
  
  # Fit the nuisance parameter
  tau.hat.ss = nuisance.tau.ss(Y = Y, X = X, Ex = Ex, W = W, mu = mu.hat, train.index = train.index.ss)

  v = sum((tau[test.index.ss]-tau.hat.ss$tau[test.index.ss])**2)/sum((tau[test.index.ss]-mean(tau[test.index.ss]))**2)
  cat("R^2 of SSRT:", 1 - v, "\n")
  bb = c(bb, 1- v)
  # Randomization tests
  test.stats.value = test.stats.group(Y = Y[test.index.ss], W = W[test.index.ss], Group = Group[test.index.ss], stats = test.stats.method, Ex = Ex[test.index.ss], 
                                      mu0.hat = tau.hat.ss$mu0[test.index.ss], mu1.hat = tau.hat.ss$mu1[test.index.ss], mu.hat = mu.hat[test.index.ss], tau.hat = tau.hat.ss$tau[test.index.ss])
  
  test.stats.ref = t(replicate(M, test.stats.group(Y = Y[test.index.ss], W = W[test.index.ss], Group = Group[test.index.ss], stats = test.stats.method, Ex = Ex[test.index.ss], 
                                                   mu0.hat = tau.hat.ss$mu0[test.index.ss], mu1.hat = tau.hat.ss$mu1[test.index.ss], mu.hat = mu.hat[test.index.ss], 
                                                   tau.hat = tau.hat.ss$tau[test.index.ss], test.index = test.index.ss))) 
  
  
  # Calculate p-values for each group
  simulation_results_ssrt[j,] = sapply(seq(1, length(test.stats.value)), function(x) { sum(test.stats.ref[, x]>=test.stats.value[x])/(M+1) + 1/(M+1)})
  
  

  cat(" ART:",round(simulation_results_art[j,] ,3), "\n")
  cat("SSRT:",round(simulation_results_ssrt[j,],3), "\n")
  cat("\n")
  

}


cat(" ART:", format(round(colMeans(simulation_results_art[1:j,]),2),2), "\n")
cat("SSRT:", format(round(colMeans(simulation_results_ssrt[1:j,]),2),2), "\n")
cat("\n")




















# Generate the Boxplots
pdf(paste0("Comparison_setting_", setting, ".pdf"), width = 13, height = 6)

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
        main = paste0("Comparison of ART and SSRT (Maximum Inference Fold: ", (1-proportion)*100, "%)"),
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











