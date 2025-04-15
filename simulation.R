# simulation
source("~/Desktop/Research/Yao/HTE inference/code/Panning/helper Yao.R")

# source("/Users/yaozhang/Documents/GitHub/Panning/helper Yao.R")

set.seed(318)

setting = 2
n = 500 # sample size: 300, 400, 500, 600, 700
num_features = 5
Group.number = 5 # total number of groups
delta = 0.5 # effect size; 0.5
sigma = 0.2 # std of noise: 0.1, 0.2, 0.4
M = 1000 # number of permutations
proportion = 0.5 # proportion of randomness in the nuisance fold
test.stats.method ="AIPW" # test statistics
n_trial = 20 # number of runs; 20
verbose = F # if true, print out intermediate results
q = 0.2 # FWER level
marginalize = T # if marginalize

# Store the experimental results
# Subgroup p-values
simulation_results_art <- matrix(0, nrow = n_trial, ncol = Group.number)
colnames(simulation_results_art) <- paste0("p", seq(1, Group.number))
simulation_results_ssrt = simulation_results_rt = simulation_results_art

# FWER
FWER = list(); FWER$ART = FWER$SSRT = FWER$RT = c()

# Number of subgroup rejections
R = list(); R$ART = R$SSRT = R$RT = c()

# Quality of nuisance estimator
R2 = list(); R2$ART = R2$SSRT = c()

# Data for evaluating the nuisance estimation performance
n.val = 10^4
X.val = matrix(runif(n.val*num_features, 0, 1), nrow = n.val, ncol = num_features) 
quantiles_1 <- quantile(c(-Inf,X.val[,1],Inf), probs = seq(0, 1, length.out = Group.number[1] + 1))
S.val = as.numeric(cut(X.val[,1], quantiles_1))
Group.val = S.val - 1
G.val = model.matrix(~ factor(Group.val))[, -1]
Ex.val = rep(0.5,n.val)
X.val[,2] = (X.val[,2]>0.75)
mu.val = X.val %*%rnorm(num_features, 1, 1)
tau.val = X.val[,setting] * delta
# tau.val = X.val %*% rep(delta, num_features)

mu0.val <- mu.val - Ex.val * tau.val
mu1.val <- mu0.val + tau.val

W.val <- sapply(Ex.val, function(p) rbinom(1, size = 1, prob = p))
Y.val = mu0.val + W.val*tau.val + rnorm(n.val, 0, sigma)

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
  tau = X[,setting] * delta
  # tau = X %*% rep(delta, num_features)
    
  mu0 <- mu - Ex * tau
  mu1 <- mu0 + tau
  
  W <- sapply(Ex, function(p) rbinom(1, size = 1, prob = p))
  Y = mu0 + W*tau + rnorm(n, 0, sigma)
  
  mu.hat = nuisance.mu(Y,X) 
  
  # Randomization test without model 
  # Randomization tests
  test.stats.value.rt = test.stats.group(Y = Y, W = W, Group = Group, stats = test.stats.method, Ex = Ex, 
                                      mu0.hat = rep(0, n), mu1.hat = rep(0, n), mu.hat = rep(0, n), tau.hat = rep(0, n))
  
  test.stats.ref.rt = t(replicate(M, test.stats.group(Y = Y, W = W, Group = Group, stats = test.stats.method, Ex = Ex, 
                                                   mu0.hat = rep(0, n), mu1.hat = rep(0, n), mu.hat = rep(0, n), tau.hat = rep(0, n), test.index = seq(1, n)))) 
  
  simulation_results_rt[j,] = sapply(seq(1, length(test.stats.value.rt)), function(x) { sum(test.stats.ref.rt[, x]>=test.stats.value.rt[x])/(M+1) + 1/(M+1)})
  
  # FWER
  adjusted_p_values.rt <- p.adjust(simulation_results_rt[j,], method = "holm")
  R$rt[j] = sum((adjusted_p_values.rt <= q))
  
  
  # Randomization test with sample splitting
  train.index.ss = select.train(Group, proportion)
  test.index.ss = setdiff(1:n, train.index.ss)
  # Nuisance function estimation
  tau.hat.ss = nuisance.tau.ss(Y = Y, X = X, Ex = Ex, W = W, mu = mu.hat, train.index = train.index.ss, marginalize = marginalize)
  
  v = sum((tau.val - cbind(1,X.val) %*% tau.hat.ss$beta)**2)/sum((tau.val - mean(tau.val))**2)
  R2$SSRT = c(R2$SSRT, 1-v)  
  if(verbose){cat("R^2 of SSRT:", 1 - v, "\n")}
  
  # Randomization tests
  test.stats.value.ss = test.stats.group(Y = Y[test.index.ss], W = W[test.index.ss], Group = Group[test.index.ss], stats = test.stats.method, Ex = Ex[test.index.ss], 
                                      mu0.hat = tau.hat.ss$mu0[test.index.ss], mu1.hat = tau.hat.ss$mu1[test.index.ss], mu.hat = mu.hat[test.index.ss], tau.hat = tau.hat.ss$tau[test.index.ss])
  
  test.stats.ref.ss = t(replicate(M, test.stats.group(Y = Y[test.index.ss], W = W[test.index.ss], Group = Group[test.index.ss], stats = test.stats.method, Ex = Ex[test.index.ss], 
                                                   mu0.hat = tau.hat.ss$mu0[test.index.ss], mu1.hat = tau.hat.ss$mu1[test.index.ss], mu.hat = mu.hat[test.index.ss], 
                                                   tau.hat = tau.hat.ss$tau[test.index.ss], test.index = test.index.ss))) 
  
  simulation_results_ssrt[j,] = sapply(seq(1, length(test.stats.value.ss)), function(x) { sum(test.stats.ref.ss[, x]>=test.stats.value.ss[x])/(M+1) + 1/(M+1)})
  
  # FWER
  adjusted_p_values.ssrt <- p.adjust(simulation_results_ssrt[j,], method = "holm")
  R$ssrt[j] = sum((adjusted_p_values.ssrt <= q))
  
  
  # Randomization test with adaptive splitting
  # Nuisance function estimation
  tau.hat.adaptive = nuisance.tau.active(Y = Y, X = X, Ex = Ex, W = W, mu = mu.hat, Group = Group, proportion = proportion, marginalize = marginalize) # estimate tau using active learning
  train.index.adaptive = tau.hat.adaptive$train.index # the data points fitted to the nuisance
  test.index.adaptive = setdiff(1:n,train.index.adaptive) # the data points used for inference
  
  v = sum((tau.val - cbind(1,X.val) %*% tau.hat.adaptive$beta)**2)/sum((tau.val - mean(tau.val))**2)
  R2$ART = c(R2$ART, 1-v)
  if(verbose){cat("R^2 of ART:", 1- v, "\n")}
  
  # FWER
  adjusted_p_values.ssrt <- p.adjust(simulation_results_ssrt[j,], method = "holm")
  R$ssrt[j] = sum((as.numeric(adjusted_p_values.ssrt) <= q))
  
  
  # Randomization tests
  test.stats.value.adaptive = test.stats.group(Y = Y[test.index.adaptive], W = W[test.index.adaptive], Group = Group[test.index.adaptive], stats = test.stats.method, Ex = Ex[test.index.adaptive], 
                                             mu0.hat = tau.hat.adaptive$mu0[test.index.adaptive], mu1.hat = tau.hat.adaptive$mu1[test.index.adaptive], mu.hat = mu.hat[test.index.adaptive], 
                                             tau.hat = tau.hat.adaptive$tau[test.index.adaptive])
  
  test.stats.ref.adaptive = t(replicate(M, test.stats.group(Y = Y[test.index.adaptive], W =  W[test.index.adaptive], Group = Group[test.index.adaptive], stats = test.stats.method, Ex = Ex[test.index.adaptive], 
                                                          mu0.hat = tau.hat.adaptive$mu0[test.index.adaptive], mu1.hat = tau.hat.adaptive$mu1[test.index.adaptive], mu.hat = mu.hat[test.index.adaptive], 
                                                          tau.hat = tau.hat.adaptive$tau[test.index.adaptive], test.index=test.index.adaptive[test.index.adaptive])))

  simulation_results_art[j,] = sapply(seq(1, length(test.stats.value.adaptive)), function(x) { sum(test.stats.ref.adaptive[, x]>=test.stats.value.adaptive[x])/(M+1) + 1/(M+1)})
  
  # FWER
  adjusted_p_values.art <- p.adjust(simulation_results_art[j,], method = "holm")
  R$art[j] = sum((as.numeric(adjusted_p_values.art) <= q))
  
  
  if(verbose){
    print(j)
    cat("RT:",round(simulation_results_rt[j,] ,3), "\n")
    cat("SSRT:",round(simulation_results_ssrt[j,],3), "\n")
    cat("ART:",round(simulation_results_art[j,] ,3), "\n")
    cat("\n")
  }  
  
}

# Visualization
# Generate the Boxplots
# pdf(paste0("Comparison_setting_", setting, ".pdf"), width = 13, height = 6)

p_values_df_rt <- as.data.frame(simulation_results_rt)
p_values_df_ssrt <- as.data.frame(simulation_results_ssrt)
p_values_df_art <- as.data.frame(simulation_results_art)

colnames(p_values_df_rt) <- paste("Group", 1:Group.number)
colnames(p_values_df_ssrt) <- paste("Group", 1:Group.number)
colnames(p_values_df_art) <- paste("Group", 1:Group.number)

p_values_df_rt$Method <- "RT"
p_values_df_art$Method <- "ART"
p_values_df_ssrt$Method <- "SSRT"

combined_p_values_df <- rbind(p_values_df_rt, p_values_df_ssrt, p_values_df_art)

combined_p_values_long <- reshape2::melt(combined_p_values_df, id.vars = "Method")
combined_p_values_long$Method <- factor(
  combined_p_values_long$Method,
  levels = c("RT", "SSRT", "ART") 
)

par(mar = c(5, 5, 4, 5))
group_spacing <- 4 
at_positions <- rep(1:Group.number, each = 3) * group_spacing + c(-1, 0, 1)

cols = c("#A3C1AD", "gold", "#8C1515")
boxplot(value ~ Method + variable, data = combined_p_values_long, 
        at = at_positions, col = cols, notch = TRUE, xaxt = "n",
        ylab = "P-values", xlab = "", # Remove the x-axis label
        main = paste0("Comparison of RT, ART, SSRT (Maximum Inference Fold: ", (1-proportion)*100, "%)"),
        cex.axis = 1.5, 
        cex.lab = 1.5,   
        cex.main = 1.6)   

axis(1, at = 1:Group.number * group_spacing, labels = paste("Group", 1:Group.number),
     cex.axis = 1.5)  

legend("topright", legend = c("RT", "SSRT", "ART"), fill = cols, 
       cex = 0.8, 
       pt.cex = 2,  
       bg = "white",
       horiz = T) 


# dev.off()


# results
data.frame("Number of rejected subgroups" = lapply(R, mean))

# nuisance worse?
data.frame("Nuisance estimation R2" = lapply(R2, mean), 
           "Nuisance estimation R2 sd" = lapply(R2, sd))

data.frame("RT" = round(colMeans(simulation_results_rt[1:j,]),2),
           "SSRT" = round(colMeans(simulation_results_ssrt[1:j,]),2),
           "ART" = round(colMeans(simulation_results_art[1:j,]),2))









