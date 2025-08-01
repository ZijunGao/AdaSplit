# Simulation
rm(list = ls())
source("~/Desktop/Research/Yao/HTE inference/code/Panning/helper Yao.R")
# source("/Users/yaozhang/Documents/GitHub/Panning/helper Yao.R")

directory = "~/Desktop/Research/Yao/HTE inference/code/Panning/July 2025"
  
n = 500 # sample size: 500, 1000
num_features = 5
Group.number = 5 # total number of groups
delta = 1 # effect size; 1, 0

sigma = 1 # std of noise: 1, sqrt(2)
M = 1000 # number of permutations
proportion = 0.5 # proportion of randomness in the nuisance fold
test.stats.method ="AIPW" # test statistics

n_trial = 100 # number of runs; 20; 100                                                              
verbose = F # if true, print out intermediate results
q = 0.2 # FWER level
marginalize = T # if marginalize
mu.learner = "linear"
settings = c("default",
             "larger sample size",
             "larger noise",
             "even larger noise",
             "larger sample size larger noise",
             "null",
             "no marginalization",
             "mu xgboost",
             "fewer for nuisance",
             "fewer for inference",
             "null more repeats",
             "no marginalization larger sample size",
             "no marginalization larger noise",
             "smaller sample size")
setting = "no marginalization larger noise" # "default", "larger sample size", "larger noise", "null more repeats"
if(setting == "larger sample size"){n = 1000}
if(setting == "larger noise"){sigma = sqrt(2)}
if(setting == "even larger noise"){sigma = 2}
if(setting == "larger sample size larger noise"){n = 1000; sigma = sqrt(2)}
if(setting == "null"){delta = 0}
if(setting == "no marginalization"){marginalize = F}
if(setting == "mu xgboost"){mu.learner = "xgboost"}
if(setting == "fewer for nuisance"){proportion = 0.25}
if(setting == "fewer for inference"){proportion = 0.75}
if(setting == "null more repeats"){delta = 0; n = 2000; n_trial = 200}
if(setting == "no marginalization larger sample size"){n = 1000; marginalize = F}
if(setting == "no marginalization larger noise"){sigma = sqrt(2); marginalize = F}
if(setting == "smaller sample size"){n = 300}

set.seed(318)

# Record
# Subgroup p-values
pval = list()
pval$ART <- matrix(0, nrow = n_trial, ncol = Group.number)
colnames(pval$ART) <- paste0("p", seq(1, Group.number))
pval$ART.no.weighting = pval$SSSSRT = pval$SSRT = pval$RT = pval$ART

# FWER
FWER = list(); FWER$ART = FWER$ART.no.weighting = FWER$SSRT = FWER$SSSSRT = FWER$RT = c()

# Number of subgroup rejections
R = list(); R$ART = R$ART.no.weighting = R$SSRT = R$SSSSRT = R$RT = c()

# Quality of nuisance estimator
R2 = list(); R2$ART = R2$ART.no.weighting = R2$SSRT = R2$SSSSRT = c()

# Proportion of inference fold
inference.proportion = list(); inference.proportion$ART.no.weighting = inference.proportion$ART = inference.proportion$SSRT = inference.proportion$SSSSRT = matrix(0, nrow = n_trial, ncol = Group.number)

# Average treatment effect in the inference fold
inference.ATE = list(); inference.ATE$ART.no.weighting = inference.ATE$ART = inference.ATE$SSRT = inference.ATE$SSSSRT = matrix(0, nrow = n_trial, ncol = Group.number)

# Data for evaluating the nuisance estimation performance
n.val = 10^4
X.val = matrix(runif(n.val*num_features, 0, 1), nrow = n.val, ncol = num_features) 
quantiles_1 <- quantile(c(-Inf,X.val[,1],Inf), probs = seq(0, 1, length.out = Group.number[1] + 1))
S.val = as.numeric(cut(X.val[,1], quantiles_1))
Group.val = S.val - 1
G.val = model.matrix(~ factor(Group.val))[, -1]
Ex.val = rep(0.5,n.val)
X.val[,2] = (X.val[,2]>0.75); X.val[,3] = (X.val[,3]>0.25)
mu.val = X.val %*%rnorm(num_features, 1, 1)
tau.val = (X.val - 0.5) %*% rep(delta, num_features) + 0.5 * delta
H = rep(1, Group.number) # H = 0: under the null; H = 1: under the alternative
if(setting == "null"){H = rep(0, Group.number)}
if(setting == "null more repeats"){H = rep(0, Group.number)}

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
  X[,3] = (X[,3]>0.25); tau = (X - 0.5) %*% rep(delta, num_features) + 0.5 * delta
    
  mu0 <- mu - Ex * tau
  mu1 <- mu0 + tau
  
  W <- sapply(Ex, function(p) rbinom(1, size = 1, prob = p))
  Y = mu0 + W*tau + rnorm(n, 0, sigma)
  
  if(mu.learner == "linear"){mu.hat = nuisance.mu(Y,X)
  }else if(mu.learner == "xgboost"){mu.hat = nuisance.mu.xgboost(Y,X)}
  
  # Standard/vanilla randomization test
  # Randomization tests
  test.stats.value.rt = test.stats.group(Y = Y, W = W, Group = Group, stats = test.stats.method, Ex = Ex, 
                                      mu0.hat = rep(0, n), mu1.hat = rep(0, n), mu.hat = rep(0, n), tau.hat = rep(0, n))
  
  test.stats.ref.rt = t(replicate(M, test.stats.group(Y = Y, W = W, Group = Group, stats = test.stats.method, Ex = Ex, 
                                                   mu0.hat = rep(0, n), mu1.hat = rep(0, n), mu.hat = rep(0, n), tau.hat = rep(0, n), test.index = seq(1, n))))  
  
  pval$RT[j,] = sapply(seq(1, length(test.stats.value.rt)), function(x) { sum(test.stats.ref.rt[, x]>=test.stats.value.rt[x])/(M+1) + 1/(M+1)})
  
  # FWER
  R$rt[j] = sum(closing(p_val = pval$RT[j,], q = q, global.null.test = "Fisher"))
  FWER$rt[j] = 1 - min(1, H[rank(pval$RT[j,]) <= R$rt[j]])
   
  # Randomization test with sample splitting
  train.index.ss = select.train(Group, proportion)
  test.index.ss = setdiff(1:n, train.index.ss)
  inference.proportion$SSRT[j,] = table(Group[test.index.ss])
  inference.ATE$SSRT[j,] = aggregate(tau[test.index.ss], by = list(Group[test.index.ss]), mean)$x
  
  # Nuisance function estimation
  tau.hat.ss = nuisance.tau.ss(Y = Y, X = X, Ex = Ex, W = W, mu = mu.hat, train.index = train.index.ss, marginalize = marginalize)
  
  v = sum((tau.val - cbind(1,X.val) %*% tau.hat.ss$beta)**2)/sum((tau.val - mean(tau.val))**2)
  R2$SSRT = c(R2$SSRT, 1-v)  
  if(verbose){cat("R^2 of SSRT:", 1 - v, "\n")}
  
  # Randomization test
  test.stats.value.ss = test.stats.group(Y = Y[test.index.ss], W = W[test.index.ss], Group = Group[test.index.ss], stats = test.stats.method, Ex = Ex[test.index.ss], 
                                      mu0.hat = tau.hat.ss$mu0[test.index.ss], mu1.hat = tau.hat.ss$mu1[test.index.ss], mu.hat = mu.hat[test.index.ss], tau.hat = tau.hat.ss$tau[test.index.ss])
  
  test.stats.ref.ss = t(replicate(M, test.stats.group(Y = Y[test.index.ss], W = W[test.index.ss], Group = Group[test.index.ss], stats = test.stats.method, Ex = Ex[test.index.ss], 
                                                   mu0.hat = tau.hat.ss$mu0[test.index.ss], mu1.hat = tau.hat.ss$mu1[test.index.ss], mu.hat = mu.hat[test.index.ss], 
                                                   tau.hat = tau.hat.ss$tau[test.index.ss], test.index = test.index.ss))) 
  
  pval$SSRT[j,] = sapply(seq(1, length(test.stats.value.ss)), function(x) { sum(test.stats.ref.ss[, x]>=test.stats.value.ss[x])/(M+1) + 1/(M+1)})

  # FWER
  R$ssrt[j] = sum(closing(p_val = pval$SSRT[j,], q = q, global.null.test = "Fisher"))
  FWER$ssrt[j] = 1 - min(1, H[rank(pval$SSRT[j,]) <= R$ssrt[j]])
  
  
  # Randomization test with subgroup selection
  test.index.ss.ss = intersect(test.index.ss, which(tau.hat.ss$tau >= 0)) 
  test.stats.value.ss.ss = test.stats.group(Y = Y[test.index.ss.ss], W = W[test.index.ss.ss], Group = Group[test.index.ss.ss], stats = test.stats.method, Ex = Ex[test.index.ss.ss], 
                                         mu0.hat = tau.hat.ss$mu0[test.index.ss.ss], mu1.hat = tau.hat.ss$mu1[test.index.ss.ss], mu.hat = mu.hat[test.index.ss.ss], tau.hat = tau.hat.ss$tau[test.index.ss.ss])
  
  test.stats.ref.ss.ss = t(replicate(M, test.stats.group(Y = Y[test.index.ss.ss], W = W[test.index.ss.ss], Group = Group[test.index.ss.ss], stats = test.stats.method, Ex = Ex[test.index.ss.ss], 
                                                      mu0.hat = tau.hat.ss$mu0[test.index.ss.ss], mu1.hat = tau.hat.ss$mu1[test.index.ss.ss], mu.hat = mu.hat[test.index.ss.ss], 
                                                      tau.hat = tau.hat.ss$tau[test.index.ss.ss], test.index = test.index.ss.ss))) 
  
  pval$SSSSRT[j, ] = 1
  pval$SSSSRT[j, 1 + sort(unique(Group[test.index.ss.ss]))] = sapply(seq(1, length(test.stats.value.ss.ss)), function(x) { sum(test.stats.ref.ss.ss[, x]>=test.stats.value.ss.ss[x])/(M+1) + 1/(M+1)}) # If a subgroup has no units left for inference, we set the associated p-value to one.
  
  
  # FWER
  R$ssssrt[j] = sum(closing(p_val = pval$SSSSRT[j,], q = q, global.null.test = "Fisher"))
  FWER$ssssrt[j] = 1 - min(1, H[rank(pval$SSSSRT[j,]) <= R$ssssrt[j]])

  # Randomization test with adaptive splitting
  # Nuisance function estimation
  tau.hat.adaptive = nuisance.tau.active(Y = Y, X = X, Ex = Ex, W = W, mu = mu.hat, Group = Group, proportion = proportion, marginalize = marginalize) # estimate tau using active learning
  train.index.adaptive = tau.hat.adaptive$train.index # the data points fitted to the nuisance
  test.index.adaptive = setdiff(1:n,train.index.adaptive) # the data points used for inference
  inference.proportion$ART[j,] = table(Group[test.index.adaptive])
  inference.ATE$ART[j,] = aggregate(tau[test.index.adaptive], by = list(Group[test.index.adaptive]), mean)$x
  
  v = sum((tau.val - cbind(1,X.val) %*% tau.hat.adaptive$beta)**2)/sum((tau.val - mean(tau.val))**2)
  R2$ART = c(R2$ART, 1-v)
  if(verbose){cat("R^2 of ART:", 1- v, "\n")}
  
  # Randomization tests
  test.stats.value.adaptive = test.stats.group(Y = Y[test.index.adaptive], W = W[test.index.adaptive], Group = Group[test.index.adaptive], stats = test.stats.method, Ex = Ex[test.index.adaptive], 
                                             mu0.hat = tau.hat.adaptive$mu0[test.index.adaptive], mu1.hat = tau.hat.adaptive$mu1[test.index.adaptive], mu.hat = mu.hat[test.index.adaptive], 
                                             tau.hat = tau.hat.adaptive$tau[test.index.adaptive])
  
  test.stats.ref.adaptive = t(replicate(M, test.stats.group(Y = Y[test.index.adaptive], W =  W[test.index.adaptive], Group = Group[test.index.adaptive], stats = test.stats.method, Ex = Ex[test.index.adaptive], 
                                                          mu0.hat = tau.hat.adaptive$mu0[test.index.adaptive], mu1.hat = tau.hat.adaptive$mu1[test.index.adaptive], mu.hat = mu.hat[test.index.adaptive], 
                                                          tau.hat = tau.hat.adaptive$tau[test.index.adaptive], test.index=test.index.adaptive[test.index.adaptive])))

  pval$ART[j,] = sapply(seq(1, length(test.stats.value.adaptive)), function(x) { sum(test.stats.ref.adaptive[, x]>=test.stats.value.adaptive[x])/(M+1) + 1/(M+1)})
  # FWER
  R$art[j] = sum(closing(p_val = pval$ART[j,], q = q, global.null.test = "Fisher"))
  FWER$art[j] = 1 - min(1, H[rank(pval$ART[j,]) <= R$art[j]])
  
  
  # AdaSplit without weighting
  # Randomization test with adaptive splitting
  # Nuisance function estimation
  tau.hat.adaptive = nuisance.tau.active(Y = Y, X = X, Ex = Ex, W = W, mu = mu.hat, Group = Group, proportion = proportion, marginalize = marginalize, weighting = F) # estimate tau using active learning
  train.index.adaptive = tau.hat.adaptive$train.index # the data points fitted to the nuisance
  test.index.adaptive = setdiff(1:n,train.index.adaptive) # the data points used for inference
  inference.proportion$ART.no.weighting[j,] = table(Group[test.index.adaptive])
  inference.ATE$ART.no.weighting[j,] = aggregate(tau[test.index.adaptive], by = list(Group[test.index.adaptive]), mean)$x
  
  v = sum((tau.val - cbind(1,X.val) %*% tau.hat.adaptive$beta)**2)/sum((tau.val - mean(tau.val))**2)
  R2$ART.no.weighting = c(R2$ART.no.weighting, 1-v)
  if(verbose){cat("R^2 of ART.no.weighting:", 1- v, "\n")}
  
  # Randomization tests
  test.stats.value.adaptive = test.stats.group(Y = Y[test.index.adaptive], W = W[test.index.adaptive], Group = Group[test.index.adaptive], stats = test.stats.method, Ex = Ex[test.index.adaptive], 
                                               mu0.hat = tau.hat.adaptive$mu0[test.index.adaptive], mu1.hat = tau.hat.adaptive$mu1[test.index.adaptive], mu.hat = mu.hat[test.index.adaptive], 
                                               tau.hat = tau.hat.adaptive$tau[test.index.adaptive])
  
  test.stats.ref.adaptive = t(replicate(M, test.stats.group(Y = Y[test.index.adaptive], W =  W[test.index.adaptive], Group = Group[test.index.adaptive], stats = test.stats.method, Ex = Ex[test.index.adaptive], 
                                                            mu0.hat = tau.hat.adaptive$mu0[test.index.adaptive], mu1.hat = tau.hat.adaptive$mu1[test.index.adaptive], mu.hat = mu.hat[test.index.adaptive], 
                                                            tau.hat = tau.hat.adaptive$tau[test.index.adaptive], test.index=test.index.adaptive[test.index.adaptive])))
  
  pval$ART.no.weighting[j,] = sapply(seq(1, length(test.stats.value.adaptive)), function(x) { sum(test.stats.ref.adaptive[, x]>=test.stats.value.adaptive[x])/(M+1) + 1/(M+1)})
  
  # FWER
  R$art.no.weighting[j] = sum(closing(p_val = pval$ART.no.weighting[j,], q = q, global.null.test = "Fisher"))
  FWER$art.no.weighting[j] = 1 - min(1, H[rank(pval$ART.no.weighting[j,]) <= R$art.no.weighting[j]])
  
  print(j)
  if(verbose){
    cat("RT:",round(pval$RT[j,] ,3), "\n")
    cat("SSRT:",round(pval$SSRT[j,],3), "\n")
    cat("ART:",round(pval$ART[j,] ,3), "\n")
    cat("ART no weighting:",round(pval$ART.no.weighting[j,] ,3), "\n")
    cat("\n")
  }  
  
}

# Visualization
# Generate the Boxplots

p_values_df_rt <- as.data.frame(pval$RT)
p_values_df_ssrt <- as.data.frame(pval$SSRT)
p_values_df_ssssrt <- as.data.frame(pval$SSSSRT)
p_values_df_art <- as.data.frame(pval$ART)
p_values_df_art_no_weighting <- as.data.frame(pval$ART.no.weighting)

colnames(p_values_df_rt) <- paste("Group", 1:Group.number)
colnames(p_values_df_ssrt) <- paste("Group", 1:Group.number)
colnames(p_values_df_ssssrt) <- paste("Group", 1:Group.number)
colnames(p_values_df_art) <- paste("Group", 1:Group.number)
colnames(p_values_df_art_no_weighting) <- paste("Group", 1:Group.number)

p_values_df_rt$Method <- "RT"
p_values_df_art$Method <- "ART"
p_values_df_art_no_weighting$Method <- "ART no weighting"
p_values_df_ssrt$Method <- "SSRT"
p_values_df_ssssrt$Method <- "SSSSRT"

combined_p_values_df <- rbind(p_values_df_rt, p_values_df_ssrt, p_values_df_ssssrt, p_values_df_art, p_values_df_art_no_weighting)

combined_p_values_long <- reshape2::melt(combined_p_values_df, id.vars = "Method")
combined_p_values_long$Method <- factor(
  combined_p_values_long$Method,
  levels = c("RT", "SSRT", "SSSSRT", "ART", "ART no weighting") 
)

par(mar = c(5, 5, 4, 5))
group_spacing <- 6
at_positions <- rep(1:Group.number, each = 5) * group_spacing + c(-1, 0, 1, 2, 3) # each = number of methods

cols = c("#A3C1AD", "gold", "coral",  "light blue", "#8C1515")
boxplot(value ~ Method + variable, data = combined_p_values_long, 
        at = at_positions, col = cols, notch = F, xaxt = "n",
        ylab = "P-values", xlab = "", # Remove the x-axis label
        main = paste0("Comparison of RT, ART, SSRT (Maximum Inference Fold: ", (1-proportion)*100, "%)"),
        cex.axis = 1.5, 
        cex.lab = 1.5,   
        cex.main = 1.6)   

axis(1, at = 1:Group.number * group_spacing, labels = paste("Group", 1:Group.number),
     cex.axis = 1.5)  

legend("topright", legend = c("RT", "SSRT", "SSSSRT", "ART", "ART no weighting"), fill = cols, 
       cex = 0.8, 
       pt.cex = 2,  
       bg = "white",
       horiz = T) 

# results
data.frame("Number of rejected subgroups" = lapply(R, mean))

data.frame("Proportion of inference fold" = lapply(inference.proportion, function(x){apply(x, 2, mean)}))

data.frame("ATE of inference fold" = lapply(inference.ATE, function(x){apply(x, 2, mean)}))

# nuisance worse?
data.frame("Nuisance estimation R2" = lapply(R2, mean), 
           "Nuisance estimation R2 sd" = lapply(R2, sd))

data.frame("RT" = round(colMeans(pval$RT[1:j,]),2),
           "SSRT" = round(colMeans(pval$SSRT[1:j,]),2),
           "SSSSRT" = round(colMeans(pval$SSSSRT[1:j,]),2),
           "ART" = round(colMeans(pval$ART[1:j,]),2),
           "ART no weighting" = round(colMeans(pval$ART.no.weighting[1:j,]),2))


# save data
record = list()
record$pval = pval; record$R = R; record$R2 = R2; record$FWER = FWER
record$inference.ATE = inference.ATE; record$inference.proportion = lapply(inference.proportion, function(x){x / (n / Group.number)})
# saveRDS(record, file = paste(file.path(directory, setting), ".rds", sep = ""))
# record2 = record





