source("~/Desktop/Research/Yao/HTE inference/code/Panning/helper Yao.R")

n = 1000 #sample size
Group.number = 8 # total number of groups
delta = 0.1 # effect size
sigma = 0.2 # std of noise
M = 1000 # number of permutations
n_trial = 30 # number of runs
m_max <- 100 # maximum number of knockoffs
proportion = 0.2 # proportion of randomness in the inference fold
test.stats.method = "AIPW" #Test statistics
nuisance.learner.method = "gradient boosting" #Nuisance estimator



Entropy_Z_given_Z_bar <- function(m, p) {
  
  # p: propensity score e(X)
  # m: number of knockoffs
  # Return the conditional entropy H(Z|Z_bar)
  
  # All possible values of Z_bar
  z_bar <- 0:(m + 1)
  
  # Binomial probabilities P(Z_bar = z)
  P_Z_bar <- dbinom(z_bar, size = m + 1, prob = p)
  
  # Conditional probabilities P(Z = 1 | Z_bar = z_bar)
  prob_Z_given_Z_bar <- z_bar / (m + 1)
  
  # Calculate entropy terms, handle 0 and 1 cases to avoid log(0)
  entropy <- ifelse(prob_Z_given_Z_bar > 0, - prob_Z_given_Z_bar * log(prob_Z_given_Z_bar), 0) +
    ifelse(prob_Z_given_Z_bar < 1, -(1 - prob_Z_given_Z_bar) * log(1 - prob_Z_given_Z_bar), 0)
  
  # Marginalize out Z_bar to compute the conditional entropy H(Z|Z_bar)
  H_Z_given_Z_bar <- sum(P_Z_bar * entropy)
  
  return(H_Z_given_Z_bar)
}

Percentage_entropy_left <- function(m, p) {
  
  # p: propensity score e(X)
  # m: number of knockoffs
  # Return the remaining entropy (%) for inference
  
  # Compute the total entropy H(Z)
  H_Z <- -p * log(p) - (1 - p) * log(1 - p)
  
  # Calculate H(Z | Z_bar)
  H_Z_given_Z_bar <- Entropy_Z_given_Z_bar(m, p)
  
  # Compute percentage of remaining entropy
  r <- (H_Z - H_Z_given_Z_bar) / H_Z
  return(r)
}



Entropy_matrix <- function(Ex, m_max) {
  
  # Ex: propensity scores for all individuals
  # m_max: the largest number of knockoffs
  # Return the remaining entropy (%) for inference for all e_x and number of knockoffs m = 0,1,...,m_max
  # PS: Every row of the output matrix is a decreasing vector from 1 to nearly 0.
  
  n <- length(Ex)
  
  # Initialize the matrix to store the entropy (%)
  entropy_matrix <- matrix(0, nrow = n, ncol = m_max + 1)
  
  # Loop over each value of p in Ex
  for (i in 1:n) {
    p <- Ex[i]
    
    # Loop over m from 0 to m_max
    for (m in 0:m_max) {
      entropy_matrix[i, m + 1] <- Percentage_entropy_left(m, p)
    }
  }
  
  return(entropy_matrix)
}

Find_interpolation_indices <- function(entropy_row, proportion) {
  
  # entropy_row: a row in the entropy matrix above
  # proportion: user-specified proportion of entropy for inference
  # Return the positions of the entropy right below and after the desire entropy (%)
  
  idx_before <- max(which(entropy_row >= proportion))
  idx_after <- idx_before + 1
  
  return(list(m_j1 = idx_before, m_j2 = idx_after))
}


Compute_interpolation_entropy <- function(value1, value2, proportion) {
  
  # value1, value2 are the entropy at positive m_j1 and m_j2
  # Return a probability for randomization to achieve the entropy proportion
  # PS: this randomization scheme itself has the lowest variance by consider two adjacent values of m.
  
  r_j <- (proportion - value2) / (value1 - value2)
  return(r_j)
}


Compute_policy_matrix <- function(entropy_matrix, proportion) {
  
  # Return the policy to randomly choose m for every row in entropy_matrix
  
  n <- nrow(entropy_matrix)
  m_max <- ncol(entropy_matrix) - 1
  
  # Initialize the zero matrix G
  Policy_matrix <- matrix(0, nrow = n, ncol = m_max + 1)
  
  # Loop through each row
  for (j in 1:n) {
    # Find indices before and after the proportion
    indices <- Find_interpolation_indices(entropy_matrix[j, ], proportion)
    m_j1 <- indices$m_j1
    m_j2 <- indices$m_j2
    
    # Extract values before and after the proportion
    value1 <- entropy_matrix[j, m_j1]
    value2 <- entropy_matrix[j, m_j2]

    # Compute r_j for the interpolation
    r_j <- Compute_interpolation_entropy(value1, value2, proportion)

    # Update Policy_matrix
    Policy_matrix[j, m_j1] <- r_j
    Policy_matrix[j, m_j2] <- 1 - r_j
  }
  
  return(Policy_matrix)
}



Sample_from_policy <- function(G_matrix) {
  # Randomly choose the number of knockoffs following the policy G
  
  samples <- apply(G_matrix, 1, function(probabilities) {
    sample(0: m_max, size = 1, prob = probabilities) 
  })
  
  return(samples)
}



# Store the experimental results
simulation_results_art <- matrix(0, nrow = n_trial, ncol = Group.number)
colnames(simulation_results_art) <- paste0("p", seq(1, Group.number))

simulation_results_ssrt <- matrix(0, nrow = n_trial, ncol = Group.number)
colnames(simulation_results_ssrt) <- paste0("p", seq(1, Group.number))

for (j in 1:n_trial){
  
  # Data generation
  X = matrix(runif(n, -3, 3), nrow = n, ncol = 1) 
  quantiles_1 <- quantile(c(-Inf,X,Inf), probs = seq(0, 1, length.out = Group.number[1] + 1))
  S = as.numeric(cut(X, quantiles_1))
  
  Group = S - 1
  G = model.matrix(~ factor(Group))[, -1]
  
  Ex = exp(1.3*X)/(1+exp(1.3*X))
  
  W <- sapply(Ex, function(p) rbinom(1, size = 1, prob = p))
  mu0 =  1 + X
  Y0 = mu0 + rnorm(n, 0, sigma)
  tau = delta*(abs(X)+1)
  Y1 = Y0 + tau 
  Y = Y1 * W + Y0 * (1 - W) 
  mu1 = mu0 + tau
  mu = mu0 * (1 - Ex) + mu1 * Ex

  
  #Run our method to decide the number of knockoffs
  entropy_matrix <- Entropy_matrix(Ex, m_max)
  policy_matrix <- Compute_policy_matrix(entropy_matrix, proportion)
  m_samples <- Sample_from_policy(policy_matrix)
  
  
  #Fit the nuisance parameter
  W.knockoff <- mapply(function(p, n) rbinom(n, size = 1, prob = p), Ex, m_samples, SIMPLIFY = FALSE) 
  W.aug <- mapply(function(w, knockoff) c(w, knockoff), W, W.knockoff, SIMPLIFY = FALSE)
  W.tilde <- sapply(W.aug, mean)
  nuisance.hat = nuisance.learner(Y = Y, X = X, prop = Ex, G = G, W =W.tilde, method = nuisance.learner.method, train.index = seq(1, n), test.index = seq(1, n))
  
  #plot(X,nuisance.hat$tau.hat, col = 'blue',pch = 17, main='ART-learner')
  #points(X,tau, col = 'red', pch = 17)
  
  #Randomization tests
  test.stats.value = test.stats.group(Y = Y, W = W, X = X, Group = Group, G = G, stats = test.stats.method, prop = Ex, mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, 
                                      mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat, W.tilde=W.tilde)
  
  test.stats.ref = t(replicate(M, test.stats.group(Y = Y, W = sapply(W.aug, function(x) sample(x, 1)), X = X, G = G, 
                                                   Group = Group, stats = test.stats.method, prop = Ex, mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, 
                                                   mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat, W.tilde=W.tilde))) 
  
  # Calculate p-values for each group
  simulation_results_art[j,] = sapply(seq(1, length(test.stats.value)), function(x) {
    permutation.p.value(stats = test.stats.value[x], stats.ref = test.stats.ref[, x])
  })
  
  
  # Sample splitting
  nuisance.index = c()
  nuisance.index = sample(n, n * (1-proportion)) 
  train.index = nuisance.index 
  test.index = setdiff(seq(1, n), nuisance.index)
  
  # Fit the nuisance parameter
  nuisance.hat = nuisance.learner(Y = Y, X = X, prop = Ex, G = G, W = W, method = nuisance.learner.method, train.index = train.index, test.index = test.index)
  
  
  #plot(X[test.index],nuisance.hat$tau.hat, col = 'blue',pch = 17,,main='SSRT-learner')
  #points(X[test.index],tau[test.index], col = 'red', pch = 17)
  
  
  # Randomization tests
  test.stats.value = test.stats.group(Y = Y[test.index], W = W[test.index], X = X[test.index,], Group = Group[test.index], G = G[test.index,], stats = test.stats.method, prop = Ex[test.index], 
                                      mu0.hat = nuisance.hat$mu0.hat, mu1.hat = nuisance.hat$mu1.hat, mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat)
  
  test.stats.ref = t(replicate(M, test.stats.group(Y = Y[test.index], W =  sapply(Ex[test.index], function(p) rbinom(1, size = 1, prob = p)), X = X[test.index,], G = G[test.index,], 
                                                   Group = Group[test.index], stats = test.stats.method, prop = Ex[test.index], mu0.hat = nuisance.hat$mu0.hat, 
                                                   mu1.hat = nuisance.hat$mu1.hat, mu.hat = nuisance.hat$mu.hat, tau.hat = nuisance.hat$tau.hat))) 
  
  # Calculate p-values for each group
  simulation_results_ssrt[j,] = sapply(seq(1, length(test.stats.value)), function(x) {
    permutation.p.value(stats = test.stats.value[x], stats.ref = test.stats.ref[, x])
  })
  
  
  if (j>=2){
    cat(" ART:", format(round(colMeans(simulation_results_art[1:j,]),2),2), "\n")
    cat("SSRT:", format(round(colMeans(simulation_results_ssrt[1:j,]),2),2), "\n")
    cat("\n")
  }


}




# Generate the Boxplots
pdf(paste0("~/Desktop/Research/Yao/HTE inference/code/Panning/0630/Comparison_plot_", proportion, ".pdf"), width = 13, height = 6)

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


