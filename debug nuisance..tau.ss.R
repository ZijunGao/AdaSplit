# debug nuisance.tau.ss
# source file
source("~/Desktop/Research/Yao/HTE inference/code/Panning/helper Yao.R")


# data generation mechanism
n = 100 # sample size; 40; 800
d = 5 # number of covariates
beta0 = 1; beta = rep(1, d) # coefficients for outcome model under control: mu0(x) = beta0 + X %*% beta
theta = rep(0, d); theta[1] = 2; theta0 = 0.5 # propensity score model: sigmoid(X %*% theta) * ((2 * theta0)
delta = 0.3 # magnitude of average treatment effect (ATE); e.g., 0.5 or 
delta0 = 1 # magnitude of treatment effect heterogeneity
sigma = 1 # standard deviation of errors in Y(0)

nuisance.proportion = 0.5 # proportion of data used for nuisance function training
mu.learner = "linear" # outcome regression method; options: "linear", "xgboost"


# data generation
set.seed(318)
X = matrix(rnorm(n * d, 0, 1), nrow = n, ncol = d) # covariates
p = 2 * theta0 * exp(X %*% theta) / (1 + exp(X %*% theta))
W = rbinom(n, 1, p) # treatment assignment
# potential outcomes
mu0 = (beta0 + X %*% beta) / 1 # mu0 is linear in X and S: (beta0 + X %*% beta + 2* X[, 1]^2) / 1
Y0 = mu0 + rnorm(n, 0, sigma)
tau =  delta - (X[,1] + X[,2] + X[,3]) * delta0 + rnorm(n, 0, 0) # (1 + X[, 2]^2); delta + 2 * X[,2]
mu1 = mu0 + tau
mu = mu0 * (1 - p) + mu1 * p
Y1 = Y0 + tau
Y = Y1 * W + Y0 * (1 - W) # observed outcome


# estimated nuisance functions
train.index.ss = sample(1:n, nuisance.proportion * n, replace = F)
if(mu.learner == "linear"){mu.hat = nuisance.mu(Y,X)
}else if(mu.learner == "xgboost"){mu.hat = nuisance.mu.xgboost(Y,X)}
tau.hat.ss = nuisance.tau.ss(Y = Y, X = X, Ex = p, W = W, mu = mu.hat, train.index = train.index.ss, marginalize = F)


# results
par(mfrow = c(1, 2))
plot(mu, mu.hat, main = "mu"); abline(0, 1, col = "red")
plot(tau, tau.hat.ss$tau, main = "tau"); abline(0, 1, col = "red")
data.frame(true.coef = c(delta, -delta0 * c(1,1,1,rep(0, d - 3))), est.coef = tau.hat.ss$beta)

