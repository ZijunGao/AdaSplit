# Power comparison of BH and Seq-Step

# Alternative distribution. A mixture of Beta(a, b) and uniform.
F1 = function(x, a, b) {
  (1-delta) * pbeta(x, shape1 = a, shape2 = b) + delta * pbeta(x, shape1 = 1, shape2 = 1)
}
f1 = function(x, a, b) {
  (1-delta) * dbeta(x, shape1 = a, shape2 = b) + delta * dbeta(x, shape1 = 1, shape2 = 1)
}
# Null distribution. Uniform.
F0 = function(x) {
  pbeta(x, 1, 1)
}
f0 = function(x) {
  dbeta(x, 1, 1)
}
# Marginal distribution.
FF = function(t) {
  pi0 * F0(t) + (1 - pi0) * F1(t, a, b)
}
f = function(t) {
  pi0 * f0(t) + (1 - pi0) * f1(t, a, b)
}

# Parameters
a = 1; b = 10
delta = 0.2
pi0 = 0.3 # Proportion of null hypotheses
q = 0.2  # FDR level
lambda = 0.5 
C = 2

# BH
# Compute Storey's null estimator
pi0_lambda = (1 - FF(lambda)) / (1 - lambda)


# Compute FDR(t)
FDR = function(t) {
  pi0_lambda * F0(t) / FF(t)
}


# Compute tau_BH by solving FDR(t) = q
tau_BH = uniroot(function(t) FDR(t) - q, lower = 0.001, upper = 1)$root

# Compute Power_BH
Power_BH = F1(tau_BH, a, b)

# Seq-Step
# Compute mu
mu = C * (1 - F1(1 - 1 / C, a, b))

# Define the function g(t)
g = function(t) {
  (1 - pi0) * f1(t, a, b) / (pi0 * f0(t) + (1 - pi0) * f1(t, a, b))
}

# Compute tau_Seq-Step
value = (1 - q) / (1 - mu)
tau_Seq_Step = ifelse(
  value >= g(0), 
  0,
  ifelse(
    g(1) < value & value < g(0), 
    uniroot(function(t) g(t) - value, lower = 0, upper = 1)$root, 
    1
  )
)

# Compute Power_Seq-Step
Power_Seq_Step = tau_Seq_Step * g(tau_Seq_Step) / g(1)

# Print results
cat("Power_BH:", Power_BH, "\n")
cat("Power_Seq_Step:", Power_Seq_Step, "\n")

# Sanity check
plot(g(seq(1, 100)/100) * seq(1, 100)) # should be increasing


