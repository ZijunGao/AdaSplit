# Partial conjunction null
# Randomized p-values
  # Seems not very helpful: adding randomness loses information?
n = 100
pi0 = 0.5
n0 = floor(n * pi0); n1 = n - n0
n0Null = floor(n * (pi0 + 0.1)) # number of nulls in the hypothesis 

delta = 2/3

m = 100
record = rep(0, m)
# set.seed(318)

# binary p-values
pval = c(rbinom(n0, 1, delta), rbinom(n1, 1, 0.05))
pval.combine = pbinom(sum(pval == 1), n0Null, delta)

pval.randomize = (1 - delta + runif(n, 0, 1) * delta) * (pval == 1) + (delta + runif(n, 0, 1) * (1 - delta)) * (pval == 0)
pval.randomize.combine = 1 - pchisq(- 2 * sum(log(sort(pval.randomize, decreasing = T)[1 : n0Null])), 2 * n0Null)

# results
print("binomial p-value"); print(pval.combine)
print("randomized p-value"); print(pval.randomize.combine)





