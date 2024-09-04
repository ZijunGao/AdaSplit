# Partial conjunction null
# Extend the binomial distribution test

n = 100
pi0 = 0.5
n0 = floor(n * pi0); n1 = n - n0
n0Null = floor(n * (pi0 + 0.1))
n.eta.seq = seq(-30, 21, by = 1)

delta = 2 / 3
m = 400
record = matrix(0, nrow = m, ncol = length(n.eta.seq))
set.seed(318)
for(j in 1 : m){
  pval.0 = c(rbinom(n0 / 2, 1, delta) * delta + (1 - delta), 
           rbinom(n1 / 2, 1, 0.05) * delta + (1 - delta))
  pval.1 = c(rbinom(n0 / 2, 1, 1 - delta) * (1 - delta) + delta, 
           rbinom(n1 / 2, 1, 0.05) * (1 - delta) + delta)
  for(k in 1:length(n.eta.seq)){
    pval.combine.0 = pbinom(sum(pval.0 == 1), n0Null / 2 + n.eta.seq[k], delta)
    pval.combine.1 = pbinom(sum(pval.1 == 1), n0Null / 2 - n.eta.seq[k], 1 - delta)
    record[j, k] = max(pval.combine.0, pval.combine.1) # pigeonhole principle, under the partial conjunction null of the entire set of p-values, at least one of the nulls for pval.combine.0 and pval.combine.1 is true
  }
}

# visualization
par(mfrow = c(1,2))
hist(apply(record, 1, min), breaks = 50, xlim = c(0,1))

