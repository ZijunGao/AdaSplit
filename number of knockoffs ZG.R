# Use FRT to test subgroup treatment effects
# Number of knockoffs
rm(list = ls())
source("~/Desktop/Research/Yao/HTE inference/code/Panning/helper.R")

n = 400 # sample size
d = 5 # number of covariates
p = 0.5 # propensity score
Group.level.number = c(3, 3) # number of levels per group
Group.number = prod(Group.level.number) # total number of groups
beta0 = 1; beta = rep(1, d) * 1; theta = rep(1, 2)
delta = 1
sigma = 1 # error magnitude in generating Y(0)

test.stats.method = "denoise" # "denoise", "ATE", "denoise + ATE", "AIPW", "ITE"; test statistic
B.seq = seq(1, 10) # number of knockoffs: seq(1, 10); c(1, 4, 10)
M = 400 # number of permutations
q = 0.2 # FDR level

setting = "ITE" # "default", "small sample", "large sample"

if(setting == "large sample"){n = 800}
if(setting == "small sample"){n = 200}
if(setting == "ITE"){test.stats.method = "ITE"} # TODO

start.time = proc.time()
m = 400 # number of trials
record = list()
record$pValue = list()
record$pValue$ORT = matrix(0, nrow = m, ncol = Group.number) 
record$pValue$ART = array(0, dim = c(m, Group.number, length(B.seq))) # ORT: oracle RT; RT (baseline): standard RT; SSRT: sample-splitting RT; DDRT: double-dipping RT; ART: augmented RT; 
record$R = list(); record$R$ORT = list(); record$R$ART = list(); for(j in seq(1, length(B.seq))){record$R$ART[[j]] = list()}
record$FDP = list(); record$FDP$ORT = rep(0, m)  
record$FDP$ART = matrix(0, m, length(B.seq)) 
record$power = record$FDP

set.seed(318)
for(i in 1:m){
  # generate data
  S = cbind(apply(rmultinom(n, 1, rep(1, Group.level.number[1])), 2, which.max),
            apply(rmultinom(n, 1, rep(1, Group.level.number[2])), 2, which.max)) # S used to define subgroups
  Group = (S[, 1] - 1) * Group.level.number[2] + S[, 2]; 
  G = model.matrix(~ factor(Group))[, -1] # one-hot encoding of group membership
  X = matrix(rnorm(n * d, 0, 1), nrow = n, ncol = d) # covariates
  W = rbinom(n, 1, p) # treatment assignment
  # potential outcomes
  # mu0 is linear in X and S
  Y0 = beta0 + X %*% beta + S %*% theta + rnorm(n, 0, sigma)
  # tau is linear in S and independent of X
  tau = delta * (S[, 1] >= Group.level.number[1]) * (S[, 2] >= Group.level.number[2]) #  (S[, 1] >= 3) * (S[, 2] >= 3)
  tau.group = sapply(seq(1, Group.number), function(x) {
    mean(tau[Group == x])
  }) # average treatment effect in each group.
  Y1 = Y0 + tau
  Y = Y1 * W + Y0 * (1 - W) # observed outcome
  # nuisance functions
  mu0 =  beta0 + X %*% beta + S %*% theta
  mu1 =  beta0 + X %*% beta + S %*% theta + tau
  mu = mu0 * (1 - p) + mu1 * p
  
  # inference
  # ORT (oracle): RT with true nuisance functions
  record$pValue$ORT[i,] = ORT(Y = Y, X = X, G = G, Group = Group, prop = p, test.stats.method = test.stats.method, mu0 = mu0, mu1 = mu1, mu = mu, tau = tau, M = M)$pval
  record$R$ORT[[i]] = which(record$pValue$ORT[i,] <= BH.threshold(pval = record$pValue$ORT[i,], q = q))
  record$FDP$ORT[i] = sum(tau.group[record$R$ORT[[i]]] == 0) / max(1, length(record$R$ORT[[i]]))
  record$power$ORT[i] = sum(tau.group[record$R$ORT[[i]]] != 0) / max(1, sum(tau.group != 0))
  
  # ART: augmented RT
  for(j in seq(1, length(B.seq))){
    B = B.seq[j]
    record$pValue$ART[i,,j] = ART(Y = Y, X = X, G = G, Group = Group, prop = p, M = M, test.stats.method = test.stats.method, B = B)$pval
    record$R$ART[[j]][[i]] = which(record$pValue$ART[i,,j] <= BH.threshold(pval = record$pValue$ART[i,,j], q = q))
    record$FDP$ART[i,j] = sum(tau.group[record$R$ART[[j]][[i]]] == 0) / max(1, length(record$R$ART[[j]][[i]]))
    record$power$ART[i,j] = sum(tau.group[record$R$ART[[j]][[i]]] != 0) / max(1, sum(tau.group != 0))
  }
  
  if(i %% 10 == 0){print(i)}
}

# results
result = rbind(c(mean(record$FDP$ORT), apply(record$FDP$ART, 2, mean)),
               c(sd(record$FDP$ORT), apply(record$FDP$ART, 2, sd)),
               c(mean(record$power$ORT), apply(record$power$ART, 2, mean)),        c(sd(record$power$ORT), apply(record$power$ART, 2, sd)))
rownames(result) = c("FDR", "FDR sd", "power", "power sd")
colnames(result) = c("ORT", paste("B = ", B.seq, sep = ""))
print(result)
end.time = proc.time()
print(end.time[3] - start.time[3])

# saveRDS(record, file.path("~/Desktop/Research/Yao/HTE inference/code/Panning/0630", paste("NumberKnockoff", setting, ".rds", sep= "")))


path = "~/Desktop/Research/Yao/HTE inference/code/Panning/0630"
plotDirectory = file.path("~/Desktop/Research/Yao/HTE inference/code/Panning/0630")
setting.seq = c("NumberKnockoffDenoise", "NumberKnockoffITE")
library(ggplot2)


for(setting in setting.seq){
   record = readRDS(file.path(path, paste(setting, "rds", sep = ".")))
  
  plot.data = data.frame(
    B = B.seq,
    FDR = apply(record$FDP$ART, 2, mean),
    power =  apply(record$power$ART, 2, mean)
  )
  # pdf(file = paste(plotDirectory, "/", setting, ".pdf", sep = ""), width = 3.5, height = 14)
  
  ggplot(plot.data, aes(x = B, y = power)) +
    geom_line(color = "blue") +  # Plot the curve
    geom_hline(yintercept = mean(record$power$ORT), linetype = "dashed", color = "red") +  # Add horizontal reference line
    labs(title = "Power ",
         x = "Number of knockoffs",
         y = "Power") +
    theme_minimal()
  
  # dev.off()
}
