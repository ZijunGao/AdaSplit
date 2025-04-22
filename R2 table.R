directory = "~/Desktop/Research/Yao/HTE inference/code/Panning" 

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
             "null more repeats")

setting1 = "default"
setting2 = "no marginalization"
record = readRDS(file.path(directory, "April 2025", paste(setting1, ".rds", sep = "")))
R2 = unlist(lapply(record$R2, function(x){mean(x)}))
R2.sd = unlist(lapply(record$R2, function(x){sd(x)}))
record = readRDS(file.path(directory, "April 2025", paste(setting2, ".rds", sep = "")))
R2 = rbind(R2, unlist(lapply(record$R2, function(x){mean(x)})))[,2,drop=F]
R2.sd = rbind(R2.sd, unlist(lapply(record$R2, function(x){sd(x)})))[,2,drop=F]
rownames(R2.sd) = rownames(R2) = c("BaR", "R")
colnames(R2.sd) = colnames(R2) = setting1

R2.total = R2
R2.total.sd = R2.sd

setting1 = "larger sample size"
setting2 = "no marginalization larger sample size"
record = readRDS(file.path(directory, "April 2025", paste(setting1, ".rds", sep = "")))
R2 = unlist(lapply(record$R2, function(x){mean(x)}))
R2.sd = unlist(lapply(record$R2, function(x){sd(x)}))
record = readRDS(file.path(directory, "April 2025", paste(setting2, ".rds", sep = "")))
R2 = rbind(R2, unlist(lapply(record$R2, function(x){mean(x)})))[,2,drop=F]
R2.sd = rbind(R2.sd, unlist(lapply(record$R2, function(x){sd(x)})))[,2,drop=F]
rownames(R2.sd) = rownames(R2) = c("BaR", "R")
colnames(R2.sd) = colnames(R2) = setting1

R2.total = cbind(R2.total, R2)
R2.total.sd = cbind(R2.total.sd, R2.sd)


setting1 = "larger noise"
setting2 = "no marginalization larger noise"
record = readRDS(file.path(directory, "April 2025", paste(setting1, ".rds", sep = "")))
R2 = unlist(lapply(record$R2, function(x){mean(x)}))
R2.sd = unlist(lapply(record$R2, function(x){sd(x)}))
record = readRDS(file.path(directory, "April 2025", paste(setting2, ".rds", sep = "")))
R2 = rbind(R2, unlist(lapply(record$R2, function(x){mean(x)})))[,2,drop=F]
R2.sd = rbind(R2.sd, unlist(lapply(record$R2, function(x){sd(x)})))[,2,drop=F]
rownames(R2.sd) = rownames(R2) = c("BaR", "R")
colnames(R2.sd) = colnames(R2) = setting1

R2.total = cbind(R2.total, R2)
R2.total.sd = cbind(R2.total.sd, R2.sd)


R2.total
R2.total.sd


# write.csv(R2, file = file.path(directory, "R2 table.csv"))
