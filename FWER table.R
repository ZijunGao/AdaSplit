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

setting = "null more repeats"
record = readRDS(file.path(directory, "April 2025", paste(setting, ".rds", sep = "")))
# Type I error
alpha = 0.1
type.I.error = as.data.frame(lapply(record$pval, function(x){apply(x, 2, function(y){mean(y < alpha)})}))

# FWER
FWER = unlist(lapply(record$FWER, function(x){mean(x)}))
FWER

# write.csv(type.I.error, file = file.path(directory, "type I error table.csv"))
# write.csv(FWER, file = file.path(directory, "FWER table.csv"))



