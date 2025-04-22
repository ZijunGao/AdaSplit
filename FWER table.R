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
FWER = unlist(lapply(record$FWER, function(x){mean(x)}))
FWER

# write.csv(FWER, file = file.path(directory, "FWER table.csv"))


