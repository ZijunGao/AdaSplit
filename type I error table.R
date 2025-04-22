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
record = readRDS(file.path(directory, paste(setting, ".rds", sep = "")))
alpha = 0.05
typeIError = unlist(lapply(record$pval, function(x){apply(x < alpha, 2, mean)}))
typeIError = matrix(typeIError, ncol = 3) # each column corresponds to a method
rownames(typeIError) = paste("G", seq(1, 5), sep = "")
colnames(typeIError) = c("RT", "RRT", "ART")

typeIError

# write.csv(typeIError, file = file.path(directory, "type I error 005 table.csv"))
