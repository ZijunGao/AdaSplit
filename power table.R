directory = "~/Desktop/Research/Yao/HTE inference/code/Panning/April 2025" 

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

n.non.nulls = 5
power.sd = power = list()
for(setting in settings[c(1, 2, 3)]){
  # setting = settings[1]
  record = readRDS(file.path(directory, paste(setting, ".rds", sep = "")))
  power[[setting]] = unlist(lapply(record$R, function(x){mean(x)/n.non.nulls}))
  power.sd[[setting]] = unlist(lapply(record$R, function(x){sd(x)/n.non.nulls}))
}
power = as.data.frame(power)
power.sd = as.data.frame(power.sd)

power
power.sd

# write.csv(power, file = file.path(directory, "power table.csv"))
# write.csv(power.sd, file = file.path(directory, "power sd table.csv"))