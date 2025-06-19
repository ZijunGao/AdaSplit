# plot for BaR learner consistency
path = "~/Desktop/Research/Yao/HTE inference/code/Panning"
plotDirectory = file.path("~/Desktop/Research/Yao/HTE inference/code/Panning")

library(ggplot2)
library(dplyr)
library(tidyr)

# Load data
record = readRDS("~/Desktop/Research/Yao/HTE inference/code/Panning/April 2025/BaR_learner_consistency.rds")
sample_sizes <- record$n.seq


# Plot for R2
pdf(file = paste(plotDirectory, "/", "BaRLearnerConsistencyR2", ".pdf", sep = ""), width = 3.8, height = 3.5)
data <- record$R2
colnames(data) <- paste0("n", sample_sizes)
data = as.data.frame(data)

# Convert to long format
data_long <- data %>%
  pivot_longer(cols = everything(), names_to = "sample_size", values_to = "R2") %>%
  mutate(sample_size = as.integer(sub("n", "", sample_size)))

# Compute mean and sd of R2
summary_df <- data_long %>%
  group_by(sample_size) %>%
  summarise(
    mean = mean(R2),
    sd = sd(R2),
    upper = mean + sd,
    lower = mean - sd
  )

# Plot: solid line for mean, dashed for mean ± sd
g = ggplot(summary_df, aes(x = sample_size)) +
  geom_line(aes(y = mean), color = "black", size = 1) +
  geom_line(aes(y = upper), linetype = "dashed", color = "gray40") +
  geom_line(aes(y = lower), linetype = "dashed", color = "gray40") +
  geom_point(aes(y = mean), color = "black") +
  labs(
    x = "Sample Size",
    y = expression(R^2),
    # title = ""
  ) +
  ylim(0, 1) +
  theme_minimal()

print(g)
dev.off()


# plot for beta error
pdf(file = paste(plotDirectory, "/", "BaRLearnerConsistencyBetaError", ".pdf", sep = ""), width = 3.8, height = 3.5)
data <- record$betaError
colnames(data) <- paste0("n", sample_sizes)
data = as.data.frame(data)

# Convert to long format
data_long <- data %>%
  pivot_longer(cols = everything(), names_to = "sample_size", values_to = "betaError") %>%
  mutate(sample_size = as.integer(sub("n", "", sample_size)))

# Compute mean and sd of R2
summary_df <- data_long %>%
  group_by(sample_size) %>%
  summarise(
    mean = mean(betaError),
    sd = sd(betaError),
    upper = mean + sd,
    lower = mean - sd
  )

# Plot: solid line for mean, dashed for mean ± sd
g = ggplot(summary_df, aes(x = sample_size)) +
  geom_line(aes(y = mean), color = "black", size = 1) +
  geom_line(aes(y = upper), linetype = "dashed", color = "gray40") +
  geom_line(aes(y = lower), linetype = "dashed", color = "gray40") +
  geom_point(aes(y = mean), color = "black") +
  labs(
    x = "Sample Size",
    y = expression("||" * hat(beta) - beta * "||"[2]^2),
    # title = ""
  ) +
  # ylim(0, 1) +
  theme_minimal()

print(g)
dev.off()


