library(tidyverse)

set.seed(123)

# Parameters
n <- 100
geo_mean_pre <- 1e6
log_sd <- 1

# Logistic model parameters
base_effect <- 1.1         # adjusts average loss (higher = more loss)
correlation_strength <- 0.5  # strength of effect of pre-bleaching values
noise_sd <- 0.5             # how much randomness is in the loss

# Simulate data
df <- tibble(
  coral_id = 1:n,
  pre_bleaching = rlnorm(n, meanlog = log(geo_mean_pre), sdlog = log_sd),
  log_pre = log(pre_bleaching),
  noise = rnorm(n, mean = 0, sd = noise_sd)
) %>%
  mutate(
    linear_combo = base_effect + correlation_strength * scale(log_pre)[,1] + noise,
    loss_proportion = 1 / (1 + exp(-linear_combo)),  # logistic transform
    post_bleaching = pre_bleaching * (1 - loss_proportion)
  )

# Preview
df

ggplot(df, aes(x = pre_bleaching, y = post_bleaching)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Pre-bleaching Symbiont:Host Ratio",
       y = "Post-bleaching Symbiont:Host Ratio",
       title = "Symbiont Loss After Bleaching (Logistic Model)") +
  theme_minimal()

ggplot(df, aes(x = pre_bleaching, y = post_bleaching/pre_bleaching)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "log Pre-bleaching Symbiont:Host Ratio",
       y = "log (Post/Pre)-bleaching Symbiont:Host Ratio",
       title = "Symbiont Loss After Bleaching (Logistic Model)") +
  theme_minimal()

hist(log10(post_bleaching/pre_bleaching))
