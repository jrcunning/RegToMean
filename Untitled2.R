set.seed(123)

n <- 100
geo_mean_pre <- 1e6
log_sd <- 1
mu_pre <- log(geo_mean_pre)

# Step 1: Simulate pre-bleaching ratios
pre_bleaching <- rlnorm(n, meanlog = mu_pre, sdlog = log_sd)
log_pre <- log(pre_bleaching)
hist(log_pre)

# Step 2: Scenario 1 - Random loss (log-scale shift of mean by log(0.25))
mu_post <- mu_pre + log(0.25)  # shift down by 75%
post_bleaching_random <- rlnorm(n, meanlog = mu_post, sdlog = log_sd)
hist(post_bleaching_random)
plot(log(pre_bleaching), log(post_bleaching_random))
abline(0, 1)

# Step 3: Scenario 2 - Same random variation as Scenario 1 + systematic loss based on initial values

# Extract the log-scale noise from scenario 1
log_noise <- log(post_bleaching_random) - mu_post  # mean zero noise

# Add systematic loss component to log values based on initial log ratio
b <- 0.2  # strength of relationship between initial log and loss
log_post_bleaching_correlated <- mu_post + log_noise - b * (log_pre - mean(log_pre))

post_bleaching_correlated <- exp(log_post_bleaching_correlated)
plot(log(pre_bleaching), log(post_bleaching_correlated))
abline(0, 1)

# Combine into a dataframe
df <- data.frame(
  coral_id = 1:n,
  pre_bleaching = pre_bleaching,
  post_bleaching_random = post_bleaching_random,
  post_bleaching_correlated = post_bleaching_correlated
)

# Optional: visualize
library(tidyr)
library(ggplot2)

df_long <- df %>%
  pivot_longer(cols = starts_with("post_bleaching"),
               names_to = "scenario",
               values_to = "post_bleaching") %>%
  mutate(scenario = factor(scenario,
                           levels = c("post_bleaching_random", "post_bleaching_correlated"),
                           labels = c("Random Loss", "Correlated Loss")))

ggplot(df_long, aes(x = pre_bleaching, y = post_bleaching/pre_bleaching)) +
  geom_point(alpha = 0.6) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~scenario) +
  labs(x = "Pre-bleaching Symbiont:Host Ratio",
       y = "log (Post/Pre)-bleaching Symbiont:Host Ratio",
       title = "Symbiont:Host Ratio Before and After Bleaching") +
  theme_minimal()

ggplot(df_long, aes(x = pre_bleaching, y = post_bleaching)) +
  geom_point(alpha = 0.6) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~scenario) +
  labs(x = "Pre-bleaching Symbiont:Host Ratio",
       y = "Post-bleaching Symbiont:Host Ratio",
       title = "Symbiont:Host Ratio Before and After Bleaching") +
  theme_minimal()



df <- tibble(df)
df  
modR <- lm(post_bleaching_random ~ pre_bleaching, data = df)
summary(modR)
