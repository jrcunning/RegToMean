library(tidyverse)
set.seed(123)  # for reproducibility

# Simulate pre-bleaching symbiont to host cell ratios (log-normal with geometric mean = 1 million)
n <- 100  # number of corals
geo_mean <- 1e6
log_sd <- 1  # standard deviation on the log scale

# Log-normal distribution: meanlog = log(geo_mean), sdlog = log_sd
pre_bleaching <- rlnorm(n, meanlog = log(geo_mean), sdlog = log_sd)
hist(pre_bleaching)

# Shared random component: noise added to a mean of 75% loss
base_loss <- 0.90
noise_sd <- 0.1
loss_noise_shared <- rnorm(n, mean = 0, sd = noise_sd)

# Scenario 1: Only random variation
loss_proportion_random <- base_loss + loss_noise_shared
# Cap only at a maximum of 100% loss (i.e., min post = 0), allow negative loss (i.e., gains)
loss_proportion_random <- pmin(loss_proportion_random, 0.999)
post_bleaching_random <- pre_bleaching * (1 - loss_proportion_random)

# Scenario 2: Same random variation + systematic effect from pre-bleaching ratio
a <- 0  # strength of systematic effect
log_pre <- log10(pre_bleaching)
systematic_component <- a * (log_pre - mean(log_pre))  # centered around 0

loss_proportion_correlated <- base_loss + loss_noise_shared + systematic_component
loss_proportion_correlated <- pmin(loss_proportion_correlated, 1)
post_bleaching_correlated <- pre_bleaching * (1 - loss_proportion_correlated)

# Combine into a dataframe
df <- tibble(
  coral_id = 1:n,
  pre_bleaching = pre_bleaching,
  logpre = log_pre,
  post_bleaching_random = post_bleaching_random,
  post_bleaching_correlated = post_bleaching_correlated,
  logchangeR = log10(post_bleaching_random/pre_bleaching),
  logchangeC = log10(post_bleaching_correlated/pre_bleaching)
)

df
# Visualize
ggplot(df, aes(x = log10(pre_bleaching), y = log10(post_bleaching_random))) +
  geom_point() +
  geom_abline(aes(slope = 1, intercept = 0))

# Create log ratios
df <- df %>%
  mutate(logpre = log10(pre_bleaching),
         )

ggplot(df, aes(x = logpre, y = logchangeR)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

ggplot(df, aes(x = logpre, y = logchangeC)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)














library(ggplot2)

# Long format for plotting
df_long <- df %>%
  pivot_longer(cols = starts_with("post_bleaching"), 
               names_to = "scenario", 
               values_to = "post_bleaching") %>%
  mutate(scenario = factor(scenario, 
                           levels = c("post_bleaching_random", "post_bleaching_correlated"),
                           labels = c("Random Loss", "Correlated Loss")))

# Plot
ggplot(df_long, aes(x = pre_bleaching, y = post_bleaching)) +
  geom_point(alpha = 0.6) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~scenario) +
  labs(x = "Pre-bleaching Symbiont:Host Ratio",
       y = "Post-bleaching Symbiont:Host Ratio",
       title = "Symbiont Loss Before and After Coral Bleaching") +
  theme_minimal()

