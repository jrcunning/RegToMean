library(tidyverse)


z_test_slope_less_than <- function(model, value) {
  # Extract coefficients
  coefs <- summary(model)$coefficients
  
  # Extract slope estimate and standard error
  slope_est <- coefs[2, "Estimate"]
  slope_se <- coefs[2, "Std. Error"]
  
  # Compute z-statistic and one-sided p-value
  z <- (slope_est - value) / slope_se
  p <- pnorm(z)  # One-sided: P(β < value)
  
  return(p)
}



set.seed(123)

# Parameters
n <- 100
geo_mean_pre <- 1e6
log_sd_pre <- 1

# Simulate pre-bleaching values (lognormal)
df <- tibble(
  coral_id = 1:n,
  pre_bleaching = rlnorm(n, meanlog = log(geo_mean_pre), sdlog = log_sd_pre),
  log_pre = log(pre_bleaching)
)

# Now simulate log(post/pre) ~ Normal(mu, sigma)
mu_logloss <- log(0.25)  # corresponds to 75% geometric mean loss
sd_logloss <- 0.5        # amount of variation in log-loss

# Optional: add a systematic relationship between pre and log-loss
b <- 0.3  # strength of effect
noise <- rnorm(n, mean = 0, sd = sd_logloss)

df <- df %>%
  mutate(log_loss = mu_logloss - b * scale(log_pre)[,1] + noise,  # correlation + noise
         post_bleaching = pre_bleaching * exp(log_loss),
         log_post = log(post_bleaching))

df %>%
  ggplot(aes(x = log_loss)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.8) +
  labs(title = "Distribution of log(post/pre)",
       x = "log(post / pre)", y = "Count") +
  theme_minimal()


# Post vs. Pre
ggplot(df, aes(x = log_pre, y = log_post)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  coord_fixed()
mod <- lm(log_post ~ log_pre, data = df)
z_test_slope_less_than(mod, 1)


# Logloss vs. Pre
ggplot(df, aes(x = log_pre, y = log(post_bleaching/pre_bleaching))) +
  geom_point() +
  coord_fixed()
mod <- lm(log_loss ~ log_pre, data = df)
z_test_slope_less_than(mod, 0)





