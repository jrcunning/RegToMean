set.seed(123)

n <- 100
pre_bleaching <- rlnorm(n, meanlog = log(1e6), sdlog = 1)
#hist(pre_bleaching)

# Linear predictor: base effect + correlation with pre-bleaching + noise
base_effect <- 1.1  # controls the average loss level
correlation_strength <- 0.2  # how much pre values affect loss
noise <- rnorm(n, mean = 0, sd = 0.5)

linear_combination <- base_effect + correlation_strength * scale(log(pre_bleaching)) + noise
#hist(linear_combination)

# Logistic transformation: smooth mapping to (0, 1)
loss_proportion <- 1 / (1 + exp(-linear_combination))
#hist(loss_proportion)

# Post-bleaching values
post_bleaching <- pre_bleaching * (1 - loss_proportion)
#hist(post_bleaching)

# Check the range of losses
summary(loss_proportion)

plot(log(pre_bleaching), log(post_bleaching))
abline(0, 1)

mod <- lm(log(post_bleaching) ~ log(pre_bleaching))
summary_model <- summary(mod)
slope <- summary_model$coefficients[2, "Estimate"]
se_slope <- summary_model$coefficients[2, "Std. Error"]
z <- (slope - 1) / se_slope
p_value <- pnorm(z)  # one-sided test (less than)
p_value
abline(mod, lty = 2)

hist(log(post_bleaching/pre_bleaching))
