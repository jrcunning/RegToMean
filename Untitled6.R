# Load tidyverse for plotting and clean data handling
library(tidyverse)

set.seed(123)

# Sample size
n <- 1000

# STEP 1: Simulate true pre and post values (X and Y)
# Let's simulate on the log scale for realism in biological ratios (e.g., SHR)

log_X <- rnorm(n, mean = -2.2, sd = 0.63)  # true pre values (log scale)
D <- rnorm(n, mean = log(0.25), sd = 0.3)  # true change (log scale), e.g., 75% reduction
log_Y <- log_X + D                         # true post = pre + change on log scale

# STEP 2: Add measurement error
error_X <- rnorm(n, mean = 0, sd = 0.4)  # error in baseline
error_Y <- rnorm(n, mean = 0, sd = 0.4)  # error in follow-up

log_x <- log_X + error_X  # observed pre
log_y <- log_Y + error_Y  # observed post

# Observed change (log scale): log(y) - log(x) = log(y/x)
log_loss <- log_y - log_x

# STEP 3: Fit models
mod1 <- lm(log_y ~ log_x)         # post vs pre
mod2 <- lm(log_loss ~ log_x)      # log ratio vs pre

# STEP 4: Report coefficients
summary(mod1)  # slope should be < 1 due to regression to the mean
summary(mod2)  # slope should be < 0 due to mathematical coupling

# STEP 5: Optional plot for visual diagnostics
tibble(log_x, log_y, log_loss) %>%
  pivot_longer(cols = c(log_y, log_loss), names_to = "response", values_to = "value") %>%
  ggplot(aes(x = log_x, y = value)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  facet_wrap(~response, scales = "free_y") +
  labs(
    x = "Observed log(pre)",
    y = "Response",
    title = "Regression to the Mean and Mathematical Coupling Due to Measurement Error"
  ) +
  theme_minimal()

