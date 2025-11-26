# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Define the functions (pitman.test and rttm.adj)
pitman.test <- function(x1, x2){
  s1 <- sd(x1)
  s2 <- sd(x2)
  r <- cor(x1, x2)
  n <- length(x1)
  stat <- sqrt(n-2)*(s1/s2 - s2/s1)/(2*sqrt(1-r^2))
  2*pt(-abs(stat), n-2)
}

rttm.adj<-function(m1, m2){
  raw.plast<-m2-m1
  vart<-var.test(m1,m2,paired = T) ## variances equal?
  vpv<-vart$p.value # var.test p value
  m1m2cor<-cor.test(m1, m2) # test correlation between m1 and m2
  rho<-m1m2cor$estimate # correlation coefficient between m1 and m2
  m1sd<-sd(m1) # m1 sd
  m2sd<-sd(m2) # m2 sd
  m1v<-var(m1) # m1 var
  m2v<-var(m2) # m2 var
  m1m<-mean(m1) # m1 mean
  m2m<-mean(m2) # m2 mean
  pm<-mean(raw.plast)
  rho2<-(2*rho*m1sd*m2sd)/(m1v+m2v) # adjusted correlation coefficient used if variances are equal
  rhof<-ifelse(vpv <= 0.05, rho, rho2) # which rho is used for dstar calculation is based on variance comparison
  dstar<-(rhof*(m1-m1m)-(m2-m2m))*-1 # adjustment values. Multiply by -1 to flip sign because Kelly and Price based on plasticity as m1-m2, not m2-m1 as in most thermal tolerance estimates
  adj.plast <- pm+dstar # corrected plasticity.
  out<-as.data.frame(cbind(raw.plast, dstar, adj.plast))
  return(out)
}

nsim <- 200
# Simulation function
run_simulation <- function(b_val, var_eps_val, var_het_val) {
  # Fixed parameters
  n <- 200
  a <- 0
  var.t <- 1
  nsim <- nsim
  
  # Add small epsilon to avoid perfect correlations when var_eps = 0
  var_eps_use <- ifelse(var_eps_val == 0, 1e-6, var_eps_val)
  
  # Storage vectors
  obs.bs <- numeric(nsim)
  blom.bs <- numeric(nsim)
  adjregs <- numeric(nsim)
  
  for(i in 1:nsim){
    tryCatch({
      # Generate data
      t1 <- rnorm(n, 0, sqrt(var.t))
      e1 <- rnorm(n, 0, sqrt(var_eps_use))
      edelt <- rnorm(n, 0, sqrt(var_het_val))
      delt <- a + b_val*t1 + edelt
      e2 <- rnorm(n, 0, sqrt(var_eps_use))
      x1 <- t1 + e1
      x2 <- t1 + delt + e2
      obs.chg <- x2 - x1
      
      # Compute estimates with error handling
      lm_fit <- lm(obs.chg ~ x1)
      if(length(coef(lm_fit)) >= 2 && !is.na(coef(lm_fit)[2])) {
        obs.bs[i] <- coef(lm_fit)[2]
      } else {
        obs.bs[i] <- NA
      }
      
      # Blomqvist correction - handle division by zero
      k <- var_eps_val / var(x1)
      if(abs(1 - k) > 1e-10 && !is.na(obs.bs[i])) {
        blom.bs[i] <- (obs.bs[i] + k)/(1 - k)
      } else {
        blom.bs[i] <- obs.bs[i]  # fallback when correction isn't applicable
      }
      
      # Kelly-Price adjusted regression with error handling
      adj <- rttm.adj(x1, x2)
      if(length(adj$dstar) > 0 && !all(is.na(adj$dstar))) {
        fit <- lm(adj$dstar ~ x1)
        if(length(coef(fit)) >= 2 && !is.na(coef(fit)[2])) {
          adjregs[i] <- coef(fit)[2]
        } else {
          adjregs[i] <- NA
        }
      } else {
        adjregs[i] <- NA
      }
    }, error = function(e) {
      # If any error occurs, set values to NA
      obs.bs[i] <<- NA
      blom.bs[i] <<- NA
      adjregs[i] <<- NA
    })
  }
  
  # Return results (remove NAs for mean calculation)
  data.frame(
    b = b_val,
    var_eps = var_eps_val,
    var_het = var_het_val,
    obs_bias = mean(obs.bs, na.rm = TRUE) - b_val,
    blom_bias = mean(blom.bs, na.rm = TRUE) - b_val,
    kp_bias = mean(adjregs, na.rm = TRUE) - b_val,
    obs_sd = sd(obs.bs, na.rm = TRUE),
    blom_sd = sd(blom.bs, na.rm = TRUE),
    kp_sd = sd(adjregs, na.rm = TRUE),
    n_valid_obs = sum(!is.na(obs.bs)),
    n_valid_blom = sum(!is.na(blom.bs)),
    n_valid_kp = sum(!is.na(adjregs))
  )
}

# Define parameter grid
b_values <- c(-.9, -0.5, -0.2, 0, 0.2, 0.5, .9)
var_eps_values <- c(.01, 0.3, 0.6)
var_het_values <- c(0.01, 0.5, 2)

# Create all combinations
param_grid <- expand.grid(
  b = b_values,
  var_eps = var_eps_values,
  var_het = var_het_values
)

# Run simulations for all parameter combinations
cat("Running simulations across parameter grid...\n")
cat("Total combinations:", nrow(param_grid), "\n")

results <- data.frame()
for(i in 1:nrow(param_grid)) {
  cat("Running combination", i, "of", nrow(param_grid), "\n")
  result <- run_simulation(
    param_grid$b[i], 
    param_grid$var_eps[i], 
    param_grid$var_het[i]
  )
  results <- rbind(results, result)
}

# Reshape data for plotting
plot_data <- results %>%
  pivot_longer(cols = c(obs_bias, blom_bias, kp_bias),
               names_to = "method", 
               values_to = "bias") %>%
  # Add corresponding standard deviations
  mutate(
    sd = case_when(
      method == "obs_bias" ~ obs_sd,
      method == "blom_bias" ~ blom_sd,
      method == "kp_bias" ~ kp_sd
    ),
    method = case_when(
      method == "obs_bias" ~ "Observed",
      method == "blom_bias" ~ "Blomqvist",
      method == "kp_bias" ~ "Kelly-Price"
    ),
    var_eps_label = paste("var.eps =", var_eps),
    var_het_label = paste("var.het =", var_het)
  )

# Create the grid of plots
p <- ggplot(plot_data, aes(x = b, y = bias, color = method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = bias - sd, ymax = bias + sd), 
                width = 0.05, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_grid(var_het_label ~ var_eps_label, scales = "free_y") +
  labs(
    title = "Bias in Slope Estimates Across Parameter Combinations",
    subtitle = "Mean Estimate - True Value (b) by Method (?1 SD error bars)",
    x = "True Slope (b)",
    y = "Bias (Mean Estimate - b)",
    color = "Method"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10),
    legend.position = "bottom"
  ) +
  scale_color_manual(values = c("Observed" = "#E31A1C", 
                                "Blomqvist" = "#1F78B4", 
                                "Kelly-Price" = "#33A02C"))

# Display the plot
print(p)

# Print summary statistics
cat("\nSummary of results:\n")
summary_stats <- results %>%
  group_by(var_eps, var_het) %>%
  summarise(
    mean_obs_bias = mean(abs(obs_bias), na.rm = TRUE),
    mean_blom_bias = mean(abs(blom_bias), na.rm = TRUE),
    mean_kp_bias = mean(abs(kp_bias), na.rm = TRUE),
    min_valid_sims = min(c(n_valid_obs, n_valid_blom, n_valid_kp)),
    .groups = "drop"
  )

print(summary_stats)

# Check for any problematic parameter combinations
problematic <- results[results$n_valid_obs < nsim * 0.8 | 
                         results$n_valid_blom < nsim * 0.8 | 
                         results$n_valid_kp < nsim * 0.8, ]

if(nrow(problematic) > 0) {
  cat("\nParameter combinations with <80% valid simulations:\n")
  print(problematic[, c("b", "var_eps", "var_het", "n_valid_obs", "n_valid_blom", "n_valid_kp")])
}