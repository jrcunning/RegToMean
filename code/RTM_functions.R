#R functions for regression to the mean with two-timepoint data.



#Function to compute Blomqvist-adjusted estimates of the slope 
#of initial value on change with two-timepoint data. 
#Required input is two numeric vectors of equal length > 1
#that represent measurements of a trait at two timepoints and either
#a repeatability (rep)---an assumed fraction of the measured variance that 
#reflects true variance---or k, which is one minus the repeatability.
#Internally, if rep is provided, it is converted to k.
#(Note: if k varies between timepoints, the estimate of 
#b remains unbiased if k is accurate at timepoint 1. The estimate
#of the heterogeneity variance, though, estimates a function of both
#the heterogeneity variance and the difference in measurement error
#between timepoints.)
#Output includes the naive and Blomqvist-adjusted slope estimates
#and standard errors, as well as an estimate of the
#variance associated with heterogeneity of change in true values.
#If requested, bootstrap standard errors are also estimated 
#along with "pivotal" bootstrap confidence intervals at the
#requested confidence level. Note: all Blomqvist standard errors
#are computed assuming that k at time 1 is known exactly. If this fraction 
#is estimated with error, then the standard errors and CIs 
#reported here are expected to understate uncertainty.
Blomqvist.est <- function(x1, x2, k = NULL, rep = NULL, bootstrap = FALSE, 
                          n_boot = 1000, conf_level = 0.95, seed = NULL) {
  # Check that exactly one  of k or rep is provided
  if (is.null(rep) && is.null(k)) {
    stop("Either 'rep' or 'k' must be provided")
  }
  if (!is.null(rep) && !is.null(k)) {
    stop("Only one of 'rep' or 'k' should be provided, not both")
  }
  # Convert rep to k if rep was provided
  if (!is.null(rep)) {
    k <- 1 - rep
  }
  if (!is.numeric(x1) || !is.numeric(x2)) {
    stop("Both x1 and x2 must be numeric vectors")
  }
  if (length(x1) != length(x2)) {
    stop("x1 and x2 must have the same length")
  }
  if (length(x1) < 3) {
    stop("Vectors must have at least 3 observations for regression")
  }
  if (any(is.na(x1)) || any(is.na(x2))) {
    stop("Missing values (NA) are not allowed")
  }
  if (!is.numeric(k) || length(k) != 1 || k < 0 || k > 1) {
    stop("k must be a number in [0,1]")
  }
  if (bootstrap && (!is.numeric(n_boot) || n_boot < 1 || n_boot != floor(n_boot))) {
    stop("n_boot must be a positive integer")
  }
  if (!is.numeric(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be between 0 and 1")
  }

  # Core estimation function
  estimate_parameters <- function(x1_vec, x2_vec, k) {
    obs.chg <- x2_vec - x1_vec
    fit.naive <- lm(obs.chg ~ x1_vec)
    obs.b <- summary(fit.naive)$coefficients[2, 1]
    se.obs.b <- summary(fit.naive)$coefficients[2, 2]
    var_err <- k*var(x1_vec)
    b <- (obs.b + k) / (1 - k)
    sigma2.het <- var(x2_vec) - var_err - (var(x1_vec) - var_err) * (1 - b)^2
    se.b <- se.obs.b / (1 - k)
    
    return(list(
      estimates = c(`Blomqvist-adjusted slope` = b, 
                    `Naive slope` = obs.b, 
                    `Heterogeneity variance` = sigma2.het),
      closed_form_se = c(`Blomqvist-adjusted slope` = se.b, 
                         `Naive slope` = se.obs.b)
    ))
  }
  
  # Point estimates
  point_est <- estimate_parameters(x1, x2, k)
  estimates <- point_est$estimates
  closed_form_se <- point_est$closed_form_se
  # Bootstrap if requested
  if (bootstrap) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    n <- length(x1)
    boot_estimates <- matrix(NA, nrow = n_boot, ncol = 3)
    colnames(boot_estimates) <- names(estimates)
    
    for (i in 1:n_boot) {
      indices <- sample(1:n, n, replace = TRUE)
      tryCatch({
        boot_result <- estimate_parameters(x1[indices], x2[indices], k)
        boot_estimates[i, ] <- boot_result$estimates
      }, error = function(e) {
        # Skip failed bootstrap samples
      })
    }
  
    # Calculate bootstrap standard errors
    boot_std_errors <- apply(boot_estimates, 2, sd, na.rm = TRUE)
    names(boot_std_errors) <- names(estimates)
    
    # Calculate confidence intervals
    alpha <- 1 - conf_level
    ci_lower <- 2*estimates - apply(boot_estimates, 2, quantile, probs = 1-alpha/2, na.rm = TRUE)
    ci_upper <- 2*estimates - apply(boot_estimates, 2, quantile, probs = alpha/2, na.rm = TRUE)
    
    # Create results object
    results <- list(
      estimates = estimates,
      closed_form_se = closed_form_se,
      bootstrap_se = boot_std_errors,
      conf_intervals = data.frame(
        Parameter = names(estimates),
        Lower = ci_lower,
        Upper = ci_upper,
        row.names = NULL
      ),
      bootstrap_samples = boot_estimates,
      n_boot = n_boot,
      conf_level = conf_level,
      n_obs = n,
      k = k
    )
  } else {
    # Results without bootstrap
    results <- list(
      estimates = estimates,
      closed_form_se = closed_form_se,
      n_obs = length(x1),
      k = k
    )
  }
  
  class(results) <- "blomqvist"
  return(results)
}

# Print method for clean output
print.blomqvist <- function(x, digits = 4, ...) {
  cat("\nBlomqvist-adjusted Parameter Estimates\n")
  cat("======================================\n")
  cat("Number of observations:", x$n_obs, "\n")
  cat("k or one minus repeatability (specified):", round(x$k, digits), "\n\n")
  
  if (!is.null(x$bootstrap_se)) {
    # With bootstrap
    cat("Bootstrap standard errors (", x$n_boot, " replications)\n", sep = "")
    cat("Confidence level:", x$conf_level * 100, "%\n\n")
    
    result_table <- data.frame(
      Parameter = names(x$estimates),
      Estimate = round(x$estimates, digits),
      `Closed-form SE` = c(round(x$closed_form_se, digits), NA),
      `Bootstrap SE` = round(x$bootstrap_se, digits),
      check.names = FALSE
    )
    
    print(result_table, row.names = FALSE)
    
    cat("\nBootstrap Confidence Intervals:\n")
    ci_table <- x$conf_intervals
    ci_table$Lower <- round(ci_table$Lower, digits)
    ci_table$Upper <- round(ci_table$Upper, digits)
    print(ci_table, row.names = FALSE)
  } else {
    # Without bootstrap
    result_table <- data.frame(
      Parameter = names(x$estimates),
      Estimate = round(x$estimates, digits),
      `Std. Error` = c(round(x$closed_form_se, digits), NA),
      check.names = FALSE
    )
    print(result_table, row.names = FALSE)
    
    cat("\nNote: Std. Error not available for Heterogeneity variance (use bootstrap)\n")
  }
  
  cat("\n")
  invisible(x)
}

# Example usage:
# Basic usage
# result <- Blomqvist.est(x1, x2, k = 0.5)
# print(result)

# With bootstrap
# result_boot <- Blomqvist.est(x1, x2, k = 0.5, bootstrap = TRUE, n_boot = 1000)
# print(result_boot)

# Access specific components
# result_boot$estimates
# result_boot$closed_form_se
# result_boot$bootstrap_se
# result_boot$conf_intervals





#Function to compute method-of-momentsestimates of 
#slope, true-score variance, and
#measurement error variance assuming no heterogeneity in true-value 
#change (beyond that which can be explained as a linear function
#of the initial value). Input is two numeric vectors of equal length > 1
#that represent measurements of a trait at two timepoints, and, if desired,
#parameters for bootstrap estimation of the standard error.
#Output includes the estimates and, if desired, bootstrap standard error 
#estimates and "pivotal" bootstrap confidence intervals at the
#requested confidence level.
prepost.moment.est.nohet <- function(x1, x2, bootstrap = FALSE, n_boot = 1000, 
                            conf_level = 0.95, seed = NULL) {
  # Error handling
  if (!is.numeric(x1) || !is.numeric(x2)) {
    stop("Both x1 and x2 must be numeric vectors")
  }
  if (length(x1) != length(x2)) {
    stop("x1 and x2 must have the same length")
  }
  if (length(x1) < 2) {
    stop("Vectors must have at least 2 observations")
  }
  if (any(is.na(x1)) || any(is.na(x2))) {
    stop("Missing values (NA) are not allowed")
  }
  if (bootstrap && (!is.numeric(n_boot) || n_boot < 1 || n_boot != floor(n_boot))) {
    stop("n_boot must be a positive integer")
  }
  if (!is.numeric(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be between 0 and 1")
  }
  # Core estimation function
  estimate_parameters <- function(x1_vec, x2_vec) {
    z1 <- var(x2_vec) - var(x1_vec)
    z2 <- cov(x2_vec - x1_vec, x1_vec) + var(x1_vec)
    w <- z1 / z2
    
    if (cov(x1_vec, x2_vec) > 0) {
      b <- (-(2 - w) + sqrt((2 - w)^2 + 4 * w)) / 2
    } else if (cov(x1_vec, x2_vec) < 0) {
      b <- (-(2 - w) - sqrt((2 - w)^2 + 4 * w)) / 2
    } else {
      b <- -1
    }
    
    if (cov(x1_vec, x2_vec) != 0) {
      sigma2t <- cov(x1_vec, x2_vec) / (1 + b)
      sigma2eps <- var(x1_vec) - sigma2t
    } else {
      sigma2eps <- var(x2_vec)
      sigma2t <- var(x1_vec) - sigma2eps
    }
    
    return(c(Slope = b, `True-value variance` = sigma2t, 
             `Error variance` = sigma2eps))
  }
  # Point estimates
  estimates <- estimate_parameters(x1, x2)
  # Bootstrap if requested
  if (bootstrap) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    n <- length(x1)
    boot_estimates <- matrix(NA, nrow = n_boot, ncol = 3)
    colnames(boot_estimates) <- names(estimates)
    for (i in 1:n_boot) {
      indices <- sample(1:n, n, replace = TRUE)
      boot_estimates[i, ] <- estimate_parameters(x1[indices], x2[indices])
    }
    # Calculate standard errors
    std_errors <- apply(boot_estimates, 2, sd, na.rm = TRUE)
    # Calculate confidence intervals
    alpha <- 1 - conf_level
    ci_lower <- 2*estimates - apply(boot_estimates, 2, quantile, probs = 1-alpha/2, na.rm = TRUE)
    ci_upper <- 2*estimates - apply(boot_estimates, 2, quantile, probs = alpha/2, na.rm = TRUE)
    
    # Create results object
    results <- list(
      estimates = estimates,
      std_errors = std_errors,
      conf_intervals = data.frame(
        Parameter = names(estimates),
        Lower = ci_lower,
        Upper = ci_upper,
        row.names = NULL
      ),
      bootstrap_samples = boot_estimates,
      n_boot = n_boot,
      conf_level = conf_level,
      n_obs = n
    )
  } else {
    # Results without bootstrap
    results <- list(
      estimates = estimates,
      n_obs = length(x1)
    )
  }
  class(results) <- "momentest"
  return(results)
}

# Print method for clean output
print.momentest <- function(x, digits = 4, ...) {
  cat("\nMoment-based Parameter Estimates\n")
  cat("================================\n")
  cat("Number of observations:", x$n_obs, "\n\n")
  
  if (!is.null(x$std_errors)) {
    # With bootstrap
    cat("Bootstrap standard errors (", x$n_boot, " replications)\n", sep = "")
    cat("Confidence level:", x$conf_level * 100, "%\n\n")
    
    result_table <- data.frame(
      Estimate = round(x$estimates, digits),
      `Std. Error` = round(x$std_errors, digits),
      check.names = FALSE
    )
    
    print(result_table)
    
    cat("\nConfidence Intervals:\n")
    ci_table <- x$conf_intervals
    ci_table$Lower <- round(ci_table$Lower, digits)
    ci_table$Upper <- round(ci_table$Upper, digits)
    print(ci_table, row.names = FALSE)
  } else {
    # Without bootstrap
    result_table <- data.frame(
      Parameter = names(x$estimates),
      Estimate = round(x$estimates, digits),
      row.names = NULL
    )
    print(result_table, row.names = FALSE)
  }
  
  cat("\n")
  invisible(x)
}

# Example usage:
# Basic usage
# result <- prepost.moment.est.nohet(x1, x2)
# print(result)

# With bootstrap
# result_boot <- prepost.moment.est.nohet(x1, x2, bootstrap = TRUE, n_boot = 1000, conf_level = 0.95)
# print(result_boot)

# Access specific components
# result_boot$estimates
# result_boot$std_errors
# result_boot$conf_intervals


