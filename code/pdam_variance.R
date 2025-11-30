library(tidyverse)

# ------------------------------------------------------------------------
# Calculate broad measurement error including short-term fluctuations using
# time-series where each colony sampled 3x (Feb., Apr., Jun.)
# ------------------------------------------------------------------------

# Load warming dataset
pdam_warm <- read_csv("data/PdamRwarming.csv") %>%
  mutate(colony = factor(colony),
         logtotal = log(total))
# Get only colonies included in bleaching study
pdam_warm <- pdam_warm %>% filter(colony %in% final$colony)


# Repeatability
library(rptR)
res <- rptGaussian(logtotal ~ time + (time|colony), grname = c("colony", "Fixed", "Residual"), 
                   data = pdam_warm,
                   adjusted = TRUE, ratio = FALSE)
res

V_between <- summary(res)$rpt[[1]]$R

V_within <- summary(res)$rpt[[2]]$R

V_total <- V_between + V_within

R_rpt <- V_between / V_total      # repeatability
k_rpt <- V_within / V_total       # Blomqvist-equivalent broad-k  (~ 1 − R)
k_rpt    # 0.497 -- higher than k from fixed effects model (0.37) because random colony slopes are shrunk toward global mean 
         # and so more of the variance ends up in residuals rather than 
         # soaked up by colony-specific slopes in the fixed effects model

