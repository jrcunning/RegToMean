library(tidyverse)
library(rptR)

# ------------------------------------------------------------------------------
# Calculate broad measurement error based on repeatability of short-term 
# fluctuations from time-series where each colony sampled 3x (Feb., Apr., Jun.)
# ------------------------------------------------------------------------------

# Load time-series dataset
pdam_warm <- read_csv("data/PdamRwarming.csv") %>%
  mutate(colony = factor(colony),
         logtotal = log(total))
# Get only colonies included in bleaching study
pdam_warm <- pdam_warm %>% filter(colony %in% final$colony)

# Calculate Repeatability (R) and k (1 - R)
res <- rptGaussian(logtotal ~ time + (time|colony), grname = c("colony", "Fixed", "Residual"), 
                   data = pdam_warm,
                   adjusted = TRUE, ratio = TRUE)
1 - summary(res)$rpt[[1]]$R     # 1 - R = Blomqvist's k = 0.4972658

# Confirm: Recalculate based on variances
res <- rptGaussian(logtotal ~ time + (time|colony), grname = c("colony", "Fixed", "Residual"), 
                   data = pdam_warm,
                   adjusted = TRUE, ratio = FALSE)

V_between <- summary(res)$rpt[[1]]$R
V_within <- summary(res)$rpt[[2]]$R
R_rpt <- V_between / (V_between + V_within)      # repeatability
k_rpt <- V_within / (V_between + V_within)       # Blomqvist k  (= 1 − R)
k_rpt    # 0.4972658

