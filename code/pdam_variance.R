library(tidyverse)
library(rptR)

# ------------------------------------------------------------------------------
# Calculate broad measurement error based on repeatability of short-term 
# fluctuations from time-series where each colony sampled 3x (Feb., Apr., Jun.)
# ------------------------------------------------------------------------------

# Load time-series dataset
pdam_warm <- read_csv("data/PdamRwarming.csv") %>%
  mutate(colony = factor(colony),
         time = match(date, c("feb", "apr", "jun")),
         logtotal = log(total)) %>%
  select(colony, sym, date, time, total, logtotal)

# Calculate colony Repeatability (R) and Blomqvist's k (1 - R)
set.seed(1234)
res <- rptGaussian(logtotal ~ time + (time|colony), grname = c("colony", "Fixed", "Residual"), 
                   data = pdam_warm, adjusted = TRUE, ratio = TRUE, nboot = 20000,
                   parallel = TRUE, ncores = 10)
summary(res)
# 1 - Colony repeatability (adjusting for effects of time)
1-summary(res)$boot[[1]]$Median
# Median residual variance across bootstraps
summary(res)$boot[[2]]$Median   # 1 - R = Blomqvist's k = 0.4648801

