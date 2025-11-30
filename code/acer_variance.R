library(tidyverse)

# Acer TLE data
dat2 <- read.csv(file = "data/Acer3DMorphologyData.csv", 
                 stringsAsFactors = TRUE, header = TRUE)

# Log-transform and calculate measurement error variance

# Based on differences in measurements
millsumm <- as.tibble(dat2) %>%
  select(Genotype, FragID, T0_TLE, T0_FieldTLE) %>%
  unite("FragID", c(Genotype, FragID)) %>%
  mutate(across(c(T0_TLE, T0_FieldTLE), log),
         diff = T0_TLE - T0_FieldTLE)

var(millsumm$diff) / 2   # 0.01301945
# use as numerator to compute Blomqvist k
# denominator mill_var_x_abs is variance of initial measurements from growth study
(var(millsumm$diff) / 2) / mill_var_x_abs      # k = 0.0166


