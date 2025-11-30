library(tidyverse)

# Acer TLE data
dat2<-read.csv(file="data/Acer3DMorphologyData.csv",stringsAsFactors=TRUE,header=TRUE)

# Log-transform and calculate measurement error variance

# Based on differences in measurements
millsumm <- as.tibble(dat2) %>%
  select(T0_TLE, T0_FieldTLE) %>%
  mutate(across(everything(), log),
         diff = T0_TLE - T0_FieldTLE,
         subj = factor(row_number()))

var(millsumm$diff) / 2   # 0.01301945
# use as numerator to compute Blomqvist k

