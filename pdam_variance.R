library(tidyverse)

# ------------------------------------------------------------------------
# Calculate measurement error as instrument error based on qPCR replicates
# ------------------------------------------------------------------------

# Import all CT data
ct <- readxl::read_xlsx("pdam_qpcr_data_20250514.xlsx") %>%
  janitor::clean_names() %>%
  select(colony = sample_name, date, target_name, ct) %>%
  mutate(ct = as.numeric(ct),
         colony = factor(colony)) %>%
  filter(date %in% c("0609", "81309", "081309")) %>%
  mutate(date = case_when(date == "081309" ~ "81309", TRUE ~ date))
n_distinct(ct$colony)

# Import list of colonies used in publication (some colonies excluded?)
final <- read_csv("PdamRbleaching.csv") %>%
  mutate(colony = factor(colony),
         sym = factor(sym))
n_distinct(final$colony)

# Get CT data just for colonies used in final study
ct <- ct %>% right_join(final, by = "colony")
n_distinct(ct$colony)

# Adjust Ct values
ct <- ct %>%
  mutate(adj_ct = case_when(
    target_name == "Pd" ~ ct - 4.478,
    target_name == "C"  ~ ct - 2.6789,
    target_name == "D"  ~ ct
  ))

# Calculate SH from ONE technical replicate on each sample
summ <- ct %>%
  select(colony, sym, date, target_name, adj_ct) %>%
  group_by(colony, date, target_name) %>%
  mutate(rep = row_number()) %>%
  group_by(colony, date, target_name, rep) %>%
  pivot_wider(names_from = target_name, values_from = adj_ct) %>%
  mutate(D.Pd = (2^(Pd - D))*2*3*1.21,
         C.Pd = (2^(Pd - C))/9*2*3*1.21,
         SH = sum(D.Pd, C.Pd, na.rm = T)) %>%
  group_by(colony, sym, date, rep) %>%
  summarize(totSH = sum(SH, na.rm = T)) %>%
  mutate(logSH = log(totSH)) %>%
  filter(is.finite(logSH))

# estimate variance of measurement differences on same sample calculated from tech rep 1 vs. tech rep 2
var_diffs_1rep <- summ %>%
  #filter(date == "0609") %>%
  group_by(colony, date) %>%
  summarize(diff = diff(logSH)) %>%
  ungroup() %>%
  summarize(var = var(diff)) %>%
  pull(var)

# Divide this value by 2 because this is the variance of difference between two measurements
# Divide by 2 again because these measurements are based on only 1 technical replicate, whereas the values we actually use are based on the average of two measurements, which halves the variance
v_err <- var_diffs_1rep / 2 / 2
v_err


# ------------------------------------------------------------------------
# Calculate measurement error including short-term fluctuations using
# timeseries where each colony sampled 3x (Feb., Apr., Jun.)
# ------------------------------------------------------------------------

# Load warming dataset
pdam_warm <- read_csv("PdamRwarming.csv") %>%
  mutate(colony = factor(colony),
         logtotal = log(total))
# Get only colonies included in bleaching study
pdam_warm <- pdam_warm %>% filter(colony %in% final$colony)

# Fit colony-specific linear trends over time, allowing each colony to have
# true heterogeneity in change (i.e., it's own linear trend over time), and then
# calculate residual variance
# (Time is numeric: 1 = Feb, 2 = Apr, 3 = Jun)

ggplot(pdam_warm, aes(x = time, y = logtotal)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~colony)

# Fit model
mod <- lm(logtotal ~ time * colony, data = pdam_warm)
res <- augment(mod)
res %>% print(n = nrow(.))

# Get variance of residuals
var(residuals(mod))

#hist(residuals(mod))








