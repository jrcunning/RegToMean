library(tidyverse)

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

# Calculate SH from AVERAGED technical replicates (how the analysis is actually done)
summ2 <- ct %>%
  group_by(colony, sym, date, target_name) %>%
  summarize(mean_ct = mean(adj_ct)) %>%
  pivot_wider(names_from = target_name, values_from = mean_ct) %>%
  mutate(D.Pd = (2^(Pd - D))*2*3*1.21,
         C.Pd = (2^(Pd - C))/9*2*3*1.21,
         SH = sum(D.Pd, C.Pd, na.rm = T)) %>%
  mutate(across(c(C.Pd, D.Pd), ~replace_na(., 0)),
         dom = if_else(C.Pd > D.Pd, "C", "D")) %>%
  group_by(colony, sym, date) %>%
  summarize(totSH = sum(SH, na.rm = T)) %>%
  mutate(logSH = log(totSH)) %>%
  filter(is.finite(logSH)) %>%
  ungroup()

# estimate total trait variance v_tot (for initial time point, 0609). (# # have to control for symbiont first...)
summ3 <- summ2 %>%
  ungroup() %>%
  filter(date == "0609") %>%
  group_by(sym) %>%
  mutate(mean_logSH = mean(logSH),
         resid = logSH - mean(logSH)) %>%
  ungroup() %>%
  mutate(adj_logSH = mean(logSH) + resid)

v_tot <- var(summ3$adj_logSH)

# estimate v_pop as v_tot - v_err
v_pop <- v_tot - v_err

# v_err relative to v_pop
v_err_rel <- v_err / v_pop





library(boot)

# Your difference vector from technical replicates
diffs <- summ %>%
  group_by(colony, date) %>%
  summarize(diff = diff(logSH)) %>%
  pull(diff)

# Bootstrap function: calculate v_err from diff vector
boot_fun <- function(data, indices) {
  v_diff <- var(data[indices])
  return(v_diff / 4)  # divide by 4 as you do in your code
}

# Bootstrap 1000 replicates
set.seed(123)
boot_out <- boot(diffs, statistic = boot_fun, R = 1000)

# Get percentile CI
ci <- boot.ci(boot_out, type = "perc")
v_err_lower <- ci$percent[4]
v_err_upper <- ci$percent[5]
v_err_lower_rel <- v_err_lower / (v_tot - v_err_lower)
v_err_upper_rel <- v_err_upper / (v_tot - v_err_upper)
