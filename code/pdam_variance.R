library(tidyverse)
library(rptR)

# ------------------------------------------------------------------------------
# Calculate broad measurement error based on repeatability of short-term 
# fluctuations from time-series where each colony sampled 3x (Feb., Apr., Jun.)
# ------------------------------------------------------------------------------

# Load time-series dataset
pdam_warm <- read_csv("data/PdamRwarming.csv") %>%
  mutate(colony = factor(colony),
         logtotal = log(total)) %>%
  select(colony, sym, time, date, total, logtotal)
write_csv(pdam_warm, file = "data/pdam_warm.csv")
# # Get only colonies included in bleaching study
pdam_bleach <- read_csv("data/PdamRbleaching.csv") %>%
  mutate(colony = factor(colony)) %>%
  select(colony, sym, jun = juntotal, aug = augtotal) %>%
  pivot_longer(3:4, names_to = "date", values_to = "total") %>%
  mutate(logtotal = log(total))
write_csv(pdam_bleach, file = "data/pdam_bleach.csv")
# pdam_warm <- pdam_warm %>% filter(colony %in% pdam_bleach$colony)

# Calculate Repeatability (R) and k (1 - R)
res <- rptGaussian(logtotal ~ time + (time|colony), grname = c("colony", "Fixed", "Residual"), 
                   data = pdam_warm,
                   adjusted = TRUE, ratio = TRUE)
1 - summary(res)$rpt[[1]]$R     # 1 - R = Blomqvist's k = 0.4972658

# Confirm: Recalculate based on variances
res <- rptGaussian(logtotal ~ time + (time|colony), grname = c("colony", "Fixed", "Residual"), 
                   data = pdam_warm,
                   adjusted = TRUE, ratio = FALSE)
res

V_between <- summary(res)$rpt[[1]]$R
V_within <- summary(res)$rpt[[2]]$R
R_rpt <- V_between / (V_between + V_within)      # repeatability
k_rpt <- V_within / (V_between + V_within)       # Blomqvist k  (= 1 − R)
k_rpt    # 0.4972658


pdam_warm
ggplot(pdam_warm, aes(x = sym, y = logtotal, color = sym)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5) +
  facet_wrap(~ time)

ggplot(pdam_warm, aes(x = time, y = logtotal, group = colony, color = sym)) +
  geom_point() +
  geom_line() +
  facet_wrap(~colony)

# Sym-Normalize?
jungrandmean <- pdam_warm %>%
  filter(date == "jun") %>%
  summarize(jungrandmean = mean(logtotal)) %>%
  pull(jungrandmean)
junsymmeans <- pdam_warm %>%
  filter(date == "jun") %>%
  group_by(sym) %>%
  summarize(mean = mean(logtotal)) %>%
  mutate(adj = mean - jungrandmean)

pdam_warm <- pdam_warm %>%
  left_join(junsymmeans) %>%
  mutate(logtotal_adj = logtotal - adj)

res <- rptGaussian(logtotal_adj ~ time + (time|colony), grname = c("colony", "Fixed", "Residual"), 
                   data = pdam_warm,
                   adjusted = TRUE, ratio = FALSE)
res
1 - summary(res)$rpt[[1]]$R     # 1 - R = Blomqvist's k = 0.4972658






####
warm <- read_csv("data/pdam_warm.csv") %>% 
  mutate(colony = factor(colony),
         log_raw = log(total))
bleach <- read_csv("data/pdam_bleach.csv") %>%
  mutate(colony = factor(colony),
         log_raw = log(total))

grand_mean <- mean(warm$log_raw)

warm_norm <- warm %>% 
  group_by(sym) %>% 
  mutate(sym_resid = log_raw - mean(log_raw),
         log_norm = grand_mean + sym_resid) %>%
  ungroup()
ggplot(warm, aes(x = factor(time), y = log_raw, color = sym, group = colony)) +
  geom_point() +
  geom_line()
ggplot(warm_norm, aes(x = factor(time), y = log_norm, color = sym, group = colony)) +
  geom_point() +
  geom_line()

mod_warm_norm <- lm(log_norm ~ time + colony, data = warm_norm)

sigma2_np <- sum(residuals(mod_warm_norm)^2) / df.residual(mod_warm_norm)
sigma2_np

bleach_norm <- bleach %>%
  group_by(sym) %>%
  mutate(sym_resid = log_raw - mean(log_raw),
         log_init_obs = grand_mean + sym_resid) %>%
  ungroup()

sigma2_init <- var(bleach_norm$log_init_obs)
sigma2_init

pdam_k_consistent <- sigma2_np / sigma2_init
pdam_k_consistent


ggplot(bleach, aes(x = date, y = log_raw, color = sym, group = colony)) +
  geom_point() +
  geom_line()
ggplot(bleach_norm, aes(x = date, y = log_init_obs, color = sym, group = colony)) +
  geom_point() +
  geom_line()

bleach %>%
  pivot_wider(id_cols = c(colony, sym), names_from = date, values_from = log_raw) %>%
  mutate(change = aug - jun) %>%
  ggplot(aes(x = jun, y = change, color = sym)) +
  geom_point()

bleach_norm %>%
  pivot_wider(id_cols = c(colony, sym), names_from = date, values_from = log_init_obs) %>%
  mutate(change = aug - jun) %>%
  ggplot(aes(x = jun, y = change, color = sym)) +
  geom_point()


mod <- bleach %>%
  pivot_wider(id_cols = c(colony, sym), names_from = date, values_from = log_raw) %>%
  mutate(change = aug - jun) %>%
  lm(change ~ jun * sym, data = .)
summary(mod)

