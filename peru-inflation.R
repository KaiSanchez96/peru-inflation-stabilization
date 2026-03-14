############################################################
# PERU INFLATION & MACRO STABILIZATION (1980–2022)
# Regime-Aware Data Analyst Portfolio Project
# Source: World Bank (wbstats)
#
# Project framing:
# This is a descriptive macro analytics case study showing
# how hyperinflation, exchange-rate depreciation, and money
# growth relate to inflation across distinct regimes in Peru.
#
# Key improvement over naive OLS:
# - explicit regime flags
# - influence diagnostics
# - robust / HAC inference
# - multiple model specifications
# - coefficient stability comparisons
############################################################

# ==========================================================
# 0) Libraries
# ==========================================================
suppressPackageStartupMessages({
  library(wbstats)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(broom)
  library(tibble)
  library(purrr)
  library(stringr)
  library(forcats)
  library(sandwich)
  library(lmtest)
  library(car)
  library(ggrepel)
})

# ==========================================================
# 1) Settings
# ==========================================================
COUNTRY  <- "PER"
YEAR_MIN <- 1980
YEAR_MAX <- 2022

OUT_DIR  <- "output"
OUT_FIG  <- file.path(OUT_DIR, "figures")
OUT_TAB  <- file.path(OUT_DIR, "tables")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

INDICATORS <- c(
  "FM.LBL.BMNY.ZG",       # broad money growth
  "NE.CON.GOVT.ZS",       # gov consumption % GDP
  "TT.PRI.MRCH.XD.WD",    # terms of trade index
  "PA.NUS.FCRF",          # official exchange rate
  "FP.CPI.TOTL.ZG"        # CPI inflation
)

STABILIZATION_YEAR <- 1995
HYPER_THRESHOLD    <- 100
HIGH_THRESHOLD     <- 20

# ==========================================================
# 2) Helpers
# ==========================================================
save_plot <- function(plot, filename, w = 8, h = 5, dpi = 300) {
  ggsave(
    filename = file.path(OUT_FIG, filename),
    plot = plot,
    width = w,
    height = h,
    dpi = dpi
  )
}

roll_sd <- function(x, window = 5) {
  out <- rep(NA_real_, length(x))
  for (i in seq_along(x)) {
    if (i >= window) {
      out[i] <- sd(x[(i - window + 1):i], na.rm = TRUE)
    }
  }
  out
}

winsorize <- function(x, p = 0.05) {
  qs <- quantile(x, probs = c(p, 1 - p), na.rm = TRUE)
  pmin(pmax(x, qs[1]), qs[2])
}

get_wb_country_indicators <- function(indicators, country = "PER") {
  raw <- wb_data(
    indicator   = indicators,
    country     = country,
    return_wide = TRUE
  )
  
  id_cols  <- c("iso2c", "iso3c", "country", "date")
  ind_cols <- setdiff(names(raw), id_cols)
  
  long <- raw %>%
    mutate(year = as.integer(date)) %>%
    select(iso3c, country, year, all_of(ind_cols)) %>%
    pivot_longer(
      cols      = all_of(ind_cols),
      names_to  = "indicator_code",
      values_to = "value"
    ) %>%
    arrange(indicator_code, year)
  
  meta <- wb_indicators() %>%
    filter(indicator_id %in% indicators) %>%
    select(indicator_id, indicator) %>%
    rename(
      indicator_code = indicator_id,
      indicator_name = indicator
    )
  
  long %>%
    left_join(meta, by = "indicator_code") %>%
    relocate(indicator_name, .after = indicator_code)
}

extract_model_results <- function(mod, model_name) {
  ols <- tidy(mod) %>%
    mutate(std_error_type = "OLS", model = model_name)
  
  hc1 <- coeftest(mod, vcov = vcovHC(mod, type = "HC1")) %>%
    tidy() %>%
    mutate(std_error_type = "HC1", model = model_name)
  
  nw <- coeftest(mod, vcov = NeweyWest(mod, lag = 1, prewhite = FALSE)) %>%
    tidy() %>%
    mutate(std_error_type = "Newey-West", model = model_name)
  
  bind_rows(ols, hc1, nw)
}

model_fit_stats <- function(mod, model_name, sample_name, outcome_name) {
  res <- mod[["residuals"]]
  fit <- mod[["fitted.values"]]
  smy <- summary(mod)
  
  tibble::tibble(
    model   = model_name,
    sample  = sample_name,
    outcome = outcome_name,
    n       = length(fit),
    r2      = smy$r.squared,
    adj_r2  = smy$adj.r.squared,
    rmse    = sqrt(mean(res^2, na.rm = TRUE)),
    aic     = AIC(mod),
    bic     = BIC(mod)
  )
}

# ==========================================================
# 3) Download data
# ==========================================================
peru_wb <- get_wb_country_indicators(INDICATORS, country = COUNTRY)

coverage <- peru_wb %>%
  group_by(indicator_code, indicator_name) %>%
  summarise(
    min_year = min(year[!is.na(value)]),
    max_year = max(year[!is.na(value)]),
    n_obs    = sum(!is.na(value)),
    .groups  = "drop"
  ) %>%
  arrange(indicator_code)

write.csv(coverage, file.path(OUT_TAB, "indicator_coverage.csv"), row.names = FALSE)

# ==========================================================
# 4) Build panel
# ==========================================================
panel <- peru_wb %>%
  filter(year >= YEAR_MIN, year <= YEAR_MAX) %>%
  select(year, indicator_code, value) %>%
  pivot_wider(names_from = indicator_code, values_from = value) %>%
  arrange(year) %>%
  rename(
    broad_money_growth = FM.LBL.BMNY.ZG,
    gov_cons_gdp       = NE.CON.GOVT.ZS,
    tot_index          = TT.PRI.MRCH.XD.WD,
    exrate_lcu_per_usd = PA.NUS.FCRF,
    inflation_cpi      = FP.CPI.TOTL.ZG
  )

# ==========================================================
# 5) Feature engineering
# ==========================================================
panel_feat <- panel %>%
  arrange(year) %>%
  mutate(
    fx_depr    = 100 * (log(exrate_lcu_per_usd) - log(lag(exrate_lcu_per_usd))),
    tot_growth = 100 * (log(tot_index) - log(lag(tot_index))),
    
    money_growth_l1 = lag(broad_money_growth),
    fx_depr_l1      = lag(fx_depr),
    tot_growth_l1   = lag(tot_growth),
    gov_cons_l1     = lag(gov_cons_gdp),
    
    post_stabilization = year >= STABILIZATION_YEAR,
    regime = case_when(
      inflation_cpi > HYPER_THRESHOLD ~ "Hyperinflation",
      year >= STABILIZATION_YEAR      ~ "Post-stabilization",
      TRUE                            ~ "Pre-stabilization"
    ),
    
    high_infl_20  = inflation_cpi > HIGH_THRESHOLD,
    high_infl_100 = inflation_cpi > HYPER_THRESHOLD,
    
    infl_vol_5y = roll_sd(inflation_cpi, window = 5),
    fx_vol_5y   = roll_sd(fx_depr, window = 5),
    
    log1p_inflation = if_else(inflation_cpi >= 0, log1p(inflation_cpi), NA_real_),
    
    inflation_cpi_w = winsorize(inflation_cpi, p = 0.05),
    fx_depr_w       = winsorize(fx_depr, p = 0.05),
    broad_money_w   = winsorize(broad_money_growth, p = 0.05)
  )

# Defensive checks
stopifnot(all(panel_feat$post_stabilization %in% c(TRUE, FALSE)))
stopifnot(all(
  panel_feat$regime %in% c("Hyperinflation", "Pre-stabilization", "Post-stabilization") |
    is.na(panel_feat$regime)
))
stopifnot(all(panel_feat$inflation_cpi >= 0 | is.na(panel_feat$inflation_cpi)))

write.csv(panel_feat, file.path(OUT_TAB, "panel_feat.csv"), row.names = FALSE)

# ==========================================================
# 6) Data dictionary + summary stats
# ==========================================================
data_dictionary <- tribble(
  ~variable, ~definition, ~unit,
  "inflation_cpi",      "Inflation, consumer prices (annual %)", "percent",
  "exrate_lcu_per_usd", "Official exchange rate (LCU per USD)", "LCU/USD",
  "broad_money_growth", "Broad money growth (annual %)", "percent",
  "tot_index",          "Net barter terms of trade index (2015=100)", "index",
  "gov_cons_gdp",       "Government consumption (% GDP)", "percent",
  "fx_depr",            "Exchange-rate depreciation = 100 x diff(log FX)", "percent",
  "tot_growth",         "Terms-of-trade growth = 100 x diff(log ToT)", "percent",
  "money_growth_l1",    "Lagged broad money growth", "percent",
  "fx_depr_l1",         "Lagged FX depreciation", "percent",
  "tot_growth_l1",      "Lagged ToT growth", "percent",
  "gov_cons_l1",        "Lagged government consumption share", "percent",
  "infl_vol_5y",        "Rolling 5-year SD of inflation", "percent",
  "fx_vol_5y",          "Rolling 5-year SD of FX depreciation", "percent",
  "post_stabilization", "Indicator for 1995 onward", "boolean",
  "regime",             "Hyperinflation / pre-stabilization / post-stabilization", "category",
  "log1p_inflation",    "log(1 + inflation)", "log points",
  "inflation_l1",       "Lagged inflation", "percent"
)

write.csv(data_dictionary, file.path(OUT_TAB, "data_dictionary.csv"), row.names = FALSE)

summary_stats_all <- panel_feat %>%
  summarise(
    n_years      = n(),
    mean_infl    = mean(inflation_cpi, na.rm = TRUE),
    sd_infl      = sd(inflation_cpi, na.rm = TRUE),
    median_infl  = median(inflation_cpi, na.rm = TRUE),
    max_infl     = max(inflation_cpi, na.rm = TRUE),
    mean_money   = mean(broad_money_growth, na.rm = TRUE),
    mean_fx_depr = mean(fx_depr, na.rm = TRUE),
    mean_tot_g   = mean(tot_growth, na.rm = TRUE)
  )

summary_stats_regime <- panel_feat %>%
  group_by(regime) %>%
  summarise(
    n = n(),
    mean_infl = mean(inflation_cpi, na.rm = TRUE),
    median_infl = median(inflation_cpi, na.rm = TRUE),
    max_infl = max(inflation_cpi, na.rm = TRUE),
    mean_money = mean(broad_money_growth, na.rm = TRUE),
    mean_fx_depr = mean(fx_depr, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_stats_all, file.path(OUT_TAB, "summary_stats_all.csv"), row.names = FALSE)
write.csv(summary_stats_regime, file.path(OUT_TAB, "summary_stats_by_regime.csv"), row.names = FALSE)

# ==========================================================
# 7) Descriptive outputs
# ==========================================================

# 7A) Episode table
episodes <- panel_feat %>%
  filter(high_infl_20) %>%
  select(year, regime, inflation_cpi, fx_depr, broad_money_growth, tot_growth, gov_cons_gdp) %>%
  arrange(desc(inflation_cpi))

write.csv(episodes, file.path(OUT_TAB, "inflation_episodes.csv"), row.names = FALSE)

# 7B) SUPPORTING / APPENDIX: time-series dashboard
panel_long <- panel_feat %>%
  select(year, inflation_cpi, broad_money_growth, fx_depr, tot_growth, gov_cons_gdp) %>%
  pivot_longer(-year, names_to = "series", values_to = "value") %>%
  filter(!is.na(value))

p_dashboard <- ggplot(panel_long, aes(x = year, y = value)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ series, scales = "free_y", ncol = 1) +
  labs(
    title = "Peru macro dashboard, 1980–2022",
    x = NULL,
    y = NULL
  ) +
  theme_minimal()

save_plot(p_dashboard, "01_dashboard_timeseries.png", w = 8, h = 11)

# 7C) CORE FIGURE: inflation regime chart
regime_rects <- tibble(
  xmin = c(1988, 1995),
  xmax = c(1991, 2022),
  ymin = c(-Inf, -Inf),
  ymax = c(Inf, Inf),
  phase = c("Hyperinflation episode", "Post-stabilization era")
)

p_regime <- ggplot(panel_feat, aes(x = year, y = inflation_cpi)) +
  geom_rect(
    data = regime_rects %>% filter(phase == "Hyperinflation episode"),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    alpha = 0.18
  ) +
  geom_rect(
    data = regime_rects %>% filter(phase == "Post-stabilization era"),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    alpha = 0.07
  ) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(
    data = panel_feat %>% filter(inflation_cpi > HYPER_THRESHOLD),
    size = 2.5
  ) +
  annotate(
    "text",
    x = 1989.5,
    y = max(panel_feat$inflation_cpi, na.rm = TRUE) * 0.95,
    label = "Hyperinflation\n1988–1991",
    size = 4
  ) +
  annotate(
    "text",
    x = 2008,
    y = max(panel_feat$inflation_cpi, na.rm = TRUE) * 0.95,
    label = "Post-stabilization\n1995–2022",
    size = 4
  ) +
  labs(
    title = "Inflation is dominated by a crisis regime",
    subtitle = "A narrow hyperinflation episode drives the scale of the full 1980–2022 series",
    x = NULL,
    y = "Inflation (%)"
  ) +
  theme_minimal()

save_plot(p_regime, "02_inflation_regime_chart.png", w = 9, h = 5)

# 7D) SUPPORTING / APPENDIX: distribution charts
dist_long <- panel_feat %>%
  select(inflation_cpi, broad_money_growth, fx_depr) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  filter(is.finite(value))

p_dist <- ggplot(dist_long, aes(x = value)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ variable, scales = "free", ncol = 1) +
  labs(
    title = "Raw distributions are highly skewed",
    x = NULL,
    y = "Count"
  ) +
  theme_minimal()

save_plot(p_dist, "03_distribution_check.png", w = 8, h = 9)

# 7E) SUPPORTING / APPENDIX: regime-colored scatterplots
scatter_df <- panel_feat %>%
  filter(complete.cases(inflation_cpi, broad_money_growth, fx_depr, tot_growth, gov_cons_gdp))

p_money_regime <- ggplot(scatter_df, aes(x = broad_money_growth, y = inflation_cpi)) +
  geom_point(aes(shape = regime), size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
  geom_text_repel(
    data = scatter_df %>% filter(inflation_cpi > HYPER_THRESHOLD),
    aes(label = year),
    max.overlaps = 10,
    size = 3
  ) +
  labs(
    title = "Inflation vs broad money growth",
    subtitle = "Relationship is heavily influenced by crisis years",
    x = "Broad money growth (%)",
    y = "Inflation (%)"
  ) +
  theme_minimal()

p_fx_regime <- ggplot(scatter_df, aes(x = fx_depr, y = inflation_cpi)) +
  geom_point(aes(shape = regime), size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
  geom_text_repel(
    data = scatter_df %>% filter(inflation_cpi > HYPER_THRESHOLD),
    aes(label = year),
    max.overlaps = 10,
    size = 3
  ) +
  labs(
    title = "Inflation vs FX depreciation",
    subtitle = "Full-sample slope is dominated by high-inflation years",
    x = "FX depreciation (%)",
    y = "Inflation (%)"
  ) +
  theme_minimal()

save_plot(p_money_regime, "04_scatter_money_regime.png", w = 8, h = 5)
save_plot(p_fx_regime, "05_scatter_fx_regime.png", w = 8, h = 5)

# 7F) SUPPORTING / APPENDIX: rolling volatility
p_vol_infl <- panel_feat %>%
  filter(!is.na(infl_vol_5y)) %>%
  ggplot(aes(x = year, y = infl_vol_5y)) +
  geom_line(linewidth = 0.6) +
  labs(
    title = "Rolling inflation volatility (5-year SD)",
    x = NULL,
    y = "SD of inflation (%)"
  ) +
  theme_minimal()

p_vol_fx <- panel_feat %>%
  filter(!is.na(fx_vol_5y)) %>%
  ggplot(aes(x = year, y = fx_vol_5y)) +
  geom_line(linewidth = 0.6) +
  labs(
    title = "Rolling FX depreciation volatility (5-year SD)",
    x = NULL,
    y = "SD of FX depreciation (%)"
  ) +
  theme_minimal()

save_plot(p_vol_infl, "06_rolling_inflation_volatility.png")
save_plot(p_vol_fx, "07_rolling_fx_volatility.png")

# ==========================================================
# 8) Modeling samples
# ==========================================================
sample_full <- panel_feat %>%
  filter(complete.cases(inflation_cpi, broad_money_growth, fx_depr, tot_growth, gov_cons_gdp))

sample_no_hyper <- sample_full %>%
  filter(inflation_cpi < HYPER_THRESHOLD)

sample_pre <- sample_full %>%
  filter(year < STABILIZATION_YEAR)

sample_post <- sample_full %>%
  filter(year >= STABILIZATION_YEAR)

sample_log <- panel_feat %>%
  filter(complete.cases(log1p_inflation, broad_money_growth, fx_depr, tot_growth, gov_cons_gdp))

sample_lagged <- panel_feat %>%
  arrange(year) %>%
  mutate(
    inflation_l1 = lag(inflation_cpi)
  ) %>%
  filter(complete.cases(
    inflation_cpi, inflation_l1,
    broad_money_growth, fx_depr, tot_growth, gov_cons_gdp
  ))

write.csv(sample_full, file.path(OUT_TAB, "sample_full.csv"), row.names = FALSE)
write.csv(sample_no_hyper, file.path(OUT_TAB, "sample_no_hyper.csv"), row.names = FALSE)
write.csv(sample_lagged, file.path(OUT_TAB, "sample_lagged.csv"), row.names = FALSE)

# ==========================================================
# 9) Models
# ==========================================================
m_full <- lm(
  inflation_cpi ~ broad_money_growth + fx_depr + tot_growth + gov_cons_gdp,
  data = sample_full
)

m_no_hyper <- lm(
  inflation_cpi ~ broad_money_growth + fx_depr + tot_growth + gov_cons_gdp,
  data = sample_no_hyper
)

m_pre <- lm(
  inflation_cpi ~ broad_money_growth + fx_depr + tot_growth + gov_cons_gdp,
  data = sample_pre
)

m_post <- lm(
  inflation_cpi ~ broad_money_growth + fx_depr + tot_growth + gov_cons_gdp,
  data = sample_post
)

m_interact <- lm(
  inflation_cpi ~ broad_money_growth * post_stabilization +
    fx_depr * post_stabilization +
    tot_growth + gov_cons_gdp,
  data = sample_full
)

m_log <- lm(
  log1p_inflation ~ broad_money_growth + fx_depr + tot_growth + gov_cons_gdp,
  data = sample_log
)

m_winsor <- lm(
  inflation_cpi_w ~ broad_money_w + fx_depr_w + tot_growth + gov_cons_gdp,
  data = sample_full
)

m_lagged <- lm(
  inflation_cpi ~ inflation_l1 + broad_money_growth + fx_depr + tot_growth + gov_cons_gdp,
  data = sample_lagged
)

model_list <- list(
  full_levels = m_full,
  no_hyper    = m_no_hyper,
  pre_1995    = m_pre,
  post_1995   = m_post,
  interaction = m_interact,
  log_outcome = m_log,
  winsorized  = m_winsor,
  lagged_aux  = m_lagged
)

# ==========================================================
# 10) Model outputs: coefficients + robust inference
# ==========================================================
all_coefs <- imap_dfr(model_list, extract_model_results)
write.csv(all_coefs, file.path(OUT_TAB, "all_model_coefficients.csv"), row.names = FALSE)

fit_stats <- bind_rows(
  model_fit_stats(m_full, "full_levels", "full", "inflation_cpi"),
  model_fit_stats(m_no_hyper, "no_hyper", "inflation<100", "inflation_cpi"),
  model_fit_stats(m_pre, "pre_1995", "1980-1994", "inflation_cpi"),
  model_fit_stats(m_post, "post_1995", "1995-2022", "inflation_cpi"),
  model_fit_stats(m_interact, "interaction", "full", "inflation_cpi"),
  model_fit_stats(m_log, "log_outcome", "full_nonnegative", "log1p_inflation"),
  model_fit_stats(m_winsor, "winsorized", "full", "inflation_cpi_w"),
  model_fit_stats(m_lagged, "lagged_aux", "lagged_sample", "inflation_cpi")
)

write.csv(fit_stats, file.path(OUT_TAB, "model_fit_stats.csv"), row.names = FALSE)

# ==========================================================
# 11) CORE DIAGNOSTICS: influence diagnostics
# ==========================================================
diag_full <- augment(m_full, data = sample_full) %>%
  mutate(
    cooks_d   = cooks.distance(m_full),
    leverage  = hatvalues(m_full),
    std_resid = rstandard(m_full),
    dffits    = dffits(m_full)
  )

top_influence <- diag_full %>%
  arrange(desc(cooks_d)) %>%
  select(
    year, inflation_cpi, broad_money_growth, fx_depr,
    .fitted, .resid, cooks_d, leverage, std_resid, dffits
  ) %>%
  slice_head(n = 10)

write.csv(diag_full, file.path(OUT_TAB, "m_full_diagnostics.csv"), row.names = FALSE)
write.csv(top_influence, file.path(OUT_TAB, "top_influential_years.csv"), row.names = FALSE)

# CORE FIGURE: Top 10 influential years
p_cooks_top10 <- top_influence %>%
  arrange(cooks_d) %>%
  mutate(year = factor(year, levels = year)) %>%
  ggplot(aes(x = year, y = cooks_d)) +
  geom_col() +
  geom_text(aes(label = round(cooks_d, 2)), hjust = -0.1, size = 3) +
  coord_flip() +
  labs(
    title = "Most influential years in the baseline inflation model",
    subtitle = "1990 dominates, but 1988–1989 also materially affect the full-sample regression",
    x = NULL,
    y = "Cook's distance"
  ) +
  theme_minimal()

save_plot(p_cooks_top10, "08_cooks_distance_top10.png", w = 8, h = 5)

# SUPPORTING / APPENDIX: full sample with log scale
p_cooks_log <- ggplot(diag_full, aes(x = year, y = cooks_d, label = year)) +
  geom_col() +
  scale_y_log10() +
  geom_text_repel(
    data = diag_full %>% arrange(desc(cooks_d)) %>% slice_head(n = 6),
    max.overlaps = 20,
    size = 3
  ) +
  labs(
    title = "Influence diagnostics across all years (log scale)",
    subtitle = "A log scale reveals secondary influential years beneath the 1990 outlier",
    x = NULL,
    y = "Cook's distance (log scale)"
  ) +
  theme_minimal()

save_plot(p_cooks_log, "08b_cooks_distance_logscale.png", w = 9, h = 5)

# SUPPORTING / APPENDIX: leverage vs residuals
p_lev_resid <- ggplot(diag_full, aes(x = leverage, y = std_resid, label = year)) +
  geom_point(size = 2) +
  geom_text_repel(
    data = diag_full %>%
      filter(abs(std_resid) > 2 | cooks_d > quantile(cooks_d, 0.90, na.rm = TRUE)),
    max.overlaps = 20,
    size = 3
  ) +
  labs(
    title = "Leverage vs standardized residuals",
    subtitle = "Crisis years combine unusually high leverage with large residual influence",
    x = "Leverage",
    y = "Standardized residual"
  ) +
  theme_minimal()

save_plot(p_lev_resid, "09_leverage_vs_stdresid.png", w = 8, h = 5)

# ==========================================================
# 12) SUPPORTING / APPENDIX: fitted values and residual checks
# ==========================================================
fitted_full <- sample_full %>%
  mutate(
    fitted_full = predict(m_full, newdata = sample_full),
    resid_full = inflation_cpi - fitted_full
  )

write.csv(fitted_full, file.path(OUT_TAB, "fitted_full.csv"), row.names = FALSE)

p_fit_time <- ggplot(fitted_full, aes(x = year)) +
  geom_line(aes(y = inflation_cpi, linetype = "Actual"), linewidth = 0.7) +
  geom_line(aes(y = fitted_full, linetype = "Fitted"), linewidth = 0.7) +
  labs(
    title = "Actual vs fitted inflation: full-sample baseline",
    subtitle = "The model tracks crisis spikes better than calmer years",
    x = NULL,
    y = "Inflation (%)",
    linetype = NULL
  ) +
  theme_minimal()

save_plot(p_fit_time, "10_actual_vs_fitted_full.png", w = 9, h = 5)

p_resid_time <- ggplot(fitted_full, aes(x = year, y = resid_full)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_col() +
  labs(
    title = "Residuals over time",
    subtitle = "Misspecification is concentrated around regime shifts",
    x = NULL,
    y = "Residual"
  ) +
  theme_minimal()

save_plot(p_resid_time, "11_residuals_over_time.png", w = 9, h = 5)

# ==========================================================
# 12B) SUPPORTING / APPENDIX: autocorrelation checks
# ==========================================================
dw_full <- car::durbinWatsonTest(m_full)
bg_full <- lmtest::bgtest(m_full, order = 1)

dw_post <- car::durbinWatsonTest(m_post)
bg_post <- lmtest::bgtest(m_post, order = 1)

autocorr_checks <- tibble::tibble(
  model = c("full_levels", "full_levels", "post_1995", "post_1995"),
  test  = c("Durbin-Watson", "Breusch-Godfrey(1)", "Durbin-Watson", "Breusch-Godfrey(1)"),
  statistic = c(
    unname(dw_full$dw),
    unname(bg_full$statistic),
    unname(dw_post$dw),
    unname(bg_post$statistic)
  ),
  p_value = c(
    unname(dw_full$p),
    unname(bg_full$p.value),
    unname(dw_post$p),
    unname(bg_post$p.value)
  )
)

write.csv(
  autocorr_checks,
  file.path(OUT_TAB, "autocorrelation_checks.csv"),
  row.names = FALSE
)

acf_obj_full <- acf(residuals(m_full), plot = FALSE)

acf_df_full <- tibble::tibble(
  lag = as.numeric(acf_obj_full$lag),
  acf = as.numeric(acf_obj_full$acf)
) %>%
  filter(lag > 0)

p_acf_full <- ggplot(acf_df_full, aes(x = lag, y = acf)) +
  geom_col() +
  labs(
    title = "Residual autocorrelation: full-sample baseline",
    x = "Lag",
    y = "ACF"
  ) +
  theme_minimal()

save_plot(p_acf_full, "13_residual_acf_full.png", w = 8, h = 4)

# ==========================================================
# 13) CORE FIGURE: coefficient stability comparison
# ==========================================================
coef_compare <- all_coefs %>%
  filter(
    std_error_type == "Newey-West",
    model %in% c("full_levels", "no_hyper", "pre_1995", "post_1995"),
    term %in% c("broad_money_growth", "fx_depr", "tot_growth", "gov_cons_gdp")
  ) %>%
  mutate(
    conf_low = estimate - 1.96 * std.error,
    conf_high = estimate + 1.96 * std.error,
    term = case_when(
      term == "broad_money_growth" ~ "Broad money growth",
      term == "fx_depr" ~ "FX depreciation",
      term == "tot_growth" ~ "Terms-of-trade growth",
      term == "gov_cons_gdp" ~ "Gov. consumption share",
      TRUE ~ term
    ),
    model = factor(
      model,
      levels = c("full_levels", "no_hyper", "pre_1995", "post_1995"),
      labels = c("Full sample", "No hyperinflation", "1980–1994", "1995–2022")
    )
  )

write.csv(
  coef_compare,
  file.path(OUT_TAB, "coefficient_stability_comparison.csv"),
  row.names = FALSE
)

p_coef_compare <- ggplot(coef_compare, aes(x = estimate, y = model)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point() +
  geom_errorbar(
    aes(xmin = conf_low, xmax = conf_high),
    width = 0.15,
    orientation = "y"
  ) +
  facet_wrap(~ term, scales = "free_x") +
  labs(
    title = "Coefficient stability across specifications",
    subtitle = "Point estimates shift materially once crisis years are handled differently",
    x = "Estimated coefficient (Newey-West SEs)",
    y = NULL
  ) +
  theme_minimal()

save_plot(p_coef_compare, "12_coefficient_stability.png", w = 10, h = 6)

# ==========================================================
# 14) Interaction model interpretation table
# ==========================================================
interaction_terms <- extract_model_results(m_interact, "interaction") %>%
  filter(std_error_type == "Newey-West") %>%
  mutate(
    conf_low = estimate - 1.96 * std.error,
    conf_high = estimate + 1.96 * std.error
  )

write.csv(
  interaction_terms,
  file.path(OUT_TAB, "interaction_model_neweywest.csv"),
  row.names = FALSE
)

# ==========================================================
# 14B) Supporting table: log-outcome model
# ==========================================================
log_terms <- extract_model_results(m_log, "log_outcome") %>%
  filter(std_error_type == "Newey-West") %>%
  mutate(
    conf_low = estimate - 1.96 * std.error,
    conf_high = estimate + 1.96 * std.error
  )

write.csv(
  log_terms,
  file.path(OUT_TAB, "log_outcome_model_neweywest.csv"),
  row.names = FALSE
)

# SUPPORTING / APPENDIX: coefficient plot for log-outcome model
log_coef_plot <- log_terms %>%
  filter(term %in% c("broad_money_growth", "fx_depr", "tot_growth", "gov_cons_gdp")) %>%
  mutate(
    term = case_when(
      term == "broad_money_growth" ~ "Broad money growth",
      term == "fx_depr" ~ "FX depreciation",
      term == "tot_growth" ~ "Terms-of-trade growth",
      term == "gov_cons_gdp" ~ "Gov. consumption share",
      TRUE ~ term
    ),
    term = factor(term, levels = rev(c(
      "Broad money growth",
      "FX depreciation",
      "Terms-of-trade growth",
      "Gov. consumption share"
    )))
  )

p_log_coef <- ggplot(log_coef_plot, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point() +
  geom_errorbar(
    aes(xmin = conf_low, xmax = conf_high),
    width = 0.15,
    orientation = "y"
  ) +
  labs(
    title = "Log-outcome model coefficients",
    subtitle = "Newey-West confidence intervals for the transformed inflation specification",
    x = "Estimated coefficient",
    y = NULL
  ) +
  theme_minimal()

save_plot(p_log_coef, "14_log_outcome_coefficients.png", w = 8, h = 5)

# ==========================================================
# 15) DEMOTED / LEGACY: correlation matrix
# ==========================================================
corr_df <- panel_feat %>%
  select(inflation_cpi, broad_money_growth, fx_depr, tot_growth, gov_cons_gdp) %>%
  na.omit()

corr_mat <- cor(corr_df)
write.csv(corr_mat, file.path(OUT_TAB, "correlation_matrix.csv"))

# ==========================================================
# 16) Portfolio-ready findings summary
# ==========================================================
portfolio_summary <- tibble(
  finding_order = 1:6,
  finding = c(
    "Inflation dynamics are dominated by a small hyperinflation regime rather than a stable full-sample relationship.",
    "Naive full-sample OLS is highly sensitive to influential crisis years.",
    "Coefficient sizes change materially when hyperinflation years are excluded or when the sample is split pre/post stabilization.",
    "Exchange-rate depreciation and money growth appear strongly associated with inflation in crisis years, but much less uniformly afterward.",
    "Robust and HAC inference is more appropriate than plain OLS standard errors in this annual macro setting.",
    "The project is descriptive and regime-aware, not causal."
  )
)

write.csv(
  portfolio_summary,
  file.path(OUT_TAB, "portfolio_summary_findings.csv"),
  row.names = FALSE
)

# ==========================================================
# 17) Project metadata / README notes
# ==========================================================
readme_notes <- tibble(
  section = c(
    "project_title",
    "project_type",
    "analytical_goal",
    "key_caution",
    "best_figures",
    "best_tables"
  ),
  content = c(
    "Peru Inflation and Macroeconomic Stabilization (1980–2022)",
    "Descriptive macro analytics case study",
    "Show how regime change and outliers alter inflation modeling conclusions",
    "Results are not causal and are sensitive to crisis-era leverage points",
    "02_inflation_regime_chart.png; 08_cooks_distance_top10.png; 12_coefficient_stability.png",
    "top_influential_years.csv; model_fit_stats.csv; coefficient_stability_comparison.csv"
  )
)

write.csv(readme_notes, file.path(OUT_TAB, "project_metadata.csv"), row.names = FALSE)

# ==========================================================
# 18) SUPPORTING / APPENDIX: simple lagged annual model
# ==========================================================
lagged_results <- extract_model_results(m_lagged, "lagged_aux")
write.csv(lagged_results, file.path(OUT_TAB, "lagged_model_results.csv"), row.names = FALSE)

lagged_fit <- model_fit_stats(
  m_lagged, "lagged_aux", "lagged_sample", "inflation_cpi"
)
write.csv(lagged_fit, file.path(OUT_TAB, "lagged_model_fit.csv"), row.names = FALSE)

############################################################
# END
############################################################