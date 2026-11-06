#### RACCOON SEASONAL ANALYSIS ####

library(tidyverse)
library(lubridate)
library(mgcv)
library(patchwork)


#### SETUP ####

out <- "revision_equalarea"

for (x in c("data", "models", "tables", "figures", "diagnostics")) {
  dir.create(file.path(out, x), recursive = TRUE, showWarnings = FALSE)
}

n_threads <- max(1L, parallel::detectCores(logical = FALSE) - 1L, na.rm = TRUE)

n_sim <- 2000L

#### DATA PREP ####

#saved timestamps and equal-area grid IDs
d_diel <- readRDS(file.path(out, "data", "raw_obs_diel.rds"))

d_temporal <- d_diel %>%
  mutate(id = as.character(id),
         grid_factor = factor(grid_id),
         year_factor = factor(year(observed_date)),
         month_num = month(observed_date),
         latitude_band = floor(latitude / 5) * 5)

#existing latitude restriction; retain cells with fewer than three records
d_season <- d_temporal %>%
  filter(latitude >= 25, latitude < 50) %>%
  droplevels()

saveRDS(d_season, file.path(out, "data", "d_season_pre_temperature.rds"))

#sample accounting
season_sample <- d_season %>%
  summarise(dataset = "Seasonal analysis",
            n = n(),
            n_arboreal = sum(arboreal),
            n_cells = n_distinct(grid_id),
            n_years = n_distinct(year_factor),
            percent_diel_retained = 100 * n() / nrow(d_diel),
            min_latitude = min(latitude),
            max_latitude = max(latitude),
            first_date = min(observed_date),
            last_date = max(observed_date))

season_sample

write_csv(season_sample, file.path(out, "tables", "Table_S1_seasonal_sample.csv"))


#### SAMPLING COVERAGE ####

#include empty month-band combinations
seasonal_coverage <- d_season %>%
  count(latitude_band, month_num, name = "n") %>%
  complete(latitude_band = seq(25, 45, 5), month_num = 1:12,
           fill = list(n = 0L))

coverage_summary <- seasonal_coverage %>%
  group_by(latitude_band) %>%
  summarise(n = sum(n),
            months_present = sum(n > 0),
            min_month_n = min(n),
            max_month_n = max(n),
            .groups = "drop")

coverage_summary

write_csv(seasonal_coverage, file.path(out, "tables", "seasonal_sampling_coverage.csv"))
write_csv(coverage_summary, file.path(out, "tables", "Table_S1_seasonal_coverage.csv"))

year_counts <- d_season %>% count(year = year(observed_date), name = "n")

write_csv(year_counts, file.path(out, "tables", "seasonal_year_counts.csv"))

p_coverage <- ggplot(seasonal_coverage, aes(month_num, latitude_band, fill = n)) +
  geom_tile() +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_y_continuous(breaks = seq(25, 45, 5),
                     labels = c("25–<30", "30–<35", "35–<40",
                                "40–<45", "45–<50")) +
  scale_fill_viridis_c(trans = "sqrt", labels = scales::comma,
                       guide = guide_colorbar(barwidth = unit(7, "cm"))) +
  labs(x = "month", y = "latitude band (°N)", fill = "observations") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

p_year <- ggplot(year_counts, aes(year, n)) +
  geom_col(fill = "gray50") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "year", y = "observations") +
  theme_classic(base_size = 11)

p_sampling <- (p_year / p_coverage) + plot_annotation(tag_levels = "a")

ggsave(file.path(out, "figures", "Figure_S2_seasonal_coverage.jpg"),
       p_sampling, width = 5, height = 7, dpi = 350, bg = "white")


#### MODELS ####

#model formulas
f_season_base <- arboreal ~
  s(solar_hour, bs = "cc", k = 20) +
  s(annual_phase, bs = "cc", k = 12) +
  s(latitude, longitude, bs = "sos", m = 2, k = 200) +
  s(grid_factor, bs = "re") +
  s(year_factor, bs = "re")

f_season_lat <- update(f_season_base, 
                       . ~ . + ti(annual_phase, latitude, bs = c("cc", "tp"), k = c(8, 5)))

#fitting settings
fit_season <- function(formula, method = "fREML") {
  set.seed(1120)
  
  mgcv::bam(formula,
            family = binomial(link = "logit"),
            data = d_season,
            method = method,
            discrete = FALSE,
            nthreads = n_threads,
            knots = list(solar_hour = c(0, 24), annual_phase = c(0, 1)),
            na.action = na.fail)
}

m_season_base <- fit_season(f_season_base)
saveRDS(m_season_base, file.path(out, "models", "m_season_base.rds"))

m_season_lat <- fit_season(f_season_lat)
saveRDS(m_season_lat, file.path(out, "models", "m_season_lat.rds"))

summary(m_season_base)
summary(m_season_lat)

#### MODEL OUTPUTS ####

#export helper function
export_gam <- function(fit, label) {
  sm <- summary(fit)
  
  write_csv(as.data.frame(sm$s.table) %>% rownames_to_column("term"),
            file.path(out, "tables", paste0("Table_S4_", label, "_smooths.csv")))
  
  write_csv(as.data.frame(sm$p.table) %>% rownames_to_column("term"),
            file.path(out, "tables", paste0("Table_S4_", label, "_parametric.csv")))
  
  write_csv(tibble(model = label, n = nobs(fit),
                   deviance_explained = sm$dev.expl,
                   converged = isTRUE(fit$converged)),
            file.path(out, "tables", paste0("Table_S4_", label, "_fit.csv")))
  
  set.seed(1120)
  k_result <- mgcv::k.check(fit, subsample = 10000, n.rep = 1000)
  
  write_csv(as.data.frame(k_result) %>% rownames_to_column("term"),
            file.path(out, "tables", paste0("Table_S4_", label, "_kcheck.csv")))
  
  conc <- mgcv::concurvity(fit, full = TRUE)
  
  write_csv(as.data.frame(conc) %>% rownames_to_column("measure"),
            file.path(out, "tables", paste0("Table_S4_", label, "_concurvity.csv")))
  
  saveRDS(mgcv::concurvity(fit, full = FALSE),
          file.path(out, "diagnostics", paste0(label, "_pairwise_concurvity.rds")))
  
  capture.output(summary(fit), print(k_result), print(round(conc, 3)),
                 file = file.path(out, "diagnostics", paste0(label, "_summary.txt")))
}

export_gam(m_season_base, "season_base")
export_gam(m_season_lat, "season_lat")

#save helper for temp analysis
dump("export_gam", file = file.path(out, "data", "export_gam.R"))


#### RESIDUAL CHECKS ####

#descriptive checks from final seasonal model
residual_data <- d_season %>%
  select(id, grid_id, arboreal, solar_hour, annual_phase, latitude, longitude) %>%
  mutate(fitted = as.numeric(fitted(m_season_lat)),
         residual = as.numeric(residuals(m_season_lat, type = "pearson")))

write_csv(residual_data, file.path(out, "diagnostics", "seasonal_residuals.csv"))

calibration <- residual_data %>%
  mutate(bin = ntile(fitted, 20)) %>%
  group_by(bin) %>%
  summarise(fitted = mean(fitted), observed = mean(arboreal),
            n = n(), .groups = "drop")

residual_bins <- bind_rows(
  residual_data %>% mutate(variable = "solar hour", bin = floor(solar_hour)),
  residual_data %>% mutate(variable = "annual phase", bin = floor(annual_phase * 12) + 1),
  residual_data %>% mutate(variable = "latitude", bin = floor(latitude))) %>%
  group_by(variable, bin) %>%
  summarise(residual = mean(residual), n = n(), .groups = "drop")

p_calibration <- ggplot(calibration, aes(fitted, observed)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "gray50") +
  geom_line() +
  geom_point() +
  labs(x = "mean fitted probability", y = "observed proportion") +
  theme_classic(base_size = 11)

p_residuals <- ggplot(residual_bins, aes(bin, residual)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
  geom_line() +
  geom_point() +
  facet_wrap(~variable, scales = "free_x", nrow = 1) +
  labs(x = "predictor bin", y = "mean Pearson residual") +
  theme_classic(base_size = 11)

ggsave(file.path(out, "figures", "Figure_S3_seasonal_diagnostics.jpg"),
       p_calibration / p_residuals,
       width = 11, height = 7, dpi = 350, bg = "white")


#### PARTIAL-EFFECT PREDICTIONS ####

#common reference values for prediction frames
reference_data <- function(nd) {
  defaults <- d_season[1, c("solar_hour", "annual_phase", "latitude",
                            "longitude", "grid_factor", "year_factor")]
  
  defaults$solar_hour <- 12
  defaults$annual_phase <- 0.5
  defaults$latitude <- median(d_season$latitude)
  defaults$longitude <- median(d_season$longitude)
  
  for (v in setdiff(names(defaults), names(nd))) {
    nd[[v]] <- defaults[[v]][rep(1, nrow(nd))]
  }
  
  nd
}

#retain intercept and selected smooth coefficients only
partial_X <- function(nd, terms) {
  X <- predict(m_season_lat, newdata = nd, type = "lpmatrix")
  keep <- which(colnames(X) == "(Intercept)")
  
  for (sm in m_season_lat$smooth) {
    if (sm$label %in% terms) {
      keep <- c(keep, seq.int(sm$first.para, sm$last.para))
    }
  }
  
  X[, unique(keep), drop = FALSE]
}

#uncertainty for curves & simulations
partial_predict <- function(nd, terms) {
  X <- partial_X(nd, terms)
  i <- match(colnames(X), names(coef(m_season_lat)))
  V <- m_season_lat$Vp[i, i, drop = FALSE]
  
  eta <- drop(X %*% coef(m_season_lat)[i])
  se <- sqrt(pmax(0, rowSums((X %*% V) * X)))
  
  nd %>%
    mutate(fit = plogis(eta),
           lower = plogis(eta - 1.96 * se),
           upper = plogis(eta + 1.96 * se))
}

season_terms <- c("s(annual_phase)", "ti(annual_phase,latitude)")

lat_vals <- as.numeric(quantile(d_season$latitude, c(0.1, 0.5, 0.9)))
lat_labels <- paste0(round(lat_vals, 1), "°N")

plot_dates <- seq.Date(as.Date("2021-01-01"), as.Date("2021-12-31"), by = "day")

#three reference latitudes
season_pred <- crossing(plot_date = plot_dates, latitude = lat_vals) %>%
  mutate(annual_phase = (yday(plot_date) - 1) / 365,
         latitude_plot = factor(latitude, levels = lat_vals, labels = lat_labels)) %>%
  reference_data() %>%
  partial_predict(season_terms)

#diel curve from the same final model
diel_pred <- tibble(solar_hour = seq(0, 24, length.out = 241)) %>%
  reference_data() %>%
  partial_predict("s(solar_hour)")

write_csv(season_pred, file.path(out, "tables", "Figure_3_seasonal_curves.csv"))
write_csv(diel_pred, file.path(out, "tables", "Figure_3_diel_curve.csv"))


#### MAIN FIGURE ####

p_diel <- ggplot(diel_pred, aes(solar_hour, fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24),
                     labels = c("solar\nmidnight", "06:00", "solar\nnoon",
                                "18:00", "solar\nmidnight"),
    expand = expansion(mult = 0)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "solar hour", y = "diel partial-effect probability") +
  theme_classic(base_size = 11) +
  theme(plot.margin = margin(6, 30, 6, 6))

p_season <- ggplot(season_pred, 
                   aes(x = plot_date, y = fit, 
                       color = as.factor(latitude_plot), 
                       fill = as.factor(latitude_plot), 
                       group = as.factor(latitude_plot))) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = 0) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c("#0072B2", "#009E73", "#E69F00"), name = "latitude (°)") +
  scale_fill_manual(values = c("#0072B2", "#009E73", "#E69F00"), name = "latitude (°)") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b", expand = expansion(mult = 0)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "month", y = "seasonal partial-effect probability") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

p_fig3 <- (p_diel / p_season) + plot_annotation(tag_levels = "a")

p_fig3

ggsave(file.path(out, "figures", "Figure_3_diel_seasonal.jpg"),
       p_fig3, width = 11, height = 8, dpi = 350, bg = "white")


#### SEASONAL AMPLITUDE & TIMING ####

#fitted curves
season_metrics <- season_pred %>%
  group_by(latitude, latitude_plot) %>%
  summarise(
    min_probability = min(fit),
    max_probability = max(fit),
    fitted_range = max(fit) - min(fit),
    fitted_range_pp = 100 * fitted_range,
    peak_date = plot_date[which.max(fit)],
    trough_date = plot_date[which.min(fit)],
    .groups = "drop")

#simulate from the same covariance used for prediction bands
X <- partial_X(season_pred, season_terms)
i <- match(colnames(X), names(coef(m_season_lat)))

set.seed(1120)
B <- MASS::mvrnorm(
  n = n_sim,
  mu = coef(m_season_lat)[i],
  Sigma = m_season_lat$Vp[i, i, drop = FALSE]
)

P <- plogis(X %*% t(B))

sim_metrics <- map_dfr(seq_along(lat_vals), function(j) {
  ii <- which(season_pred$latitude == lat_vals[j])
  M <- P[ii, , drop = FALSE]
  
  tibble(
    latitude = lat_vals[j],
    latitude_plot = lat_labels[j],
    simulation = seq_len(n_sim),
    seasonal_range = apply(M, 2, function(x) max(x) - min(x)),
    peak_date = season_pred$plot_date[ii][apply(M, 2, which.max)],
    trough_date = season_pred$plot_date[ii][apply(M, 2, which.min)]
  )
})

range_summary <- sim_metrics %>%
  group_by(latitude, latitude_plot) %>%
  summarise(
    simulated_range_median_pp = 100 * median(seasonal_range),
    simulated_range_lower_pp = 100 * quantile(seasonal_range, 0.025),
    simulated_range_upper_pp = 100 * quantile(seasonal_range, 0.975),
    .groups = "drop"
  )

season_table <- season_metrics %>%
  mutate(latitude_plot = as.character(latitude_plot)) %>%
  left_join(range_summary, by = c("latitude", "latitude_plot"))

season_table

write_csv(season_table, file.path(out, "tables", "Table_S5_seasonal_amplitude.csv"))

#peak and trough month frequencies
timing_summary <- sim_metrics %>%
  transmute(latitude, latitude_plot, simulation,
            Peak = month(peak_date), Trough = month(trough_date)) %>%
  pivot_longer(c(Peak, Trough), names_to = "event", values_to = "month_num") %>%
  count(latitude, latitude_plot, event, month_num, name = "n") %>%
  complete(nesting(latitude, latitude_plot),
           event = c("Peak", "Trough"), month_num = 1:12,
           fill = list(n = 0L)) %>%
  group_by(latitude, latitude_plot, event) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

write_csv(timing_summary, file.path(out, "tables", "Table_S5_seasonal_timing.csv"))

saveRDS(sim_metrics, file.path(out, "data", "seasonal_simulation_metrics.rds"))


#### SEASONAL SURFACE ####

#support defined by occupied one-degree latitude/month bins
surface_support <- d_season %>%
  mutate(latitude_bin = floor(latitude)) %>%
  count(latitude_bin, month_num, name = "n_support")

surface_pred <- crossing(
  plot_date = seq.Date(as.Date("2021-01-01"), as.Date("2021-12-31"), by = "3 days"),
  latitude = seq(min(d_season$latitude), max(d_season$latitude), length.out = 60)) %>%
  mutate(annual_phase = (yday(plot_date) - 1) / 365,
         month_num = month(plot_date),
         latitude_bin = floor(latitude)) %>%
  reference_data() %>%
  partial_predict(season_terms) %>%
  left_join(surface_support, by = c("latitude_bin", "month_num")) %>%
  mutate(n_support = replace_na(n_support, 0L),
         fit = if_else(n_support > 0, fit, NA_real_))

write_csv(surface_pred, file.path(out, "tables", "Figure_S4_seasonal_surface.csv"))

p_surface <- ggplot(surface_pred, aes(plot_date, latitude, fill = fit)) +
  geom_raster() +
  scale_x_date(date_breaks = "3 months", date_labels = "%b",
               expand = expansion(mult = 0)) +
  scale_fill_viridis_c(labels = scales::percent, na.value = "white") +
  labs(x = "month", y = "latitude (°N)",
       fill = "seasonal partial-effect\nprobability") +
  theme_classic(base_size = 11)

p_timing <- timing_summary %>%
  mutate(latitude_plot = factor(latitude_plot, levels = lat_labels)) %>%
  ggplot(aes(month_num, percent, linetype = event)) +
  geom_line() +
  geom_point(size = 1.3) +
  facet_wrap(~latitude_plot, nrow = 1) +
  scale_x_continuous(breaks = c(1, 4, 7, 10),
                     labels = month.abb[c(1, 4, 7, 10)]) +
  labs(x = "month", y = "simulations (%)", linetype = NULL) +
  theme_classic(base_size = 11)

p_uncertainty <- (p_surface / p_timing) +
  plot_annotation(tag_levels = "a")

ggsave(file.path(out, "figures", "Figure_S4_seasonal_uncertainty.jpg"),
       p_uncertainty, width = 11, height = 8, dpi = 350, bg = "white")


#### SAVE ####

#templates for same-sample comparisons in coon_temp.R
saveRDS(list(season_pred = season_pred, surface_pred = surface_pred, lat_vals = lat_vals),
        file.path(out, "data", "seasonal_prediction_templates.rds"))

capture.output(sessionInfo(), file = file.path(out, "diagnostics", "session_seasonal.txt"))


#### SENSITIVITY CHECK ####
season_template <- season_pred
set.seed(1120)

m_season_geo400 <- mgcv::bam(
  arboreal ~
    s(solar_hour, bs = "cc", k = 20) +
    s(annual_phase, bs = "cc", k = 12) +
    s(latitude, longitude, bs = "sos", m = 2, k = 400) +
    s(grid_factor, bs = "re") +
    s(year_factor, bs = "re") +
    ti(annual_phase, latitude, bs = c("cc", "tp"), k = c(8, 5)),
  family = binomial(link = "logit"),
  data = d_season,
  method = "fREML",
  discrete = FALSE,
  nthreads = n_threads,
  knots = list(solar_hour = c(0, 24), annual_phase = c(0, 1)),
  na.action = na.fail)

saveRDS(m_season_geo400, file.path(out, "models", "m_season_geo400.rds"))

summary(m_season_geo400)
export_gam(m_season_geo400, "season_geo400")

#seasonal partial effects 
get_season_curve <- function(model, label) {
  smooth_labels <- vapply(model$smooth, function(x) x$label, character(1))
  season_terms <- c("s(annual_phase)", "ti(annual_phase,latitude)")
  
  pr <- predict(model,
                newdata = season_template,
                type = "link",
                se.fit = TRUE,
                exclude = setdiff(smooth_labels, season_terms),
                unconditional = FALSE)
  
  season_template %>%
    transmute(plot_date, latitude, latitude_plot,
              model = label,
              fit = as.numeric(plogis(pr$fit)),
              lower = as.numeric(plogis(pr$fit - 1.96 * pr$se.fit)),
              upper = as.numeric(plogis(pr$fit + 1.96 * pr$se.fit)))
}

check_curves <- bind_rows(get_season_curve(m_season_lat, "k200"),
                          get_season_curve(m_season_geo400, "k400"))

write_csv(check_curves, file.path(out, "tables", "seasonal_geographic_k_curves.csv"))

#changes in fitted probabilities and annual ranges
check_summary <- check_curves %>%
  select(plot_date, latitude_plot, model, fit) %>%
  pivot_wider(names_from = model, values_from = fit) %>%
  group_by(latitude_plot) %>%
  summarise(max_absolute_change_pp = 100 * max(abs(k400 - k200)),
            mean_absolute_change_pp = 100 * mean(abs(k400 - k200)),
            range_k200_pp = 100 * (max(k200) - min(k200)),
            range_k400_pp = 100 * (max(k400) - min(k400)),
            .groups = "drop")

check_summary

write_csv(check_summary, file.path(out, "tables", "seasonal_geographic_k_summary.csv"))

#diagnostic comparison figure
p_geo_check <- ggplot(check_curves, aes(plot_date, fit, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~latitude_plot, nrow = 1) +
  scale_color_manual(values = c(k200 = "gray30", k400 = "steelblue")) +
  scale_fill_manual(values = c(k200 = "gray30", k400 = "steelblue")) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "month", y = "seasonal partial-effect probability",
       color = "geographic basis", fill = "geographic basis") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

p_geo_check

ggsave(file.path(out, "figures", "seasonal_geographic_k_comparison.jpg"),
       p_geo_check, width = 11, height = 4, dpi = 350, bg = "white")


ggsave(file.path(out, "figures", "Figure_3_diel_seasonal.jpg"),
       p_fig3, width = 8, height = 8, dpi = 350, bg = "white")
