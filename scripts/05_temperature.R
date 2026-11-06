#### RACCOON TEMPERATURE ANALYSIS ####

library(tidyverse)
library(lubridate)
library(mgcv)
library(patchwork)

#### SETUP ####
out <- "revision_equalarea"


cache_file <- file.path(out, "data", "era5_temperature_cache.rds")

for (x in c("data", "models", "tables", "figures", "diagnostics")) {
  dir.create(file.path(out, x), recursive = TRUE, showWarnings = FALSE)
}

n_threads <- min(6L, max(1L, parallel::detectCores(logical = FALSE) - 1L, na.rm = TRUE))

d_season <- readRDS(file.path(out, "data", "d_season_pre_temperature.rds"))

#keep revised grid IDs and valid UTC timestamps
d_season <- d_season %>%
  mutate(id = as.character(id), grid_factor = factor(grid_id),
         year_factor = factor(year(observed_date)),
         era_datetime_utc = round_date(with_tz(observed_datetime_utc, "UTC"), "hour"),
         era_datetime_str = format(era_datetime_utc, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
         era_offset_min = as.numeric(difftime(
           era_datetime_utc, observed_datetime_utc, units = "mins")))

d_request <- d_season %>%
  select(id, latitude, longitude, era_datetime_str) %>%
  mutate(clim_key = format(d_season$era_datetime_utc, "%m_%H", tz = "UTC"))

#normalize cached times only
utc_text <- function(x) {
  if (!inherits(x, "POSIXt")) {
    x <- if (is.numeric(x)) as.POSIXct(x, origin = "1970-01-01", tz = "UTC") else
      ymd_hms(x, tz = "UTC", quiet = TRUE)
  }
  format(x, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

temp_cols <- c("temp_obs_k", "temp_clim_k", "temp_obs_k_fb15", "temp_clim_k_fb15")


#### REUSE COMPLETED EXTRACTION ####
temp_cache <- readRDS(file.path("data", "derived", "era5_temperature_cache.rds"))
#working copy; any additional extraction writes to the ignored output folder
saveRDS(temp_cache, cache_file)

temp_cache <- temp_cache %>%
  mutate(id = as.character(id), era_datetime_str = utc_text(era_datetime_str),
         era_matched_iso = utc_text(era_matched_iso),
         across(all_of(temp_cols), as.numeric))

#verify ID, coordinates and rounded UTC hour
cache_matches <- function(cache) {
  d_request %>%
    inner_join(cache %>%
                 transmute(id, cache_lat = latitude, cache_lon = longitude,
                           cache_hour = era_datetime_str), by = "id") %>%
    filter(abs(latitude - cache_lat) < 1e-7, abs(longitude - cache_lon) < 1e-7,
           era_datetime_str == cache_hour) %>%
    pull(id)
}

reused_ids <- cache_matches(temp_cache)

pending <- d_request %>% filter(!id %in% reused_ids)

tibble(n_requested = nrow(d_request), n_reused = length(reused_ids),
       n_new_requests = nrow(pending)) %>% print()


#### EXTRACT UNMATCHED RECORDS ONLY ####
if (nrow(pending) > 0) {
  if (!nzchar(Sys.getenv("GEE_PROJECT"))) {
    stop("Temperature cache does not cover these IDs, coordinates and hours. Check the input snapshot; fresh extraction requires GEE_PROJECT and rgee.")
  }
  Sys.setenv(RETICULATE_USE_UV = "FALSE")
  reticulate::use_virtualenv("rgee_env", required = TRUE)
  library(rgee)
  rgee::ee_Initialize(project = Sys.getenv("GEE_PROJECT"), drive = FALSE)
  era5_hourly <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY")$select("temperature_2m")
  clim_source <- ee$ImageCollection("ECMWF/ERA5_LAND/MONTHLY_BY_HOUR")$
    filterDate("1991-01-01", "2021-01-01")$select("temperature_2m")
  #1991-2020 climatology for each month x UTC hour
  clim_images <- vector("list", 288)
  i <- 1L
  for (m in 1:12) {
    for (h in 0:23) {
      clim_images[[i]] <- clim_source$
        filter(ee$Filter$calendarRange(m, m, "month"))$
        filter(ee$Filter$eq("hour", h))$mean()$rename("temp_clim_k")$
        set("clim_key", sprintf("%02d_%02d", m, h))
      i <- i + 1L
    }
  }
  
  era5_clim <- ee$ImageCollection$fromImages(clim_images)
  extract_temp <- function(feature) {
    f <- ee$Feature(feature)
    dt <- ee$Date(f$get("era_datetime_str"))
    hourly <- era5_hourly$filterDate(dt, dt$advance(1, "hour"))
    empty <- ee$Image$constant(0)$rename("temperature_2m")$selfMask()$
      set("system:time_start", dt$millis())
    obs <- ee$Image(ee$Algorithms$If(hourly$size()$gt(0), hourly$first(), empty))$
      rename("temp_obs_k")
    clim <- ee$Image(era5_clim$
                       filter(ee$Filter$eq("clim_key", f$get("clim_key")))$first())
    image <- obs$addBands(clim)
    point <- image$reduceRegion(reducer = ee$Reducer$first(),
                                geometry = f$geometry(), scale = 11132, tileScale = 4)
    buffer <- image$reduceRegion(reducer = ee$Reducer$mean(),
                                 geometry = f$geometry()$buffer(15000),
                                 scale = 11132, maxPixels = 1000, tileScale = 4)
    matched <- ee$Algorithms$If(hourly$size()$gt(0),
                                ee$Date(obs$get("system:time_start"))$format("YYYY-MM-dd'T'HH:mm:ss", "UTC"), "")
    ee$Feature(NULL, list(id = f$get("id"), era_matched_iso = matched,
                          temp_obs_k = point$get("temp_obs_k", -9999),
                          temp_clim_k = point$get("temp_clim_k", -9999),
                          temp_obs_k_fb15 = buffer$get("temp_obs_k", -9999),
                          temp_clim_k_fb15 = buffer$get("temp_clim_k", -9999)))
  }
  
  get_temp <- function(p, x) if (is.null(p[[x]])) NA_real_ else as.numeric(p[[x]])
  batch_size <- 50L
  for (b in seq_len(ceiling(nrow(pending) / batch_size))) {
    message("ERA5 batch ", b, "/", ceiling(nrow(pending) / batch_size))
    x <- pending %>% slice(seq.int((b - 1L) * batch_size + 1L,
                                   min(b * batch_size, nrow(pending))))
    x_sf <- sf::st_as_sf(x, coords = c("longitude", "latitude"),
                         crs = 4326, remove = FALSE)
    x_ee <- rgee::sf_as_ee(x_sf, via = "getInfo")
    info <- ee$FeatureCollection(x_ee$map(extract_temp))$getInfo()
    values <- map_dfr(info$features, function(z) {
      p <- z$properties
      tibble(id = as.character(p$id), era_matched_iso = as.character(p$era_matched_iso),
             temp_obs_k = get_temp(p, "temp_obs_k"),
             temp_clim_k = get_temp(p, "temp_clim_k"),
             temp_obs_k_fb15 = get_temp(p, "temp_obs_k_fb15"),
             temp_clim_k_fb15 = get_temp(p, "temp_clim_k_fb15"))
    }) %>%
      mutate(across(all_of(temp_cols), ~ na_if(.x, -9999)),
             era_matched_iso = utc_text(era_matched_iso))
    stopifnot(!anyDuplicated(values$id), setequal(values$id, x$id))
    add <- x %>% select(-clim_key) %>% left_join(values, by = "id")
    temp_cache <- bind_rows(anti_join(temp_cache, add, by = "id"), add)
    saveRDS(temp_cache, cache_file)
  }
}


#### TEMPERATURE DATA ####
d_temp_all <- d_season %>%
  left_join(temp_cache %>% select(id, era_datetime_str, era_matched_iso, all_of(temp_cols)),
            by = c("id", "era_datetime_str")) %>%
  mutate(exact_ok = is.finite(temp_obs_k) & is.finite(temp_clim_k),
         fallback_ok = is.finite(temp_obs_k_fb15) & is.finite(temp_clim_k_fb15),
         temp_source = case_when(
           exact_ok ~ "Exact ERA5-Land pixel",
           fallback_ok ~ "15-km ERA5-Land buffer mean",
           TRUE ~ "Missing"),
         temp_obs_c = case_when(
           exact_ok ~ temp_obs_k - 273.15,
           fallback_ok ~ temp_obs_k_fb15 - 273.15, TRUE ~ NA_real_),
         temp_clim_c = case_when(
           exact_ok ~ temp_clim_k - 273.15,
           fallback_ok ~ temp_clim_k_fb15 - 273.15, TRUE ~ NA_real_),
         temp_anom_c = temp_obs_c - temp_clim_c,
         has_temp = is.finite(temp_obs_c) & is.finite(temp_clim_c),
         era_matched_utc = ymd_hms(era_matched_iso, tz = "UTC", quiet = TRUE))

time_errors <- d_temp_all %>%
  filter(has_temp & (is.na(era_matched_utc) | era_matched_utc != era_datetime_utc)) %>%
  select(id, era_datetime_utc, era_matched_utc, temp_source)

write_csv(time_errors, file.path(out, "diagnostics", "temperature_time_errors.csv"))

stopifnot(nrow(time_errors) == 0, nrow(d_temp_all) == nrow(d_season))

d_temp_model <- d_temp_all %>% filter(has_temp) %>% droplevels()
d_temp_exact <- d_temp_all %>% filter(exact_ok, has_temp) %>% droplevels()

saveRDS(d_temp_all, file.path(out, "data", "raccoon_temperature_all_records.rds"))
saveRDS(d_temp_model, file.path(out, "data", "raccoon_temperature_model_data.rds"))
saveRDS(d_temp_exact, file.path(out, "data", "raccoon_temperature_exact_pixel_data.rds"))
saveRDS(list(hourly_collection = "ECMWF/ERA5_LAND/HOURLY",
             climatology_collection = "ECMWF/ERA5_LAND/MONTHLY_BY_HOUR",
             climatology_years = 1991:2020, climatology_groups = "month x UTC hour",
             scale_m = 11132, fallback_radius_m = 15000),
        file.path(out, "data", "temperature_extraction_settings.rds"))

#### SAMPLE & EXTRACTION QC ####
temp_sample <- imap_dfr(list(seasonal_input = d_season, temperature_models = d_temp_model,
                             exact_pixel_sensitivity = d_temp_exact), function(x, nm) {
                               tibble(dataset = nm, n = nrow(x), n_arboreal = sum(x$arboreal),
                                      n_cells = n_distinct(x$grid_id), n_years = n_distinct(x$year_factor),
                                      first_date = min(x$observed_date), last_date = max(x$observed_date))
                             })

write_csv(temp_sample, file.path(out, "tables", "Table_S1_temperature_sample.csv"))

median_finite <- function(x) if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_

temp_sources <- d_temp_all %>%
  group_by(temp_source) %>%
  summarise(n = n(), n_arboreal = sum(arboreal), prop_arboreal = mean(arboreal),
            median_latitude = median(latitude),
            median_temperature_c = median_finite(temp_obs_c),
            median_anomaly_c = median_finite(temp_anom_c), .groups = "drop") %>%
  mutate(percent_all = 100 * n / nrow(d_temp_all),
         percent_retained = if_else(temp_source == "Missing", NA_real_,
                                    100 * n / nrow(d_temp_model)))

write_csv(temp_sources, file.path(out, "tables", "Table_S1_temperature_sources.csv"))

d_temp_all %>%
  count(arboreal, has_temp, name = "n") %>%
  group_by(arboreal) %>% mutate(percent = 100 * n / sum(n)) %>% ungroup() %>%
  write_csv(file.path(out, "tables", "temperature_missingness_by_behavior.csv"))

tibble(n = nrow(d_season), n_reused = length(reused_ids), n_new_requests = nrow(pending),
       missing_timestamp = sum(is.na(d_season$observed_datetime_utc)),
       maximum_offset_min = max(abs(d_season$era_offset_min)),
       mean_offset_min = mean(abs(d_season$era_offset_min)),
       mismatched_usable_hours = nrow(time_errors)) %>%
  write_csv(file.path(out, "tables", "temperature_timestamp_QC.csv"))

temp_distribution <- map_dfr(c("temp_obs_c", "temp_anom_c"), function(v) {
  x <- d_temp_model[[v]]
  tibble(variable = v, n = length(x), minimum = min(x), q01 = quantile(x, .01),
         q025 = quantile(x, .025), median = median(x), q975 = quantile(x, .975),
         q99 = quantile(x, .99), maximum = max(x))
})

write_csv(temp_distribution, file.path(out, "tables", "temperature_distribution.csv"))

d_temp_all %>% filter(!has_temp) %>%
  select(id, observed_datetime_utc, latitude, longitude, arboreal, any_of("url")) %>%
  write_csv(file.path(out, "tables", "temperature_exclusions.csv"))

temp_sample
temp_sources
temp_distribution

#### MODELS ####
f_temp_base <- arboreal ~
  s(solar_hour, bs = "cc", k = 20) + s(annual_phase, bs = "cc", k = 12) +
  s(latitude, longitude, bs = "sos", m = 2, k = 200) +
  s(grid_factor, bs = "re") + s(year_factor, bs = "re") +
  ti(annual_phase, latitude, bs = c("cc", "tp"), k = c(8, 5))

f_temp_anom <- update(f_temp_base, . ~ . + s(temp_anom_c, bs = "ts", k = 10))

f_temp_abs <- update(f_temp_base, . ~ . + s(temp_obs_c, bs = "ts", k = 12))

fit_temp <- function(formula, data, filename) {
  set.seed(1120)
  message("Fitting ", filename)
  m <- mgcv::bam(formula, family = binomial(link = "logit"), data = data,
                 method = "fREML", discrete = FALSE, nthreads = n_threads,
                 knots = list(solar_hour = c(0, 24), annual_phase = c(0, 1)),
                 na.action = na.fail)
  saveRDS(m, file.path(out, "models", filename))
  m
}

#baseline, anomaly and absolute models
m_temp_base <- fit_temp(f_temp_base, d_temp_model, "m_temp_baseline.rds")

m_temp_anom <- fit_temp(f_temp_anom, d_temp_model, "m_temp_anomaly.rds")

m_temp_absolute <- fit_temp(f_temp_abs, d_temp_model, "m_temp_absolute.rds")

m_temp_anom_exact <- fit_temp(f_temp_anom, d_temp_exact, "m_temp_anomaly_exact_pixels.rds")

models <- list(baseline = m_temp_base, anomaly = m_temp_anom,
               absolute = m_temp_absolute, exact_pixels = m_temp_anom_exact)

model_data <- list(baseline = d_temp_model, anomaly = d_temp_model,
                   absolute = d_temp_model, exact_pixels = d_temp_exact)

#### MODEL EXPORTS ####
model_fit <- imap_dfr(models, function(m, nm) {
  sm <- summary(m)
  capture.output(sm, file = file.path(out, "diagnostics", paste0("temp_", nm, "_summary.txt")))
  tibble(model = nm, n = nobs(m), deviance_explained = sm$dev.expl,
         total_edf = sum(m$edf), converged = m$converged)
})

smooth_table <- imap_dfr(models, ~ as.data.frame(summary(.x)$s.table) %>%
                           rownames_to_column("term") %>% mutate(model = .y, .before = 1))

parametric_table <- imap_dfr(models, ~ as.data.frame(summary(.x)$p.table) %>%
                               rownames_to_column("term") %>% mutate(model = .y, .before = 1))

k_table <- imap_dfr(models, function(m, nm) {
  set.seed(1120)
  as.data.frame(mgcv::k.check(m, subsample = 10000, n.rep = 1000)) %>%
    rownames_to_column("term") %>% mutate(model = nm, .before = 1)
})

concurvity_table <- imap_dfr(models, ~ as.data.frame(mgcv::concurvity(.x, full = TRUE)) %>%
                               rownames_to_column("measure") %>% mutate(model = .y, .before = 1))

write_csv(model_fit, file.path(out, "tables", "Table_S6_temperature_fit.csv"))
write_csv(smooth_table, file.path(out, "tables", "Table_S6_temperature_smooths.csv"))
write_csv(parametric_table, file.path(out, "tables", "Table_S6_temperature_parametric.csv"))
write_csv(k_table, file.path(out, "tables", "Table_S6_temperature_kcheck.csv"))
write_csv(concurvity_table, file.path(out, "tables", "Table_S6_temperature_concurvity.csv"))

saveRDS(map(models, ~ mgcv::concurvity(.x, full = FALSE)),
        file.path(out, "diagnostics", "temperature_pairwise_concurvity.rds"))

model_fit
smooth_table
k_table

#### TEMPERATURE EFFECTS ####
smooth_columns <- function(model, terms) {
  keep <- vapply(model$smooth, function(x) x$label %in% terms, logical(1))
  unlist(lapply(model$smooth[keep], function(x) seq.int(x$first.para, x$last.para)),
         use.names = FALSE)
}

#pointwise 95% intervals conditional on smoothing parameters
or_curve <- function(model, data, variable, x, reference, label) {
  nd <- data[rep(1L, length(x)), ]
  nd[[variable]] <- x
  nd0 <- nd
  nd0[[variable]] <- reference
  j <- smooth_columns(model, paste0("s(", variable, ")"))
  D <- predict(model, nd, type = "lpmatrix")[, j, drop = FALSE] -
    predict(model, nd0, type = "lpmatrix")[, j, drop = FALSE]
  log_or <- drop(D %*% coef(model)[j])
  se <- sqrt(pmax(0, rowSums((D %*% model$Vp[j, j, drop = FALSE]) * D)))
  #calculate before tibble creates the variable column
  observed_range <- range(data[[variable]], na.rm = TRUE)
  in_range <- x >= observed_range[1] & x <= observed_range[2]
  tibble(model = label, variable, x = as.numeric(x), reference,
         odds_ratio = exp(log_or), lower = exp(log_or - 1.96 * se),
         upper = exp(log_or + 1.96 * se),
         in_observed_range = in_range)
}

anom_x <- sort(unique(c(seq(quantile(d_temp_model$temp_anom_c, .01),
                            quantile(d_temp_model$temp_anom_c, .99), length.out = 200), 0)))

temp_reference <- median(d_temp_model$temp_obs_c)

absolute_x <- sort(unique(c(seq(quantile(d_temp_model$temp_obs_c, .01),
                                quantile(d_temp_model$temp_obs_c, .99), length.out = 200),
                            temp_reference)))

anom_effect <- or_curve(m_temp_anom, d_temp_model, "temp_anom_c", anom_x, 0, "anomaly")

absolute_effect <- or_curve(m_temp_absolute, d_temp_model, "temp_obs_c",
                            absolute_x, temp_reference, "absolute")

write_csv(bind_rows(anom_effect, absolute_effect),
          file.path(out, "tables", "Figure_4_temperature_curves.csv"))

#evaluate contrasts at their exact values
or_report <- bind_rows(
  or_curve(m_temp_anom, d_temp_model, "temp_anom_c", c(-5, 5), 0, "anomaly"),
  or_curve(m_temp_anom_exact, d_temp_exact, "temp_anom_c", c(-5, 5), 0, "exact_pixels"),
  or_curve(m_temp_absolute, d_temp_model, "temp_obs_c",
           temp_reference + c(-5, 5), temp_reference, "absolute")) %>%
  mutate(percent_change_in_odds = 100 * (odds_ratio - 1))

write_csv(or_report, file.path(out, "tables", "Table_S7_temperature_OR_contrasts.csv"))

#compare extraction methods over common central support
common_limits <- c(max(quantile(d_temp_model$temp_anom_c, .01),
                       quantile(d_temp_exact$temp_anom_c, .01)),
                   min(quantile(d_temp_model$temp_anom_c, .99),
                       quantile(d_temp_exact$temp_anom_c, .99)))

common_x <- sort(unique(c(seq(common_limits[1], common_limits[2], length.out = 200), 0)))

sensitivity <- bind_rows(
  or_curve(m_temp_anom, d_temp_model, "temp_anom_c", common_x, 0, "Exact + 15-km fallback"),
  or_curve(m_temp_anom_exact, d_temp_exact, "temp_anom_c", common_x, 0, "Exact pixels only"))

write_csv(sensitivity, file.path(out, "tables", "Figure_S6_temperature_sensitivity.csv"))

or_report

#### MAIN FIGURE ####
set.seed(1120)

rug_data <- d_temp_model %>% slice_sample(n = min(5000L, nrow(d_temp_model)))

p_anom <- ggplot(anom_effect, aes(x, odds_ratio)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray50", alpha = .2) +
  geom_line(linewidth = .9) +
  geom_hline(yintercept = 1, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) +
  geom_rug(data = rug_data, aes(x = temp_anom_c), inherit.aes = FALSE,
           sides = "b", alpha = .06) +
  coord_cartesian(xlim = range(anom_effect$x)) + scale_y_log10() +
  labs(x = "temperature anomaly (°C)",
       y = "adjusted odds ratio\n(relative to 0 °C anomaly)") +
  theme_classic(base_size = 11)

p_absolute <- ggplot(absolute_effect, aes(x, odds_ratio)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray50", alpha = .2) +
  geom_line(linewidth = .9) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_vline(xintercept = temp_reference, linetype = 2) +
  geom_rug(data = rug_data, aes(x = temp_obs_c), inherit.aes = FALSE,
           sides = "b", alpha = .06) +
  coord_cartesian(xlim = range(absolute_effect$x)) + scale_y_log10() +
  labs(x = "ERA5-Land 2-m air temperature (°C)",
       y = paste0("adjusted odds ratio\n(relative to ", round(temp_reference, 1), " °C)")) +
  theme_classic(base_size = 11)

p_fig4 <- (p_anom / p_absolute) + plot_annotation(tag_levels = "a")

p_fig4

ggsave(file.path(out, "figures", "Figure_4_temperature.jpg"),
       p_fig4, width = 6, height = 7, dpi = 350, bg = "white")

#### SUPPLEMENTARY FIGURES ####
temp_long <- d_temp_model %>%
  select(temp_obs_c, temp_anom_c) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "temperature") %>%
  mutate(variable = factor(variable, levels = c("temp_obs_c", "temp_anom_c"),
                           labels = c("2-m air temperature (°C)", "temperature anomaly (°C)")))

p_distribution <- ggplot(temp_long, aes(temperature)) +
  geom_histogram(bins = 40, fill = "gray50", color = "white", linewidth = .15) +
  facet_wrap(~ variable, ncol = 1, scales = "free") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "observations") + 
  theme_classic(base_size = 11)

ggsave(file.path(out, "figures", "Figure_S5_temperature_distributions.jpg"),
       p_distribution, width = 6, height = 6, dpi = 350, bg = "white")

p_sensitivity <- ggplot(sensitivity, aes(x, odds_ratio, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .15, color = NA) +
  geom_line(aes(linetype = model), linewidth = .9) +
  geom_hline(yintercept = 1, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) +
  scale_y_log10() +
  scale_color_manual(values = c("gray25", "#0072B2")) +
  scale_fill_manual(values = c("gray25", "#0072B2")) +
  labs(x = "temperature anomaly (°C)", y = "adjusted odds ratio\n(relative to 0 °C anomaly)",
       color = NULL, fill = NULL, linetype = NULL) +
  theme_classic(base_size = 11) + theme(legend.position = "bottom")

ggsave(file.path(out, "figures", "Figure_S6_temperature_sensitivity.jpg"),
       p_sensitivity, width = 6, height = 4.5, dpi = 350, bg = "white")

#### SEASONAL PATTERNS ON IDENTICAL OBSERVATIONS ####
lat_vals <- as.numeric(quantile(d_season$latitude, c(.1, .5, .9)))
season_template <- crossing(plot_date = seq(as.Date("2021-01-01"), as.Date("2021-12-31"),
                                            by = "day"), latitude = lat_vals) %>%
  mutate(annual_phase = (yday(plot_date) - 1) / 365,
         latitude_plot = factor(latitude, levels = lat_vals,
                                labels = paste0(round(lat_vals, 1), "°N")))

#retain intercept and seasonal terms only, as in the seasonal script
season_curve <- function(model, label) {
  nd <- d_temp_model[rep(1L, nrow(season_template)), ]
  nd$latitude <- season_template$latitude
  nd$annual_phase <- season_template$annual_phase
  j <- c(match("(Intercept)", names(coef(model))),
         smooth_columns(model, c("s(annual_phase)", "ti(annual_phase,latitude)")))
  X <- predict(model, nd, type = "lpmatrix")[, j, drop = FALSE]
  eta <- drop(X %*% coef(model)[j])
  se <- sqrt(pmax(0, rowSums((X %*% model$Vp[j, j, drop = FALSE]) * X)))
  season_template %>%
    mutate(model = label, eta, fit = plogis(eta),
           lower = plogis(eta - 1.96 * se), upper = plogis(eta + 1.96 * se))
}

season_compare <- imap_dfr(models[c("baseline", "anomaly", "absolute")], season_curve) %>%
  mutate(model = factor(model, levels = c("baseline", "anomaly", "absolute")))

season_ranges <- season_compare %>%
  group_by(model, latitude, latitude_plot) %>%
  summarise(fitted_range_pp = 100 * (max(fit) - min(fit)),
            fitted_range_log_odds = max(eta) - min(eta),
            fitted_peak_date = plot_date[which.max(fit)],
            fitted_trough_date = plot_date[which.min(fit)], .groups = "drop")

write_csv(season_compare, file.path(out, "tables", "Figure_S7_temperature_seasonal_curves.csv"))
write_csv(season_ranges, file.path(out, "tables", "Table_S7_temperature_seasonal_ranges.csv"))

p_season_temp <- ggplot(season_compare, aes(plot_date, fit, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .12, color = NA) +
  geom_line(linewidth = .9) + facet_wrap(~ latitude_plot, ncol = 1) +
  scale_color_manual(values = c(baseline = "gray25", anomaly = "#0072B2", absolute = "#E69F00")) +
  scale_fill_manual(values = c(baseline = "gray25", anomaly = "#0072B2", absolute = "#E69F00")) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b", expand = expansion(mult = 0)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "month", y = "seasonal partial-effect probability", color = NULL, fill = NULL) +
  theme_classic(base_size = 11) + theme(legend.position = "bottom")

ggsave(file.path(out, "figures", "Figure_S7_temperature_seasonal_comparison.jpg"),
       p_season_temp, width = 6, height = 8, dpi = 350, bg = "white")

#### RESIDUAL DIAGNOSTICS ####
temp_residuals <- imap_dfr(models, function(m, nm) {
  model_data[[nm]] %>%
    select(id, grid_id, year_factor, latitude, longitude, solar_hour, annual_phase,
           temp_obs_c, temp_anom_c) %>%
    mutate(model = nm, predicted = as.numeric(fitted(m)),
           pearson = as.numeric(residuals(m, type = "pearson")))
})

write_csv(temp_residuals, file.path(out, "diagnostics", "temperature_residuals.csv"))

calibration <- temp_residuals %>%
  left_join(d_temp_all %>% select(id, arboreal), by = "id") %>%
  group_by(model) %>% mutate(bin = ntile(predicted, 20)) %>%
  group_by(model, bin) %>%
  summarise(n = n(), predicted = mean(predicted), observed = mean(arboreal), .groups = "drop")

residual_bins <- temp_residuals %>%
  mutate(anomaly_bin = floor(temp_anom_c) + .5) %>%
  group_by(model, anomaly_bin) %>%
  summarise(n = n(), mean_pearson = mean(pearson), .groups = "drop")

write_csv(calibration, file.path(out, "diagnostics", "temperature_calibration.csv"))
write_csv(residual_bins, file.path(out, "diagnostics", "temperature_residual_bins.csv"))

p_calibration <- ggplot(filter(calibration, model == "anomaly"), aes(predicted, observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray50") +
  geom_line() + geom_point(size = 1.8) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "mean fitted probability", y = "observed arboreal proportion") +
  theme_classic(base_size = 11)

p_residual <- ggplot(filter(residual_bins, model == "anomaly", n >= 30),
                     aes(anomaly_bin, mean_pearson)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
  geom_line() + geom_point(size = 1.8) +
  labs(x = "temperature anomaly (°C)", y = "mean Pearson residual") +
  theme_classic(base_size = 11)

p_diagnostics <- (p_calibration / p_residual) + plot_annotation(tag_levels = "a")

ggsave(file.path(out, "diagnostics", "temperature_diagnostics.jpg"),
       p_diagnostics, width = 6, height = 7, dpi = 350, bg = "white")

#### SAVE ####
capture.output(sessionInfo(), file = file.path(out, "diagnostics", "session_temperature.txt"))
