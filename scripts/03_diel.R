#### RACCOON TIME OF DAY ANALYSIS ####

library(tidyverse)
library(lubridate)
library(suncalc)
library(patchwork)


#### SETUP ####

out <- "revision_equalarea"

for (x in c("data", "tables", "figures", "diagnostics")) {
  dir.create(file.path(out, x), recursive = TRUE, showWarnings = FALSE)
}

#observation-level data with revised grid IDs
d <- readRDS(file.path(out, "data", "observations_prepared.rds"))


#### TIMESTAMPS ####

#parse timestamps
d_time <- d %>%
  mutate(id = as.character(id),
         time_observed_clean = as.character(time_observed_at) %>%
           str_replace_all("\u00A0", " ") %>%
           str_squish() %>%
           na_if(""),
         observed_date = ymd(observed_on, quiet = TRUE),
         observed_datetime_utc = ymd_hms(time_observed_clean, tz = "UTC", quiet = TRUE),
         has_time = !is.na(observed_datetime_utc),
         hour_utc = hour(observed_datetime_utc),
         minute_utc = minute(observed_datetime_utc),
         second_utc = second(observed_datetime_utc),
         exact_midnight_utc = has_time & hour_utc == 0 & minute_utc == 0 & second_utc == 0)

#timestamp audit
timestamp_QC <- d_time %>%
  summarise(
    n_total = n(),
    missing_time_raw = sum(is.na(time_observed_clean)),
    failed_time_parse = sum(!is.na(time_observed_clean) & !has_time),
    missing_or_invalid_date = sum(is.na(observed_date)),
    usable_timestamp = sum(has_time),
    percent_usable_timestamp = 100 * mean(has_time),
    exact_midnight_utc = sum(exact_midnight_utc))

timestamp_QC

write_csv(timestamp_QC, file.path(out, "tables", "Table_S1_timestamp_QC.csv"))

#records requiring timestamp review
d_time %>%
  filter(is.na(observed_date) | !has_time) %>%
  select(id, observed_on, time_observed_at, time_observed_clean,
         observed_date, observed_datetime_utc, arboreal,
         latitude, longitude, url) %>%
  write_csv(file.path(out, "tables", "timestamp_exclusions.csv"))

#flag midnight values
d_time %>%
  filter(exact_midnight_utc) %>%
  select(id, observed_on, time_observed_at, observed_datetime_utc,
         arboreal, latitude, longitude, url) %>%
  write_csv(file.path(out, "tables", "midnight_timestamp_review.csv"))

#timestamp availability by behavioral classification
timestamp_missingness <- d_time %>%
  count(arboreal, has_time, name = "n") %>%
  group_by(arboreal) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

timestamp_missingness

write_csv(timestamp_missingness, file.path(out, "tables", "timestamp_missingness_by_behavior.csv"))

#timestamp rounding
timestamp_rounding <- d_time %>%
  filter(has_time) %>%
  summarise(
    n = n(),
    exact_midnight_utc = sum(exact_midnight_utc),
    exact_hour = sum(minute_utc == 0 & second_utc == 0),
    percent_exact_hour = 100 * mean(minute_utc == 0 & second_utc == 0),
    exact_half_hour = sum(minute_utc %in% c(0, 30) & second_utc == 0),
    percent_exact_half_hour = 100 * mean(minute_utc %in% c(0, 30) & second_utc == 0))

write_csv(timestamp_rounding, file.path(out, "tables", "timestamp_rounding.csv"))

#all eligible observations
d_diel <- d_time %>% filter(has_time, !is.na(observed_date))


#### SOLAR POSITION ####

#sun elevation at the actual observation time
sun_position <- suncalc::getSunlightPosition(
  data = data.frame(date = d_diel$observed_datetime_utc,
                    lat = d_diel$latitude,
                    lon = d_diel$longitude))

d_diel$sun_elevation <- sun_position$altitude * 180 / pi

#existing solar-period definitions
d_diel <- d_diel %>%
  mutate(solar_period = case_when(
    sun_elevation >= 0 ~ "Day",
    sun_elevation >= -6 ~ "Twilight",
    sun_elevation < -6 ~ "Night",
    TRUE ~ NA_character_),
    solar_period = factor(solar_period, levels = c("Night", "Twilight", "Day")))


#### SOLAR TIME & ANNUAL PHASE ####

#solar noon for recorded observation date
solar_noon <- suncalc::getSunlightTimes(
  data = data.frame(date = d_diel$observed_date,
                    lat = d_diel$latitude,
                    lon = d_diel$longitude),
  keep = "solarNoon", tz = "UTC")

d_diel$solar_noon_utc <- solar_noon$solarNoon

#calculate annual phase 
d_diel <- d_diel %>%
  mutate(hours_from_solar_noon = as.numeric(
    difftime(observed_datetime_utc, solar_noon_utc, units = "hours")),
    solar_hour = (hours_from_solar_noon + 12) %% 24,
    day_of_year = yday(observed_date),
    annual_phase = (day_of_year - 1) / if_else(leap_year(observed_date), 366, 365))

#check variables
solar_QC <- d_diel %>%
  summarise(n = n(),
            missing_elevation = sum(!is.finite(sun_elevation)),
            missing_solar_noon = sum(is.na(solar_noon_utc)),
            missing_solar_hour = sum(!is.finite(solar_hour)),
            missing_annual_phase = sum(!is.finite(annual_phase)))

solar_QC

write_csv(solar_QC, file.path(out, "tables", "solar_QC.csv"))


#### RAW DIEL SUMMARIES ####

#unadjusted proportions
diel_summary <- d_diel %>%
  group_by(solar_period) %>%
  summarise(n_total = n(),
            n_arboreal = sum(arboreal == 1),
            prop_arboreal = mean(arboreal),
            .groups = "drop") %>%
  rowwise() %>%
  mutate(ci = list(binom.test(n_arboreal, n_total)$conf.int),
         ci_low = ci[[1]],
         ci_high = ci[[2]]) %>%
  ungroup() %>%
  select(-ci)

diel_summary

write_csv(diel_summary, file.path(out, "tables", "Table_S1_raw_diel_proportions.csv"))

#counts and proportions by solar hour
d_hourly <- d_diel %>%
  mutate(solar_hour_bin = floor(solar_hour)) %>%
  group_by(solar_hour_bin) %>%
  summarise(n = n(),
            n_arboreal = sum(arboreal == 1),
            prop_arboreal = mean(arboreal),
            .groups = "drop") %>%
  complete(solar_hour_bin = 0:23, fill = list(n = 0L, n_arboreal = 0L))

write_csv(d_hourly, file.path(out, "tables", "solar_hour_summary.csv"))


#### DESCRIPTIVE PLOTS ####

#raw proportions by solar period
p_period <- ggplot(diel_summary, aes(solar_period, prop_arboreal)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15) +
  geom_point(size = 2.5) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "solar period", y = "raw arboreal proportion") +
  theme_classic(base_size = 11)

#raw hourly proportions
p_hour_prop <- ggplot(d_hourly, aes(solar_hour_bin + 0.5, prop_arboreal)) +
  geom_line(color = "gray40") +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = seq(0, 24, 6), limits = c(0, 24)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.45)) +
  labs(x = "", y = "raw arboreal proportion") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_blank())

#sampling distribution
p_hour_n <- ggplot(d_hourly, aes(solar_hour_bin + 0.5, n)) +
  geom_col(width = 1, fill = "gray50") +
  scale_x_continuous(breaks = seq(0, 24, 6), limits = c(-1, 25),
                     expand = expansion(mult = 0)) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "solar hour", y = "observations") +
  theme_classic(base_size = 11)

p_diel_raw <- (p_hour_prop / p_hour_n) + plot_annotation(tag_levels = "a")

p_diel_raw

ggsave(file.path(out, "figures", "diel_descriptive.jpg"),
       p_diel_raw, width = 6, height = 7, dpi = 350, bg = "white")


#### SAVE ####

#RDS preserves UTC datetime
saveRDS(d_diel, file.path(out, "data", "raw_obs_diel.rds"))

#sample accounting
d_diel %>%
  summarise(dataset = "Diel preparation",
            n = n(),
            n_arboreal = sum(arboreal),
            n_cells = n_distinct(grid_id),
            n_years = n_distinct(year(observed_date)),
            first_date = min(observed_date),
            last_date = max(observed_date)) %>%
  write_csv(file.path(out, "tables", "Table_S1_diel_sample.csv"))

capture.output(sessionInfo(), file = file.path(out, "diagnostics", "session_diel.txt"))
