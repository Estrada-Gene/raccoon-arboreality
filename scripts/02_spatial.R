#### RACCOON MODELING ####


library(tidyverse)
library(spaMM)
library(spdep)
library(DHARMa)
library(patchwork)

#### SETUP ####

out <- "revision_equalarea"

for (x in c("data", "models", "tables", "figures", "diagnostics")) {
  dir.create(file.path(out, x), recursive = TRUE, showWarnings = FALSE)
}

set.seed(1120)
n_threads <- max(1L, parallel::detectCores(logical = FALSE) - 1L, na.rm = TRUE)

# CHANGED: revised equal-area data
df_grid_model <- readRDS(file.path(out, "data", "raccoon_grid_model.rds"))


#### DATA PREP ####

#identical cells for both effort models
required <- c("n_arboreal", "n_total", "latitude", "human_mod", "tree_cover",
              "ndvi", "ann_precip", "precip_seas", "effortRG", "effortRG_mammals",
              "X_coord", "Y_coord")

df_model_prep <- df_grid_model %>%
  filter(n_total >= 3,
         if_all(all_of(required), ~is.finite(.x))) %>%
  arrange(grid_id) %>%
  mutate(prop_arboreal = n_arboreal / n_total,
         log_effortRG = log1p(effortRG),
         log_effortMam = log1p(effortRG_mammals))

#save scaling constants
scale_vars <- c("latitude", "human_mod", "tree_cover", "ndvi",
                "ann_precip", "precip_seas", "log_effortRG", "log_effortMam")

scaling <- map_dfr(scale_vars, function(v) {
  tibble(variable = v,
         center = mean(df_model_prep[[v]]),
         SD = sd(df_model_prep[[v]]))
})

#center and scale
for (v in scale_vars) {
  i <- match(v, scaling$variable)
  df_model_prep[[paste0(v, ".sc")]] <-
    (df_model_prep[[v]] - scaling$center[i]) / scaling$SD[i]
}

write_csv(scaling, file.path(out, "tables", "Table_S2_scaling.csv"))
saveRDS(df_model_prep, file.path(out, "data", "spatial_analysis_frame.rds"))

#sample accounting
df_model_prep %>%
  summarise(
    dataset = "Spatial models",
    n_cells = n(),
    n_observations = sum(n_total),
    n_arboreal = sum(n_arboreal),
    min_latitude = min(latitude),
    max_latitude = max(latitude)) %>%
  write_csv(file.path(out, "tables", "Table_S1_spatial_sample.csv"))


#### PREDICTOR CHECKS ####
library(corrplot)

#correlations using analysis cells
candidate_vars <- c("log_effortRG", "log_effortMam", "human_mod", 
                    "tree_cover", "latitude", "ndvi", "ann_precip", "precip_seas")

cor_mat <- cor(df_model_prep[, candidate_vars], use = "pairwise.complete.obs", method = "spearman")

round(cor_mat, 2)

write_csv(as.data.frame(cor_mat) %>% rownames_to_column("variable"),
          file.path(out, "tables", "Table_S2_correlations.csv"))

clean_names <- c(
  "effort RG (log)", 
  "effort mammal (log)", 
  "human modification", 
  "tree cover %", 
  "latitude", 
  "NDVI", 
  "annual precip.", 
  "precip. seasonality")

rownames(cor_mat) <- clean_names
colnames(cor_mat) <- clean_names

jpeg(file.path(out, "figures", "Figure_s1_correlations.jpg"), 
     width = 8, height = 7, units = "in", res = 300)

corrplot(cor_mat, 
         method = "circle",     
         type = "lower",        
         addCoef.col = "black", 
         number.cex = 0.7,      
         tl.col = "black",      
         tl.srt = 45,           
         diag = FALSE)
dev.off()

#number of complete pairs underlying each correlation
pair_n <- crossprod(1L * !is.na(as.matrix(df_model_prep[, candidate_vars])))

write_csv(as.data.frame(pair_n) %>% rownames_to_column("variable"),
          file.path(out, "tables", "Table_S2_correlation_sample_sizes.csv"))

#final predictors
fixed_vars <- c("log_effortRG", "human_mod", "tree_cover", "latitude",
                "ndvi", "ann_precip", "precip_seas")

fixed_terms <- paste0(fixed_vars, ".sc")


#### SPATIAL MODELS ####

#common scaling and ML fitting for model comparison
f_allRG <- cbind(n_arboreal, n_total - n_arboreal) ~
  log_effortRG.sc +
  human_mod.sc +
  tree_cover.sc +
  latitude.sc +
  ndvi.sc +
  ann_precip.sc +
  precip_seas.sc +
  Matern(1 | X_coord + Y_coord)

f_mammals <- update(f_allRG, . ~ . - log_effortRG.sc + log_effortMam.sc)

m_allRG <- fitme(
  f_allRG,
  family = binomial(link = "logit"),
  data = df_model_prep,
  method = "ML",
  control.HLfit = list(NbThreads = n_threads))

saveRDS(m_allRG, file.path(out, "models", "m_allRG.rds"))
summary(m_allRG)

m_mammals <- fitme(
  f_mammals,
  family = binomial(link = "logit"),
  data = df_model_prep,
  method = "ML",
  control.HLfit = list(NbThreads = n_threads))

saveRDS(m_mammals, file.path(out, "models", "m_mammals.rds"))
summary(m_mammals)

models <- list(allRG = m_allRG, mammals = m_mammals)

#archive full summaries
capture.output(summary(m_allRG, details = TRUE),
               summary(m_mammals, details = TRUE),
               file = file.path(out, "diagnostics", "spatial_model_summaries.txt"))


#### MODEL COMPARISON ####

ic_table <- imap_dfr(models, function(fit, label) {
  ic <- spaMM::get_any_IC(fit, short.names = TRUE, verbose = FALSE)
  tibble(model = label, criterion = names(ic), value = as.numeric(ic)) }) %>%
  group_by(criterion) %>%
  mutate(delta = value - min(value)) %>%
  ungroup()

ic_table <- ic_table %>%
  mutate(delta = if_else(criterion == "GoFdf", NA_real_, delta))

ic_table

write_csv(ic_table, file.path(out, "tables", "Table_S3_information_criteria.csv"))

#recording-intensity correlation on the final model cells
tibble(comparison = "All taxa vs mammals",
       pearson_r = cor(df_model_prep$effortRG, df_model_prep$effortRG_mammals)) %>%
  write_csv(file.path(out, "tables", "Table_S3_effort_correlation.csv"))


#### COEFFICIENT TABLES ####

#export both models with consistent units
extract_coef <- function(fit, label) {
  as.data.frame(summary(fit, details = TRUE)$beta_table) %>%
    rownames_to_column("term") %>%
    rename(SE = `Cond. SE`) %>%
    mutate(model = label,
           OR = exp(Estimate),
           lower = exp(Estimate - 1.96 * SE),
           upper = exp(Estimate + 1.96 * SE),
           unit = if_else(term == "(Intercept)", "Intercept odds",
                          "One SD of grid-level predictor"))
}

coef_table <- imap_dfr(models, extract_coef)

write_csv(coef_table, file.path(out, "tables", "Table_S3_coefficients.csv"))

#save spatial variance and correlation parameters
spatial_parameters <- imap_dfr(models, function(fit, label) {
  pars <- spaMM::get_fittedPars(fit)
  pars <- unlist(pars[intersect(c("lambda", "corrPars"), names(pars))])
  tibble(model = label, parameter = names(pars), estimate = as.numeric(pars))
})

write_csv(spatial_parameters, file.path(out, "tables", "Table_S3_spatial_parameters.csv"))


#### MODEL DIAGNOSTICS ####

#VIF for model fixed predictors
vif_table <- imap_dfr(models, function(fit, label) {
  X <- model.matrix(fit)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  
  tibble(model = label, term = colnames(X),
         VIF = diag(solve(cor(X))))
})

vif_table

write_csv(vif_table, file.path(out, "tables", "Table_S2_VIF.csv"))

#condition simulations on the fitted spatial field
sim_res <- DHARMa::simulateResiduals(m_allRG, n = 1000, seed = 1120, 
                                    simulateREs = "user-specified",  re.form = NULL)

saveRDS(sim_res, file.path(out, "diagnostics", "spatial_DHARMa.rds"))

#raw and residual Moran tests
coords <- as.matrix(df_model_prep[, c("X_coord", "Y_coord")])

knn <- spdep::knearneigh(coords, k = 8)
nb <- spdep::knn2nb(knn)

#symmetrize the row-standardized nearest-neighbor weights
weights <- spdep::listw2U(spdep::nb2listw(nb, style = "W"))

saveRDS(list(grid_id = df_model_prep$grid_id, coords = coords, weights = weights),
        file.path(out, "diagnostics", "spatial_neighbors.rds"))

raw_moran <- spdep::moran.test(df_model_prep$prop_arboreal, weights, alternative = "two.sided")
residual_moran <- spdep::moran.test(sim_res$scaledResiduals, weights, alternative = "two.sided")

#distance-based DHARMa check
distance_test <- DHARMa::testSpatialAutocorrelation(sim_res, x = coords[, 1], 
                                                    y = coords[, 2], plot = FALSE)

diagnostic_tests <- list(
  uniformity = DHARMa::testUniformity(sim_res, plot = FALSE),
  dispersion = DHARMa::testDispersion(sim_res, plot = FALSE),
  raw_moran = raw_moran,
  residual_moran = residual_moran,
  distance_spatial_test = distance_test)

diagnostic_tests

saveRDS(diagnostic_tests, file.path(out, "diagnostics", "spatial_diagnostic_tests.rds"))

capture.output(diagnostic_tests, file = file.path(out, "diagnostics", "spatial_diagnostic_tests.txt"))

#export Moran statistics for supplementary reporting
bind_rows(tibble(response = "Raw arboreal proportion",
                 Moran_I = unname(raw_moran$estimate[1]),
                 statistic = unname(raw_moran$statistic),
                 p_value = raw_moran$p.value),
          tibble(response = "Conditional DHARMa residual",
                 Moran_I = unname(residual_moran$estimate[1]),
                 statistic = unname(residual_moran$statistic),
                 p_value = residual_moran$p.value)) %>%
  write_csv(file.path(out, "tables", "Table_S3_Moran_tests.csv"))

jpeg(file.path(out, "figures", "Figure_S3_spatial_diagnostics.jpg"), 
     width = 10, height = 5, units = "in", res = 300)
plot(sim_res)
dev.off()

#residual source data
df_model_prep %>%
  select(grid_id, X_coord, Y_coord, n_total, prop_arboreal) %>%
  mutate(residual = sim_res$scaledResiduals) %>%
  write_csv(file.path(out, "diagnostics", "spatial_residuals.csv"))


#### COEFFICIENT PLOT ####

labels <- c("latitude.sc" = "latitude", "human_mod.sc" = "human modification",
            "tree_cover.sc" = "tree cover %", "ndvi.sc" = "NDVI",
            "ann_precip.sc" = "annual precip.", "precip_seas.sc" = "precip. seasonality",
            "log_effortRG.sc" = "recording effort")

coef_df <- coef_table %>%
  filter(model == "allRG", term != "(Intercept)") %>%
  mutate(label = unname(labels[term]),
         lower_50 = exp(Estimate - qnorm(0.75) * SE),
         upper_50 = exp(Estimate + qnorm(0.75) * SE))

#plot
p_coef <- ggplot(coef_df, aes(reorder(label, OR), OR)) +
  geom_hline(yintercept = 1, linetype = 2, color = "gray50") +
  geom_linerange(aes(ymin = lower, ymax = upper), linewidth = 0.6) +
  geom_linerange(aes(ymin = lower_50, ymax = upper_50), linewidth = 1.3) +
  geom_point(size = 2.5) +
  scale_y_log10() +
  coord_flip() +
  labs(x = NULL, y = "odds ratio per one SD") +
  theme_classic(base_size = 11)


#### CONDITIONAL EFFECTS ####

#fixed effects only
make_curve <- function(variable) {
  i <- match(variable, scaling$variable)
  
  x <- seq(min(df_model_prep[[variable]]),
           max(df_model_prep[[variable]]), length.out = 200)
  
  nd <- as.data.frame(matrix(
    0, nrow = length(x), ncol = length(fixed_terms),
    dimnames = list(NULL, fixed_terms)
  ))
  
  nd[[paste0(variable, ".sc")]] <-
    (x - scaling$center[i]) / scaling$SD[i]
  
  nd$X_coord <- median(df_model_prep$X_coord)
  nd$Y_coord <- median(df_model_prep$Y_coord)
  
  #predict logits; uncertainty in fixed effects
  pr <- predict(
    m_allRG, newdata = nd, type = "link", re.form = NA,
    variances = list(fixefVar = TRUE)
  )
  
  variance <- attr(pr, "fixefVar")
  stopifnot(!is.null(variance), length(variance) == length(x))
  
  eta <- as.numeric(pr)
  se <- sqrt(pmax(0, as.numeric(variance)))
  
  tibble(
    variable = variable, x = x,
    fit = plogis(eta),
    lower = plogis(eta - 1.96 * se),
    upper = plogis(eta + 1.96 * se)
  )
}

curve_data <- bind_rows(make_curve("latitude"), make_curve("human_mod"), make_curve("ann_precip"))

write_csv(curve_data, file.path(out, "tables", "Figure_2_conditional_curves.csv"))

#shared probability axis
prob_limits <- range(curve_data$lower, curve_data$upper)

#### 1. UPDATE PLOT CURVE FUNCTION WITH TITLES & SELECTIVE TAGS ####

plot_curve <- function(variable, x_label, color, show_y_axis = TRUE, plot_title = NULL, tag_label = NULL) {
  d <- curve_data[curve_data$variable == variable, ]
  rug <- tibble(x = df_model_prep[[variable]])
  
  p <- ggplot(d, aes(x, fit)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                fill = color, alpha = 0.2) +
    geom_line(color = color, linewidth = 0.9) +
    geom_rug(data = rug, aes(x = x), inherit.aes = FALSE,
             alpha = 0.12) +
    coord_cartesian(ylim = prob_limits) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = x_label, 
         y = if (show_y_axis) "arboreal probability" else NULL,
         title = plot_title,
         tag = tag_label) +
    theme_classic(base_size = 11) +
    theme(plot.tag.position = c(0, 1.02),
          plot.tag = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 10, face = "bold", margin = margin(b = 6)))
  
  if (variable == "ann_precip") {
    p <- p + scale_x_continuous(labels = scales::label_comma())
  }
  
  if (show_y_axis) {
    p <- p + theme(axis.title.y = element_text(margin = margin(r = 2)))
  } else {
    p <- p + theme(axis.text.y = element_blank(),
                   axis.title.y = element_blank())
  }
  
  return(p)
}


#top coefficient plot (narrowed via spacers) with tag 'a'
p_top_main <- p_coef + 
  labs(tag = "a") +
  theme(plot.tag.position = c(0, 1.02),
        plot.tag = element_text(face = "bold"))

p_top_row <- plot_spacer() + p_top_main + plot_spacer() + 
  plot_layout(widths = c(0.1, 1, 0.1))

#bottom row conditional effects subplots 
p_lat <- plot_curve("latitude", "latitude (°N)", "grey30", show_y_axis = TRUE, 
                    plot_title = "latitude", tag_label = "b")

p_hm  <- plot_curve("human_mod", "human mod. index", "forestgreen", show_y_axis = FALSE, 
                    plot_title = "human modification", tag_label = NULL)

p_ap  <- plot_curve("ann_precip", "precipitation (mm)", "steelblue", show_y_axis = FALSE, 
                    plot_title = "annual precip.", tag_label = NULL)

p_bottom <- (p_lat + p_hm + p_ap) + 
  plot_layout(nrow = 1, widths = c(1, 1, 1))


#multiplot
p_fig2 <- (p_top_row / free(p_bottom)) +
  plot_layout(heights = c(1, 1))

p_fig2

ggsave(file.path(out, "figures", "Figure_2_spatial_associations.jpg"),
       p_fig2, width = 10, height = 9, dpi = 350, bg = "white")

write_csv(coef_df, file.path(out, "tables", "Figure_2_coefficient_data.csv"))

capture.output(sessionInfo(), file = file.path(out, "diagnostics", "session_spatial.txt"))
