########## DATA PREP (same as before) ##########

response_data_raw <- read.csv("data/flies_data.csv", sep = ";")
key               <- read.csv("data/flies_info.csv", sep = ";")

treatment_order  <- c("naive", "PBS", "heatkill", "Pentomophila")
my_color_palette <- stats::setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7"),
  treatment_order
)

MINUTES_PER_HOUR <- 60
ZOOM_MINUTES     <- 3 * 24 * 60

remove_trailing_zeros <- function(df) {
  df |>
    dplyr::group_by(FlyID) |>
    dplyr::mutate(
      last_nonzero = ifelse(
        any(Activity != 0),
        max(which(Activity != 0), na.rm = TRUE),
        0L
      )
    ) |>
    dplyr::filter(dplyr::row_number() <= last_nonzero) |>
    dplyr::select(-last_nonzero) |>
    dplyr::ungroup()
}

response_data <- response_data_raw |>
  tidyr::pivot_longer(
    cols      = -c(Date, Time, Light1_Dark0),
    names_to  = "FlyID",
    values_to = "Activity"
  ) |>
  dplyr::mutate(FlyID = stringr::str_remove(FlyID, "^.")) |>
  dplyr::group_by(FlyID) |>
  dplyr::mutate(Minute = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::left_join(
    key |> dplyr::select(Monitor_TubeLoc, Fly_Sex, Treatment, Status_dead0_alive1),
    by = dplyr::join_by(FlyID == Monitor_TubeLoc)
  ) |>
  dplyr::filter(Treatment != "BLANK", Minute <= ZOOM_MINUTES) |>
  remove_trailing_zeros() |>
  dplyr::mutate(
    Treatment = factor(Treatment, levels = treatment_order),
    Fly_Sex   = dplyr::recode(Fly_Sex, "F" = "Female", "M" = "Male"),
    Fly_Sex   = factor(Fly_Sex, levels = c("Female", "Male")),
    Hour      = (Minute - 1) / MINUTES_PER_HOUR  # Convert to hours
  )

########## COSINOR MODEL FITTING ##########

fit_cosinor <- function(hours, activity) {
  # Remove NAs
  valid_idx <- !is.na(activity) & !is.na(hours)
  hours     <- hours[valid_idx]
  activity  <- activity[valid_idx]
  
  if (length(activity) < 10) return(NULL)  # Need enough data points
  
  tryCatch({
    # Fit: Activity = Mesor + Amplitude * cos(2π * t / 24 + Phase)
    # Rearrange: Activity = Mesor + A_cos * cos(2πt/24) + A_sin * sin(2πt/24)
    # This is linear in Mesor, A_cos, A_sin
    
    t_rad <- 2 * pi * hours / 24
    
    # Linear regression: Activity ~ Mesor + A_cos*cos(t_rad) + A_sin*sin(t_rad)
    fit <- lm(activity ~ cos(t_rad) + sin(t_rad))
    
    mesor      <- coef(fit)[1]
    a_cos      <- coef(fit)[2]
    a_sin      <- coef(fit)[3]
    amplitude  <- sqrt(a_cos^2 + a_sin^2)
    phase_rad  <- atan2(a_sin, a_cos)
    phase_h    <- (phase_rad / (2 * pi)) * 24
    
    # Calculate R²
    r_squared  <- summary(fit)$r.squared
    
    # Predict fitted values
    fitted <- mesor + amplitude * cos(t_rad - phase_rad)
    residuals_val <- activity - fitted
    rmse <- sqrt(mean(residuals_val^2, na.rm = TRUE))
    
    return(tibble(
      mesor     = mesor,
      amplitude = amplitude,
      phase_h   = phase_h,
      r_squared = r_squared,
      rmse      = rmse,
      n_points  = length(activity)
    ))
  }, error = function(e) {
    return(NULL)
  })
}

########## FIT COSINOR FOR EACH FLY ##########

cosinor_results <- response_data |>
  dplyr::group_by(FlyID, Fly_Sex, Treatment) |>
  dplyr::summarise(
    cosinor = list(fit_cosinor(Hour, Activity)),
    .groups = "drop"
  ) |>
  tidyr::unnest(cosinor) |>
  dplyr::filter(!is.na(mesor))

cat("Cosinor model fitted for", nrow(cosinor_results), "flies\n")
cat("Flies excluded (insufficient data):", 
    dplyr::n_distinct(response_data$FlyID) - nrow(cosinor_results), "\n\n")

########## SUMMARY BY TREATMENT & SEX ##########

cat("=== COSINOR PARAMETERS BY TREATMENT AND SEX ===\n\n")
cosinor_summary <- cosinor_results |>
  dplyr::group_by(Treatment, Fly_Sex) |>
  dplyr::summarise(
    n_flies           = dplyr::n(),
    mean_mesor        = mean(mesor, na.rm = TRUE),
    sd_mesor          = sd(mesor, na.rm = TRUE),
    mean_amplitude    = mean(amplitude, na.rm = TRUE),
    sd_amplitude      = sd(amplitude, na.rm = TRUE),
    mean_phase_h      = mean(phase_h, na.rm = TRUE),
    sd_phase_h        = sd(phase_h, na.rm = TRUE),
    median_r_squared  = median(r_squared, na.rm = TRUE),
    .groups = "drop"
  )

print(cosinor_summary)

########## PLOT 1: AMPLITUDE BY TREATMENT ##########

p_amplitude <- ggplot(cosinor_results, 
                      aes(x = Treatment, y = amplitude, fill = Treatment)) +
  geom_boxplot(alpha = 0.6, outlier.size = 2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  facet_wrap(~ Fly_Sex, ncol = 2) +
  labs(
    x        = "",
    y        = "Circadian Amplitude (activity units)",
    title    = "Circadian Rhythm Strength Across Treatments",
    caption  = "Amplitude = strength of 24h oscillation"
  ) +
  scale_fill_manual(values = my_color_palette, guide = "none") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    strip.text      = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

########## PLOT 2: MESOR BY TREATMENT ##########

p_mesor <- ggplot(cosinor_results, 
                  aes(x = Treatment, y = mesor, fill = Treatment)) +
  geom_boxplot(alpha = 0.6, outlier.size = 2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  facet_wrap(~ Fly_Sex, ncol = 2) +
  labs(
    x        = "",
    y        = "Mesor (Rhythm-adjusted Mean, activity units)",
    title    = "Baseline Activity Level (Mesor)",
    caption  = "Mesor = average activity after removing 24h rhythm"
  ) +
  scale_fill_manual(values = my_color_palette, guide = "none") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    strip.text      = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

########## PLOT 3: PHASE BY TREATMENT ##########

p_phase <- ggplot(cosinor_results, 
                  aes(x = Treatment, y = phase_h, fill = Treatment)) +
  geom_boxplot(alpha = 0.6, outlier.size = 2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  facet_wrap(~ Fly_Sex, ncol = 2) +
  labs(
    x        = "",
    y        = "Phase Shift (hours)",
    title    = "Circadian Phase Across Treatments",
    caption  = "Phase = time of peak activity within 24h cycle"
  ) +
  scale_fill_manual(values = my_color_palette, guide = "none") +
  scale_y_continuous(limits = c(-12, 12), breaks = seq(-12, 12, 6)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    strip.text      = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

########## PLOT 4: FIT QUALITY ##########

p_fit_quality <- ggplot(cosinor_results, 
                        aes(x = Treatment, y = r_squared, fill = Treatment)) +
  geom_boxplot(alpha = 0.6, outlier.size = 2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  facet_wrap(~ Fly_Sex, ncol = 2) +
  labs(
    x        = "",
    y        = "R² (Model Fit Quality)",
    title    = "Cosinor Model Fit Quality",
    caption  = "R² = variance explained by 24h cosine model (higher = stronger rhythm)"
  ) +
  scale_fill_manual(values = my_color_palette, guide = "none") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    strip.text      = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

########## ASSEMBLE PLOTS ##########

Fig_cosinor <- cowplot::plot_grid(
  p_amplitude,
  p_mesor,
  p_phase,
  p_fit_quality,
  labels     = c("A", "B", "C", "D"),
  label_size = 14,
  ncol       = 2
)

print(Fig_cosinor)
cowplot::save_plot("outputs/Figure_cosinor_analysis.png", Fig_cosinor,
                   base_width = 14, base_height = 12)

########## PLOT 5: FITTED CURVES OVERLAID ON RAW DATA ##########

# Select 2 representative flies per treatment+sex combo for visualization
representative_flies <- cosinor_results |>
  dplyr::group_by(Treatment, Fly_Sex) |>
  dplyr::filter(amplitude == median(amplitude)) |>  # Median amplitude fly
  dplyr::slice(1) |>
  dplyr::pull(FlyID)

# Get fitted curves for these flies
response_with_cosinor <- response_data |>
  dplyr::filter(FlyID %in% representative_flies) |>
  dplyr::left_join(
    cosinor_results |> dplyr::select(FlyID, mesor, amplitude, phase_h),
    by = "FlyID"
  ) |>
  dplyr::mutate(
    phase_rad = (phase_h / 24) * 2 * pi,
    t_rad     = (2 * pi * Hour / 24),
    fitted    = mesor + amplitude * cos(t_rad - phase_rad)
  )

p_fitted_curves <- ggplot(response_with_cosinor, aes(x = Hour)) +
  geom_line(aes(y = Activity), alpha = 0.4, size = 0.5, color = "gray50") +
  geom_line(aes(y = fitted, color = Treatment), size = 1.2) +
  facet_grid(Fly_Sex ~ Treatment, scales = "free_y") +
  labs(
    x     = "Time (hours)",
    y     = "Activity (counts/min)",
    title = "Cosinor Model Fit: Representative Flies (Median Amplitude per Group)",
    caption = "Gray line = raw activity | Colored line = fitted 24h cosine curve"
  ) +
  scale_color_manual(values = my_color_palette, guide = "none") +
  scale_x_continuous(breaks = c(0, 24, 48, 72)) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    strip.text  = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )

print(p_fitted_curves)
cowplot::save_plot("outputs/Figure_cosinor_fits.png", p_fitted_curves,
                   base_width = 14, base_height = 8)

########## STATISTICAL TESTS ##########

cat("\n=== STATISTICAL TESTS ===\n\n")

# ANOVA on Amplitude by Treatment (within each sex)
cat("Amplitude differences across treatments:\n")
for (sex in c("Female", "Male")) {
  data_sex <- cosinor_results |> dplyr::filter(Fly_Sex == sex)
  fit_aov  <- aov(amplitude ~ Treatment, data = data_sex)
  p_val    <- summary(fit_aov)[[1]]$`Pr(>F)`[1]
  cat(sprintf("  %s: F-test p = %.4f\n", sex, p_val))
}

cat("\nMesor differences across treatments:\n")
for (sex in c("Female", "Male")) {
  data_sex <- cosinor_results |> dplyr::filter(Fly_Sex == sex)
  fit_aov  <- aov(mesor ~ Treatment, data = data_sex)
  p_val    <- summary(fit_aov)[[1]]$`Pr(>F)`[1]
  cat(sprintf("  %s: F-test p = %.4f\n", sex, p_val))
}

cat("\nPhase differences across treatments:\n")
for (sex in c("Female", "Male")) {
  data_sex <- cosinor_results |> dplyr::filter(Fly_Sex == sex)
  fit_aov  <- aov(phase_h ~ Treatment, data = data_sex)
  p_val    <- summary(fit_aov)[[1]]$`Pr(>F)`[1]
  cat(sprintf("  %s: F-test p = %.4f\n", sex, p_val))
}

########## EXPORT RESULTS ##########

write.csv(cosinor_results, "outputs/cosinor_results_per_fly.csv", row.names = FALSE)
write.csv(cosinor_summary, "outputs/cosinor_summary_by_treatment.csv", row.names = FALSE)

cat("\nResults saved to:\n")
cat("  - outputs/cosinor_results_per_fly.csv\n")
cat("  - outputs/cosinor_summary_by_treatment.csv\n")
cat("  - outputs/Figure_cosinor_analysis.png\n")
cat("  - outputs/Figure_cosinor_fits.png\n")
