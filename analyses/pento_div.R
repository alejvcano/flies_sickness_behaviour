########## DATA INPUT & PREP ##########

response_data_raw <- read.csv("data/flies_data.csv", sep = ";")
key               <- read.csv("data/flies_info.csv", sep = ";")

treatment_order  <- c("naive", "infected (acute)", "infected (chronic)")
my_color_palette <- stats::setNames(
  c("#0072B2", "#CC79A7", "#E69F00"),
  treatment_order
)

MINUTES_PER_HOUR <- 60
WINDOW_HOURS     <- 6
WINDOW_MINUTES   <- WINDOW_HOURS * MINUTES_PER_HOUR
ZOOM_MINUTES     <- 3 * 24 * 60
PENTO_THRESHOLD  <- 12000

# ── Compute Last_mov ──
last_mov_raw <- response_data_raw |>
  tidyr::pivot_longer(
    cols      = -c(Date, Time, Light1_Dark0),
    names_to  = "FlyID",
    values_to = "Activity"
  ) |>
  dplyr::mutate(FlyID = stringr::str_remove(FlyID, "^.")) |>
  dplyr::group_by(FlyID) |>
  dplyr::summarise(
    Last_mov = ifelse(any(Activity != 0), max(which(Activity != 0)), 0L),
    .groups  = "drop"
  )

# ── Main data prep ──
response_data_pre <- response_data_raw |>
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
  dplyr::left_join(last_mov_raw, by = "FlyID") |>
  dplyr::mutate(
    Treatment = dplyr::case_when(
      Treatment == "Pentomophila" & Last_mov <  PENTO_THRESHOLD ~ "infected (acute)",
      Treatment == "Pentomophila" & Last_mov >= PENTO_THRESHOLD ~ "infected (chronic)",
      TRUE ~ Treatment
    )
  ) |>
  dplyr::filter(Treatment %in% treatment_order) |>
  remove_trailing_zeros() |>
  dplyr::mutate(
    Treatment    = factor(Treatment, levels = treatment_order),
    Fly_Sex      = dplyr::recode(Fly_Sex, "F" = "Female", "M" = "Male"),
    Fly_Sex      = factor(Fly_Sex, levels = c("Female", "Male")),
    Window       = ((Minute - 1) %/% WINDOW_MINUTES),
    Window_label = sprintf("h%02d-%02d", Window * WINDOW_HOURS + 1,
                           (Window + 1) * WINDOW_HOURS)
  )

cat("Final counts:\n")
print(dplyr::distinct(response_data_pre, FlyID, Treatment, Fly_Sex) |>
        dplyr::count(Treatment, Fly_Sex))

window_levels <- unique(response_data_pre$Window_label[order(response_data_pre$Window)])
response_data_pre <- response_data_pre |>
  dplyr::mutate(Window_label = factor(Window_label, levels = window_levels))

########## PENTO TRAJECTORIES (Hourly, Day 1 only) ##########

pento_trajectories <- response_data_pre |>
  dplyr::filter(Treatment %in% c("infected (acute)", "infected (chronic)")) |>
  dplyr::filter(Minute <= 24 * 60) |>
  dplyr::mutate(Hour = ((Minute - 1) %/% 60) + 1) |>
  dplyr::group_by(FlyID, Treatment, Fly_Sex, Hour) |>
  dplyr::summarise(
    Activity_mean = mean(Activity, na.rm = TRUE),
    n_minutes     = dplyr::n(),
    .groups = "drop"
  )

########## PENTO B BASELINE + DIVERGENCE SCORES ##########

pentoB_baseline <- pento_trajectories |>
  dplyr::filter(Treatment == "infected (chronic)") |>
  dplyr::group_by(Hour) |>
  dplyr::summarise(
    mu    = mean(Activity_mean, na.rm = TRUE),
    sigma = sd(Activity_mean,   na.rm = TRUE),
    n_B   = dplyr::n(),
    .groups = "drop"
  )

divergence_scores <- pento_trajectories |>
  dplyr::filter(Treatment == "infected (acute)") |>
  dplyr::left_join(pentoB_baseline, by = "Hour") |>
  dplyr::mutate(
    z_score    = (Activity_mean - mu) / sigma,
    divergence = abs(z_score)
  ) |>
  dplyr::group_by(FlyID, Fly_Sex) |>
  dplyr::arrange(Hour) |>
  dplyr::mutate(
    cum_div   = cumsum(ifelse(is.na(divergence), 0, divergence)),
    first_1sd = {
      idx <- which(abs(z_score) > 1)
      if (length(idx) > 0) Hour[min(idx)] else NA_real_
    }
  ) |>
  dplyr::ungroup()

########## PER-FLY SUMMARY — ALL infected (acute) flies ##########

fly_cum_div <- divergence_scores |>
  dplyr::group_by(FlyID, Fly_Sex) |>
  dplyr::summarise(
    total_cum_div = max(cum_div, na.rm = TRUE),
    mean_div_hour = mean(divergence, na.rm = TRUE),
    n_hours       = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::right_join(
    response_data_pre |>
      dplyr::filter(Treatment == "infected (acute)") |>
      dplyr::distinct(FlyID, Fly_Sex) |>
      dplyr::left_join(last_mov_raw, by = "FlyID"),
    by = c("FlyID", "Fly_Sex")
  )

########## DIVERGENCE SUMMARY ##########

div_summary <- divergence_scores |>
  dplyr::group_by(Fly_Sex) |>
  dplyr::summarise(
    n_A_flies       = dplyr::n_distinct(FlyID),
    n_diverged_1sd  = sum(!is.na(first_1sd)),
    pct_diverged    = round(100 * mean(!is.na(first_1sd)), 1),
    mean_first_hour = mean(first_1sd, na.rm = TRUE),
    mean_z          = mean(z_score, na.rm = TRUE),
    max_cum_div     = max(cum_div, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== PENTO A → B DIVERGENCE SUMMARY (Day 1, Hourly) ===\n")
print(div_summary)

########## THEIL-SEN TREND ##########

ts_stats <- fly_cum_div |>
  dplyr::filter(!is.na(mean_div_hour), !is.na(Last_mov)) |>
  dplyr::group_by(Fly_Sex) |>
  dplyr::group_modify(~ fit_theilsen(.x)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    label = sprintf("%s: slope = %.1f, p = %.3f", Fly_Sex, ts_slope, p_val)
  )

cat("\n=== THEIL-SEN TREND: DIVERGENCE vs LAST_MOV ===\n")
print(ts_stats)

########## SEX COLORS ##########

sex_colors <- c("Female" = "#CC79A7", "Male" = "#0072B2")

########## SUPPLEMENTARY FIGURE: panels A + B ##########

p_supp_A <- ggplot2::ggplot(
  pento_trajectories,
  ggplot2::aes(x = Hour, y = Activity_mean, color = Treatment)
) +
  ggplot2::geom_smooth(method = "loess", se = TRUE, alpha = 0.3, linewidth = 1.2) +
  ggplot2::geom_point(alpha = 0.4, size = 0.8) +
  ggplot2::facet_wrap(~ Fly_Sex) +
  ggplot2::scale_x_continuous(breaks = seq(1, 24, by = 2)) +
  ggplot2::labs(
    x = "Hour post-infection",
    y = "Mean Activity (per hour)"           # caption removed
  ) +
  ggplot2::scale_color_manual(
    values = my_color_palette[c("infected (acute)", "infected (chronic)")]
  ) +
  ggplot2::theme_minimal(base_size = 16) +   # increased from 12
  ggplot2::theme(
    legend.position  = "bottom",
    legend.title     = ggplot2::element_blank(),
    legend.text      = ggplot2::element_text(size = 14),
    axis.text        = ggplot2::element_text(size = 13),
    axis.title       = ggplot2::element_text(size = 15),
    strip.text       = ggplot2::element_text(face = "bold", size = 15),
    panel.grid.minor = ggplot2::element_blank()
  )

p_supp_B <- ggplot2::ggplot(
  divergence_scores |> dplyr::filter(!is.na(z_score)),
  ggplot2::aes(x = Hour, y = z_score)
) +
  ggplot2::geom_hline(yintercept = 0,        linetype = "solid",  color = "gray50") +
  ggplot2::geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "#E69F00",
                      alpha = 0.8) +
  ggplot2::geom_line(ggplot2::aes(group = FlyID), alpha = 0.3, color = "#CC79A7") +
  ggplot2::geom_smooth(method = "loess", se = TRUE, color = "#CC79A7",
                       fill = "#CC79A7", alpha = 0.3) +
  ggplot2::facet_wrap(~ Fly_Sex) +
  ggplot2::scale_x_continuous(breaks = seq(1, 24, by = 2)) +
  ggplot2::labs(                             # subtitle removed
    x = "Hour post-infection",
    y = "Z-score from chronic mean"
  ) +
  ggplot2::theme_minimal(base_size = 16) +   # increased from 12
  ggplot2::theme(
    legend.position  = "none",
    axis.text        = ggplot2::element_text(size = 13),
    axis.title       = ggplot2::element_text(size = 15),
    strip.text       = ggplot2::element_text(face = "bold", size = 15),
    panel.grid.minor = ggplot2::element_blank()
  )

Fig_supplementary <- cowplot::plot_grid(
  p_supp_A, p_supp_B,
  labels = c("A", "B"), label_size = 18, ncol = 2  # label_size increased from 14
)

print(Fig_supplementary)
cowplot::save_plot("outputs/Figure_supp_PentoAB_trajectories.png", Fig_supplementary,
                   base_width = 14, base_height = 6)

########## MAIN FIGURE ##########

cum_div_summary <- divergence_scores |>
  dplyr::group_by(Fly_Sex, Hour) |>
  dplyr::summarise(mean_cum_div = mean(cum_div, na.rm = TRUE), .groups = "drop")

p_main_A <- ggplot2::ggplot(
  cum_div_summary,
  ggplot2::aes(x = Hour, y = mean_cum_div,
               color = Fly_Sex, linetype = Fly_Sex, group = Fly_Sex)
) +
  ggplot2::geom_line(linewidth = 1.2) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_color_manual(values = sex_colors, name = "Sex") +
  ggplot2::scale_linetype_manual(
    values = c("Female" = "solid", "Male" = "dashed"),
    name   = "Sex"
  ) +
  ggplot2::scale_x_continuous(breaks = seq(1, 24, by = 2)) +
  ggplot2::labs(
    x = "Hour post-infection",
    y = "Mean Cumulative |Z-score|"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position  = "bottom",
    legend.title     = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank()
  )

max_lastmov_h <- ceiling(max(fly_cum_div$Last_mov, na.rm = TRUE) / 60 / 12) * 12
y_breaks_h    <- seq(0, max_lastmov_h, by = 12)

ts_lines <- ts_stats |>
  dplyr::mutate(
    slope_h     = ts_slope / 60,
    intercept_h = ts_intercept / 60
  )

p_main_B <- ggplot2::ggplot(
  fly_cum_div |> dplyr::filter(!is.na(mean_div_hour)),
  ggplot2::aes(x = mean_div_hour, y = Last_mov / 60,
               color = Fly_Sex, shape = Fly_Sex)
) +
  ggplot2::geom_point(size = 3, alpha = 0.8) +
  ggplot2::geom_abline(
    data = ts_lines,
    ggplot2::aes(slope = slope_h, intercept = intercept_h, color = Fly_Sex),
    linewidth = 1, alpha = 0.9
  ) +
  ggplot2::geom_text(
    data = ts_stats,
    ggplot2::aes(x = -Inf, y = Inf, label = label, color = Fly_Sex),
    hjust = -0.05, vjust = c(1.5, 3.0),
    size = 3.5, fontface = "italic", inherit.aes = FALSE
  ) +
  ggplot2::scale_color_manual(values = sex_colors, name = "Sex") +
  ggplot2::scale_shape_manual(
    values = c("Female" = 16, "Male" = 17),
    name   = "Sex"
  ) +
  ggplot2::scale_y_continuous(
    breaks = y_breaks_h,
    labels = y_breaks_h
  ) +
  ggplot2::labs(
    x = "Mean Hourly Divergence |Z-score| (Day 1)",
    y = "Last movement (h)"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position  = "bottom",
    legend.title     = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank()
  )

Fig_main <- cowplot::plot_grid(
  p_main_A, p_main_B,
  labels = c("A", "B"), label_size = 14, ncol = 2
)

print(Fig_main)
cowplot::save_plot("outputs/Figure_PentoA_divergence_vs_lastmov.png", Fig_main,
                   base_width = 12, base_height = 6)

########## EXPORT ##########

write.csv(divergence_scores, "outputs/pentoA_divergence_scores_day1.csv", row.names = FALSE)
write.csv(fly_cum_div,       "outputs/pentoA_cumulative_div_per_fly.csv", row.names = FALSE)
write.csv(div_summary,       "outputs/pento_divergence_summary_day1.csv", row.names = FALSE)
write.csv(ts_stats,          "outputs/pento_theilsen_trend.csv",          row.names = FALSE)

cat("\nDivergence analysis complete!\n")

########## STATISTICAL ANALYSIS: ACUTE vs CHRONIC ACTIVITY DAY 1 ##########

hourly_wilcox <- pento_trajectories |>
  dplyr::group_by(Fly_Sex, Hour) |>
  dplyr::summarise(
    n_acute   = sum(Treatment == "infected (acute)"),
    n_chronic = sum(Treatment == "infected (chronic)"),
    p_value   = tryCatch(
      wilcox.test(
        Activity_mean[Treatment == "infected (acute)"],
        Activity_mean[Treatment == "infected (chronic)"],
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    ),
    median_acute   = median(Activity_mean[Treatment == "infected (acute)"],   na.rm = TRUE),
    median_chronic = median(Activity_mean[Treatment == "infected (chronic)"], na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    p_adj       = p.adjust(p_value, method = "BH"),
    acute_lower = median_acute < median_chronic
  )

cat("\n=== [S1] HOURLY WILCOXON: ACUTE vs CHRONIC ACTIVITY (Day 1, BH-corrected) ===\n")
print(hourly_wilcox)

cat("\nHours where acute < chronic (median), per sex:\n")
hourly_wilcox |>
  dplyr::group_by(Fly_Sex) |>
  dplyr::summarise(
    n_hours_acute_lower   = sum(acute_lower, na.rm = TRUE),
    n_hours_sig_p05       = sum(p_value < 0.05, na.rm = TRUE),
    n_hours_sig_adj_p05   = sum(p_adj   < 0.05, na.rm = TRUE),
    pct_hours_acute_lower = round(100 * mean(acute_lower, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  print()

cat("\n=== [S2] LMM: ACTIVITY ~ TREATMENT * HOUR + (1|FlyID), Day 1 ===\n")

lmm_day1 <- pento_trajectories |>
  dplyr::mutate(Hour = as.numeric(Hour)) |>
  dplyr::group_by(Fly_Sex) |>
  dplyr::group_map(~ {
    fit <- lmerTest::lmer(
      Activity_mean ~ Treatment * Hour + (1 | FlyID),
      data    = .x,
      REML    = TRUE,
      control = lme4::lmerControl(optimizer = "bobyqa")
    )
    anova_res <- anova(fit, type = "III") |>
      as.data.frame() |>
      tibble::rownames_to_column("Term") |>
      dplyr::mutate(Fly_Sex = .y$Fly_Sex)
    anova_res
  }, .keep = TRUE) |>
  dplyr::bind_rows()

print(lmm_day1)

cat("\n=== [S3] ONE-SAMPLE WILCOXON: Z-SCORES < 0 ACROSS DAY 1 ===\n")

z_below_zero <- divergence_scores |>
  dplyr::filter(!is.na(z_score)) |>
  dplyr::group_by(Fly_Sex) |>
  dplyr::summarise(
    n_obs       = dplyr::n(),
    median_z    = median(z_score, na.rm = TRUE),
    pct_below_0 = round(100 * mean(z_score < 0, na.rm = TRUE), 1),
    p_value     = wilcox.test(z_score, mu = 0, alternative = "less", exact = FALSE)$p.value,
    .groups = "drop"
  )

cat("\nMedian Z-score and one-sided test (H1: Z < 0):\n")
print(z_below_zero)

cat("\nPer-hour mean Z-score (acute vs chronic baseline):\n")
divergence_scores |>
  dplyr::filter(!is.na(z_score)) |>
  dplyr::group_by(Fly_Sex, Hour) |>
  dplyr::summarise(
    mean_z      = mean(z_score, na.rm = TRUE),
    median_z    = median(z_score, na.rm = TRUE),
    pct_below_0 = round(100 * mean(z_score < 0, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  dplyr::group_by(Fly_Sex) |>
  dplyr::mutate(n_hours_below_0 = sum(mean_z < 0)) |>
  print(n = Inf)
