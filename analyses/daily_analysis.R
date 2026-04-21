########## DATA INPUT & PREP - DAILY BINNING (FULL EXPERIMENT) ##########

response_data_raw <- read.csv("data/flies_data.csv", sep = ";")
key               <- read.csv("data/flies_info.csv", sep = ";")

treatment_order <- c("naive", "injured", "immune stimulated", "infected")

my_color_palette <- stats::setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7"),
  treatment_order
)

treatment_remap <- c(
  "naive"        = "naive",
  "PBS"          = "injured",
  "heatkill"     = "immune stimulated",
  "Pentomophila" = "infected"
)

MINUTES_PER_DAY <- 1440
endtime <- 40000

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
  dplyr::filter(Treatment != "BLANK", Minute <= endtime) |>
  dplyr::mutate(Treatment = dplyr::recode(Treatment, !!!treatment_remap)) |>
  remove_trailing_zeros() |>
  dplyr::mutate(
    Treatment = factor(Treatment, levels = treatment_order),
    Fly_Sex   = dplyr::recode(Fly_Sex, "F" = "Female", "M" = "Male"),
    Fly_Sex   = factor(Fly_Sex, levels = c("Female", "Male"))
  )

daily_activity <- response_data_pre |>
  dplyr::mutate(
    Day = ((Minute - 1) %/% MINUTES_PER_DAY) + 1
  ) |>
  dplyr::group_by(FlyID, Fly_Sex, Treatment, Day) |>
  dplyr::summarise(
    total_activity   = sum(Activity, na.rm = TRUE),
    minutes_recorded = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    total_activity      = pmax(total_activity, 1e-6),
    activity_per_minute = total_activity / minutes_recorded
  )

max_day <- max(daily_activity$Day, na.rm = TRUE)
cat("Experiment runs for", max_day, "days\n")

daily_means <- daily_activity |>
  dplyr::group_by(Treatment, Fly_Sex, Day) |>
  dplyr::filter(dplyr::n() >= 2) |>
  dplyr::summarise(
    mean_activity = mean(total_activity, na.rm = TRUE),
    n_flies       = dplyr::n(),
    se_activity   = stats::sd(total_activity, na.rm = TRUE) / sqrt(n_flies),
    mean_minutes  = mean(minutes_recorded),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    ci_upper = mean_activity + qt(0.975, df = n_flies - 1) * se_activity,
    ci_lower = pmax(mean_activity - qt(0.975, df = n_flies - 1) * se_activity, 0)
  )

cat("Days with data per treatment and sex:\n")
print(table(daily_means$Treatment, daily_means$Fly_Sex, daily_means$Day))

########## PLOT - FIRST 3 DAYS MEAN ACTIVITY SUMMARY ##########

daily_means_3d <- daily_means |>
  dplyr::filter(Day <= 3)

Fig_early_activity <- ggplot2::ggplot(
  daily_means_3d,
  ggplot2::aes(x = Day, y = mean_activity, color = Treatment, fill = Treatment)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
    alpha = 0.2, color = NA
  ) +
  ggplot2::geom_line(linewidth = 1.2, alpha = 0.9) +
  ggplot2::geom_point(size = 3.0, alpha = 0.9) +
  ggplot2::facet_wrap(
    ~ Fly_Sex,
    ncol           = 2,
    strip.position = "top"
  ) +
  ggplot2::labs(
    x = "Day post-infection",
    y = "Mean total activity per day (a.u.)"
  ) +
  ggplot2::scale_color_manual(values = my_color_palette) +
  ggplot2::scale_fill_manual(values = my_color_palette) +
  ggplot2::scale_x_continuous(breaks = 1:3) +
  ggplot2::scale_y_continuous(labels = scales::comma_format()) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    legend.position  = "bottom",
    legend.title     = ggplot2::element_blank(),
    legend.text      = ggplot2::element_text(size = 12),
    legend.key.width = ggplot2::unit(1.5, "cm"),
    axis.text        = ggplot2::element_text(size = 11),
    axis.title       = ggplot2::element_text(size = 13),
    strip.text       = ggplot2::element_text(face = "bold", size = 13),
    panel.grid.minor = ggplot2::element_blank()
  )

print(Fig_early_activity)
cowplot::save_plot("outputs/Figure_early_activity_day1_3.png", Fig_early_activity,
                   base_width = 10, base_height = 5)

########## PLOT - TWO PANELS (FEMALE TOP / MALE BOTTOM) ##########

Fig_daily_activity <- ggplot2::ggplot(
  daily_means,
  ggplot2::aes(x = Day, y = mean_activity, color = Treatment, fill = Treatment)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
    alpha = 0.2, color = NA
  ) +
  ggplot2::geom_line(size = 1.2, alpha = 0.9) +
  ggplot2::geom_point(size = 2.5, alpha = 0.9) +
  ggplot2::facet_wrap(
    ~ Fly_Sex,
    ncol           = 1,
    strip.position = "top"
  ) +
  ggplot2::labs(
    x = "Days",
    y = "Mean total activity per day (a.u.)"
  ) +
  ggplot2::scale_color_manual(values = my_color_palette) +
  ggplot2::scale_fill_manual(values = my_color_palette) +
  ggplot2::scale_x_continuous(breaks = seq(1, max_day, by = 2)) +
  ggplot2::scale_y_log10(labels = scales::scientific_format()) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    legend.position  = "bottom",
    legend.title     = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(size = 10),
    axis.title       = ggplot2::element_text(size = 12),
    strip.text       = ggplot2::element_text(face = "bold", size = 12),
    panel.grid.minor = ggplot2::element_blank()
  )

print(Fig_daily_activity)
cowplot::save_plot("outputs/Figure_daily_activity.png", Fig_daily_activity,
                   base_width = 12, base_height = 10)

########## SUMMARY STATISTICS ##########

cat("\n=== FULL EXPERIMENT SUMMARY BY SEX ===\n")
daily_summary <- daily_activity |>
  dplyr::group_by(Treatment, Fly_Sex) |>
  dplyr::summarise(
    n_flies         = dplyr::n_distinct(FlyID),
    total_days      = max(Day),
    mean_activity   = mean(total_activity),
    median_activity = median(total_activity),
    .groups = "drop"
  )
print(daily_summary)

cat("\nFlies per treatment and sex:\n")
print(table(daily_activity$Treatment, daily_activity$Fly_Sex))

########## DATA QUALITY CHECK ##########

cat("\n=== DATA QUALITY ===\n")
cat("Days with observations per treatment and sex:\n")
daily_counts <- daily_activity |>
  dplyr::group_by(Treatment, Fly_Sex, Day) |>
  dplyr::summarise(n_flies = dplyr::n(), .groups = "drop")
print(daily_counts |> head(20))

low_coverage <- daily_counts |> dplyr::filter(n_flies < 3)
if (nrow(low_coverage) > 0) {
  cat("\nWARNING: Days with <3 flies (sparse data):\n")
  print(low_coverage)
}




########## PLOT - MEAN ACTIVITY PER FLY (FULL EXPERIMENT, BOXPLOT) ##########

# One value per fly: mean daily activity across all days the fly was alive
fly_mean_activity <- daily_activity |>
  dplyr::group_by(FlyID, Fly_Sex, Treatment) |>
  dplyr::summarise(
    mean_daily_activity = mean(total_activity, na.rm = TRUE),
    n_days              = dplyr::n(),
    .groups = "drop"
  )

Fig_activity_boxplot <- ggplot2::ggplot(
  fly_mean_activity,
  ggplot2::aes(x = Treatment, y = mean_daily_activity,
               color = Treatment, fill = Treatment)
) +
  ggplot2::geom_boxplot(
    alpha    = 0.3,
    width    = 0.5,
    outlier.shape = NA    # outliers shown via jitter instead
  ) +
  ggplot2::geom_jitter(
    width = 0.15, size = 2, alpha = 0.7
  ) +
  ggplot2::facet_wrap(
    ~ Fly_Sex,
    ncol           = 2,
    strip.position = "top"
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Mean daily activity per fly (a.u.)"
  ) +
  ggplot2::scale_color_manual(values = my_color_palette) +
  ggplot2::scale_fill_manual(values  = my_color_palette) +
  ggplot2::scale_x_discrete(
    labels = scales::label_wrap(12)   # wrap long treatment names
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    legend.position  = "none",
    axis.text.x      = ggplot2::element_text(size = 11),
    axis.text.y      = ggplot2::element_text(size = 11),
    axis.title.y     = ggplot2::element_text(size = 13),
    strip.text       = ggplot2::element_text(face = "bold", size = 13),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank()
  )

print(Fig_activity_boxplot)
cowplot::save_plot("outputs/Figure_activity_boxplot_full.png", Fig_activity_boxplot,
                   base_width = 10, base_height = 5)





fly_to_plot <- "12T001"   # change to the FlyID you want

one_fly <- response_data_pre |>
  dplyr::filter(FlyID == fly_to_plot) |>
  dplyr::mutate(
    Day = ((Minute - 1) %/% 1440) + 1,
    Minute_in_day = ((Minute - 1) %% 1440) + 1
  ) |>
  dplyr::group_by(FlyID, Treatment, Fly_Sex, Day) |>
  dplyr::arrange(Minute, .by_group = TRUE) |>
  dplyr::mutate(
    cum_activity = cumsum(tidyr::replace_na(Activity, 0)),
    day_total    = sum(Activity, na.rm = TRUE)
  ) |>
  dplyr::ungroup()

check_table <- one_fly |>
  dplyr::group_by(FlyID, Treatment, Fly_Sex, Day) |>
  dplyr::summarise(
    cum_end   = dplyr::last(cum_activity),
    day_total = dplyr::first(day_total),
    diff      = cum_end - day_total,
    .groups   = "drop"
  )

print(check_table)

one_fly_day1 <- response_data_pre |>
  dplyr::filter(FlyID == fly_to_plot) |>
  dplyr::mutate(
    Day = ((Minute - 1) %/% 1440) + 1
  ) |>
  dplyr::filter(Day == 1) |>
  dplyr::arrange(Minute)

# Plot: activity per minute on day 1
p_day1 <- ggplot2::ggplot(
  data = one_fly_day1,
  ggplot2::aes(x = Minute, y = Activity)
) +
  ggplot2::geom_line(alpha = 0.8, linewidth = 0.7) +
  ggplot2::labs(
    x = "Minute (day 1)",
    y = "Activity (a.u.)",
    title = paste("Activity per minute for Fly", fly_to_plot, "- Day 1")
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1440),
    expand = ggplot2::expansion(add = c(0, 0))
  ) +
  ggplot2::theme_minimal(base_size = 13)

print(p_day1)
