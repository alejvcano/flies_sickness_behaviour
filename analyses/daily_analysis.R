########## HELPER FUNCTIONS ##########

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

########## DATA INPUT & PREP - DAILY BINNING (FULL EXPERIMENT) ##########

response_data_raw <- read.csv("data/flies_data.csv", sep = ";")
key               <- read.csv("data/flies_info.csv", sep = ";")

treatment_order <- c("naive", "PBS", "heatkill", "Pentomophila")
my_color_palette <- stats::setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7"),
  treatment_order
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
  remove_trailing_zeros() |>
  dplyr::mutate(
    Treatment = factor(Treatment, levels = treatment_order),
    # Recode sex labels for clean panel titles — adjust if your codes differ
    Fly_Sex   = dplyr::recode(Fly_Sex, "F" = "Female", "M" = "Male"),
    Fly_Sex   = factor(Fly_Sex, levels = c("Female", "Male"))  # Female first = top panel
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

# Daily means by Treatment AND Sex
daily_means <- daily_activity |>
  dplyr::group_by(Treatment, Fly_Sex, Day) |>          # ← added Fly_Sex
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
    ncol = 1,                   # stacked vertically: Female top, Male bottom
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
                   base_width = 12, base_height = 10)   # taller for 2 stacked panels

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
