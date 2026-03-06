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

########## DATA INPUT & PREP - HOURLY BINNING (FIRST 3 DAYS) ##########

response_data_raw <- read.csv("data/flies_data.csv", sep = ";")
key               <- read.csv("data/flies_info.csv", sep = ";")

treatment_order <- c("naive", "PBS", "heatkill", "Pentomophila")
my_color_palette <- stats::setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7"),
  treatment_order
)

MINUTES_PER_HOUR <- 60
ZOOM_MINUTES     <- 3 * 24 * 60   # First 3 days = 4320 minutes

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
  dplyr::filter(Treatment != "BLANK", Minute <= ZOOM_MINUTES) |>  # ← zoom to 3 days
  remove_trailing_zeros() |>
  dplyr::mutate(
    Treatment = factor(Treatment, levels = treatment_order),
    Fly_Sex   = dplyr::recode(Fly_Sex, "F" = "Female", "M" = "Male"),
    Fly_Sex   = factor(Fly_Sex, levels = c("Female", "Male"))
  )

# Bin by HOUR instead of day
hourly_activity <- response_data_pre |>
  dplyr::mutate(
    Hour = ((Minute - 1) %/% MINUTES_PER_HOUR) + 1   # Hour 1, 2, ... 72
  ) |>
  dplyr::group_by(FlyID, Fly_Sex, Treatment, Hour) |>
  dplyr::summarise(
    total_activity   = sum(Activity, na.rm = TRUE),
    minutes_recorded = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    total_activity      = pmax(total_activity, 1e-6),
    activity_per_minute = total_activity / minutes_recorded
  )

max_hour <- max(hourly_activity$Hour, na.rm = TRUE)
cat("Zoom window:", max_hour, "hours (", max_hour / 24, "days)\n")

# Hourly means by Treatment and Sex
hourly_means <- hourly_activity |>
  dplyr::group_by(Treatment, Fly_Sex, Hour) |>
  dplyr::filter(dplyr::n() >= 2) |>
  dplyr::summarise(
    mean_activity = mean(total_activity, na.rm = TRUE),
    n_flies       = dplyr::n(),
    se_activity   = stats::sd(total_activity, na.rm = TRUE) / sqrt(n_flies),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    ci_upper = mean_activity + qt(0.975, df = n_flies - 1) * se_activity,
    ci_lower = pmax(mean_activity - qt(0.975, df = n_flies - 1) * se_activity, 0)
  )

########## PLOT - HOURLY ZOOM (FEMALE TOP / MALE BOTTOM) ##########

# Day boundary lines at hours 24 and 48
day_breaks <- c(24, 48)

Fig_hourly_activity <- ggplot2::ggplot(
  hourly_means,
  ggplot2::aes(x = Hour, y = mean_activity, color = Treatment, fill = Treatment)
) +
  # Day boundary vertical lines
  ggplot2::geom_vline(
    xintercept = day_breaks,
    linetype   = "dashed",
    color      = "gray60",
    linewidth  = 0.6
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
    alpha = 0.2, color = NA
  ) +
  ggplot2::geom_line(size = 1.0, alpha = 0.9) +
  ggplot2::geom_point(size = 1.8, alpha = 0.9) +
  ggplot2::facet_wrap(
    ~ Fly_Sex,
    ncol           = 1,
    strip.position = "top"
  ) +
  ggplot2::labs(
    x = "Hour (0–72 h)",
    y = "Mean total activity per hour (a.u.)"
  ) +
  ggplot2::scale_color_manual(values = my_color_palette) +
  ggplot2::scale_fill_manual(values = my_color_palette) +
  # X-axis: every 12 hours, label as day+hour for readability
  ggplot2::scale_x_continuous(
    breaks = seq(0, 72, by = 12),
    labels = c("0h", "12h", "24h\n(Day 2)", "36h", "48h\n(Day 3)", "60h", "72h")
  ) +
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

print(Fig_hourly_activity)
cowplot::save_plot("outputs/Figure_hourly_activity_3days.png", Fig_hourly_activity,
                   base_width = 12, base_height = 10)

########## SUMMARY STATISTICS ##########

cat("\n=== HOURLY ZOOM SUMMARY BY SEX ===\n")
hourly_summary <- hourly_activity |>
  dplyr::group_by(Treatment, Fly_Sex) |>
  dplyr::summarise(
    n_flies         = dplyr::n_distinct(FlyID),
    mean_activity   = mean(total_activity),
    median_activity = median(total_activity),
    .groups = "drop"
  )
print(hourly_summary)

########## DATA QUALITY CHECK ##########

cat("\n=== DATA QUALITY ===\n")
hourly_counts <- hourly_activity |>
  dplyr::group_by(Treatment, Fly_Sex, Hour) |>
  dplyr::summarise(n_flies = dplyr::n(), .groups = "drop")
print(hourly_counts |> head(20))

low_coverage <- hourly_counts |> dplyr::filter(n_flies < 3)
if (nrow(low_coverage) > 0) {
  cat("\nWARNING: Hours with <3 flies (sparse data):\n")
  print(low_coverage)
}



############### SIX HOURS BINS ########################

########## DATA INPUT & PREP - 6-HOUR BINNING (FIRST 3 DAYS) ##########

response_data_raw <- read.csv("data/flies_data.csv", sep = ";")
key               <- read.csv("data/flies_info.csv", sep = ";")

treatment_order <- c("naive", "PBS", "heatkill", "Pentomophila")
my_color_palette <- stats::setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7"),
  treatment_order
)

MINUTES_PER_HOUR <- 60
BIN_HOURS        <- 6
MINUTES_PER_BIN  <- BIN_HOURS * MINUTES_PER_HOUR   # 360 minutes
ZOOM_MINUTES     <- 3 * 24 * 60                     # First 3 days = 4320 minutes
N_BINS           <- ZOOM_MINUTES / MINUTES_PER_BIN  # 12 bins

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
  dplyr::filter(Treatment != "BLANK", Minute <= ZOOM_MINUTES) |>
  remove_trailing_zeros() |>
  dplyr::mutate(
    Treatment = factor(Treatment, levels = treatment_order),
    Fly_Sex   = dplyr::recode(Fly_Sex, "F" = "Female", "M" = "Male"),
    Fly_Sex   = factor(Fly_Sex, levels = c("Female", "Male"))
  )

# Bin by 6-HOUR windows
sixhour_activity <- response_data_pre |>
  dplyr::mutate(
    Bin = ((Minute - 1) %/% MINUTES_PER_BIN) + 1   # Bin 1, 2, ... 12
  ) |>
  dplyr::group_by(FlyID, Fly_Sex, Treatment, Bin) |>
  dplyr::summarise(
    total_activity   = sum(Activity, na.rm = TRUE),
    minutes_recorded = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    total_activity      = pmax(total_activity, 1e-6),
    activity_per_minute = total_activity / minutes_recorded
  )

max_bin <- max(sixhour_activity$Bin, na.rm = TRUE)
cat("Zoom window:", max_bin, "bins of", BIN_HOURS, "h (", max_bin * BIN_HOURS, "h /",
    max_bin * BIN_HOURS / 24, "days)\n")

# 6-hour means by Treatment and Sex
sixhour_means <- sixhour_activity |>
  dplyr::group_by(Treatment, Fly_Sex, Bin) |>
  dplyr::filter(dplyr::n() >= 2) |>
  dplyr::summarise(
    mean_activity = mean(total_activity, na.rm = TRUE),
    n_flies       = dplyr::n(),
    se_activity   = stats::sd(total_activity, na.rm = TRUE) / sqrt(n_flies),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    ci_upper = mean_activity + qt(0.975, df = n_flies - 1) * se_activity,
    ci_lower = pmax(mean_activity - qt(0.975, df = n_flies - 1) * se_activity, 0)
  )

########## PLOT - 6-HOUR BINS (FEMALE TOP / MALE BOTTOM) ##########

# Day boundaries fall at bins 4 and 8 (every 4 bins = 24 h)
day_breaks_bins <- c(4, 8)

# X-axis: one label per bin (bins 1–12 → midpoints at 3h, 9h, 15h, ...)
bin_labels <- paste0(seq(0, (N_BINS - 1) * BIN_HOURS, by = BIN_HOURS), "h")
bin_labels[c(1, 5, 9)] <- paste0(bin_labels[c(1, 5, 9)],
                                  c("", "\n(Day 2)", "\n(Day 3)"))

Fig_6h_activity <- ggplot2::ggplot(
  sixhour_means,
  ggplot2::aes(x = Bin, y = mean_activity, color = Treatment, fill = Treatment)
) +
  # Day boundary vertical lines
  ggplot2::geom_vline(
    xintercept = day_breaks_bins + 0.5,   # place between bins
    linetype   = "dashed",
    color      = "gray60",
    linewidth  = 0.6
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
    alpha = 0.2, color = NA
  ) +
  ggplot2::geom_line(linewidth = 1.0, alpha = 0.9) +
  ggplot2::geom_point(size = 2.2, alpha = 0.9) +
  ggplot2::facet_wrap(
    ~ Fly_Sex,
    ncol           = 1,
    strip.position = "top"
  ) +
  ggplot2::labs(
    x = "Time post infection (6-hour bins)",
    y = "Mean total activity per 6-hour bin (a.u.)"
  ) +
  ggplot2::scale_color_manual(values = my_color_palette) +
  ggplot2::scale_fill_manual(values = my_color_palette) +
  ggplot2::scale_x_continuous(
    breaks = 1:N_BINS,
    labels = bin_labels
  ) +
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

print(Fig_6h_activity)
cowplot::save_plot("outputs/Figure_6h_activity_3days.png", Fig_6h_activity,
                   base_width = 12, base_height = 10)

########## SUMMARY STATISTICS ##########

cat("\n=== 6-HOUR BIN SUMMARY BY SEX ===\n")
sixhour_summary <- sixhour_activity |>
  dplyr::group_by(Treatment, Fly_Sex) |>
  dplyr::summarise(
    n_flies         = dplyr::n_distinct(FlyID),
    mean_activity   = mean(total_activity),
    median_activity = median(total_activity),
    .groups = "drop"
  )
print(sixhour_summary)

########## DATA QUALITY CHECK ##########

cat("\n=== DATA QUALITY ===\n")
sixhour_counts <- sixhour_activity |>
  dplyr::group_by(Treatment, Fly_Sex, Bin) |>
  dplyr::summarise(n_flies = dplyr::n(), .groups = "drop")
print(sixhour_counts |> head(20))

low_coverage <- sixhour_counts |> dplyr::filter(n_flies < 3)
if (nrow(low_coverage) > 0) {
  cat("\nWARNING: Bins with <3 flies (sparse data):\n")
  print(low_coverage)
}
