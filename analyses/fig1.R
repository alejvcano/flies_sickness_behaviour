########## TREATMENT LABEL REMAP ##########

treatment_remap <- c(
  "naive"        = "Naive",
  "PBS"          = "Injured",
  "heatkill"     = "Immune stimulated",
  "Pentomophila" = "Infected"
)

treatment_order_new <- c("Naive", "Injured", "Immune stimulated", "Infected")

my_color_palette <- stats::setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7"),
  treatment_order_new
)

########## DATA INPUT & PREP ##########

response_data_raw <- read.csv("data/flies_data.csv", sep = ";")
key               <- read.csv("data/flies_info.csv", sep = ";")

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
  dplyr::filter(Treatment != "BLANK") |>
  dplyr::mutate(Treatment = dplyr::recode(Treatment, !!!treatment_remap)) |>
  remove_trailing_zeros() |>
  dplyr::group_by(FlyID) |>
  dplyr::mutate(Last_mov = length(Activity)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    Treatment = factor(Treatment, levels = treatment_order_new),
    class     = ifelse(Status_dead0_alive1 == 0, "dead", "alive")
  )

########## PANEL A – ACTIVITY TIME SERIES  ##########

fly_data_hourly <- response_data_pre |>
  dplyr::rename(id = "FlyID") |>
  dplyr::mutate(hour = floor(Minute / 60)) |>
  dplyr::group_by(id, Treatment, hour) |>
  dplyr::summarise(Activity = mean(Activity, na.rm = TRUE), .groups = "drop")

mean_data_hourly <- fly_data_hourly |>
  dplyr::group_by(Treatment, hour) |>
  dplyr::summarise(meanActivity = mean(Activity, na.rm = TRUE), .groups = "drop")

n_labels <- response_data_pre |>
  dplyr::distinct(FlyID, Treatment, Fly_Sex) |>
  dplyr::count(Treatment, Fly_Sex, name = "n") |>
  tidyr::pivot_wider(names_from = Fly_Sex, values_from = n, values_fill = 0) |>
  dplyr::mutate(
    label     = paste0("F = ", F, "   M = ", M),
    Treatment = factor(Treatment, levels = treatment_order_new)
  )

activitytimeseries <- ggplot2::ggplot() +
  ggplot2::geom_line(
    data  = fly_data_hourly,
    ggplot2::aes(x = hour, y = Activity, group = id, color = Treatment),
    alpha = 0.08
  ) +
  ggplot2::geom_line(
    data      = mean_data_hourly,
    ggplot2::aes(x = hour, y = meanActivity, color = Treatment),
    linewidth = 1.4
  ) +
  ggplot2::geom_text(
    data     = n_labels,
    ggplot2::aes(label = label, color = Treatment),
    x        = Inf,
    y        = Inf,
    hjust    = 1.05,
    vjust    = 1.8,
    size     = 4.0,
    fontface = "italic"
  ) +
  ggplot2::facet_wrap(~ Treatment, ncol = 2) +   
    expand = ggplot2::expansion(add = c(0, 10)),
    breaks = seq(0, max(fly_data_hourly$hour, na.rm = TRUE), by = 100)
  ) +
  ggplot2::scale_y_continuous(limits = c(0, 20)) +
  ggplot2::scale_color_manual(values = my_color_palette) +
  ggplot2::labs(
    x = "Time post-pricking (hours)",  
    y = "Mean activity (a.u.)"
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "none",
    axis.text       = ggplot2::element_text(size = 11),
    axis.title      = ggplot2::element_text(size = 14),
    strip.text      = ggplot2::element_text(face = "bold", size = 15),
    panel.spacing   = ggplot2::unit(1.0, "lines"),
    plot.margin     = ggplot2::margin(5, 15, 5, 5, "pt")
  )

cowplot::save_plot(
  "outputs/FigA_activity_timeseries.png",
  activitytimeseries,
  base_width  = 8.4,    # <-- swapped: narrow width (left half of slide)
  base_height = 6.6     # <-- swapped: tall height for 4 stacked panels
)
