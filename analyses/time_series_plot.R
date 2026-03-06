########## DATA INPUT & PREP ##########

response_data_raw <- read.csv("data/flies_data.csv", sep = ";")
key               <- read.csv("data/flies_info.csv", sep = ";")

treatment_order <- c("naive", "PBS", "heatkill", "Pentomophila")

my_color_palette <- stats::setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7"),
  treatment_order
)

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
  remove_trailing_zeros() |>
  dplyr::group_by(FlyID) |>
  dplyr::mutate(Last_mov = length(Activity)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    Treatment = factor(Treatment, levels = treatment_order),
    class     = ifelse(Status_dead0_alive1 == 0, "dead", "alive")
  )

########## PANEL A – ACTIVITY TIME SERIES (2x2 grid) ##########

fly_data_hourly <- response_data_pre |>
  dplyr::rename(id = "FlyID") |>
  dplyr::mutate(hour = floor(Minute / 60)) |>
  dplyr::group_by(id, Treatment, hour) |>
  dplyr::summarise(Activity = mean(Activity, na.rm = TRUE), .groups = "drop")

mean_data_hourly <- fly_data_hourly |>
  dplyr::group_by(Treatment, hour) |>
  dplyr::summarise(meanActivity = mean(Activity, na.rm = TRUE), .groups = "drop")

# Count flies per treatment per sex, then build "F=x  M=y" label
n_labels <- response_data_pre |>
  dplyr::distinct(FlyID, Treatment, Fly_Sex) |>
  dplyr::count(Treatment, Fly_Sex, name = "n") |>
  tidyr::pivot_wider(names_from = Fly_Sex, values_from = n, values_fill = 0) |>
  dplyr::mutate(
    label     = paste0("F = ", F, "   M = ", M),
    Treatment = factor(Treatment, levels = treatment_order)
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
    size     = 3.8,
    fontface = "italic"
  ) +
  ggplot2::facet_wrap(~ Treatment, ncol = 2) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(add = c(0, 10)),
    breaks = seq(0, max(fly_data_hourly$hour, na.rm = TRUE), by = 100)
  ) +
  ggplot2::scale_y_continuous(limits = c(0, 20)) +
  ggplot2::scale_color_manual(values = my_color_palette) +
  ggplot2::labs(x = "Hours post infection", y = "Mean activity (a.u.)") +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "none",
    axis.text       = ggplot2::element_text(size = 11),
    axis.title      = ggplot2::element_text(size = 13),
    strip.text      = ggplot2::element_text(face = "bold", size = 12),
    panel.spacing   = ggplot2::unit(1.2, "lines"),
    plot.margin     = ggplot2::margin(5, 15, 5, 5, "pt")
  )

########## PANEL B – SURVIVAL CURVE split by sex ##########

total_per_treatment_sex <- response_data_pre |>
  dplyr::distinct(FlyID, Treatment, Fly_Sex) |>
  dplyr::count(Treatment, Fly_Sex, name = "total")

max_hour <- floor(max(response_data_pre$Last_mov, na.rm = TRUE) / 60)

fly_status <- response_data_pre |>
  dplyr::distinct(FlyID, Treatment, Fly_Sex, class, Last_mov) |>
  dplyr::mutate(
    death_hour = dplyr::if_else(class == "dead", floor(Last_mov / 60), Inf)
  )

survival_curve <- tidyr::expand_grid(
  Treatment = unique(fly_status$Treatment),
  Fly_Sex   = unique(fly_status$Fly_Sex),
  hour      = 0:max_hour
) |>
  dplyr::left_join(fly_status, by = c("Treatment", "Fly_Sex"),
                   relationship = "many-to-many") |>
  dplyr::group_by(Treatment, Fly_Sex, hour) |>
  dplyr::summarise(
    n_alive = sum(death_hour > hour, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::left_join(total_per_treatment_sex, by = c("Treatment", "Fly_Sex")) |>
  dplyr::mutate(
    fraction_alive = n_alive / total,
    Treatment      = factor(Treatment, levels = treatment_order)
  )

survivalplot <- ggplot2::ggplot(
  survival_curve,
  ggplot2::aes(x = hour, y = fraction_alive, color = Treatment)
) +
  ggplot2::geom_step(linewidth = 1.2) +
  ggplot2::facet_wrap(~ Fly_Sex, ncol = 1) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(add = c(0, 10)),
    breaks = seq(0, max_hour, by = 100)
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0.5, 1),
    breaks = seq(0.5, 1, by = 0.1),
    oob    = scales::squish
  ) +
  ggplot2::scale_color_manual(values = my_color_palette) +
  ggplot2::labs(
    x = "Hours post infection",
    y = "Fraction alive"
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "none",
    axis.text       = ggplot2::element_text(size = 11),
    axis.title      = ggplot2::element_text(size = 13),
    strip.text      = ggplot2::element_text(face = "bold", size = 12),
    panel.spacing   = ggplot2::unit(1.2, "lines"),
    plot.margin     = ggplot2::margin(5, 15, 5, 5, "pt")
  )

########## ASSEMBLE FIGURE 1 (A = 2/3, B = 1/3) ##########

Fig1 <- cowplot::plot_grid(
  activitytimeseries,
  survivalplot,
  labels      = c("A", "B"),
  label_size  = 18,
  ncol        = 2,
  rel_widths  = c(2, 1),
  align       = "hv",
  axis        = "tb"
)

cowplot::save_plot(
  "outputs/Figure1_immune.png",
  Fig1,
  base_width  = 12,
  base_height = 7                                        # <-- increased for B's two panels
)
