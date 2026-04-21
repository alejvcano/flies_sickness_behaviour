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
endtime         <- 40000

########## SHARED DATA PREP ##########

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
  dplyr::group_by(FlyID) |>
  dplyr::mutate(Last_mov = length(Activity)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    Treatment = factor(Treatment, levels = treatment_order),
    Fly_Sex   = dplyr::recode(Fly_Sex, "F" = "Female", "M" = "Male"),
    Fly_Sex   = factor(Fly_Sex, levels = c("Female", "Male")),
    class     = ifelse(Status_dead0_alive1 == 0, "dead", "alive")
  )

########## PANEL A – DAILY ACTIVITY (FULL EXPERIMENT) ##########

daily_activity <- response_data_pre |>
  dplyr::mutate(Day = ((Minute - 1) %/% MINUTES_PER_DAY) + 1) |>
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

daily_means <- daily_activity |>
  dplyr::group_by(Treatment, Fly_Sex, Day) |>
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

panel_A <- ggplot2::ggplot(
  daily_means,
  ggplot2::aes(x = Day, y = mean_activity, color = Treatment, fill = Treatment)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
    alpha = 0.2, color = NA
  ) +
  ggplot2::geom_line(linewidth = 1.2, alpha = 0.9) +
  ggplot2::geom_point(size = 2.5, alpha = 0.9) +
  ggplot2::facet_wrap(~ Fly_Sex, ncol = 1, strip.position = "top") +
  ggplot2::labs(
    x = "Days post-pricking",
    y = "Mean total activity per day (a.u.)"
  ) +
  ggplot2::scale_color_manual(values = my_color_palette) +
  ggplot2::scale_fill_manual(values = my_color_palette) +
  ggplot2::scale_x_continuous(breaks = seq(1, max_day, by = 2)) +
  ggplot2::scale_y_log10(labels = scales::scientific_format()) +
  ggplot2::theme_minimal(base_size = 15) +   # base +2
  ggplot2::theme(
    legend.position    = "none",
    axis.text.x        = ggplot2::element_text(size = 12),
    axis.text.y        = ggplot2::element_text(size = 12),
    axis.title         = ggplot2::element_text(size = 14),
    strip.text         = ggplot2::element_text(face = "bold", size = 14),
    panel.grid.minor   = ggplot2::element_blank(),
    plot.margin        = ggplot2::margin(5, 15, 5, 5, "pt")
  )

########## PANEL B – SURVIVAL CURVES ##########

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

panel_B <- ggplot2::ggplot(
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
    x = "Hours post-pricking",
    y = "Fraction alive"
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_minimal(base_size = 15) +   # base +2
  ggplot2::theme(
    legend.position  = "none",
    axis.text        = ggplot2::element_text(size = 12),
    axis.title       = ggplot2::element_text(size = 14),
    strip.text       = ggplot2::element_text(face = "bold", size = 14),
    panel.spacing    = ggplot2::unit(1.2, "lines"),
    plot.margin      = ggplot2::margin(5, 15, 5, 5, "pt")
  )

########## SHARED LEGEND — built as a bare grob, no panel/line ##########

legend_grob <- ggplot2::ggplot(
  data = data.frame(
    Treatment = factor(treatment_order, levels = treatment_order),
    x = 1, y = 1
  ),
  ggplot2::aes(x = x, y = y, color = Treatment)
) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = my_color_palette) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      nrow      = 1,
      keywidth  = ggplot2::unit(1.8, "cm"),
      keyheight = ggplot2::unit(0.4, "cm")
    )
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    legend.position       = "bottom",
    legend.title          = ggplot2::element_blank(),
    legend.text           = ggplot2::element_text(size = 14),   # +2
    legend.background     = ggplot2::element_rect(fill = NA, color = NA),
    legend.box.background = ggplot2::element_rect(fill = NA, color = NA),
    legend.box.margin     = ggplot2::margin(0, 0, 0, 0),
    legend.margin         = ggplot2::margin(0, 0, 0, 0),
    plot.background       = ggplot2::element_rect(fill = NA, color = NA),
    panel.background      = ggplot2::element_rect(fill = NA, color = NA)
  )

shared_legend <- cowplot::get_legend(legend_grob)

########## ASSEMBLE ##########

top_row <- cowplot::plot_grid(
  panel_A,
  panel_B,
  labels     = c("A", "B"),
  label_size = 20,
  ncol       = 2,
  rel_widths = c(2, 1),
  align      = "hv",
  axis       = "tb"
)

Fig_final <- cowplot::plot_grid(
  top_row,
  shared_legend,
  ncol        = 1,
  rel_heights = c(1, 0.05)
)

cowplot::save_plot(
  "outputs/Figure_activity_survival.png",
  Fig_final,
  base_width  = 12,
  base_height = 8
)
