########## DATA PREP ##########

response_data_raw <- read.csv("data/flies_data.csv", sep = ";")
key               <- read.csv("data/flies_info.csv", sep = ";")

treatment_order  <- c("naive", "injured", "immune stimulated",
                      "infected (acute)", "infected (chronic)")
my_color_palette <- stats::setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"),
  treatment_order
)

MINUTES_PER_HOUR <- 60
DAY1_MINUTES     <- 24 * 60
PENTO_THRESHOLD  <- 12000

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

response_data_day1 <- response_data_raw |>
  tidyr::pivot_longer(
    cols      = -c(Date, Time, Light1_Dark0),
    names_to  = "FlyID",
    values_to = "Activity"
  ) |>
  dplyr::mutate(FlyID = stringr::str_remove(FlyID, "^.")) |>
  dplyr::group_by(FlyID) |>
  dplyr::mutate(Minute = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::filter(Minute <= DAY1_MINUTES) |>
  dplyr::left_join(
    key |> dplyr::select(Monitor_TubeLoc, Fly_Sex, Treatment, Status_dead0_alive1),
    by = dplyr::join_by(FlyID == Monitor_TubeLoc)
  ) |>
  dplyr::left_join(last_mov_raw, by = "FlyID") |>
  dplyr::filter(Treatment != "BLANK") |>
  remove_trailing_zeros() |>
  dplyr::mutate(
    Treatment = dplyr::case_when(
      Treatment == "Pentomophila" & Last_mov <  PENTO_THRESHOLD ~ "infected (acute)",
      Treatment == "Pentomophila" & Last_mov >= PENTO_THRESHOLD ~ "infected (chronic)",
      Treatment == "PBS"                                         ~ "injured",
      Treatment == "heatkill"                                    ~ "immune stimulated",
      TRUE ~ Treatment
    ),
    Treatment = factor(Treatment, levels = treatment_order),
    Fly_Sex   = dplyr::recode(Fly_Sex, "F" = "Female", "M" = "Male"),
    Fly_Sex   = factor(Fly_Sex, levels = c("Female", "Male")),
    State     = ifelse(Activity > 0, "Active", "Inactive")
  )

cat("Day 1 sample sizes:\n")
print(
  response_data_day1 |>
    dplyr::distinct(FlyID, Treatment, Fly_Sex) |>
    dplyr::count(Treatment, Fly_Sex)
)

########## TRANSITION MATRICES ##########

compute_transition_matrix <- function(states_vector) {
  if (length(states_vector) < 2) return(NULL)
  states <- c("Active", "Inactive")
  mat    <- matrix(0, nrow = 2, ncol = 2, dimnames = list(states, states))
  for (i in 1:(length(states_vector) - 1)) {
    mat[states_vector[i], states_vector[i + 1]] <-
      mat[states_vector[i], states_vector[i + 1]] + 1
  }
  mat_prob <- mat / rowSums(mat)
  list(counts = mat, probs = mat_prob, n_transitions = sum(mat))
}

transition_results <- response_data_day1 |>
  dplyr::group_by(FlyID, Treatment, Fly_Sex) |>
  dplyr::summarise(state_seq = list(State), .groups = "drop") |>
  dplyr::mutate(trans_mat = purrr::map(state_seq, compute_transition_matrix)) |>
  dplyr::filter(purrr::map_lgl(trans_mat, ~!is.null(.))) |>
  dplyr::mutate(
    counts_mat = purrr::map(trans_mat, ~.$counts),
    probs_mat  = purrr::map(trans_mat, ~.$probs),
    n_trans    = purrr::map_dbl(trans_mat, ~.$n_transitions)
  )

per_fly_probs <- transition_results |>
  dplyr::mutate(
    active_to_active     = purrr::map_dbl(probs_mat, ~.[1, 1]),
    active_to_inactive   = purrr::map_dbl(probs_mat, ~.[1, 2]),
    inactive_to_active   = purrr::map_dbl(probs_mat, ~.[2, 1]),
    inactive_to_inactive = purrr::map_dbl(probs_mat, ~.[2, 2]),
    on_off_ratio         = inactive_to_active / active_to_inactive
  ) |>
  dplyr::select(FlyID, Treatment, Fly_Sex,
                active_to_active, active_to_inactive,
                inactive_to_active, inactive_to_inactive,
                on_off_ratio)

aggregated_transitions <- transition_results |>
  dplyr::mutate(
    active_to_active     = purrr::map_dbl(counts_mat, ~.[1, 1]),
    active_to_inactive   = purrr::map_dbl(counts_mat, ~.[1, 2]),
    inactive_to_active   = purrr::map_dbl(counts_mat, ~.[2, 1]),
    inactive_to_inactive = purrr::map_dbl(counts_mat, ~.[2, 2])
  ) |>
  dplyr::group_by(Treatment, Fly_Sex) |>
  dplyr::summarise(
    dplyr::across(c(active_to_active, active_to_inactive,
                    inactive_to_active, inactive_to_inactive), sum),
    n_flies = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    total_active   = active_to_active + active_to_inactive,
    total_inactive = inactive_to_active + inactive_to_inactive,
    p_aa = active_to_active     / total_active,
    p_ai = active_to_inactive   / total_active,
    p_ia = inactive_to_active   / total_inactive,
    p_ii = inactive_to_inactive / total_inactive
  )

########## PANEL A – HEATMAP ##########

heatmap_data <- aggregated_transitions |>
  tidyr::pivot_longer(
    cols      = c(p_aa, p_ai, p_ia, p_ii),
    names_to  = "Transition",
    values_to = "Probability"
  ) |>
  dplyr::mutate(
    Transition_label = dplyr::recode(Transition,
      p_aa = "Active → Active",
      p_ai = "Active → Inactive",
      p_ia = "Inactive → Active",
      p_ii = "Inactive → Inactive"
    )
  )

p_heatmap <- ggplot2::ggplot(
  heatmap_data,
  ggplot2::aes(x = Treatment, y = Transition_label,
               fill = Treatment, alpha = Probability)
) +
  ggplot2::geom_tile(color = "white", linewidth = 1.5) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.3f", Probability), alpha = 1),
    size = 2.5, color = "black"
  ) +
  ggplot2::facet_wrap(~ Fly_Sex) +
  ggplot2::scale_fill_manual(values = my_color_palette, guide = "none") +
  ggplot2::scale_alpha_continuous(range = c(0.2, 1), guide = "none") +
  ggplot2::labs(x = NULL, y = NULL, title = NULL) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    axis.text.x      = ggplot2::element_text(angle = 35, hjust = 1,
                                             size = 9, face = "bold"),
    strip.text       = ggplot2::element_text(face = "bold", size = 11),
    panel.grid       = ggplot2::element_blank()
  )

########## PANEL B – ON/OFF RATIO ##########

on_off_scatter <- per_fly_probs |>
  dplyr::filter(!is.infinite(on_off_ratio), !is.nan(on_off_ratio),
                !is.na(on_off_ratio))

on_off_summary <- on_off_scatter |>
  dplyr::group_by(Treatment, Fly_Sex) |>
  dplyr::summarise(
    median_val = median(on_off_ratio, na.rm = TRUE),
    q25        = quantile(on_off_ratio, 0.25, na.rm = TRUE),
    q75        = quantile(on_off_ratio, 0.75, na.rm = TRUE),
    .groups    = "drop"
  )

p_on_off <- ggplot2::ggplot() +
  ggplot2::geom_jitter(
    data  = on_off_scatter,
    ggplot2::aes(x = Treatment, y = on_off_ratio, color = Treatment),
    width = 0.2, alpha = 0.4, size = 1.5
  ) +
  ggplot2::geom_crossbar(
    data  = on_off_summary,
    ggplot2::aes(x = Treatment, y = median_val,
                 ymin = q25, ymax = q75, fill = Treatment),
    width = 0.5, color = "black", linewidth = 0.5
  ) +
  ggplot2::facet_wrap(~ Fly_Sex) +
  ggplot2::scale_color_manual(values = my_color_palette) +
  ggplot2::scale_fill_manual(values = my_color_palette) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(x = NULL, y = "On/Off ratio (log scale)", title = NULL) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    axis.text.x      = ggplot2::element_text(angle = 35, hjust = 1,
                                             size = 9, face = "bold"),
    strip.text       = ggplot2::element_text(face = "bold", size = 11),
    legend.position  = "none",
    panel.grid.minor = ggplot2::element_blank()
  )


########## ASSEMBLE FIGURE (A half-width | B full | C full) ##########
# Use a nested plot_grid: inner row holds A+B, then combine with C below.
# rel_widths c(1, 2) makes A half the width of B.

Fig_markov <- cowplot::plot_grid(
  p_heatmap,
  p_on_off,
  labels      = c("A", "B"),
  label_size  = 12,
  ncol        = 2,
  rel_widths  = c(6, 4)
)

print(Fig_markov)
cowplot::save_plot("outputs/Figure_markov_day1.png", Fig_markov,
                   base_width = 8, base_height =4)

########## PANEL D – POWER LAW (separate figure) ##########

bout_activity_raw <- response_data_day1 |>
  dplyr::group_by(FlyID, Treatment, Fly_Sex) |>
  dplyr::summarise(
    bout_totals = {
      r         <- rle(Activity > 0)
      is_active <- r$values
      lengths   <- r$lengths
      vals      <- Activity
      ends      <- cumsum(lengths)
      starts    <- ends - lengths + 1
      totals    <- mapply(function(s, e, active) {
        if (active) sum(vals[s:e]) else NA_real_
      }, starts, ends, is_active)
      list(totals[!is.na(totals)])
    },
    .groups = "drop"
  ) |>
  tidyr::unnest_longer(bout_totals) |>
  dplyr::filter(bout_totals > 0)

powerlaw_counts <- bout_activity_raw |>
  dplyr::group_by(Treatment, bout_totals) |>
  dplyr::summarise(n_bouts = dplyr::n(), .groups = "drop") |>
  dplyr::filter(bout_totals >= 1)

# Fit log-log OLS per treatment — use all rows so naive is never excluded
powerlaw_fits <- powerlaw_counts |>
  dplyr::group_by(Treatment) |>
  dplyr::summarise(
    fit     = list(lm(log10(n_bouts) ~ log10(bout_totals),
                      data = dplyr::pick(dplyr::everything()))),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    slope     = purrr::map_dbl(fit, ~coef(.)[2]),
    intercept = purrr::map_dbl(fit, ~coef(.)[1]),
    r2        = purrr::map_dbl(fit, ~summary(.)$r.squared),
    label     = sprintf("\u03b1 = %.2f\nR\u00b2 = %.2f", slope, r2)
  )

cat("\nPower law fits:\n")
print(dplyr::select(powerlaw_fits, Treatment, slope, intercept, r2))

# Predicted lines
pred_lines <- powerlaw_fits |>
  dplyr::mutate(
    x_seq = purrr::map(Treatment, ~{
      sub  <- dplyr::filter(powerlaw_counts, Treatment == .x)
      seq(log10(min(sub$bout_totals)), log10(max(sub$bout_totals)),
          length.out = 100)
    }),
    y_seq = purrr::map2(x_seq, fit, ~predict(
      .y, newdata = data.frame(bout_totals = 10^.x)))
  ) |>
  tidyr::unnest(cols = c(x_seq, y_seq))

# Label positions per treatment facet
label_pos <- powerlaw_counts |>
  dplyr::group_by(Treatment) |>
  dplyr::summarise(
    x_lab = 10^(log10(max(bout_totals)) * 0.5),
    y_lab = max(n_bouts),
    .groups = "drop"
  ) |>
  dplyr::left_join(dplyr::select(powerlaw_fits, Treatment, label),
                   by = "Treatment")

p_powerlaw <- ggplot2::ggplot(
  powerlaw_counts,
  ggplot2::aes(x = bout_totals, y = n_bouts, color = Treatment)
) +
  ggplot2::geom_point(alpha = 0.5, size = 1.4) +
  ggplot2::geom_line(
    data      = pred_lines,
    ggplot2::aes(x = 10^x_seq, y = 10^y_seq, color = Treatment),
    linewidth = 1.0, linetype = "dashed"
  ) +
  ggplot2::geom_text(
    data    = label_pos,
    ggplot2::aes(x = x_lab, y = y_lab, label = label, color = Treatment),
    size    = 3.2, hjust = 0, vjust = 1.2, fontface = "italic",
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~ Treatment, ncol = 3, scales = "free_y") +
  ggplot2::scale_x_log10(labels = scales::label_comma()) +
  ggplot2::scale_y_log10(labels = scales::label_comma()) +
  ggplot2::scale_color_manual(values = my_color_palette) +
  ggplot2::labs(
    x       = "Total activity per active bout (a.u., log scale)",
    y       = "Number of bouts (log scale)",
    title   = NULL,
    caption = "Dashed line = log-log OLS fit  |  \u03b1 = slope"
  ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    legend.position  = "none",
    strip.text       = ggplot2::element_text(face = "bold", size = 11),
    panel.grid.minor = ggplot2::element_blank(),
    plot.caption     = ggplot2::element_text(size = 9, color = "gray50")
  )

print(p_powerlaw)
cowplot::save_plot("outputs/Figure_powerlaw_day1.png", p_powerlaw,
                   base_width = 14, base_height = 9)

########## EXPORT ##########

write.csv(aggregated_transitions,
          "outputs/markov_transitions_day1.csv",  row.names = FALSE)
write.csv(per_fly_probs,
          "outputs/per_fly_transition_probs_day1.csv", row.names = FALSE)
write.csv(powerlaw_counts,
          "outputs/powerlaw_bout_activity_day1.csv", row.names = FALSE)
write.csv(dplyr::select(powerlaw_fits, -fit),
          "outputs/powerlaw_fits_day1.csv", row.names = FALSE)
