########## DATA PREP ##########

response_data_raw <- read.csv("data/flies_data.csv", sep = ";")
key               <- read.csv("data/flies_info.csv", sep = ";")

treatment_order  <- c("naive", "PBS", "heatkill", "PentomophilaA", "PentomophilaB")
my_color_palette <- stats::setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"),
  treatment_order
)

MINUTES_PER_HOUR <- 60
ZOOM_MINUTES     <- 3 * 24 * 60
DAY1_MINUTES     <- 24 * 60

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

# ── Compute Last_mov from raw data ──
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

# ── Build dataset, split Pentomophila into A/B ──
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
      Treatment == "Pentomophila" & Last_mov <  12000 ~ "PentomophilaA",
      Treatment == "Pentomophila" & Last_mov >= 12000 ~ "PentomophilaB",
      TRUE ~ Treatment
    ),
    Treatment = factor(Treatment, levels = treatment_order),
    Fly_Sex   = dplyr::recode(Fly_Sex, "F" = "Female", "M" = "Male"),
    Fly_Sex   = factor(Fly_Sex, levels = c("Female", "Male")),
    State     = ifelse(Activity > 0, "Active", "Inactive")
  )

cat("Day 1 sample sizes by treatment:\n")
print(
  response_data_day1 |>
    dplyr::distinct(FlyID, Treatment, Fly_Sex) |>
    dplyr::count(Treatment, Fly_Sex)
)

########## TRANSITION MATRIX CALCULATION ##########

compute_transition_matrix <- function(states_vector) {
  if (length(states_vector) < 2) return(NULL)

  states <- c("Active", "Inactive")
  mat <- matrix(0, nrow = 2, ncol = 2, dimnames = list(states, states))

  for (i in 1:(length(states_vector) - 1)) {
    from <- states_vector[i]
    to   <- states_vector[i + 1]
    mat[from, to] <- mat[from, to] + 1
  }

  mat_prob <- mat / rowSums(mat)

  return(list(
    counts = mat,
    probs  = mat_prob,
    n_transitions = sum(mat)
  ))
}

########## AGGREGATE TRANSITION MATRICES BY TREATMENT+SEX ##########

transition_results <- response_data_day1 |>
  dplyr::group_by(FlyID, Treatment, Fly_Sex) |>
  dplyr::summarise(
    state_seq = list(State),
    n_minutes = dplyr::n(),
    .groups   = "drop"
  ) |>
  dplyr::mutate(
    trans_mat  = purrr::map(state_seq, compute_transition_matrix)
  ) |>
  dplyr::filter(!is.na(purrr::map_lgl(trans_mat, ~!is.null(.)))) |>
  dplyr::mutate(
    counts_mat = purrr::map(trans_mat, ~.$counts),
    probs_mat  = purrr::map(trans_mat, ~.$probs),
    n_trans    = purrr::map_dbl(trans_mat, ~.$n_transitions)
  )

########## EXTRACT PER-FLY PROBABILITIES & RATIOS ##########

per_fly_probs <- transition_results |>
  dplyr::mutate(
    active_to_active     = purrr::map_dbl(probs_mat, ~.[1, 1]),
    active_to_inactive   = purrr::map_dbl(probs_mat, ~.[1, 2]),
    inactive_to_active   = purrr::map_dbl(probs_mat, ~.[2, 1]),
    inactive_to_inactive = purrr::map_dbl(probs_mat, ~.[2, 2])
  ) |>
  dplyr::select(FlyID, Treatment, Fly_Sex,
                active_to_active, active_to_inactive,
                inactive_to_active, inactive_to_inactive) |>
  dplyr::mutate(
    persistence_ratio = active_to_active / inactive_to_inactive,
    on_off_ratio      = inactive_to_active / active_to_inactive
  )

########## AGGREGATE BY TREATMENT+SEX (COUNTS FOR HEATMAP & FLOW) ##########

aggregated_transitions <- transition_results |>
  dplyr::group_by(Treatment, Fly_Sex) |>
  dplyr::mutate(
    active_to_active     = purrr::map_dbl(counts_mat, ~.[1, 1]),
    active_to_inactive   = purrr::map_dbl(counts_mat, ~.[1, 2]),
    inactive_to_active   = purrr::map_dbl(counts_mat, ~.[2, 1]),
    inactive_to_inactive = purrr::map_dbl(counts_mat, ~.[2, 2])
  ) |>
  dplyr::summarise(
    n_flies                  = dplyr::n(),
    sum_active_to_active     = sum(active_to_active),
    sum_active_to_inactive   = sum(active_to_inactive),
    sum_inactive_to_active   = sum(inactive_to_active),
    sum_inactive_to_inactive = sum(inactive_to_inactive),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    p_active_to_active     = sum_active_to_active / (sum_active_to_active + sum_active_to_inactive),
    p_active_to_inactive   = sum_active_to_inactive / (sum_active_to_active + sum_active_to_inactive),
    p_inactive_to_active   = sum_inactive_to_active / (sum_inactive_to_active + sum_inactive_to_inactive),
    p_inactive_to_inactive = sum_inactive_to_inactive / (sum_inactive_to_active + sum_inactive_to_inactive),
    Treatment = factor(Treatment, levels = treatment_order)
  )

########## SUMMARY STATS FOR B & C (MEDIAN + IQR) ##########

persistence_summary <- per_fly_probs |>
  dplyr::filter(!is.infinite(persistence_ratio), !is.na(persistence_ratio)) |>
  dplyr::group_by(Treatment, Fly_Sex) |>
  dplyr::summarise(
    median_ratio = median(persistence_ratio, na.rm = TRUE),
    q1_ratio     = quantile(persistence_ratio, 0.25, na.rm = TRUE),
    q3_ratio     = quantile(persistence_ratio, 0.75, na.rm = TRUE),
    .groups      = "drop"
  ) |>
  dplyr::mutate(Treatment = factor(Treatment, levels = treatment_order))

on_off_summary <- per_fly_probs |>
  dplyr::filter(!is.infinite(on_off_ratio), !is.na(on_off_ratio)) |>
  dplyr::group_by(Treatment, Fly_Sex) |>
  dplyr::summarise(
    median_ratio = median(on_off_ratio, na.rm = TRUE),
    q1_ratio     = quantile(on_off_ratio, 0.25, na.rm = TRUE),
    q3_ratio     = quantile(on_off_ratio, 0.75, na.rm = TRUE),
    .groups      = "drop"
  ) |>
  dplyr::mutate(Treatment = factor(Treatment, levels = treatment_order))

########## PLOT A: HEATMAP ##########

heatmap_data <- aggregated_transitions |>
  dplyr::select(Treatment, Fly_Sex,
                p_active_to_active, p_active_to_inactive,
                p_inactive_to_active, p_inactive_to_inactive) |>
  tidyr::pivot_longer(
    cols      = starts_with("p_"),
    names_to  = "Transition",
    values_to = "Probability"
  ) |>
  dplyr::mutate(
    Transition_label = dplyr::recode(
      stringr::str_replace_all(Transition, "p_", ""),
      "active_to_active"     = "Active→Active",
      "active_to_inactive"   = "Active→Inactive",
      "inactive_to_active"   = "Inactive→Active",
      "inactive_to_inactive" = "Inactive→Inactive"
    ),
    Transition_label = factor(Transition_label,
                              levels = c("Active→Active", "Active→Inactive",
                                         "Inactive→Active", "Inactive→Inactive"))
  )

p_heatmap <- ggplot(heatmap_data,
                    aes(x = Treatment, y = Transition_label,
                        fill = Treatment, alpha = Probability)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = sprintf("%.3f", Probability), alpha = 1),
            size = 3.5, fontface = "bold", color = "black") +
  facet_wrap(~ Fly_Sex) +
  scale_fill_manual(values = my_color_palette, guide = "none") +
  scale_alpha_continuous(range = c(0.2, 1), guide = "none") +
  labs(
    x       = "Treatment",
    y       = "Transition",
    title   = "Markov Chain Transition Probabilities - Day 1",
    caption = "Color = treatment | Alpha = probability strength"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 10),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
    strip.text       = element_text(face = "bold", size = 11),
    panel.grid       = element_blank()
  )

########## PLOT B: PERSISTENCE RATIO - BAR + SCATTER + IQR ERRORBARS ##########

persistence_scatter <- per_fly_probs |>
  dplyr::filter(!is.infinite(persistence_ratio), !is.na(persistence_ratio)) |>
  dplyr::mutate(Treatment = factor(Treatment, levels = treatment_order))

p_persistence <- ggplot() +
  # Gray scatter (individual flies) — drawn first so bars appear on top
  geom_jitter(data = persistence_scatter,
              aes(x = Treatment, y = persistence_ratio),
              width = 0.2, size = 2, alpha = 0.4, color = "lightgray") +
  # Median bar
  geom_col(data = persistence_summary,
           aes(x = Treatment, y = median_ratio, fill = Treatment),
           color = "black", linewidth = 0.8, alpha = 0.85, width = 0.6) +
  # IQR error bars
  geom_errorbar(data = persistence_summary,
                aes(x = Treatment, ymin = q1_ratio, ymax = q3_ratio),
                width = 0.25, linewidth = 1.2, color = "black") +
  facet_wrap(~ Fly_Sex) +
  scale_fill_manual(values = my_color_palette, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x       = "Treatment",
    y       = "Active / Inactive Persistence Ratio",
    title   = "Persistence Ratio - Day 1",
    caption = "Bar = median | Error bars = IQR (Q1–Q3) | Gray dots = individual flies"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1),
    plot.title         = element_text(hjust = 0.5, face = "bold"),
    strip.text         = element_text(face = "bold", size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )

########## PLOT C: ON/OFF RATIO - BAR + SCATTER + IQR ERRORBARS ##########

on_off_scatter <- per_fly_probs |>
  dplyr::filter(!is.infinite(on_off_ratio), !is.na(on_off_ratio)) |>
  dplyr::mutate(Treatment = factor(Treatment, levels = treatment_order))

p_on_off <- ggplot() +
  # Gray scatter (individual flies) — drawn first so bars appear on top
  geom_jitter(data = on_off_scatter,
              aes(x = Treatment, y = on_off_ratio),
              width = 0.2, size = 2, alpha = 0.4, color = "lightgray") +
  # Median bar
  geom_col(data = on_off_summary,
           aes(x = Treatment, y = median_ratio, fill = Treatment),
           color = "black", linewidth = 0.8, alpha = 0.85, width = 0.6) +
  # IQR error bars
  geom_errorbar(data = on_off_summary,
                aes(x = Treatment, ymin = q1_ratio, ymax = q3_ratio),
                width = 0.25, linewidth = 1.2, color = "black") +
  facet_wrap(~ Fly_Sex) +
  scale_fill_manual(values = my_color_palette, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x       = "Treatment",
    y       = "Turn ON / Turn OFF Ratio",
    title   = "Activity Initiation vs Termination Ratio - Day 1",
    caption = "Bar = median | Error bars = IQR (Q1–Q3) | Gray dots = individual flies"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1),
    plot.title         = element_text(hjust = 0.5, face = "bold"),
    strip.text         = element_text(face = "bold", size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )

########## PLOT D: TRANSITION FLOW (NORMALIZED PER FLY) ##########

flow_normalized <- transition_results |>
  dplyr::mutate(
    active_to_active     = purrr::map_dbl(counts_mat, ~.[1, 1]),
    active_to_inactive   = purrr::map_dbl(counts_mat, ~.[1, 2]),
    inactive_to_active   = purrr::map_dbl(counts_mat, ~.[2, 1]),
    inactive_to_inactive = purrr::map_dbl(counts_mat, ~.[2, 2])
  ) |>
  dplyr::group_by(Treatment, Fly_Sex) |>
  dplyr::summarise(
    mean_active_to_active     = mean(active_to_active),
    mean_active_to_inactive   = mean(active_to_inactive),
    mean_inactive_to_active   = mean(inactive_to_active),
    mean_inactive_to_inactive = mean(inactive_to_inactive),
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(
    cols      = starts_with("mean_"),
    names_to  = "Transition",
    values_to = "Count"
  ) |>
  dplyr::mutate(
    From = dplyr::case_when(
      Transition %in% c("mean_active_to_active",
                        "mean_active_to_inactive")   ~ "Active",
      Transition %in% c("mean_inactive_to_active",
                        "mean_inactive_to_inactive") ~ "Inactive"
    ),
    To = dplyr::case_when(
      Transition %in% c("mean_active_to_active",
                        "mean_inactive_to_active")   ~ "Active",
      Transition %in% c("mean_active_to_inactive",
                        "mean_inactive_to_inactive") ~ "Inactive"
    ),
    Treatment = factor(Treatment, levels = treatment_order)
  )

p_flow <- ggplot(flow_normalized, aes(x = From, y = Count, fill = To)) +
  geom_col(position = "stack", color = "black", linewidth = 0.8) +
  facet_grid(Fly_Sex ~ Treatment) +
  scale_fill_manual(values = c("Active" = "#3498DB", "Inactive" = "#E74C3C"),
                    name = "Next State") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x       = "Current State",
    y       = "Mean Transitions per Fly",
    title   = "Transition Flow - Normalized per Fly - Day 1",
    caption = "Blue = transitions to Active | Red = transitions to Inactive"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(hjust = 0.5, face = "bold"),
    strip.text         = element_text(face = "bold", size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom"
  )

########## ASSEMBLE ALL PLOTS ##########

Fig_markov <- cowplot::plot_grid(
  p_heatmap,
  p_persistence,
  p_on_off,
  p_flow,
  labels      = c("A", "B", "C", "D"),
  label_size  = 14,
  ncol        = 2,
  rel_heights = c(1, 1)
)

print(Fig_markov)
cowplot::save_plot("outputs/Figure_markov_day1_v3.png", Fig_markov,
                   base_width = 16, base_height = 14)

########## EXPORT RESULTS ##########

write.csv(aggregated_transitions,
          "outputs/markov_transitions_day1.csv", row.names = FALSE)

write.csv(per_fly_probs,
          "outputs/markov_per_fly_ratios_day1.csv", row.names = FALSE)

cat("\nResults saved to:\n")
cat("  - outputs/markov_transitions_day1.csv\n")
cat("  - outputs/markov_per_fly_ratios_day1.csv\n")
cat("  - outputs/Figure_markov_day1_v3.png\n")
