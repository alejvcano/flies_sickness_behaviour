########## DATA INPUT & PREP ##########

response_data_raw <- read.csv("data/flies_data.csv", sep = ";")
key               <- read.csv("data/flies_info.csv", sep = ";")

treatment_order <- c("naive", "infected (acute)", "infected (chronic)")

my_color_palette <- stats::setNames(
  c("#0072B2", "#CC79A7", "#E69F00"),
  treatment_order
)

MINUTES_PER_HOUR <- 60
WINDOW_HOURS     <- 6
WINDOW_MINUTES   <- WINDOW_HOURS * MINUTES_PER_HOUR
ZOOM_MINUTES     <- 3 * 24 * 60
PENTO_THRESHOLD  <- 12000

# ── Step 1: compute Last_mov from RAW data (before any trimming) ──
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

cat("Last_mov summary:\n"); print(summary(last_mov_raw$Last_mov))
cat("infected (acute)  (< ", PENTO_THRESHOLD, "):", sum(last_mov_raw$Last_mov <  PENTO_THRESHOLD), "\n")
cat("infected (chronic)(>=", PENTO_THRESHOLD, "):", sum(last_mov_raw$Last_mov >= PENTO_THRESHOLD), "\n")

# ── Step 2: single pivot + join, recode BEFORE filter ──
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
      Treatment == "naive"                                       ~ "naive",
      TRUE ~ Treatment
    )
  ) |>
  dplyr::filter(Treatment %in% treatment_order, Minute <= ZOOM_MINUTES) |>
  remove_trailing_zeros() |>
  dplyr::mutate(
    Treatment    = factor(Treatment, levels = treatment_order),
    Fly_Sex      = dplyr::recode(Fly_Sex, "F" = "Female", "M" = "Male"),
    Fly_Sex      = factor(Fly_Sex, levels = c("Female", "Male")),
    Window       = ((Minute - 1) %/% WINDOW_MINUTES),
    Window_label = sprintf("h%02d-%02d",
                           Window * WINDOW_HOURS + 1,
                           (Window + 1) * WINDOW_HOURS)
  )

# ── Verify all three groups present ──
final_counts <- dplyr::distinct(response_data_pre, FlyID, Treatment, Fly_Sex) |>
  dplyr::count(Treatment, Fly_Sex)
cat("\nFinal treatment × sex counts:\n"); print(final_counts)
stopifnot(all(treatment_order %in% levels(droplevels(response_data_pre$Treatment))))

window_levels <- unique(response_data_pre$Window_label[order(response_data_pre$Window)])
response_data_pre <- response_data_pre |>
  dplyr::mutate(Window_label = factor(Window_label, levels = window_levels))

########## BOUT STRUCTURE PER FLY PER 6H WINDOW ##########

bout_features <- response_data_pre |>
  dplyr::group_by(FlyID, Fly_Sex, Treatment, Window, Window_label) |>
  dplyr::summarise(
    n_active_bouts     = sum(rle(Activity > 0)$values),
    n_inactive_bouts   = sum(!rle(Activity > 0)$values),
    mean_active_bout   = ifelse(
      sum(rle(Activity > 0)$values) > 0,
      mean(rle(Activity > 0)$lengths[rle(Activity > 0)$values]),
      NA_real_
    ),
    mean_inactive_bout = ifelse(
      sum(!rle(Activity > 0)$values) > 0,
      mean(rle(Activity > 0)$lengths[!rle(Activity > 0)$values]),
      NA_real_
    ),
    total_activity     = sum(Activity, na.rm = TRUE),
    .groups = "drop"
  )

########## MEDIAN SUMMARY BY TREATMENT, SEX, WINDOW ##########

make_window_summary <- function(df, metric) {
  df |>
    dplyr::group_by(Treatment, Fly_Sex, Window, Window_label) |>
    dplyr::filter(dplyr::n() >= 2, !is.na(!!rlang::sym(metric))) |>
    dplyr::summarise(
      median_val = median(!!rlang::sym(metric), na.rm = TRUE),
      n_flies    = dplyr::n(),
      spread     = IQR(!!rlang::sym(metric), na.rm = TRUE) / 2,
      .groups    = "drop"
    ) |>
    dplyr::mutate(
      ci_upper = median_val + spread,
      ci_lower = pmax(median_val - spread, 0)
    )
}

summary_active_bout   <- make_window_summary(bout_features, "mean_active_bout")
summary_inactive_bout <- make_window_summary(bout_features, "mean_inactive_bout")
summary_n_bouts       <- make_window_summary(bout_features, "n_active_bouts")

########## PLOT FUNCTION ##########

make_bout_plot <- function(summary_df, y_label) {
  all_labels  <- levels(summary_df$Window_label)
  show_labels <- ifelse(seq_along(all_labels) %% 2 == 1, all_labels, "")

  ggplot2::ggplot(
    summary_df,
    ggplot2::aes(x = Window_label, y = median_val,
                 color = Treatment, fill = Treatment, group = Treatment)
  ) +
    ggplot2::geom_vline(
      xintercept = c(4.5, 8.5),
      linetype   = "dashed",
      color      = "gray60",
      linewidth  = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      alpha       = 0.2,
      color       = NA,
      show.legend = FALSE          # removes ribbon key from legend
    ) +
    ggplot2::geom_line(linewidth = 1.2, alpha = 0.9) +
    ggplot2::geom_point(size = 2.5, alpha = 0.9) +
    ggplot2::facet_wrap(~ Fly_Sex, ncol = 1, strip.position = "top") +
    ggplot2::labs(
      x = "6-hour window",
      y = y_label
    ) +
    ggplot2::scale_color_manual(values = my_color_palette) +
    ggplot2::scale_fill_manual(values = my_color_palette) +
    ggplot2::scale_x_discrete(labels = show_labels) +
    ggplot2::theme_minimal(base_size = 18) +    # increased from 17
    ggplot2::theme(
      legend.position  = "none",
      axis.text.x      = ggplot2::element_text(size = 13, angle = 45, hjust = 1),
      axis.text.y      = ggplot2::element_text(size = 14),
      axis.title       = ggplot2::element_text(size = 16),
      strip.text       = ggplot2::element_text(face = "bold", size = 16),
      panel.grid.minor = ggplot2::element_blank()
    )
}

########## GENERATE THREE PLOTS ##########

p_n_bouts       <- make_bout_plot(summary_n_bouts,       "Median number of active bouts per 6h")
p_active_bout   <- make_bout_plot(summary_active_bout,   "Median active bout length (min)")
p_inactive_bout <- make_bout_plot(summary_inactive_bout, "Median inactive bout length (min)")

########## EXTRACT SHARED LEGEND ##########

legend_plot <- make_bout_plot(summary_n_bouts, "") +
  ggplot2::theme(
    legend.position       = "bottom",
    legend.title          = ggplot2::element_blank(),
    legend.text           = ggplot2::element_text(size = 15),
    legend.key.width      = ggplot2::unit(1.5, "cm"),
    legend.background     = ggplot2::element_rect(fill = "white", color = NA),
    legend.box.background = ggplot2::element_rect(fill = "white", color = NA)
  )

shared_legend <- cowplot::get_legend(legend_plot)

plots_row <- cowplot::plot_grid(
  p_n_bouts,
  p_active_bout,
  p_inactive_bout,
  labels     = c("A", "B", "C"),
  label_size = 18,               
  ncol       = 3
)

Fig_bout_structure <- cowplot::plot_grid(
  plots_row,
  shared_legend,
  ncol        = 1,
  rel_heights = c(1, 0.06)
)

print(Fig_bout_structure)
cowplot::save_plot("outputs/Figure_bout_structure_PentomophilaAB.png", Fig_bout_structure,
                   base_width = 18, base_height = 12)

########## SUMMARY STATISTICS ##########

cat("\n=== BOUT STRUCTURE SUMMARY (MEDIAN) BY GROUP AND SEX ===\n")
bout_summary <- bout_features |>
  dplyr::group_by(Treatment, Fly_Sex) |>
  dplyr::summarise(
    n_flies              = dplyr::n_distinct(FlyID),
    median_active_bout   = median(mean_active_bout,   na.rm = TRUE),
    median_inactive_bout = median(mean_inactive_bout, na.rm = TRUE),
    median_n_bouts       = median(n_active_bouts,     na.rm = TRUE),
    iqr_active_bout      = IQR(mean_active_bout,      na.rm = TRUE),
    iqr_inactive_bout    = IQR(mean_inactive_bout,    na.rm = TRUE),
    iqr_n_bouts          = IQR(n_active_bouts,        na.rm = TRUE),
    .groups = "drop"
  )
print(bout_summary)

########## DATA QUALITY CHECK ##########

cat("\n=== DATA QUALITY ===\n")
bout_counts <- bout_features |>
  dplyr::group_by(Treatment, Fly_Sex, Window_label) |>
  dplyr::summarise(n_flies = dplyr::n(), .groups = "drop")

low_coverage <- bout_counts |> dplyr::filter(n_flies < 3)
if (nrow(low_coverage) > 0) {
  cat("\nWARNING: Windows with <3 flies (sparse data):\n")
  print(low_coverage)
}

########## TREATMENT × WINDOW INTERACTION: ACTIVE BOUT NUMBER ##########

lmm_nbouts_time <- lmerTest::lmer(
  n_active_bouts ~ Treatment * Window * Fly_Sex + (1 | FlyID),
  data    = bout_features |> dplyr::mutate(Window = as.numeric(Window)),
  REML    = TRUE,
  control = lme4::lmerControl(optimizer = "bobyqa")
)

cat("\n===== TYPE III ANOVA: n_active_bouts ~ Treatment * Window * Sex =====\n")
print(anova(lmm_nbouts_time, type = "III"))

########## TARGETED CONTRASTS: infected (acute) vs others ##########

emm_nbouts_time <- emmeans::emmeans(
  lmm_nbouts_time,
  ~ Treatment | Fly_Sex * Window,
  at = list(Window = 0:11)
)

nbouts_targeted <- emmeans::contrast(
  emm_nbouts_time,
  method = list(
    "acute - naive"   = c(-1,  1,  0),
    "acute - chronic" = c( 0,  1, -1)
  ),
  adjust = "BH"
)

cat("\n===== TARGETED CONTRASTS: infected (acute) vs naive/chronic, BY SEX AND WINDOW =====\n")
print(nbouts_targeted)

nbouts_targeted_df <- as.data.frame(nbouts_targeted)
print(nbouts_targeted_df)
