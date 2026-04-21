########## DATA INPUT & PREP ##########

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

MINUTES_PER_HOUR <- 60
WINDOW_HOURS     <- 6
WINDOW_MINUTES   <- WINDOW_HOURS * MINUTES_PER_HOUR
ZOOM_MINUTES     <- 3 * 24 * 60

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
  dplyr::mutate(Treatment = dplyr::recode(Treatment, !!!treatment_remap)) |>
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
    ggplot2::geom_line(linewidth = 1.0, alpha = 0.9) +
    ggplot2::geom_point(size = 2.0, alpha = 0.9) +
    ggplot2::facet_wrap(~ Fly_Sex, ncol = 1, strip.position = "top") +
    ggplot2::labs(
      x = "6-hour window",
      y = y_label
    ) +
    ggplot2::scale_color_manual(values = my_color_palette) +
    ggplot2::scale_fill_manual(values = my_color_palette) +
    ggplot2::scale_x_discrete(labels = show_labels) +
    ggplot2::theme_minimal(base_size = 18) +    # increased from 15
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
    legend.position  = "bottom",
    legend.title     = ggplot2::element_blank(),
    legend.text      = ggplot2::element_text(size = 15),
    legend.key.width = ggplot2::unit(1.5, "cm")
  )

shared_legend <- cowplot::get_legend(legend_plot)

########## ASSEMBLE INTO ONE FIGURE ##########

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
cowplot::save_plot("outputs/Figure_bout_structure.png", Fig_bout_structure,
                   base_width = 18, base_height = 10)

########## SUMMARY STATISTICS ##########

cat("\n=== BOUT STRUCTURE SUMMARY (MEDIAN) BY TREATMENT AND SEX ===\n")
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

########## [PV1] SEX DIFFERENCE IN BOUT RHYTHM REGULARITY ##########

bout_cv <- bout_features |>
  dplyr::group_by(FlyID, Fly_Sex, Treatment) |>
  dplyr::summarise(
    mean_bouts = mean(n_active_bouts, na.rm = TRUE),
    sd_bouts   = sd(n_active_bouts,   na.rm = TRUE),
    cv_bouts   = sd_bouts / mean_bouts,
    .groups = "drop"
  )

cv_summary <- bout_cv |>
  dplyr::group_by(Fly_Sex, Treatment) |>
  dplyr::summarise(
    n_flies    = dplyr::n(),
    median_cv  = median(cv_bouts, na.rm = TRUE),
    iqr_cv     = IQR(cv_bouts,    na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== [PV1] BOUT RHYTHM REGULARITY (CV) BY SEX AND TREATMENT ===\n")
print(cv_summary)

wilcox_sex_cv <- wilcox.test(
  cv_bouts ~ Fly_Sex,
  data  = bout_cv,
  exact = FALSE
)
cat("\nWilcoxon test — CV of active bouts, Female vs Male (all treatments pooled):\n")
print(wilcox_sex_cv)

cat("\nPer-treatment Wilcoxon tests (Female vs Male):\n")
bout_cv |>
  dplyr::group_by(Treatment) |>
  dplyr::summarise(
    p_value = tryCatch(
      wilcox.test(cv_bouts ~ Fly_Sex, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) |>
  print()

########## [PV2] NAIVE vs PRICKING GROUPS — ACTIVE BOUT LENGTH ##########

cat("\n=== [PV2] ACTIVE BOUT LENGTH: NAIVE vs OTHER TREATMENTS ===\n")

bout_active_flat <- bout_features |>
  dplyr::filter(!is.na(mean_active_bout))

kw_results <- bout_active_flat |>
  dplyr::group_by(Fly_Sex) |>
  dplyr::summarise(
    kw_stat = kruskal.test(mean_active_bout ~ Treatment)$statistic,
    kw_p    = kruskal.test(mean_active_bout ~ Treatment)$p.value,
    .groups = "drop"
  )
cat("\nKruskal-Wallis test — active bout length by treatment:\n")
print(kw_results)

cat("\nDunn post-hoc tests (BH-corrected):\n")
bout_active_flat |>
  dplyr::group_by(Fly_Sex) |>
  dplyr::group_map(~ {
    dunn <- FSA::dunnTest(mean_active_bout ~ Treatment,
                          data   = .x,
                          method = "bh")
    dunn$res |>
      dplyr::mutate(Fly_Sex = .y$Fly_Sex) |>
      dplyr::filter(grepl("naive", Comparison))
  }) |>
  dplyr::bind_rows() |>
  dplyr::select(Fly_Sex, Comparison, Z, P.unadj, P.adj) |>
  print()

cat("\nMedian fold-difference in active bout length (naive / other group):\n")
naive_medians <- bout_features |>
  dplyr::filter(Treatment == "naive", !is.na(mean_active_bout)) |>
  dplyr::group_by(Fly_Sex) |>
  dplyr::summarise(naive_med = median(mean_active_bout), .groups = "drop")

bout_features |>
  dplyr::filter(Treatment != "naive", !is.na(mean_active_bout)) |>
  dplyr::group_by(Fly_Sex, Treatment) |>
  dplyr::summarise(group_med = median(mean_active_bout), .groups = "drop") |>
  dplyr::left_join(naive_medians, by = "Fly_Sex") |>
  dplyr::mutate(fold_diff = round(naive_med / group_med, 2)) |>
  dplyr::select(Fly_Sex, Treatment, group_med, naive_med, fold_diff) |>
  print()

########## [PF1] FORMAL TEST: TREATMENT EFFECT ON ACTIVE BOUT NUMBER ##########

lmm_nbouts <- lme4::lmer(
  n_active_bouts ~ Treatment * Fly_Sex + (1 | FlyID),
  data    = bout_features,
  REML    = TRUE,
  control = lme4::lmerControl(optimizer = "bobyqa")
)

cat("\n===== LMM: Active bout number ~ Treatment * Sex =====\n")
print(summary(lmm_nbouts))

lmm_nbouts_lmerTest <- lmerTest::lmer(
  n_active_bouts ~ Treatment * Fly_Sex + (1 | FlyID),
  data    = bout_features,
  REML    = TRUE,
  control = lme4::lmerControl(optimizer = "bobyqa")
)

cat("\n===== TYPE III ANOVA: Active bout number ~ Treatment * Sex =====\n")
print(anova(lmm_nbouts_lmerTest, type = "III"))

emm_nbouts_treatment <- emmeans::emmeans(lmm_nbouts, ~ Treatment)
cat("\n===== PAIRWISE TREATMENT CONTRASTS: Active bout number (averaged over sex) =====\n")
print(emmeans::contrast(emm_nbouts_treatment, method = "pairwise", adjust = "tukey"))

emm_nbouts_sex <- emmeans::emmeans(lmm_nbouts, ~ Fly_Sex)
cat("\n===== SEX CONTRAST: Active bout number =====\n")
print(emmeans::contrast(emm_nbouts_sex, method = "pairwise"))

emm_nbouts_interaction <- emmeans::emmeans(lmm_nbouts, ~ Treatment | Fly_Sex)
cat("\n===== TREATMENT CONTRASTS BY SEX (Treatment x Sex interaction) =====\n")
print(emmeans::contrast(emm_nbouts_interaction, method = "pairwise", adjust = "tukey"))
