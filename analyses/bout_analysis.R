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

get_letters_for_sex <- function(pairs_df, sex_level) {
  sub_pairs      <- pairs_df |> dplyr::filter(Fly_Sex == sex_level)
  contrast_names <- gsub(" - ", "-", sub_pairs$contrast)
  pvals          <- stats::setNames(sub_pairs$p.value, contrast_names)
  letters_obj    <- multcompView::multcompLetters(pvals, threshold = 0.05)
  cld_vals       <- tolower(unname(letters_obj$Letters))

  if (sex_level == "M") {
    cld_vals <- chartr("abcdefghijklmnopqrstuvwxyz",
                       "efghijklmnopqrstuvwxyzabcd", cld_vals)
  }

  tibble::tibble(
    Fly_Sex   = sex_level,
    Treatment = names(letters_obj$Letters),
    .group    = cld_vals
  )
}

add_letters <- function(model, data_for_levels) {
  emm           <- emmeans::emmeans(model, ~ Treatment | Fly_Sex)
  pairs_results <- as.data.frame(
    emmeans::contrast(emm, method = "pairwise", adjust = "tukey")
  )
  purrr::map_dfr(
    unique(pairs_results$Fly_Sex),
    ~ get_letters_for_sex(pairs_results, .x)
  ) |>
    dplyr::mutate(
      Treatment = factor(Treatment, levels = levels(data_for_levels$Treatment))
    )
}

make_bar_plot <- function(raw_df, summary_df, cld_df, y_var, y_lab,
                          scientific_y = FALSE) {

  plot_df <- raw_df |>
    dplyr::filter(!is.na(!!rlang::sym(y_var))) |>
    dplyr::select(FlyID, Fly_Sex, Treatment, !!rlang::sym(y_var)) |>
    dplyr::rename(value = !!rlang::sym(y_var))

  summary_with_cld <- summary_df |>
    dplyr::left_join(cld_df, by = c("Fly_Sex", "Treatment"))

  dodge_width  <- 0.8
  jitter_width <- 0.15
  max_y        <- max(plot_df$value, na.rm = TRUE)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Treatment, fill = Treatment)) +
    ggplot2::geom_point(
      ggplot2::aes(y = value),
      position = ggplot2::position_jitterdodge(
        jitter.width = jitter_width,
        dodge.width  = dodge_width
      ),
      size   = 2.0,
      alpha  = 0.5,
      color  = "gray60",
      stroke = 0
    ) +
    ggplot2::geom_col(
      data = summary_df,
      ggplot2::aes(y = mean_value),
      position = ggplot2::position_dodge(width = dodge_width),
      width    = 0.7,
      color    = "black",
      alpha    = 0.8
    ) +
    ggplot2::geom_errorbar(
      data = summary_df,
      ggplot2::aes(
        x = Treatment,
        ymin = pmax(mean_value - ci_value, 0),
        ymax = mean_value + ci_value
      ),
      position = ggplot2::position_dodge(width = dodge_width),
      width    = 0.2,
      color    = "black"
    ) +
    ggplot2::geom_text(
      data = summary_with_cld,
      ggplot2::aes(
        x     = Treatment,
        y     = mean_value + ci_value + 0.05 * max_y,
        label = .group
      ),
      position     = ggplot2::position_dodge(width = dodge_width),
      vjust        = 0,
      size         = 4,
      fontface     = "bold",
      inherit.aes  = FALSE
    ) +
    ggplot2::facet_wrap(~ Fly_Sex) +
    ggplot2::labs(x = NULL, y = y_lab) +
    ggplot2::scale_fill_manual(values = my_color_palette) +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y     = ggplot2::element_text(size = 10),
      axis.title.y    = ggplot2::element_text(size = 11),
      strip.text      = ggplot2::element_text(face = "bold", size = 11)
    )

  if (scientific_y) {
    p <- p + ggplot2::scale_y_continuous(labels = scales::scientific)
  }

  p
}

########## DATA INPUT & PREP ##########

response_data_raw <- read.csv("data/flies_data.csv", sep = ";")
key               <- read.csv("data/flies_info.csv", sep = ";")

treatment_order <- c("naive", "PBS", "heatkill", "Pentomophila")

my_color_palette <- stats::setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7"),
  treatment_order
)

endtime <- 40000 # endtime in minutes

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
  dplyr::filter(Treatment != "BLANK", Minute < endtime) |>
  remove_trailing_zeros() |>
  dplyr::group_by(FlyID) |>
  dplyr::mutate(Last_mov = length(Activity)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    Treatment = factor(Treatment, levels = treatment_order)
  )

fly_data_features <- response_data_pre |>
  dplyr::group_by(FlyID, Fly_Sex, Treatment) |>
  dplyr::mutate(
    state       = ifelse(Activity > 0, "active", "inactive"),
    bout_id     = data.table::rleid(state),
    bout_length = sequence(rle(state)$lengths)
  ) |>
  dplyr::ungroup()

transition_stats <- fly_data_features |>
  dplyr::group_by(FlyID, Fly_Sex, Treatment) |>
  dplyr::summarise(
    mean_active_bout   = mean(bout_length[state == "active"], na.rm = TRUE),
    mean_inactive_bout = mean(bout_length[state == "inactive"], na.rm = TRUE),
    activity_rate      = sum(Activity[state == "active"], na.rm = TRUE),
    prop_active        = sum(state == "active", na.rm = TRUE) / dplyr::n(),
    transition_prob    = {
      mc_fit <- try(markovchain::markovchainFit(data = state), silent = TRUE)
      if (inherits(mc_fit, "try-error") || is.null(mc_fit$estimate)) {
        NA_real_
      } else {
        tm <- mc_fit$estimate@transitionMatrix
        if (all(c("inactive", "active") %in% rownames(tm)) &&
            all(c("inactive", "active") %in% colnames(tm))) {
          as.numeric(tm["inactive", "active"])
        } else {
          mean(tm[row(tm) != col(tm)], na.rm = TRUE)
        }
      }
    },
    .groups = "drop"
  ) |>
  dplyr::mutate(
    activity_rate = pmax(activity_rate, 1e-6)
  )

########## EXTRA DATA FOR FIRST 24 HOURS ##########

response_data_24h <- response_data_pre |>
  dplyr::filter(Minute <= 1440)

fly_data_features_24h <- response_data_24h |>
  dplyr::group_by(FlyID, Fly_Sex, Treatment) |>
  dplyr::mutate(
    state       = ifelse(Activity > 0, "active", "inactive"),
    bout_id     = data.table::rleid(state),
    bout_length = sequence(rle(state)$lengths)
  ) |>
  dplyr::ungroup()

transition_stats_24h <- fly_data_features_24h |>
  dplyr::group_by(FlyID, Fly_Sex, Treatment) |>
  dplyr::summarise(
    activity_rate   = sum(Activity[state == "active"], na.rm = TRUE),
    prop_active     = sum(state == "active", na.rm = TRUE) / dplyr::n(),
    .groups         = "drop"
  ) |>
  dplyr::mutate(
    activity_rate = pmax(activity_rate, 1e-6)
  )

########## MODELS + SUMMARY DATA + PLOTS ##########
## Full duration

# Activity rate (Panel A) - scientific notation
glm_activity_rate <- stats::glm(
  activity_rate ~ Treatment * Fly_Sex,
  data = transition_stats, family = stats::Gamma(link = "log")
)
cld_rate <- add_letters(glm_activity_rate, transition_stats)
rate_summary_full <- transition_stats |>
  dplyr::group_by(Fly_Sex, Treatment) |>
  dplyr::summarise(
    mean_value = mean(activity_rate, na.rm = TRUE),
    ci_value   = 1.96 * stats::sd(activity_rate, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups    = "drop"
  )
plot_activity_rate <- make_bar_plot(
  raw_df        = transition_stats,
  summary_df    = rate_summary_full,
  cld_df        = cld_rate,
  y_var         = "activity_rate",
  y_lab         = "Total activity (a.u.)",
  scientific_y  = TRUE
)

# Proportion active (Panel B)
glm_prop_active <- stats::glm(
  prop_active ~ Treatment * Fly_Sex,
  data = transition_stats, family = stats::quasibinomial(link = "logit")
)
cld_prop <- add_letters(glm_prop_active, transition_stats)
prop_summary_full <- transition_stats |>
  dplyr::group_by(Fly_Sex, Treatment) |>
  dplyr::summarise(
    mean_value = mean(prop_active, na.rm = TRUE),
    ci_value   = 1.96 * stats::sd(prop_active, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups    = "drop"
  )
plot_prop_active <- make_bar_plot(
  raw_df        = transition_stats,
  summary_df    = prop_summary_full,
  cld_df        = cld_prop,
  y_var         = "prop_active",
  y_lab         = "Proportion time active",
  scientific_y  = FALSE
)

## First 24 hours (Panels C–D)

# Activity rate 24 h (Panel C) - scientific notation
glm_activity_rate_24h <- stats::glm(
  activity_rate ~ Treatment * Fly_Sex,
  data = transition_stats_24h, family = stats::Gamma(link = "log")
)
cld_rate_24h <- add_letters(glm_activity_rate_24h, transition_stats_24h)
rate_summary_24h <- transition_stats_24h |>
  dplyr::group_by(Fly_Sex, Treatment) |>
  dplyr::summarise(
    mean_value = mean(activity_rate, na.rm = TRUE),
    ci_value   = 1.96 * stats::sd(activity_rate, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups    = "drop"
  )
plot_activity_rate_24h <- make_bar_plot(
  raw_df        = transition_stats_24h,
  summary_df    = rate_summary_24h,
  cld_df        = cld_rate_24h,
  y_var         = "activity_rate",
  y_lab         = "Total activity 0–24 h (a.u.)",
  scientific_y  = TRUE
)

# Proportion active 24 h (Panel D)
glm_prop_active_24h <- stats::glm(
  prop_active ~ Treatment * Fly_Sex,
  data = transition_stats_24h, family = stats::quasibinomial(link = "logit")
)
cld_prop_24h <- add_letters(glm_prop_active_24h, transition_stats_24h)
prop_summary_24h <- transition_stats_24h |>
  dplyr::group_by(Fly_Sex, Treatment) |>
  dplyr::summarise(
    mean_value = mean(prop_active, na.rm = TRUE),
    ci_value   = 1.96 * stats::sd(prop_active, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups    = "drop"
  )
plot_prop_active_24h <- make_bar_plot(
  raw_df        = transition_stats_24h,
  summary_df    = prop_summary_24h,
  cld_df        = cld_prop_24h,
  y_var         = "prop_active",
  y_lab         = "Proportion time active 0–24 h",
  scientific_y  = FALSE
)

########## ASSEMBLE FIGURE 2 ##########

Fig2 <- cowplot::plot_grid(
  plot_activity_rate,       # A: total activity (full, scientific)
  plot_prop_active,         # B: proportion active (full)
  plot_activity_rate_24h,   # C: total activity first 24 h (scientific)
  plot_prop_active_24h,     # D: proportion active first 24 h
  labels     = c("A", "B", "C", "D"),
  label_size = 16,
  ncol       = 2,
  nrow       = 2,
  align      = "hv",
  axis       = "tblr"
)

print(Fig2)

cowplot::save_plot("outputs/Figure2_immune.png", Fig2, base_width = 12, base_height = 10)
