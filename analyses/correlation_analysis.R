########### DATA INPUT & PREP ############

response_data_raw <- read.csv("data/flies_data.csv", sep = ";")
key <- read.csv("data/flies_info.csv", sep = ";")

treatment_order <- c("naive", "PBS", "heatkill", "Pentomophila")

# Named color palette to ensure colors lock to treatments properly
my_color_palette <- stats::setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7"),
  treatment_order
)

endtime <- 12000 #endtime in minutes

response_data_pre <- response_data_raw |>
  tidyr::pivot_longer(
    cols = -c(Date, Time, Light1_Dark0),
    names_to = "FlyID",
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
  dplyr::mutate(Treatment = factor(Treatment, levels = treatment_order))

fly_data_hourly <- response_data_pre |> 
  dplyr::rename(id = "FlyID") |>
  dplyr::mutate(
    hour         = floor(Minute / 60),
    lastmov_hour = floor(Last_mov / 60)
  ) |>
  dplyr::group_by(id, Treatment, Fly_Sex, hour, lastmov_hour) |>
  dplyr::summarise(Activity = mean(Activity, na.rm = TRUE), .groups  = "drop")

########### EXTRACT PAIRWISE CORRELATIONS ############

pairwise_tau_data <- fly_data_hourly |> 
  dplyr::filter(hour < 24) |>
  dplyr::group_by(Fly_Sex, Treatment) |>
  tidyr::nest() |>
  dplyr::mutate(
    corr_matrix = purrr::map(data, function(group_df) {
      wide_df <- group_df |>
        tidyr::pivot_wider(names_from = id, values_from = Activity, id_cols = hour) |>
        dplyr::select(-hour)
      stats::cor(wide_df, method = "kendall", use = "pairwise.complete.obs")
    })
  ) |>
  # Rowwise extraction of the lower triangle for the stats model
  dplyr::rowwise() |>
  dplyr::mutate(
    tau_vals = list({
      mat <- corr_matrix
      vals <- mat[lower.tri(mat)]
      vals[!is.na(vals)] # Keep only valid pairwise correlations
    })
  ) |>
  dplyr::ungroup() |>
  dplyr::select(Fly_Sex, Treatment, tau_vals) |>
  tidyr::unnest(tau_vals) |>
  dplyr::rename(tau = tau_vals) |>
  dplyr::mutate(Treatment = factor(Treatment, levels = treatment_order))

########### STATISTICS (ANOVA & LETTERS) ############

# Kendall's tau is [-1, 1], standard lm() is robust here
glm_tau <- stats::lm(tau ~ Treatment * Fly_Sex, data = pairwise_tau_data)

print("--- ANOVA Results for Kendall's Tau ---")
print(car::Anova(glm_tau, type = "II"))

emm <- emmeans::emmeans(glm_tau, ~ Treatment | Fly_Sex)
pairs_results <- as.data.frame(emmeans::contrast(emm, method = "pairwise", adjust = "tukey"))

print(pairs_results)

# Function to generate letters: 'a,b,c' for Female, 'e,f,g' for Male
get_letters_for_sex <- function(sex_level) {
  sub_pairs <- pairs_results |> dplyr::filter(Fly_Sex == sex_level)
  contrast_names <- gsub(" - ", "-", sub_pairs$contrast)
  pvals <- stats::setNames(sub_pairs$p.value, contrast_names)
  
  letters_obj <- multcompView::multcompLetters(pvals, threshold = 0.05)
  cld_vals <- tolower(unname(letters_obj$Letters))
  
  if (sex_level == "F") {
    # Keep standard a,b,c...
    cld_vals <- chartr("abcdefghijklmnopqrstuvwxyz", "abcdefghijklmnopqrstuvwxyz", cld_vals)
  } else if (sex_level == "M") {
    # Shift a->e, b->f, c->g...
    cld_vals <- chartr("abcdefghijklmnopqrstuvwxyz", "efghijklmnopqrstuvwxyzabcd", cld_vals)
  }
  
  tibble::tibble(
    Fly_Sex   = sex_level,
    Treatment = names(letters_obj$Letters),
    .group    = cld_vals
  )
}

cld_results <- purrr::map_dfr(unique(pairs_results$Fly_Sex), get_letters_for_sex) |>
  dplyr::mutate(Treatment = factor(Treatment, levels = treatment_order))



########### SUMMARY & PLOT ############

tau_summary <- pairwise_tau_data |>
  dplyr::group_by(Fly_Sex, Treatment) |>
  dplyr::summarise(
    mean_tau = mean(tau),
    se_tau   = sd(tau) / sqrt(dplyr::n()),
    .groups  = "drop"
  ) |>
  dplyr::left_join(cld_results, by = c("Fly_Sex", "Treatment"))

tau_plot <- ggplot2::ggplot(tau_summary, ggplot2::aes(x = Treatment, y = mean_tau, fill = Treatment)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7, color = "black") +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mean_tau - se_tau, ymax = mean_tau + se_tau),
    width = 0.2,
    position = ggplot2::position_dodge(width = 0.8)
  ) +
  ggplot2::geom_text(
    ggplot2::aes(y = mean_tau + se_tau + 0.01, label = .group), # 0.01 offset above error bar
    position = ggplot2::position_dodge(width = 0.8),
    vjust    = 0,
    size     = 5,
    fontface = "plain"
  ) +
  ggplot2::facet_wrap(~ Fly_Sex) +
  ggplot2::labs(
    x = "Treatment",
    y = "Mean Pairwise Kendall's Tau"
  ) +
  ggplot2::scale_fill_manual(values = my_color_palette) +
  # Extend the y-axis dynamically to fit the letters
  ggplot2::coord_cartesian(ylim = c(0, max(tau_summary$mean_tau + tau_summary$se_tau, na.rm = TRUE) + 0.05)) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "none", # Removed legend
    axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
    strip.text      = ggplot2::element_text(face = "bold", size = 12)
  )

print(tau_plot)
