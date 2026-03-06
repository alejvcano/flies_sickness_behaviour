########### DATA INPUT & PREP ############

response_data_raw <- read.csv("data/flies_data.csv", sep = ";")
key <- read.csv("data/flies_info.csv", sep = ";")

treatment_order <- c("naive", "PBS", "heatkill", "Pentomophila")
my_color_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")

endtime <- 1440 #endtime in minutes

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
  dplyr::filter(Treatment != "BLANK",Minute < endtime) |> 
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

########### PREPARE GROUPS ############
fly_data_hourly_list <- fly_data_hourly |>
  dplyr::arrange(id, hour) |>
  split(~ id)

fly_metadata <- fly_data_hourly |>
  dplyr::distinct(id, Fly_Sex, Treatment) |>
  dplyr::filter(!is.na(id), !is.na(Fly_Sex), !is.na(Treatment))

group_list <- fly_metadata |>
  dplyr::group_by(Fly_Sex, Treatment) |>
  dplyr::group_split()

group_names <- fly_metadata |>
  dplyr::group_by(Fly_Sex, Treatment) |>
  dplyr::group_keys() |>
  dplyr::mutate(group = paste(Fly_Sex, Treatment, sep = "__")) |>
  dplyr::pull(group)

########### COMPUTE GROUPED PLV ############
grouped_plv <- group_list |>
  purrr::set_names(group_names) |>
  purrr::map(function(meta_group) {
    requested_ids  <- stats::na.omit(meta_group$id)
    valid_ids      <- intersect(requested_ids, names(fly_data_hourly_list))
    flies_in_group <- fly_data_hourly_list[valid_ids]
    compute_group_plv_matrix(flies_in_group, min_length = 10) 
  })

pairwise_plv_data <- purrr::imap_dfr(grouped_plv, function(mat, group_name) {
  n <- nrow(mat)
  if (n < 2) return(tibble::tibble(group = group_name, plv = NA_real_))
  
  plv_values <- mat[lower.tri(mat)]
  tibble::tibble(group = group_name, plv = plv_values[!is.na(plv_values)])
}) |>
  tidyr::separate(group, into = c("Fly_Sex", "Treatment"), sep = "__") |>
  dplyr::filter(!is.na(plv)) |>
  dplyr::mutate(Treatment = factor(Treatment, levels = c("naive", "PBS", "heatkill", "Pentomophila")))

########### Statistics ############
glm_plv <- stats::glm(
  plv ~ Treatment * Fly_Sex,
  data = pairwise_plv_data,
  family = stats::quasibinomial(link = "logit")
)

print("--- ANOVA Results for PLV ---")
car::Anova(glm_plv, type = "II")

emm <- emmeans::emmeans(glm_plv, ~ Treatment | Fly_Sex)
pairs_results <- as.data.frame(emmeans::contrast(emm, method = "pairwise", adjust = "tukey"))

print(pairs_results)

get_letters_for_sex <- function(sex_level) {
  sub_pairs <- pairs_results |> dplyr::filter(Fly_Sex == sex_level)
  contrast_names <- gsub(" - ", "-", sub_pairs$contrast)
  pvals <- stats::setNames(sub_pairs$p.value, contrast_names)
  
  letters_obj <- multcompView::multcompLetters(pvals, threshold = 0.05)
  cld_vals <- unname(letters_obj$Letters)
  
  if (sex_level == "F") {
    cld_vals <- tolower(cld_vals) 
    cld_vals <- chartr("abcdefghijklmn", "abcdefghijklmn", cld_vals)
  } else if (sex_level == "M") {
    cld_vals <- tolower(cld_vals)
    cld_vals <- chartr("abcdefghijklmn", "efghijklmnopqrs", cld_vals)
  }
  
  tibble::tibble(
    Fly_Sex   = sex_level,
    Treatment = names(letters_obj$Letters),
    .group    = cld_vals
  )
}

cld_results <- purrr::map_dfr(unique(pairs_results$Fly_Sex), get_letters_for_sex) |>
  dplyr::mutate(Treatment = factor(Treatment, levels = levels(pairwise_plv_data$Treatment)))

########### PLOT OUTPUT ############

# Build the summary and enforce the factor order using your global setting
plv_summary <- pairwise_plv_data |>
  dplyr::group_by(Fly_Sex, Treatment) |>
  dplyr::summarise(
    mean_PLV = mean(plv),
    se_PLV   = sd(plv) / sqrt(dplyr::n()),
    .groups  = "drop"
  ) |>
  dplyr::left_join(cld_results, by = c("Fly_Sex", "Treatment")) |>
  dplyr::mutate(Treatment = factor(Treatment, levels = treatment_order))

# Plot
plv_plot <- ggplot2::ggplot(plv_summary, ggplot2::aes(x = Treatment, y = mean_PLV, fill = Treatment)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7, color = "black") +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mean_PLV - se_PLV, ymax = mean_PLV + se_PLV),
    width = 0.2,
    position = ggplot2::position_dodge(width = 0.8)
  ) +
  ggplot2::geom_text(
    ggplot2::aes(y = mean_PLV + se_PLV + 0.005, label = .group),
    position = ggplot2::position_dodge(width = 0.8),
    vjust    = 0,
    size     = 5,
    fontface = "plain" 
  ) +
  ggplot2::facet_wrap(~ Fly_Sex) +
  ggplot2::labs(
    x = "Treatment",
    y = "Mean Pairwise PLV"
  ) +
  ggplot2::scale_fill_manual(values = my_color_palette) +
  ggplot2::coord_cartesian(ylim = c(0, max(plv_summary$mean_PLV + plv_summary$se_PLV, na.rm = TRUE) + 0.02)) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "none", 
    axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
    strip.text      = ggplot2::element_text(face = "bold", size = 12)
  )

print(plv_plot)
