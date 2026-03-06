########### data input ############

response_data_raw = read.csv("data/flies_data.csv",sep = ";")
key = read.csv("data/flies_info.csv",sep = ";")

my_color_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")

response_data_pre <- response_data_pre |>
  dplyr::mutate(Treatment = factor(Treatment, levels = c("T1", "T2", "T3", "T4")))

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
  ) |> dplyr::filter(Treatment != "BLANK") |> remove_trailing_zeros() |> 
  dplyr::group_by(FlyID) |>
  dplyr::mutate(Last_mov = length(Activity)) |>
  dplyr::ungroup() |>
  dplyr::mutate(Treatment = factor(Treatment, levels = c("naive", "heatkill", "PBS", "Pentomophila")))

######### Last movement Plots ##########

response_data_pre |>
  dplyr::group_by(FlyID) |>
  dplyr::summarise(
    last_mov  = dplyr::first(Last_mov) / 60,
    Treatment = dplyr::first(Treatment),
    Fly_Sex   = dplyr::first(Fly_Sex),
    .groups   = "drop"
  ) |>
  ggplot(aes(x = Treatment, y = last_mov, fill = Treatment)) +
  geom_boxplot(
    outlier.shape = NA, alpha = 0.8, na.rm = TRUE,
    position = position_dodge()
  ) +
  geom_point(
    na.rm = TRUE,
    position = position_jitterdodge(0.27),
    pch = 21
  ) +
  facet_grid(rows = vars(Fly_Sex)) +
  labs(
    x = "Treatment",
    y = "Last Movement Time (hours)"
  ) +
  scale_fill_manual(values = my_color_palette) +
  theme_light() +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold"),
    axis.title      = element_text(size = 12),  # axis labels
    axis.text       = element_text(size = 10)   # tick labels
  )
########### PREPARE FLY DATA LIST ###########

fly_data_hourly <- response_data_pre |> dplyr::rename(id="FlyID") |>
  dplyr::mutate(
    hour     = floor(Minute / 60),
    lastmov_hour = floor(Last_mov / 60)
  ) |>
  dplyr::group_by(id, Treatment, Fly_Sex, hour, lastmov_hour) |>
  dplyr::summarise(
    Activity = mean(Activity, na.rm = TRUE),
    .groups  = "drop"
  )

fly_data <- ecoFAST::prep_data(fly_data_hourly,id="id",time="hour",var="Activity")

########## RUN DYNFOOT #############

fly_output <- ecoFAST::run_dynfoot(
  fly_data,
  metrics.type = "all",
  detrend.type = "none",
  group        = "Treatment",
  cores        = 5,
)

fly_output <- do.call(rbind,fly_output)

fly_output <- fly_output |> dplyr::left_join(key,by = dplyr::join_by(id == Monitor_TubeLoc))

fly_output <- split(fly_output,fly_output$id)

ecoFAST::dynfoot_summary_plot(fly_output,metrics = "all",group = "Treatment")

######## DEFINE METRICS AND BUILD EWS SUMMARY TABLE ############

metrics <- c(
  "mean.full", "SD.full", "ar1.full", "CV.full",
  "mean.roll", "SD.roll", "ar1.roll", "CV.roll"
)

group <- "Treatment"

ews_df <- ecoFAST:::sort_ews(do.call(rbind, fly_output), group, metrics)

str(ews_df)

########## RANK-NORMALISE FULL METRICS + GROUPING #############

metrics.norm       <- metrics[stringr::str_detect(metrics, "full")]
ews_df[metrics.norm] <- lapply(ews_df[metrics.norm], norm_ts, type = "rank")

ews_df <- ecoFAST:::stats_grouping(ews_df, metrics, group)

names(ews_df)[names(ews_df) == group] <- "clase"

# ==============================================================================
# DF2 FOR STATS
# No block variable in this dataset (original used block as random effect)
# ==============================================================================

df2 <- ews_df |>
  dplyr::mutate(
    clase_f = factor(clase, levels = c("naive", "heatkill", "PBS", "Pentomophila")),
    FlyID   = ID
  )

# ==============================================================================
# FIGURE 2A: FULL METRICS PLOT
# ==============================================================================

metrics.full <- metrics[stringr::str_detect(metrics, "full")]

temp_full <- as.data.frame(ews_df) |>
  dplyr::add_count(clase) |>
  dplyr::mutate(clase = paste0(clase, " (", n, ")")) |>
  reshape2::melt(measure.vars = metrics.full, variable.name = "metric")

lab_lookup_full        <- setNames(
  as.character(ews_df[1, paste0("labels.", metrics.full)]),
  metrics.full
)
temp_full$labels       <- unname(lab_lookup_full[as.character(temp_full$metric)])

metric_labs_full <- c(
  "mean.full" = "Mean",
  "SD.full"   = "SD",
  "ar1.full"  = "Lag-1 AC",
  "CV.full"   = "CV"
)

full_metric_plot <- ggplot(temp_full, aes(x = metric, y = value, fill = clase)) +
  geom_boxplot(
    outlier.shape = NA, alpha = 0.8, na.rm = TRUE,
    position = position_dodge()
  ) +
  geom_point(
    na.rm = TRUE,
    position = position_jitterdodge(0.07), pch = 21
  ) +
  labs(x = "Metric", y = "Percentile rank") +
  scale_x_discrete(labels = metric_labs_full) +
  theme_light() +
  scale_fill_manual(values = my_color_palette) +
  theme(
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 16)
  ) +
  geom_text(
    aes(x = metric, y = 1.05, label = labels, fontface = "bold"),
    check_overlap = TRUE
  ) +
  labs(fill = "Treatment (# obs)")

full_metric_plot

# ==============================================================================
# FIGURE 2B: ROLLING METRICS PLOT
# ==============================================================================

metrics.roll <- metrics[stringr::str_detect(metrics, "roll")]

temp_roll <- as.data.frame(ews_df) |>
  dplyr::add_count(clase) |>
  dplyr::mutate(clase = paste0(clase, " (", n, ")")) |>
  reshape2::melt(measure.vars = metrics.roll, variable.name = "metric")

lab_lookup_roll        <- setNames(
  as.character(ews_df[1, paste0("labels.", metrics.roll)]),
  metrics.roll
)
temp_roll$labels       <- unname(lab_lookup_roll[as.character(temp_roll$metric)])

metric_labs_roll <- c(
  "mean.roll" = "Mean",
  "SD.roll"   = "SD",
  "ar1.roll"  = "Lag-1 AC",
  "CV.roll"   = "CV"
)

rolling_metric_plot <- ggplot(temp_roll, aes(x = metric, y = value, fill = clase)) +
  geom_boxplot(
    outlier.shape = NA, alpha = 0.8, na.rm = TRUE,
    position = position_dodge()
  ) +
  geom_point(
    na.rm = TRUE,
    position = position_jitterdodge(0.05), pch = 21
  ) +
  labs(x = "Metric", y = "Kendall's tau") +
  scale_x_discrete(labels = metric_labs_roll) +
  theme_light() +
  scale_fill_manual(values = my_color_palette) +
  theme(
    axis.text    = element_text(size = 14),
    axis.title   = element_text(size = 16),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 14)
  ) +
  geom_text(
    aes(x = metric, y = 1.05, label = labels, fontface = "bold"),
    check_overlap = TRUE
  ) +
  labs(fill = "Treatment (# obs)")

rolling_metric_plot

# ==============================================================================
# COMBINE FIGURE 2
# ==============================================================================

legend <- cowplot::get_legend(
  rolling_metric_plot + theme(legend.position = "right")
)

full_metric_plot_clean    <- full_metric_plot    + theme(legend.position = "none")
rolling_metric_plot_clean <- rolling_metric_plot + theme(legend.position = "none")

plots    <- cowplot::plot_grid(
  full_metric_plot_clean, rolling_metric_plot_clean,
  nrow = 2, labels = c("A", "B")
)
Figure_2 <- cowplot::plot_grid(plots, legend, ncol = 2, rel_widths = c(1, 0.2))
Figure_2

# ==============================================================================
# STATISTICS
# Wilcoxon (2 groups) → Kruskal-Wallis (4 groups) + Dunn post-hoc
# ==============================================================================

full_metrics <- c("mean.full", "SD.full", "ar1.full", "CV.full")
roll_metrics <- c("mean.roll", "SD.roll", "ar1.roll", "CV.roll")

##### Kruskal-Wallis (replaces Wilcoxon) #####
p_kruskal_full <- sapply(full_metrics, \(m) {
  stats::kruskal.test(df2[[m]] ~ df2$clase_f)$p.value
})
p_kruskal_full

p_kruskal_roll <- sapply(roll_metrics, \(m) {
  stats::kruskal.test(df2[[m]] ~ df2$clase_f)$p.value
})
p_kruskal_roll

##### Dunn post-hoc with BH correction (replaces pairwise Wilcoxon) #####
p_dunn_full <- lapply(full_metrics, \(m) {
  rstatix::dunn_test(
    data            = df2,
    formula         = stats::as.formula(paste(m, "~ clase_f")),
    p.adjust.method = "BH"
  )
})
names(p_dunn_full) <- full_metrics
p_dunn_full

p_dunn_roll <- lapply(roll_metrics, \(m) {
  rstatix::dunn_test(
    data            = df2,
    formula         = stats::as.formula(paste(m, "~ clase_f")),
    p.adjust.method = "BH"
  )
})
names(p_dunn_roll) <- roll_metrics
p_dunn_roll

##### LM per metric (replaces lmer with block; no block in this dataset) #####
# If you want to control for Sex, change to: value ~ clase_f + Fly_Sex
long <- df2 |>
  tidyr::pivot_longer(
    cols      = dplyr::all_of(metrics),
    names_to  = "metric",
    values_to = "value"
  ) |>
  dplyr::filter(!is.na(value)) |>
  dplyr::mutate(metric = factor(metric, levels = metrics))

models_by_metric <- long |>
  dplyr::group_by(metric) |>
  dplyr::group_map(~{
    stats::lm(value ~ clase_f, data = .x)
  })
names(models_by_metric) <- levels(long$metric)

get_p <- function(a, term) {
  if (!term %in% rownames(a))        return(NA_real_)
  if ("Pr(>Chisq)" %in% colnames(a)) return(a[term, "Pr(>Chisq)"])
  if ("Pr(>F)"     %in% colnames(a)) return(a[term, "Pr(>F)"])
  NA_real_
}

anova_main <- purrr::imap_dfr(models_by_metric, ~{
  a <- car::Anova(.x, type = 2)
  tibble::tibble(
    metric         = .y,
    p_treat_main   = get_p(a, "clase_f")
  )
}) |>
  dplyr::mutate(p_treat_BH = stats::p.adjust(p_treat_main, method = "BH"))

anova_main




























fly_data_hourly <- response_data_pre |>
  dplyr::mutate(hour = floor(Minute / 60)) |>
  dplyr::filter(hour < 2) |>
  dplyr::group_by(FlyID, Treatment, hour) |>
  dplyr::summarise(activity = mean(Activity, na.rm = TRUE), .groups = "drop")

mean_data_hourly = fly_data_hourly |>
  dplyr::group_by(Treatment, hour) |>
  dplyr::summarise(mean_activity = mean(activity), .groups = "drop")

ggplot() +
  geom_line(data = fly_data_hourly, aes(x = hour, y = activity, group = FlyID, color = Treatment), alpha = 0.1) +
  geom_line(data = mean_data_hourly, aes(x = hour, y = mean_activity, color = "black"), size = 1.0) +
  facet_wrap(~Treatment, ncol = 2) +  # Separate panels for each Status
  labs(x = "hour",
       y = "mean activity",
       color = "Treatment") +
  theme_minimal() +
  theme(legend.position = "none")

fly_data_response = split(response_data_pre,response_data_pre$FlyID)

