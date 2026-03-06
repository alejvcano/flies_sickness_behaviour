########## DATA INPUT & PREP ##########

response_data_raw <- read.csv("data/flies_data.csv", sep = ";")
key               <- read.csv("data/flies_info.csv", sep = ";")

treatment_order <- c("naive", "PBS", "heatkill", "Pentomophila")

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
  dplyr::mutate(Treatment = factor(Treatment, levels = treatment_order))

########## BUILD PER-FLY LIST FOR run_dynfoot() ##########
# run_dynfoot() expects a list of data frames, each with:
#   ID (fly identifier), x (time), y (activity), class (grouping)

flydata <- split(response_data_pre,response_data_pre$FlyID)

########## RUN ecoFAST ##########

flydata <- ecoFAST::prep_data(fl)

flyoutput <- ecoFAST::run_dynfoot(
  flydata,
  metrics.type = "all",
  detrend.type = "gaussian",
  group        = "Treatment",
  cores        = 4,     
)

########## TIDY EWS RESULTS ##########

metrics <- c("mean.full", "SD.full", "ar1.full", "CV.full",
             "mean.roll", "SD.roll",  "ar1.roll",  "CV.roll")

ewsdf <- ecoFAST::sort_ewsdata(
  do.call(rbind, flyoutput),
  group   = "class",
  metrics = metrics
)

########## ROLLING METRICS BOXPLOT ##########

metrics_roll <- metrics[stringr::str_detect(metrics, "roll")]

metriclabs_roll <- c(
  "mean.roll" = "Mean",
  "SD.roll"   = "SD",
  "ar1.roll"  = "Lag-1 AC",
  "CV.roll"   = "CV"
)

temproll <- as.data.frame(ewsdf) |>
  dplyr::add_count(clase) |>
  dplyr::mutate(clase = paste0(clase, " (n=", n, ")")) |>
  reshape2::melt(measure.vars = metrics_roll, variable.name = "metric")

lablookup_roll <- stats::setNames(
  as.character(ewsdf[[1]]),
  paste0("labels.", metrics_roll)
)
temproll$labels <- unname(lablookup_roll[as.character(temproll$metric)])

rollingmetricplot <- ggplot2::ggplot(
    temproll,
    ggplot2::aes(x = metric, y = value, fill = clase)
  ) +
  ggplot2::geom_boxplot(
    outlier.shape = NA, alpha = 0.8, na.rm = TRUE,
    position = ggplot2::position_dodge()
  ) +
  ggplot2::geom_point(
    na.rm = TRUE,
    position = ggplot2::position_jitterdodge(0.05),
    pch = 21
  ) +
  ggplot2::labs(x = "Metric", y = "Kendall's tau", fill = "Class obs") +
  ggplot2::scale_x_discrete(labels = metriclabs_roll) +
  ggplot2::theme_light() +
  ggplot2::theme(
    axis.text    = ggplot2::element_text(size = 14),
    axis.title   = ggplot2::element_text(size = 16),
    legend.text  = ggplot2::element_text(size = 12),
    legend.title = ggplot2::element_text(size = 14)
  ) +
  ggplot2::scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  ggplot2::geom_text(
    ggplot2::aes(x = metric, y = 1.05, label = labels),
    fontface = "bold", check_overlap = TRUE
  )

print(rollingmetricplot)
ggplot2::ggsave("Figure2_rolling_metrics.png", rollingmetricplot,
                width = 10, height = 6, dpi = 300)
