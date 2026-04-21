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

########## DATA INPUT & PREP ##########

response_data_raw <- read.csv("data/flies_data.csv", sep = ";")
key               <- read.csv("data/flies_info.csv", sep = ";")

treatment_order  <- c("PentomophilaA", "PentomophilaB")
my_color_palette <- stats::setNames(
  c("#CC79A7", "#E69F00"),
  treatment_order
)

MINUTES_PER_HOUR <- 60
WINDOW_HOURS     <- 6
WINDOW_MINUTES   <- WINDOW_HOURS * MINUTES_PER_HOUR
ZOOM_MINUTES     <- 6 * 60                      # first 3 days = 4320 min
PENTO_THRESHOLD  <- 12000                            # same as in your reference plot

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
  # ── Recode Pentomophila into A/B using raw Last_mov ──
  dplyr::mutate(
    Treatment = dplyr::case_when(
      Treatment == "Pentomophila" & Last_mov <  PENTO_THRESHOLD ~ "PentomophilaA",
      Treatment == "Pentomophila" & Last_mov >= PENTO_THRESHOLD ~ "PentomophilaB",
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
  ) |> dplyr::filter(Last_mov>1440)

unique(response_data_pre$FlyID)

fly_data_pre <- ecoFAST::prep_data(df = response_data_pre,id = "FlyID",time = "Minute",var = "Activity")

fly_data_pre$ts <- fly_data_pre$ts |>
  purrr::map(~ dplyr::left_join(
    .x,
    response_data_pre |>
      dplyr::select(FlyID, Fly_Sex, Treatment, Status_dead0_alive1) |>
      dplyr::distinct(FlyID, .keep_all = TRUE),
    by = dplyr::join_by(FlyID)
  ))

fly_data_dynfoot <- ecoFAST::run_dynfoot(fly_data_pre,detrend.type = "gaussian",group = "Treatment",cores = 4,winsize = 10,winsize_is_percentage = TRUE,stepsize = 5,min.length = 50)

test <- fly_data_dynfoot |>
  purrr::map(~ .x |>
               dplyr::rename(id = FlyID) |>
               dplyr::mutate(Treatment = factor(Treatment, levels = treatment_order)))
ecoFAST::dynfoot_summary_plot(test, metrics = "all", group = "Treatment")

ml_test <- ecoFAST::run_ml_model(test,target = "Treatment")

ecoFAST::ml_summary_plot(ml_test)

table(ml_test$output_data$Treatment,ml_test$output_data$Predictedclass)
