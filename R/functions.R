remove_trailing_zeros <- function(df) {
  df |>
    dplyr::group_by(FlyID) |>
    dplyr::mutate(
      reversed_activity = rev(Activity),
      trailing_zeros = cumsum(reversed_activity != 0) == 0,
      trailing_zeros = rev(trailing_zeros),
    ) |>
    dplyr::filter(!trailing_zeros) |>
    dplyr::select(-reversed_activity, -trailing_zeros) |>
    dplyr::ungroup()
}

align_time_series <- function(fly_list) {
  min_len <- min(sapply(fly_list, nrow))
  lapply(fly_list, function(df) df[1:min_len, , drop = FALSE])
}

compute_plv <- function(fly1_df, fly2_df) {
  # Find the overlapping time range
  common_hours <- intersect(fly1_df$hour, fly2_df$hour)
  
  # If they overlap by fewer than 2 hours, phase sync is impossible
  if (length(common_hours) < 2) return(NA_real_)
  
  # Subset both flies to the exact same hours
  d1 <- fly1_df |> dplyr::filter(hour %in% common_hours) |> dplyr::arrange(hour)
  d2 <- fly2_df |> dplyr::filter(hour %in% common_hours) |> dplyr::arrange(hour)
  
  t1 <- as.matrix(d1[, c("hour", "Activity")])
  t2 <- as.matrix(d2[, c("hour", "Activity")])
  
  # Try/Catch protects against any lingering internal errors in phase.sync
  plv_result <- tryCatch(
    synchrony::phase.sync(t1, t2, nrands = 0, mod = 1, nbreaks = 100, method = "fft", quiet = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(plv_result)) return(NA_real_)
  
  pd <- plv_result$deltaphase$phasediff
  if (length(pd) == 0 || all(is.na(pd))) return(NA_real_)
  
  abs(mean(exp(1i * pd), na.rm = TRUE))
}

is_valid_fly <- function(df, min_length = 10) {
  nrow(df) >= min_length &&
    !all(is.na(df$Activity)) &&
    length(unique(df$Activity)) > 1 &&
    !all(df$Activity == 0)
}

compute_group_plv_matrix <- function(group_flies, min_length = 10) {
  # Filter out invalid or excessively short flies
  group_flies <- Filter(function(df) is_valid_fly(df, min_length = min_length), group_flies)
  
  n <- length(group_flies)
  if (n < 2) return(matrix(NA_real_, nrow = n, ncol = n))
  
  # REMOVED: group_flies <- align_time_series(group_flies) 
  
  plv_mat <- matrix(NA_real_, nrow = n, ncol = n)
  rownames(plv_mat) <- names(group_flies)
  colnames(plv_mat) <- names(group_flies)
  
  for (i in 1:n) {
    for (j in i:n) {
      if (i == j) {
        plv_mat[i, j] <- 1
      } else {
        plv <- compute_plv(group_flies[[i]], group_flies[[j]])
        plv_mat[i, j] <- plv
        plv_mat[j, i] <- plv
      }
    }
  }
  plv_mat
}