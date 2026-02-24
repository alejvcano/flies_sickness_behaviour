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

