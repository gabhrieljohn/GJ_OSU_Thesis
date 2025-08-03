library(dplyr)
library(lubridate)

daily_bai_stats <- function(df, date_col, circumference_col, basal_area_col, bai_col) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install the 'dplyr' package.")
  }
  if (!requireNamespace("lubridate", quietly = TRUE)) {
    stop("Please install the 'lubridate' package.")
  }
  
  df_with_bai <- df %>%
    dplyr::mutate(
      date = {{ date_col }},
      circumference_mm = {{ circumference_col }},
      basal_area_mm2 = {{ basal_area_col }},
      bai_mm2 = {{ bai_col }}
    ) %>%
    dplyr::filter(!is.na(bai_mm2))
  
  daily_bai <- df_with_bai %>%
    dplyr::mutate(date_day = lubridate::as_date(date)) %>%
    dplyr::group_by(date_day) %>%
    dplyr::summarise(
      total_bai_mm2 = sum(bai_mm2, na.rm = TRUE),
      max_bai_mm2 = max(bai_mm2, na.rm = TRUE),
      daily_circ_amplitude_mm = max(
        circumference_mm, na.rm = TRUE) - min(
          circumference_mm, na.rm = TRUE),
      daily_bai_amplitude_mm2 = max(
        bai_mm2, na.rm = TRUE) - min(bai_mm2, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::arrange(date_day) %>%
    dplyr::mutate(
      cumu_bai_mm2 = cumsum(total_bai_mm2),
      cumu_bai_mm2_norm = cumu_bai_mm2 / max(cumu_bai_mm2, na.rm = TRUE),
      delta_max_bai_mm2 = c(NA, diff(max_bai_mm2)),
      DOY = lubridate::yday(date_day),
      year = as.factor(lubridate::year(date_day))
    ) %>%
    dplyr::rename(date = date_day)
  
  return(daily_bai)
}

#example usage
#daily_df  <- daily_bai_stats(df, date, circumference_mm, basal_area_mm2, bai_mm2)