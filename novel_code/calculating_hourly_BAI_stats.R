library(dplyr)
library(lubridate)

hourly_bai_stats <- function(df, date_col, circumference_col, bai_col) {
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
      bai_mm2 = {{ bai_col }}
    ) %>%
    dplyr::filter(!is.na(bai_mm2))
  
  hourly_bai <- df_with_bai %>%
    dplyr::mutate(
      date_hour = lubridate::floor_date(date, "hour")
    ) %>%
    dplyr::group_by(date_hour) %>%
    dplyr::summarise(
      total_hourly_bai_mm2 = sum(bai_mm2, na.rm = TRUE),
      max_circumference_mm = max(circumference_mm, na.rm = TRUE),
      min_circumference_mm = min(circumference_mm, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::rename(date = date_hour)
  
  return(hourly_bai)
}

#example usage
#hourly_df <- hourly_bai_stats(df, date, circumference_mm, bai_mm2)