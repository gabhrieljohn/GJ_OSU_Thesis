library(dplyr)
#circumference to BAI is inspired by the area of a circle: BA=pi*(DBH/2)^2
circ_to_bai <- function(df, date_col = NULL, radius_col, circumference_col) {
  
  # Check if dplyr is available
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("The 'dplyr' package is required. Please install it using install.packages('dplyr').")
  }
  
  # If date_col is provided, check chronological order
  if (!is.null(substitute(date_col))) {
    date_values <- df %>% dplyr::pull({{ date_col }})
    if (is.unsorted(date_values[!is.na(date_values)])) {
      warning("Date column is not in chronological order. BAI calculations may be incorrect.")
    }
  }
  
  result <- df %>%
    dplyr::mutate(
      basal_area_mm2 = pi * radius_mm^2,
      bai_mm2 = c(NA, diff(basal_area_mm2))
    )
  
  return(result)
}