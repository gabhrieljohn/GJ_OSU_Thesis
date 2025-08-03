dc3_to_circumference <- function(df, date_column = NULL, Vi_micrometers, C0, V0 = NULL, 
                                 na_method = "skip") {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("The 'dplyr' package is required. Please install it using install.packages('dplyr').")
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("The 'tibble' package is required. Please install it using install.packages('tibble').")
  }
  
  if (missing(Vi_micrometers) || missing(C0)) {
    stop("Both Vi_micrometers (dendrometer records in microns) and C0 (initial circumference in mm) are required")
  }
  if (length(Vi_micrometers) == 0) stop("Vi_micrometers must contain at least one value")
  if (C0 <= 0) stop("C0 (initial circumference in mm) must be positive")
  
  # Convert to mm
  Vi_mm <- Vi_micrometers / 1000
  
  # Handle V0
  if (is.null(V0)) {
    # Find first non-NA value
    first_valid_idx <- which(!is.na(Vi_mm))[1]
    if (is.na(first_valid_idx)) {
      stop("All Vi_micrometers values are NA")
    }
    V0 <- Vi_mm[first_valid_idx]
    message("V0 not provided. Using first non-NA Vi measurement: ", V0, " mm")
  }
  
  # Constants
  Dini <- 128
  R0 <- C0 / (2 * pi)
  
  # Calculate Cwire
  sqrt_term_cwire <- sqrt((Dini + R0 - V0)^2 - R0^2)
  acos_term_cwire <- acos(R0 / (Dini + R0 - V0))
  Cwire <- C0 + 2 * sqrt_term_cwire - 2 * R0 * acos_term_cwire
  
  # Initialize output vectors
  Ri <- numeric(length(Vi_mm))
  Ci <- numeric(length(Vi_mm))
  
  # Method 1: Linear interpolation for NA handling
  if (na_method == "interpolate") {
    # First, interpolate Vi_mm values
    Vi_mm_interp <- Vi_mm
    if (any(is.na(Vi_mm))) {
      # Simple linear interpolation
      na_indices <- which(is.na(Vi_mm))
      valid_indices <- which(!is.na(Vi_mm))
      
      if (length(valid_indices) >= 2) {
        Vi_mm_interp <- approx(valid_indices, Vi_mm[valid_indices], 
                               xout = 1:length(Vi_mm), rule = 2)$y
      } else if (length(valid_indices) == 1) {
        # If only one valid value, use it for all NAs
        Vi_mm_interp[na_indices] <- Vi_mm[valid_indices[1]]
      }
    }
    
    # Now calculate with interpolated values
    for (i in seq_along(Vi_mm_interp)) {
      Ri_prev <- if (i == 1) R0 else Ri[i - 1]
      
      sqrt_arg <- (Dini + Ri_prev - Vi_mm_interp[i])^2 - Ri_prev^2
      
      if (sqrt_arg < 0) {
        Ri[i] <- NA
        Ci[i] <- NA
        next
      }
      
      sqrt_term <- sqrt(sqrt_arg)
      acos_arg <- Ri_prev / (Dini + Ri_prev - Vi_mm_interp[i])
      acos_arg <- min(max(acos_arg, -1), 1)
      acos_term <- acos(acos_arg)
      
      Ri[i] <- (Cwire - 2 * sqrt_term + 2 * Ri_prev * acos_term) / (pi * 2)
      Ci[i] <- Ri[i] * 2 * pi
    }
    
    # Set results back to NA where original data was NA
    original_na <- is.na(Vi_mm)
    Ri[original_na] <- NA
    Ci[original_na] <- NA
  }
  
  # Method 2: Use last valid Ri for NA periods
  else if (na_method == "last_valid") {
    for (i in seq_along(Vi_mm)) {
      if (is.na(Vi_mm[i])) {
        Ri[i] <- NA
        Ci[i] <- NA
        next
      }
      
      # Find last valid Ri
      if (i == 1) {
        Ri_prev <- R0
      } else {
        # Look backwards for last valid Ri
        valid_ri_idx <- max(which(!is.na(Ri[1:(i-1)])))
        if (length(valid_ri_idx) == 0 || is.na(valid_ri_idx)) {
          Ri_prev <- R0
        } else {
          Ri_prev <- Ri[valid_ri_idx]
        }
      }
      
      sqrt_arg <- (Dini + Ri_prev - Vi_mm[i])^2 - Ri_prev^2
      
      if (sqrt_arg < 0) {
        Ri[i] <- NA
        Ci[i] <- NA
        next
      }
      
      sqrt_term <- sqrt(sqrt_arg)
      acos_arg <- Ri_prev / (Dini + Ri_prev - Vi_mm[i])
      acos_arg <- min(max(acos_arg, -1), 1)
      acos_term <- acos(acos_arg)
      
      Ri[i] <- (Cwire - 2 * sqrt_term + 2 * Ri_prev * acos_term) / (pi * 2)
      Ci[i] <- Ri[i] * 2 * pi
    }
  }
  
  # Method 3: Skip NA values but continue calculation
  else if (na_method == "skip") {
    last_valid_ri <- R0
    
    for (i in seq_along(Vi_mm)) {
      if (is.na(Vi_mm[i])) {
        Ri[i] <- NA
        Ci[i] <- NA
        next
      }
      
      sqrt_arg <- (Dini + last_valid_ri - Vi_mm[i])^2 - last_valid_ri^2
      
      if (sqrt_arg < 0) {
        Ri[i] <- NA
        Ci[i] <- NA
        next
      }
      
      sqrt_term <- sqrt(sqrt_arg)
      acos_arg <- last_valid_ri / (Dini + last_valid_ri - Vi_mm[i])
      acos_arg <- min(max(acos_arg, -1), 1)
      acos_term <- acos(acos_arg)
      
      Ri[i] <- (Cwire - 2 * sqrt_term + 2 * last_valid_ri * acos_term) / (pi * 2)
      Ci[i] <- Ri[i] * 2 * pi
      
      # Update last valid Ri
      last_valid_ri <- Ri[i]
    }
  }
  
  # Build result tibble
  if (!is.null(date_column)) {
    if (!inherits(date_column, c("POSIXct", "POSIXt", "Date"))) {
      stop("The provided 'date_column' must be of class POSIXct, POSIXt, or Date")
    }
    
    if (is.unsorted(date_column[!is.na(date_column)])) {
      warning("Date column is not sorted in ascending order. Results may be misaligned.")
    }
    
    result <- tibble::tibble(
      date = date_column,
      radius_mm = Ri,
      circumference_mm = Ci,
      change_from_c0_mm = Ci - C0
    )
  } else {
    result <- tibble::tibble(
      radius_mm = Ri,
      circumference_mm = Ci,
      change_from_c0_mm = Ci - C0
    )
  }
  
  return(result)
}