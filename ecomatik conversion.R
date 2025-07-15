
# Ecomatik DC3 Function to Convert Raw Micron Output to Circumference in mm 
# This function is adapted from the Ecomatik-provided Excel calculations
# Disclaimer: Generative AI (ChatGPT and Claude) was consulted in the creation of this function

dc3_to_circumference <- function(Vi_micrometers, C0, V0 = NULL) {
  # Input validation
  if (missing(Vi_micrometers) || missing(C0)) {
    stop("Both Vi_micrometers (dendrometer records in microns) and C0 (initial circumference in mm) are required")
  }
  if (length(Vi_micrometers) == 0) stop("Vi_micrometers must contain at least one value")
  if (C0 <= 0) stop("C0 (initial circumference in mm) must be positive")
  
  # Convert to mm "Vi", raw dendrometer records in microns, i.e., DEND_AVG 
  Vi_mm <- Vi_micrometers / 1000
  
  # Set "V0" (First valid value record of the dendrometer after installation in mm, i.e., plunger value or DEND_PLNG) if not provided
  if (is.null(V0)) {
    V0 <- Vi_mm[1]
    message("V0 not provided. Using first Vi measurement: ", V0, "mm")
  }
  
  # Constants
  Dini <- 128 #length of DC3 sensor body
  R0 <- C0 / (2 * pi)
  
  # Calculation of fixed wire circumference in mm
  sqrt_term_cwire <- sqrt((Dini + R0 - V0)^2 - R0^2)
  acos_term_cwire <- acos(R0 / (Dini + R0 - V0))
  Cwire <- C0 + 2 * sqrt_term_cwire - 2 * R0 * acos_term_cwire
  
  # Output vectors
  Ri <- numeric(length(Vi_mm))
  Ci <- numeric(length(Vi_mm))
  
  for (i in seq_along(Vi_mm)) {
    Ri_prev <- if (i == 1) R0 else Ri[i - 1]
    sqrt_term <- sqrt((Dini + Ri_prev - Vi_mm[i])^2 - Ri_prev^2)
    acos_term <- acos(Ri_prev / (Dini + Ri_prev - Vi_mm[i]))
    Ri[i] <- (Cwire - 2 * sqrt_term + 2 * Ri_prev * acos_term) / (pi * 2)
    Ci[i] <- Ri[i] * 2 * pi
  }
  
  # Output
  result <- data.frame(
    measurement = seq_along(Vi_micrometers),
    circumference_mm = Ci,
    circumference_microns = Ci * 1000,
    change_from_c0_mm = (Ci - C0),
    change_from_c0_microns = (Ci - C0) * 1000
  )
  
  return(result)
}