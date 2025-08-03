library(dplyr)
library(tidyverse)

smart_bind_rows <- function(...) {
  dfs <- list(...)
  
  # Find all column names
  all_names <- unique(unlist(map(dfs, names)))
  
  # Fill missing columns
  dfs_filled <- map(dfs, ~mutate(.x, !!!setNames(rep(list(NA), length(setdiff(all_names, names(.x)))), setdiff(all_names, names(.x)))))
  
  # Determine common type per column
  common_types <- map(all_names, function(col) {
    col_classes <- map_chr(dfs_filled, ~class(.x[[col]])[1])
    if ("character" %in% col_classes) {
      return("character")
    } else if ("numeric" %in% col_classes | "double" %in% col_classes) {
      return("numeric")
    } else if ("integer" %in% col_classes) {
      return("integer")
    } else if ("logical" %in% col_classes) {
      return("logical")
    } else {
      return("character") # fallback
    }
  })
  names(common_types) <- all_names
  
  # Coerce to common type
  dfs_coerced <- map(dfs_filled, function(df) {
    for (col in all_names) {
      if (!col %in% names(df)) next
      target_type <- common_types[[col]]
      df[[col]] <- switch(
        target_type,
        "character" = as.character(df[[col]]),
        "numeric"   = as.numeric(df[[col]]),
        "integer"   = as.integer(df[[col]]),
        "logical"   = as.logical(df[[col]]),
        df[[col]]
      )
    }
    df
  })
  
  bind_rows(dfs_coerced)
}
