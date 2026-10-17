inspect_df <- function(df) {
  cat("--- Dataset Overview ---\n")
  cat("Dimensions:", nrow(df), "rows x", ncol(df), "columns\n")
  cat("Total NA values:", sum(is.na(df)), "\n\n")
  
  cat("--- Column Breakdown ---\n")
  summary_info <- data.frame(
    Column = names(df),
    Type = sapply(df, function(x) class(x)[1]),
    Missing = sapply(df, function(x) sum(is.na(x))),
    Missing_Pct = paste0(round(sapply(df, function(x) mean(is.na(x))) * 100, 1), "%"),
    Unique = sapply(df, function(x) length(unique(x))),
    row.names = NULL
  )
  print(summary_info)
}