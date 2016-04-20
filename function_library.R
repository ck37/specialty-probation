# Review distributions of variables that match a string.
review_by_string = function(data, string, ignore_case = T) {
  cols = grep(string, names(data), ignore.case=ignore_case)
  print(names(data)[cols])
  for (col in cols) {
    cat(names(data)[col], ":\n")
    print(summary(data[, col]))
  }
  return(invisible(names(data)[cols]))
}

# Create a new dataset with just the desired fields, which are also renamed.
clean_dataset = function(data, keep_fields) {
  # Create a blank data frame with the correct number of rows.
  clean_data = data.frame(matrix(nrow=nrow(data)))
  # Remove the weird default column that was created.
  clean_data[, 1] = NULL
  dim(clean_data)

  cat("Saving fields... ")
  keep_fields
  names(keep_fields)
  for (name in names(keep_fields)) {
    cat(name, " ")
    clean_data[, keep_fields[[name]]] = data[, name]
  }
  cat("\n")
  cat("New dimensions:", paste(dim(clean_data), collapse=" "), "\n")
  cat("New variables:", paste(names(clean_data), collapse=", "))
  clean_data
}

impute_missing_values = function(data, add_indicators = T) {
  # Loop over each feature.
  missing_indicators = NULL
  for (i in 1:ncol(data)){
    nas = sum(is.na(data[, i]))
    # Nothing to impute, continue to next column.
    if (nas == 0) {
      #cat("Skipping", i, "\n")
      next
    }
    #cat("Found", nas, "missing values for", colnames(data)[i], "\n")


    missing_indicator = matrix(1*is.na(data[, i]))
    new_name = paste0("miss_", colnames(data)[i])
    colnames(missing_indicator) = new_name
    next_col = ncol(missing_indicators) + 1
    if (is.null(missing_indicators)) {
      missing_indicators = matrix(missing_indicator, nrow=nrow(data), ncol=1)
      colnames(missing_indicators) = new_name
    } else {
      missing_indicators = cbind(missing_indicators, missing_indicator)
    }
    #new_colnames = c(colnames(missing_indicators), )
    #cat("New colnames:", new_colnames, "\n")
    #colnames(missing_indicators) = new_colnames

    if (class(data[, i]) == "factor") {
      # Impute factors to the mode.
      # Choose the first mode in case of ties.
      data[is.na(data[, i]), i] = Mode(data[, i])[1]
    } else {
      # Impute numeric values to the median.
      data[is.na(data[, i]), i] = median(data[, i], na.rm=T)
    }
  }
  if (add_indicators) {
    #cat("New indicator names:\n")
    data = cbind(data, missing_indicators)
  }
  data
}

# Via http://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
Mode = function(x) {
  ux = unique(x)
  tab <- tabulate(match(x, ux))
  # This will return multiple modes if they exist.
  ux[tab == max(tab)]
}
