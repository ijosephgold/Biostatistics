# ================================================================
# Joe's BIOL 330 BIOSTATISTICS TOOLKIT
# Revision 18 - adds assignment-pattern preparation helpers and clearer clues.
# ================================================================
# This is the single master file to source from .Rprofile.
# Recommended startup command:
# source("C:/Rstudio/Biostatistics/biostats_toolkit.R")
#
# Quick starts (replace data and column names with yours):
#   inspect_df(data)
#   analysis_start(data)  # first look at any new worksheet
#   analysis_start(data, question = "Does the outcome change by treatment?")
#   analysis_start(data, outcome = "PlantHeight", predictors = "Treatment",
#                  goal = "difference", design = "independent")
#   analysis_guide(data, outcome = "RatDnsty",
#                  predictors = c("ShrbCvr", "SeedPrd"))
#   correlation_workflow(data, "RatDnsty", "ShrbCvr")
#   linear_regression_workflow(data, RatDnsty ~ ShrbCvr)
#   logistic_regression_workflow(data, Survival ~ BodySize, event = "yes")
#   anova_workflow(data, PlantHeight ~ Treatment)
#   independent_t_test(data, outcome = "PlantHeight", group = "Treatment")
#   paired_t_test(data, before = "Before", after = "After")
#   data <- add_paired_difference(data, before = "Before", after = "After")
#   data <- add_product_variable(data, first = "Length", second = "Width")
#   data <- add_ratio_variable(data, numerator = "Protein", denominator = "DryMass")
#   data <- add_tree_height(data, distance = "Distance", angle = "Angle",
#                           observer_height = "ObserverHeight")
#   diversity <- shannon_by_sample(data, "Transect", "Species", "Abundance", "Habitat")
#   dispersion_summary(data, response = "Count", group = "Site")
#   data <- add_gender_contribution_index(data, "MaleCones", "FemaleCones")
#   goodness_of_fit_workflow(c(12, 8, 10), probabilities = c(1/3, 1/3, 1/3))
#   contingency_test_workflow(data, row = "Species", column = "Habitat")
#   binomial_test_workflow(successes = 18, trials = 25, probability = 0.5)
#   repeated_measures_anova_workflow(data, "HeartRate", "Subject", "Test")
#   random_effects_anova_workflow(data, "Area", "Genet")
#   permanova_workflow(data, variables = c("Length", "Width"), group = "Site")
#   categorical_test_guide()
#   utils::View(data)    # R distinguishes upper- and lower-case letters
#
# inspect_df() can flag likely sequential IDs and suggest (but never assume)
# a paired two-condition layout. add_paired_difference() adds a signed
# after-minus-before column when you explicitly request it.
# Matrices and clearly named lists of equal-length columns can be converted
# automatically. One-column vectors are wrapped automatically for inspect_df()
# and assumption_report(). Other ambiguous inputs receive a plain-English tip.
# ================================================================


# ---------------------------------------------------------------
# 1. PACKAGE SETUP
# ---------------------------------------------------------------

load_biostats_packages <- function(install_missing = FALSE,
                                   attach = TRUE,
                                   packages = c(
                                     "tidyverse", "car", "moments", "DescTools",
                                     "agricolae", "lme4", "vegan", "readxl",
                                     "rstatix", "afex", "emmeans", "ggpubr"
                                   ),
                                   optional_packages = c(
                                     "performance", "effectsize", "parameters",
                                     "see", "insight", "ggResidpanel"
                                   )) {
  if (!is.logical(install_missing) || length(install_missing) != 1 ||
      is.na(install_missing)) {
    stop("load_biostats_packages(): `install_missing` should be TRUE to install missing packages or FALSE to leave installation alone.", call. = FALSE)
  }
  if (!is.logical(attach) || length(attach) != 1 || is.na(attach)) {
    stop("load_biostats_packages(): `attach` should be TRUE to make installed packages' functions available by name, or FALSE to only check availability.", call. = FALSE)
  }
  if (!is.character(packages) || length(packages) == 0 || anyNA(packages) ||
      any(!nzchar(packages))) {
    stop("load_biostats_packages(): `packages` should be a character vector of package names, such as c(\"tidyverse\", \"car\").", call. = FALSE)
  }
  packages <- unique(packages)
  if (!is.character(optional_packages) || anyNA(optional_packages) ||
      any(!nzchar(optional_packages))) {
    stop("load_biostats_packages(): `optional_packages` should be a character vector of package names, or character(0) if you do not want to check optional add-ons.", call. = FALSE)
  }
  all_packages <- unique(c(packages, optional_packages))
  missing_packages <- all_packages[
    !vapply(all_packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0 && install_missing) {
    install.packages(missing_packages, repos = "https://cloud.r-project.org")
    missing_packages <- missing_packages[
      !vapply(missing_packages, requireNamespace, logical(1), quietly = TRUE)
    ]
  }

  attached_packages <- character(0)
  attach_failures <- character(0)
  if (attach) {
    available_to_attach <- setdiff(packages, missing_packages)
    for (package_name in available_to_attach) {
      attached <- tryCatch({
        suppressPackageStartupMessages(
          library(package_name, character.only = TRUE)
        )
        TRUE
      }, error = function(e) {
        attach_failures <<- c(attach_failures,
                              paste0(package_name, " (", conditionMessage(e), ")"))
        FALSE
      })
      if (attached) attached_packages <- c(attached_packages, package_name)
    }
  }

  if (length(missing_packages) == 0 && length(attach_failures) == 0) {
    if (attach) {
      message("Biostatistics toolkit ready. Attached packages: ",
              paste(attached_packages, collapse = ", "), ".")
    } else {
      message("Biostatistics toolkit ready. All requested packages are installed; none were attached.")
    }
  } else {
    if (length(missing_packages) > 0) {
      message("Not installed: ", paste(missing_packages, collapse = ", "),
              ". To install these, run install_biostats_packages().")
    }
    if (length(attach_failures) > 0) {
      message("Installed but could not be attached: ",
              paste(attach_failures, collapse = "; "),
              ". Restart R and check the package installation if this repeats.")
    }
    if (length(attached_packages) > 0) {
      message("Packages attached and ready: ",
              paste(attached_packages, collapse = ", "), ".")
    }
    message("The toolkit itself is ready. A missing optional package only affects the functions that use it.")
  }

  invisible(list(installed = setdiff(all_packages, missing_packages),
                 missing = missing_packages,
                 attached = attached_packages,
                 attach_failures = attach_failures))
}

install_biostats_packages <- function() {
  load_biostats_packages(install_missing = TRUE)
}


# ---------------------------------------------------------------
# 1B. FRIENDLY INPUT CHECKS
# ---------------------------------------------------------------

# Turn common table-shaped inputs into a data frame when that is safe.
# A matrix has rows and columns, so converting it preserves its layout.
# A single vector or an arbitrary list is ambiguous, so we explain it instead
# of guessing which values are observations or which items are variables.
.toolkit_prepare_data <- function(data, function_name = "This function",
                                  allow_vector = FALSE) {
  if (is.data.frame(data)) return(data)

  if (is.matrix(data)) {
    message("I recognized your input as a matrix (a table where every column has the same basic type). ",
            "I converted it to a data frame so the function can use it.")
    return(as.data.frame(data, stringsAsFactors = FALSE))
  }

  if (allow_vector && is.null(dim(data)) && length(data) > 0 &&
      (is.numeric(data) || is.logical(data) || is.factor(data) ||
       (is.character(data) && length(data) > 1))) {
    message("I recognized your input as one column of values and wrapped it in a one-column table named `value`.")
    return(data.frame(value = data, stringsAsFactors = FALSE))
  }

  if (is.list(data) && length(data) > 0 && !is.null(names(data)) &&
      !anyNA(names(data)) && all(nzchar(names(data))) &&
      all(vapply(data, function(column) {
        is.atomic(column) && is.null(dim(column))
      }, logical(1))) && length(unique(lengths(data))) == 1) {
    message("I recognized this as a named list of equal-length columns and converted it to a data frame.")
    return(as.data.frame(data, stringsAsFactors = FALSE))
  }

  object_class <- paste(class(data), collapse = "/")
  object_type <- typeof(data)
  if (is.null(data)) {
    explanation <- "The input is empty (NULL). This often means an object was not assigned, or a function returned no result."
    suggestion <- "Load or create the dataset first, then give the function the dataset's name."
  } else if (is.function(data)) {
    explanation <- "The input is a function, not the data produced by that function."
    suggestion <- "Run the function with parentheses and pass the result, or pass the dataset object it should use."
  } else if (inherits(data, "lm") || inherits(data, "aov")) {
    explanation <- "The input is a fitted statistical model, not the original data table."
    suggestion <- "For a new analysis, pass the original dataset. For model diagnostics, use enhanced_model_report() with the fitted model."
  } else if (is.character(data) && length(data) == 1) {
    explanation <- paste0("The input is text (\"", data, "\"), not the contents of a dataset. This may be a file name or a column name.")
    suggestion <- "To import a CSV, use read.csv(file.choose()) and then pass the name you assigned to the result."
  } else if (is.atomic(data) && is.null(dim(data))) {
    explanation <- paste0("The input is one vector (one column of values), not a complete data table. Its R type is ", object_type, ".")
    suggestion <- "If you only have one measurement column, make it a one-column table with data.frame(measurement = your_vector). If the analysis needs groups or other variables, pass the full dataset instead."
  } else if (is.list(data)) {
    explanation <- "The input is a list. A list can hold columns, but it can also hold unrelated things, so I cannot safely guess how to turn it into a dataset."
    suggestion <- "If each list item is a column of equal length, use as.data.frame(your_list). Otherwise, pass the dataset you imported or built."
  } else {
    explanation <- paste0("The input has class ", object_class, " and R type ", object_type, ", which this function cannot use as a dataset.")
    suggestion <- "Pass a data frame (a spreadsheet-like object with rows and columns), or a matrix that can be converted to one."
  }

  stop(function_name, " needs a dataset table.\n\n",
       "What I received: ", explanation, "\n\n",
       "What to try: ", suggestion, call. = FALSE)
}

.toolkit_check_alpha <- function(alpha, function_name = "This function") {
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop(function_name, ": `alpha` is the cutoff used for a p-value decision. Set it to one number greater than 0 and less than 1, usually 0.05.",
         call. = FALSE)
  }
  invisible(alpha)
}

.toolkit_require_inputs <- function(function_name, inputs,
                                    descriptions = NULL) {
  missing_names <- names(inputs)[vapply(inputs, function(value) {
    is.null(value) || length(value) == 0
  }, logical(1))]
  if (length(missing_names) > 0) {
    labels <- missing_names
    if (!is.null(descriptions)) {
      described <- descriptions[missing_names]
      has_description <- !is.na(described) & nzchar(described)
      labels[has_description] <- paste0(missing_names[has_description],
                                         " (", described[has_description], ")")
    }
    stop(function_name, " cannot run yet. Please provide the important missing information: ",
         paste(labels, collapse = "; "), ".", call. = FALSE)
  }
  invisible(TRUE)
}

.toolkit_check_column <- function(data, variable, argument_name = "variable",
                                 function_name = "This function") {
  if (!is.character(variable) || length(variable) != 1 || is.na(variable) ||
      !nzchar(variable)) {
    stop(function_name, ": `", argument_name, "` should be the name of one column written in quotation marks, such as \"RatDnsty\".",
         call. = FALSE)
  }
  if (!variable %in% names(data)) {
    if (grepl("^YOUR_", variable, ignore.case = TRUE)) {
      stop(function_name, ": \"", variable, "\" is an example placeholder, not a column in your worksheet. Replace it with the exact column name shown in the table. For example, use outcome = \"", names(data)[1], "\" only if that column is truly the result your question asks about.", call. = FALSE)
    }
    available <- if (length(names(data)) == 0) "(the dataset has no named columns)" else paste(names(data), collapse = ", ")
    stop(function_name, ": I could not find a column named \"", variable, "\".\n",
         "Columns I can see: ", available, ".\n",
         "Check the spelling and capitalization, or choose one of the column names shown above.",
         call. = FALSE)
  }
  invisible(variable)
}

.toolkit_is_id_like <- function(name, x) {
  nonmissing <- x[!is.na(x)]
  clean_name <- tolower(gsub("[ .-]+", "_", name))
  name_hint <- grepl(
    "(^id$|_id$|^(subject|individual|sample|specimen|tree|plant|plot|transect|quadrat)(_?(id|number|no|code))?$)",
    clean_name
  )
  unique_values <- unique(nonmissing)
  sequence_id <- FALSE
  if (is.numeric(x) && length(nonmissing) >= 3 &&
      length(unique_values) == length(nonmissing)) {
    ordered_values <- sort(as.numeric(unique_values))
    sequence_id <- identical(ordered_values, as.numeric(seq_along(ordered_values))) ||
      identical(ordered_values, as.numeric(seq_along(ordered_values) - 1L))
  }
  many_unique_labels <- (is.character(x) || is.factor(x)) &&
    length(nonmissing) >= 5 && length(unique_values) >= 0.9 * length(nonmissing)
  name_hint || sequence_id || many_unique_labels
}

.toolkit_paired_wide_hint <- function(data) {
  id_columns <- names(data)[vapply(names(data), function(name) {
    .toolkit_is_id_like(name, data[[name]])
  }, logical(1))]
  numeric_columns <- names(data)[vapply(data, is.numeric, logical(1))]
  measurement_columns <- setdiff(numeric_columns, id_columns)
  if (nrow(data) >= 3 && length(measurement_columns) == 2 &&
      sum(is.finite(data[[measurement_columns[1]]]) &
            is.finite(data[[measurement_columns[2]]])) >= 2) {
    return(list(id = id_columns, measurements = measurement_columns))
  }
  NULL
}

add_paired_difference <- function(data = NULL,
                                  before = NULL,
                                  after = NULL,
                                  name = NULL) {
  .toolkit_require_inputs("add_paired_difference()",
                          list(data = data, before = before, after = after),
                          c(data = "the worksheet", before = "the first condition's numeric measurement column", after = "the second condition's numeric measurement column"))
  data <- .toolkit_prepare_data(data, "add_paired_difference()")
  .toolkit_check_column(data, before, "before", "add_paired_difference()")
  .toolkit_check_column(data, after, "after", "add_paired_difference()")
  if (identical(before, after)) {
    stop("add_paired_difference(): `before` and `after` must name two different columns.", call. = FALSE)
  }
  for (column in c(before, after)) {
    if (!is.numeric(data[[column]])) {
      stop("add_paired_difference(): \"", column,
           "\" is not stored as a numeric measurement. Check how the worksheet was imported before converting it.", call. = FALSE)
    }
  }
  if (is.null(name)) {
    name <- paste0("Difference_", make.names(after), "_minus_", make.names(before))
  }
  if (!is.character(name) || length(name) != 1 || is.na(name) || !nzchar(name)) {
    stop("add_paired_difference(): `name` should be one non-empty column name, such as \"Difference_After_minus_Before\".", call. = FALSE)
  }
  if (name %in% names(data)) {
    stop("add_paired_difference(): A column named \"", name,
         "\" already exists. Choose another name so I do not replace existing data.", call. = FALSE)
  }
  data[[name]] <- data[[after]] - data[[before]]
  cat("Added column:", name, "=", after, "minus", before, "\n")
  cat("A positive value means the measurement was higher in ", after,
      " than in ", before, ". Missing values stay missing.\n", sep = "")
  cat("Use this only if each row's two measurements belong to the same individual or matched pair.\n")
  data
}

add_product_variable <- function(data = NULL,
                                 first = NULL,
                                 second = NULL,
                                 name = NULL) {
  .toolkit_require_inputs("add_product_variable()",
                          list(data = data, first = first, second = second),
                          c(data = "the worksheet", first = "the first numeric column", second = "the second numeric column"))
  data <- .toolkit_prepare_data(data, "add_product_variable()")
  .toolkit_check_column(data, first, "first", "add_product_variable()")
  .toolkit_check_column(data, second, "second", "add_product_variable()")
  for (column in c(first, second)) {
    if (!is.numeric(data[[column]])) {
      stop("add_product_variable(): \"", column,
           "\" is not stored as numeric. Confirm that it contains measurements before converting it.", call. = FALSE)
    }
  }
  if (is.null(name)) {
    name <- paste0("Product_", make.names(first), "_times_", make.names(second))
  }
  if (!is.character(name) || length(name) != 1 || is.na(name) || !nzchar(name)) {
    stop("add_product_variable(): `name` should be one non-empty column name, such as \"Length_times_Width\".", call. = FALSE)
  }
  if (name %in% names(data)) {
    stop("add_product_variable(): A column named \"", name,
         "\" already exists. Choose another name so I do not replace existing data.", call. = FALSE)
  }
  data[[name]] <- data[[first]] * data[[second]]
  cat("Added column:", name, "=", first, "times", second, "\n")
  cat("Use this only when the problem defines the measured quantity as this product (for example, a formula such as area = length × width). The product's units are the units of both columns multiplied together.\n")
  cat("If the problem asks whether one predictor's effect changes with another predictor, that is a regression interaction instead; use a formula such as outcome ~ predictor1 * predictor2.\n")
  data
}

# Add a defined per-unit measure, such as protein / tissue mass.
add_ratio_variable <- function(data = NULL, numerator = NULL, denominator = NULL,
                               name = "Ratio", scale = 1) {
  .toolkit_require_inputs("add_ratio_variable()",
                          list(data = data, numerator = numerator,
                               denominator = denominator),
                          c(data = "the worksheet", numerator = "the numerator column",
                            denominator = "the denominator column"))
  data <- .toolkit_prepare_data(data, "add_ratio_variable()")
  .toolkit_check_column(data, numerator, "numerator", "add_ratio_variable()")
  .toolkit_check_column(data, denominator, "denominator", "add_ratio_variable()")
  if (!is.numeric(data[[numerator]]) || !is.numeric(data[[denominator]]))
    stop("add_ratio_variable(): Both columns must contain numeric measurements. Check the imported values before converting text to numbers.", call. = FALSE)
  if (!is.numeric(scale) || length(scale) != 1 || is.na(scale) || !is.finite(scale))
    stop("add_ratio_variable(): `scale` must be one finite number. Use 1 for a ratio or 100 for a percentage.", call. = FALSE)
  if (!is.character(name) || length(name) != 1 || is.na(name) || !nzchar(name) || name %in% names(data))
    stop("add_ratio_variable(): Choose one new, non-empty column name that is not already in the worksheet.", call. = FALSE)
  zero <- !is.na(data[[denominator]]) & data[[denominator]] == 0
  data[[name]] <- scale * data[[numerator]] / data[[denominator]]
  data[[name]][zero] <- NA_real_
  cat("Added", name, "=", scale, "×", numerator, "/", denominator, "\n")
  if (any(zero)) cat(sum(zero), "row(s) had a zero denominator; their new values are NA, because division by zero is undefined.\n")
  cat("Use this only when the problem defines this numerator/denominator. The result is per unit of the denominator; check the units.\n")
  data
}

# A common field-biology height formula; R's tan() expects radians.
add_tree_height <- function(data = NULL, distance = NULL, angle = NULL,
                            observer_height = NULL, name = "TreeHeight",
                            angle_unit = "degrees") {
  .toolkit_require_inputs("add_tree_height()", list(data = data, distance = distance,
      angle = angle, observer_height = observer_height),
      c(data = "the worksheet", distance = "horizontal distance to the tree",
        angle = "measured angle", observer_height = "observer height"))
  data <- .toolkit_prepare_data(data, "add_tree_height()")
  for (arg in c("distance", "angle", "observer_height")) {
    column <- get(arg)
    .toolkit_check_column(data, column, arg, "add_tree_height()")
    if (!is.numeric(data[[column]])) stop("add_tree_height(): Column '", column, "' must be numeric. Check the imported measurements.", call. = FALSE)
  }
  if (!is.character(angle_unit) || length(angle_unit) != 1 || is.na(angle_unit) || !angle_unit %in% c("degrees", "radians")) stop("add_tree_height(): Set angle_unit to 'degrees' or 'radians' according to how the angle was recorded.", call. = FALSE)
  if (!is.character(name) || length(name) != 1 || is.na(name) || !nzchar(name) || name %in% names(data)) stop("add_tree_height(): Choose one new column name that is not already used.", call. = FALSE)
  angle_value <- data[[angle]]
  if (angle_unit == "degrees") angle_value <- angle_value * pi / 180
  data[[name]] <- data[[distance]] * tan(angle_value) + data[[observer_height]]
  cat("Added", name, "= distance × tan(angle) + observer height. Angle was interpreted as", angle_unit, "and converted to radians for R's tan().\n")
  cat("Confirm this is the formula and measurement setup specified by your lab; use the same units for distance and observer height.\n")
  data
}

# Compute a Shannon diversity index from long-format counts/abundances.
shannon_by_sample <- function(data = NULL, sample_unit = NULL,
                              category = NULL, abundance = NULL,
                              group = NULL) {
  .toolkit_require_inputs("shannon_by_sample()", list(data = data,
      sample_unit = sample_unit, category = category, abundance = abundance),
      c(data = "the worksheet", sample_unit = "the independent sample/transect ID",
        category = "the species/category column", abundance = "the count or abundance column"))
  data <- .toolkit_prepare_data(data, "shannon_by_sample()")
  for (arg in c("sample_unit", "category", "abundance", if (!is.null(group)) "group")) {
    column <- get(arg)
    .toolkit_check_column(data, column, arg, "shannon_by_sample()")
  }
  if (!is.numeric(data[[abundance]]) || any(data[[abundance]] < 0, na.rm = TRUE))
    stop("shannon_by_sample(): Abundance must be numeric counts/amounts at least zero. Negative values cannot be used.", call. = FALSE)
  if (anyNA(data[[sample_unit]]))
    stop("shannon_by_sample(): Some rows have no sample ID. Fill or remove those rows so every category belongs to the correct independent sample.", call. = FALSE)
  if (!is.null(group)) {
    group_counts <- tapply(as.character(data[[group]]), as.character(data[[sample_unit]]), function(x) length(unique(x[!is.na(x)])))
    if (any(group_counts > 1))
      stop("shannon_by_sample(): At least one sample ID is assigned to more than one group. Check the sample IDs and group column before calculating.", call. = FALSE)
  }
  key <- interaction(data[[sample_unit]], drop = TRUE, lex.order = TRUE)
  out <- lapply(split(seq_len(nrow(data)), key), function(rows) {
    a <- data[[abundance]][rows]
    cats <- data[[category]][rows]
    valid <- !is.na(a) & !is.na(cats)
    totals <- tapply(a[valid], as.character(cats[valid]), sum)
    totals <- totals[totals > 0]
    if (!length(totals) || sum(totals) <= 0) h <- NA_real_ else {
      p <- totals / sum(totals)
      h <- -sum(p * log(p))
    }
    ans <- data.frame(sample = as.character(data[[sample_unit]][rows[1]]),
                      Shannon_H = h, stringsAsFactors = FALSE)
    if (!is.null(group)) ans[[group]] <- as.character(data[[group]][rows[1]])
    ans
  })
  result <- do.call(rbind, out)
  names(result)[1] <- sample_unit
  rownames(result) <- NULL
  cat("Calculated Shannon H for", nrow(result), "sample unit(s), using H = -sum(p × log(p)).\n")
  cat("Each output row is one sample unit. Compare groups only when these sample units are independent replicates; do not treat species rows as replicates. Confirm the abundance definition and the scale your course expects.\n")
  result
}

# Descriptive spread summaries; these do not themselves test equal variance.
dispersion_summary <- function(data = NULL, response = NULL, group = NULL) {
  .toolkit_require_inputs("dispersion_summary()", list(data = data,
      response = response, group = group),
      c(data = "the worksheet", response = "the numeric measurement", group = "the grouping column"))
  data <- .toolkit_prepare_data(data, "dispersion_summary()")
  .toolkit_check_column(data, response, "response", "dispersion_summary()")
  .toolkit_check_column(data, group, "group", "dispersion_summary()")
  if (!is.numeric(data[[response]])) stop("dispersion_summary(): The response must be numeric. Check the column meaning and import format.", call. = FALSE)
  groups <- split(data[[response]], data[[group]], drop = TRUE)
  result <- do.call(rbind, lapply(names(groups), function(g) {
    x <- groups[[g]]; x <- x[is.finite(x)]; n <- length(x)
    m <- if (n) mean(x) else NA_real_; s <- if (n > 1) sd(x) else NA_real_
    data.frame(group = g, n = n, mean = m,
      variance = if (n > 1) var(x) else NA_real_, SD = s,
      CV_percent = if (is.finite(m) && m > 0 && is.finite(s)) 100*s/m else NA_real_)
  }))
  names(result)[1] <- group; rownames(result) <- NULL
  print(result, row.names = FALSE)
  cat("Variance describes absolute spread; CV describes spread relative to a positive mean. CV is shown as NA when a group's mean is zero or negative because it is not useful there. These are summaries, not a test that group variances differ. Use plots and, if required, Bartlett (normal data) or Levene (more robust) for that question.\n")
  invisible(result)
}

# Implements the course-style relative male/female contribution index.
# Use only when this exact formula is defined by the assignment.
add_gender_contribution_index <- function(data = NULL, male_count = NULL,
                                          female_count = NULL,
                                          name = "GenderIndex") {
  .toolkit_require_inputs("add_gender_contribution_index()",
    list(data = data, male_count = male_count, female_count = female_count),
    c(data = "the worksheet", male_count = "male count per individual",
      female_count = "female count per individual"))
  data <- .toolkit_prepare_data(data, "add_gender_contribution_index()")
  for (arg in c("male_count", "female_count")) {
    column <- get(arg)
    .toolkit_check_column(data, column, arg, "add_gender_contribution_index()")
    if (!is.numeric(data[[column]]) || any(data[[column]] < 0, na.rm = TRUE))
      stop("add_gender_contribution_index(): Both columns must be nonnegative numeric counts. Check the worksheet and assignment definitions.", call. = FALSE)
  }
  if (!is.character(name) || length(name) != 1 || is.na(name) || !nzchar(name) || name %in% names(data))
    stop("add_gender_contribution_index(): Choose one new column name that is not already used.", call. = FALSE)
  M <- sum(data[[male_count]], na.rm = TRUE)
  F <- sum(data[[female_count]], na.rm = TRUE)
  if (M == 0 || F == 0) stop("add_gender_contribution_index(): The worksheet has a total of zero for one sex, so the requested index cannot be calculated.", call. = FALSE)
  male_share <- data[[male_count]] / M
  female_share <- data[[female_count]] / F
  denom <- male_share + female_share
  data[[name]] <- ifelse(denom > 0, male_share / denom, NA_real_)
  cat("Added", name, "using each individual's share of all male and all female counts.\n")
  cat("Use only if the assignment defines this exact index. A value near 0.5 represents equal relative contribution; it is not automatically a general measure of biological sex or fitness.\n")
  data
}

.toolkit_clean_numeric_vector <- function(x, function_name, argument_name = "x") {
  if (is.data.frame(x)) {
    if (ncol(x) != 1) {
      stop(function_name, ": `", argument_name,
           "` is a data frame with more than one column. Select one numeric column, such as mydata$measurement.", call. = FALSE)
    }
    x <- x[[1]]
  } else if (is.matrix(x)) {
    if (ncol(x) != 1) {
      stop(function_name, ": `", argument_name,
           "` is a matrix with more than one column. Select one numeric column.", call. = FALSE)
    }
    x <- x[, 1]
  }
  if (!is.numeric(x) || !is.null(dim(x))) {
    stop(function_name, ": `", argument_name,
         "` must be one numeric column/vector. If you have a data frame, select its numeric column with mydata$column_name. If numbers were imported as text, inspect them before converting.", call. = FALSE)
  }
  usable <- is.finite(x)
  removed <- sum(!usable)
  x <- x[usable]
  if (removed > 0) {
    message(function_name, ": I skipped ", removed,
            " missing or non-finite value(s) and used the remaining numeric values.")
  }
  x
}


# ---------------------------------------------------------------
# 2. DATASET INSPECTION
# ---------------------------------------------------------------

inspect_df <- function(data = NULL,
                       digits = 3,
                       print_report = TRUE,
                       variable_types = NULL) {
  data <- .toolkit_prepare_data(data, "inspect_df()", allow_vector = TRUE)
  if (nrow(data) == 0 || ncol(data) == 0) {
    stop("inspect_df(): The table has no rows or no columns, so there is nothing to summarize. Check that the correct dataset was imported.", call. = FALSE)
  }
  if (!is.numeric(digits) || length(digits) != 1 || is.na(digits) ||
      !is.finite(digits) || digits < 0 || digits != floor(digits)) {
    stop("inspect_df(): `digits` controls how many decimal places to show. Use one whole number, such as 2 or 3.", call. = FALSE)
  }
  if (!is.logical(print_report) || length(print_report) != 1 || is.na(print_report)) {
    stop("inspect_df(): `print_report` should be TRUE to display the report or FALSE to return the summary quietly.", call. = FALSE)
  }

  if (!is.null(variable_types)) {
    if (!is.character(variable_types) || is.null(names(variable_types)) ||
        anyNA(names(variable_types)) || any(!nzchar(names(variable_types)))) {
      stop("inspect_df(): `variable_types` is an optional way to label columns. Give it names and use text labels, for example c(treatment = \"nominal\").", call. = FALSE)
    }
    unknown_names <- setdiff(names(variable_types), names(data))
    if (length(unknown_names) > 0) {
      stop("inspect_df(): I could not find these column(s) listed in `variable_types`: ",
           paste(unknown_names, collapse = ", "), ". Check spelling and capitalization.", call. = FALSE)
    }
    unknown_types <- setdiff(variable_types, c("nominal", "ordinal"))
    if (length(unknown_types) > 0) {
      stop("inspect_df(): Column labels in `variable_types` must be \"nominal\" (categories with no order) or \"ordinal\" (categories with an order).", call. = FALSE)
    }
  }

  classify_variable <- function(variable, x) {
    if (!is.null(variable_types) && variable %in% names(variable_types)) {
      return(c("categorical", variable_types[[variable]]))
    }
    if (.toolkit_is_id_like(variable, x)) return(c("identifier", "ID-like"))
    if (is.ordered(x)) return(c("categorical", "ordinal"))
    if (is.factor(x) || is.character(x) || is.logical(x)) {
      return(c("categorical", "nominal"))
    }
    if (is.numeric(x)) return(c("quantitative", "numeric"))
    c("other", "other")
  }

  summary_table <- do.call(rbind, lapply(names(data), function(variable) {
    x <- data[[variable]]
    classification <- classify_variable(variable, x)

    if (is.numeric(x)) {
      finite_x <- x[is.finite(x)]
      n <- length(finite_x)
      mean_value <- if (n > 0) mean(finite_x) else NA_real_
      sd_value <- if (n > 1) sd(finite_x) else NA_real_
      standard_error_value <- if (n > 1) sd_value / sqrt(n) else NA_real_

      data.frame(
        variable = variable,
        category = classification[1],
        measurement_scale = classification[2],
        class = paste(class(x), collapse = "/"),
        n = n,
        missing = sum(is.na(x)),
        unique_values = length(unique(x[!is.na(x)])),
        minimum = if (n > 0) min(finite_x) else NA_real_,
        maximum = if (n > 0) max(finite_x) else NA_real_,
        mean = mean_value,
        median = if (n > 0) median(finite_x) else NA_real_,
        sd = sd_value,
        variance = if (n > 1) var(finite_x) else NA_real_,
        standard_error = standard_error_value,
        SEM = standard_error_value,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        variable = variable,
        category = classification[1],
        measurement_scale = classification[2],
        class = paste(class(x), collapse = "/"),
        n = sum(!is.na(x)),
        missing = sum(is.na(x)),
        unique_values = length(unique(x[!is.na(x)])),
        minimum = NA_real_, maximum = NA_real_, mean = NA_real_,
        median = NA_real_, sd = NA_real_, variance = NA_real_,
        standard_error = NA_real_, SEM = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }))

  numeric_names <- c(
    "minimum", "maximum", "mean", "median", "sd", "variance",
    "standard_error", "SEM"
  )
  summary_table[numeric_names] <- lapply(
    summary_table[numeric_names],
    function(x) round(x, digits)
  )

  numeric_text_columns <- names(data)[vapply(data, function(x) {
    if (!is.character(x) && !is.factor(x)) return(FALSE)
    values <- unique(as.character(na.omit(x)))
    if (length(values) == 0) return(FALSE)
    converted <- suppressWarnings(as.numeric(values))
    all(!is.na(converted))
  }, logical(1))]

  if (print_report) {
    cat("\n====================================\n")
    cat("DATASET INSPECTION\n")
    cat("====================================\n")
    cat("Rows:", nrow(data), " | Columns:", ncol(data), "\n")
    cat("Each row is an observation. Each column is one variable.\n\n")
    cat("RESPONSE VARIABLE: NOT IDENTIFIED BY THIS INSPECTION\n")
    cat("This report describes columns; it cannot know which result your question asks about. Use analysis_start(data, question = \"paste the problem question here\") to get a route and wording-based preparation clues.\n")
    cat("\n")

    cat("HOW TO READ THE TABLES\n")
    cat("- class = how R stores the column (for example, numeric, character, or factor).\n")
    cat("- n = number of non-missing observations used for that column.\n")
    cat("- missing = number of blank/NA entries.\n")
    cat("- unique_values = number of different values recorded.\n")
    cat("- Numeric summaries (minimum through SEM) are only meaningful for numeric columns.\n\n")

    show_section <- function(title, description, rows) {
      cat(title, "\n", sep = "")
      cat(description, "\n")
      if (nrow(rows) == 0) {
        cat("(none found)\n")
      } else {
        print(rows, row.names = FALSE)
      }
    }

    show_section(
      "1. CATEGORICAL VARIABLES WITH NO ORDER (NOMINAL)",
      "Labels without a ranking, such as treatment, species, or location.",
      summary_table[summary_table$measurement_scale == "nominal",
                    c("variable", "class", "n", "missing", "unique_values")]
    )

    show_section(
      "\n2. CATEGORICAL VARIABLES WITH AN ORDER (ORDINAL)",
      "Ordered labels, such as low / medium / high.",
      summary_table[summary_table$measurement_scale == "ordinal",
                    c("variable", "class", "n", "missing", "unique_values")]
    )

    show_section(
      "\n3. POSSIBLE ROW / SUBJECT IDENTIFIERS",
      "These values appear to label rows or individuals. Confirm with your notes; do not analyze an ID as a measurement.",
      summary_table[summary_table$category == "identifier",
                    c("variable", "class", "n", "missing", "unique_values")]
    )

    numeric_rows <- summary_table[summary_table$category == "quantitative", ]
    cat("\n4. NUMERIC MEASUREMENTS OR COUNTS\n")
    if (nrow(numeric_rows) == 0) {
      cat("(none found)\n")
    } else {
      show_section(
        "\nA. HOW MANY VALUES?",
        "n = usable values; missing = blank/NA values; unique_values = distinct values.",
        numeric_rows[, c("variable", "class", "n", "missing", "unique_values")]
      )
      show_section(
        "\nB. WHERE ARE THE VALUES?",
        "minimum / maximum = ends of the data; mean = average; median = middle value.",
        numeric_rows[, c("variable", "minimum", "maximum", "mean", "median")]
      )
      show_section(
        "\nC. HOW SPREAD OUT ARE THE VALUES?",
        "sd = typical spread; variance = spread squared; standard_error / SEM = uncertainty in the sample mean.",
        numeric_rows[, c("variable", "sd", "variance", "standard_error", "SEM")]
      )
    }

    paired_hint <- .toolkit_paired_wide_hint(data)
    if (!is.null(paired_hint)) {
      before_column <- paired_hint$measurements[1]
      after_column <- paired_hint$measurements[2]
      cat("\nPOSSIBLE TWO-CONDITION / PAIRED LAYOUT\n")
      cat("I found two numeric measurement columns: ", before_column, " and ",
          after_column, ". This may be one measurement per condition for the same individuals, but the worksheet cannot prove that.\n", sep = "")
      if (length(paired_hint$id) > 0) {
        cat("Possible ID column(s): ", paste(paired_hint$id, collapse = ", "),
            ". A sequence such as 1 to 15 is likely individual labels, not a measured variable; confirm that from the study description.\n", sep = "")
      }
      cat("If each row really is the same individual measured in both conditions, you can add a signed difference column (second condition minus first):\n")
      cat("  data <- add_paired_difference(data, before = \"", before_column,
          "\", after = \"", after_column, "\")\n", sep = "")
      cat("  View(data)\n")
      cat("The added column is for viewing/reporting; the paired test can use the two original columns directly:\n")
      cat("  paired_t_test(data, before = \"", before_column,
          "\", after = \"", after_column, "\")\n", sep = "")
      cat("A positive difference means the second condition's measurement was higher. If rows are not matched by individual, do not use this paired approach. inspect_df() only suggests this pattern; it does not change your data.\n")
    }

    cat("\nNote: standard_error and SEM are two names for the same calculation here.\n")
    if (length(numeric_text_columns) > 0) {
      cat("\nCHECK THESE COLUMNS\n")
      cat("These column(s) contain text that looks numeric: ",
          paste(numeric_text_columns, collapse = ", "), ".\n", sep = "")
      cat("If they are measurements, they may need conversion to numbers. If they are IDs or category codes, keep them as labels.\n")
      cat("To inspect one, use str(your_data[[\"column_name\"]]). Do not convert until you know what the values represent.\n")
      cat("For a confirmed numeric measurement, the conversion is: your_data[[\"column_name\"]] <- as.numeric(as.character(your_data[[\"column_name\"]])).\n")
    }
  }

  invisible(summary_table)
}


# ---------------------------------------------------------------
# 3. TRANSFORMATIONS
# ---------------------------------------------------------------

add_transformations <- function(data = NULL, variable = NULL, save_to = NULL) {
  .toolkit_require_inputs("add_transformations()",
                          list(data = data, variable = variable),
                          c(data = "the imported worksheet", variable = "the numeric column to transform"))
  data <- .toolkit_prepare_data(data, "add_transformations()")
  .toolkit_check_column(data, variable, function_name = "add_transformations()")
  if (!is.numeric(data[[variable]])) {
    stop("add_transformations(): The column \"", variable, "\" is stored as ",
         paste(class(data[[variable]]), collapse = "/"), ", not numeric.\n",
         "These transformations need numbers. If this column contains numbers stored as text, inspect it with str(data[[\"",
         variable, "\"]]) and convert it only after confirming every value is really numeric. If it contains group names, choose a measurement column instead.",
         call. = FALSE)
  }

  x <- data[[variable]]

  sqrt_values <- rep(NA_real_, length(x))
  nonnegative <- !is.na(x) & x >= 0
  sqrt_values[nonnegative] <- sqrt(x[nonnegative])

  log_values <- rep(NA_real_, length(x))
  positive <- !is.na(x) & x > 0
  log_values[positive] <- log(x[positive])

  log10_values <- rep(NA_real_, length(x))
  log10_values[positive] <- log10(x[positive])

  data[[paste0(variable, "_sqrt")]] <- sqrt_values
  data[[paste0(variable, "_square")]] <- x^2
  data[[paste0(variable, "_log")]] <- log_values
  data[[paste0(variable, "_log10")]] <- log10_values

  if (!is.null(save_to)) {
    if (!is.character(save_to) || length(save_to) != 1 || is.na(save_to) ||
        !nzchar(save_to)) {
      stop("add_transformations(): `save_to` should be one file path in quotation marks, for example \"transformed_data.csv\". Leave it as NULL if you do not want to save a file.", call. = FALSE)
    }
    tryCatch(
      write.csv(data, save_to, row.names = FALSE),
      error = function(e) stop("add_transformations(): I created the transformed data, but could not save it to \"",
                               save_to, "\". Check that the folder exists and that you have permission to save there. R's message was: ",
                               conditionMessage(e), call. = FALSE)
    )
    cat("Updated data saved to:", save_to, "\n")
  }

  data
}


# ---------------------------------------------------------------
# 4. Q-Q PLOTS WITH CONFIDENCE ENVELOPES
# ---------------------------------------------------------------

qq_plot_with_ci <- function(x,
                            conf = 0.95,
                            simulations = 1000,
                            main = "Q-Q Plot with Confidence Envelope") {
  if (!is.numeric(x) || !is.null(dim(x))) {
    stop("qq_plot_with_ci(): `x` needs to be one numeric column or vector. For a data frame column, pass it like mydata$measurement. A whole data frame or a text/category column cannot be plotted this way.", call. = FALSE)
  }
  if (!is.numeric(conf) || length(conf) != 1 || is.na(conf) ||
      conf <= 0 || conf >= 1) {
    stop("qq_plot_with_ci(): `conf` is the confidence level for the blue envelope. Use one number between 0 and 1, usually 0.95.", call. = FALSE)
  }
  if (!is.numeric(simulations) || length(simulations) != 1 ||
      is.na(simulations) || simulations < 100 ||
      simulations != floor(simulations)) {
    stop("qq_plot_with_ci(): `simulations` is how many simulated samples build the envelope. Use a whole number of at least 100, such as 1000.", call. = FALSE)
  }
  if (all(is.na(x))) {
    stop("qq_plot_with_ci(): The numeric input contains no usable values. Check for missing values or whether the selected column is correct.", call. = FALSE)
  }
  x <- x[is.finite(x)]
  n <- length(x)

  if (n < 3) stop("qq_plot_with_ci(): I found fewer than 3 usable numeric values. A Q-Q plot needs at least 3.", call. = FALSE)
  if (sd(x) == 0) stop("qq_plot_with_ci(): All usable values are identical, so there is no spread to compare with a normal distribution.", call. = FALSE)

  theoretical_values <- qnorm(ppoints(n))
  ordered_values <- sort(x)

  simulated_values <- matrix(
    rnorm(n * simulations, mean = mean(x), sd = sd(x)),
    nrow = simulations,
    ncol = n
  )
  simulated_values <- t(apply(simulated_values, 1, sort))

  lower_probability <- (1 - conf) / 2
  upper_probability <- 1 - lower_probability
  lower_limit <- apply(simulated_values, 2, quantile,
                       probs = lower_probability)
  upper_limit <- apply(simulated_values, 2, quantile,
                       probs = upper_probability)

  plot(theoretical_values, ordered_values,
       main = main,
       xlab = "Theoretical normal quantiles",
       ylab = "Observed quantiles",
       pch = 19)
  lines(theoretical_values, lower_limit,
        col = "blue", lty = 2, lwd = 2)
  lines(theoretical_values, upper_limit,
        col = "blue", lty = 2, lwd = 2)
  qqline(ordered_values, col = "red", lwd = 2)

  legend("topleft",
         legend = c("Reference line", paste0(conf * 100, "% envelope")),
         col = c("red", "blue"),
         lty = c(1, 2),
         lwd = 2,
         bty = "n")

  invisible(list(
    theoretical = theoretical_values,
    observed = ordered_values,
    lower = lower_limit,
    upper = upper_limit
  ))
}


# ---------------------------------------------------------------
# 5. GENERAL NORMALITY CHECK
# ---------------------------------------------------------------

assumption_report <- function(data = NULL,
                              response = NULL,
                              group = NULL,
                              alpha = 0.05,
                              plots = TRUE) {
  data <- .toolkit_prepare_data(data, "assumption_report()", allow_vector = TRUE)
  .toolkit_check_alpha(alpha, "assumption_report()")
  if (!is.logical(plots) || length(plots) != 1 || is.na(plots)) {
    stop("assumption_report(): `plots` should be TRUE to make plots or FALSE to skip them.", call. = FALSE)
  }

  if (is.null(response)) {
    response <- names(data)[vapply(data, is.numeric, logical(1))]
    if (length(response) == 0) {
      stop("assumption_report(): I could not find any numeric columns to check. The outcome should be a numeric measurement. Use str(your_data) to inspect how R imported the columns; a column of numbers may have been imported as text.", call. = FALSE)
    }
  }
  if (!is.character(response) || length(response) == 0 || anyNA(response)) {
    stop("assumption_report(): `response` should be one or more column names in quotation marks, such as \"RatDnsty\". Leave it blank to use all numeric columns.", call. = FALSE)
  }
  non_numeric <- response[response %in% names(data) &
                            !vapply(response[response %in% names(data)],
                                    function(name) is.numeric(data[[name]]),
                                    logical(1))]
  if (length(non_numeric) > 0) {
    stop("assumption_report(): These response column(s) are not stored as numbers: ",
         paste(non_numeric, collapse = ", "), ". Choose a numeric measurement column. If the values are numbers stored as text, check them before converting.", call. = FALSE)
  }
  for (response_name in response) {
    .toolkit_check_column(data, response_name, "response", "assumption_report()")
  }

  if (!is.null(group)) {
    .toolkit_check_column(data, group, "group", "assumption_report()")
  }

  results <- list()

  for (variable in response) {
    selected_columns <- c(variable, group)
    working_data <- data[, selected_columns, drop = FALSE]
    working_data <- working_data[complete.cases(working_data), , drop = FALSE]
    values <- working_data[[variable]]

    cat("\n====================================\n")
    cat("Response variable:", variable, "\n")
    cat("Complete observations:", length(values), "\n")
    cat("====================================\n")

    overall_shapiro <- NULL
    if (length(values) >= 3 && length(values) <= 5000) {
      overall_shapiro <- shapiro.test(values)
      cat("\nOverall Shapiro-Wilk test:\n")
      print(overall_shapiro)
    }

    group_shapiro <- NULL
    levene_result <- NULL

    if (!is.null(group)) {
      group_shapiro <- lapply(split(values, working_data[[group]]), function(x) {
        if (length(x) >= 3 && length(x) <= 5000) shapiro.test(x) else NULL
      })

      cat("\nShapiro-Wilk tests by group:\n")
      for (group_name in names(group_shapiro)) {
        cat("\nGroup:", group_name, "\n")
        if (is.null(group_shapiro[[group_name]])) {
          cat("Unavailable: group needs 3 to 5000 observations.\n")
        } else {
          print(group_shapiro[[group_name]])
        }
      }

      if (requireNamespace("car", quietly = TRUE)) {
        levene_data <- data.frame(
          response = values,
          group = as.factor(working_data[[group]])
        )
        levene_result <- car::leveneTest(response ~ group, data = levene_data)
        cat("\nLevene's test for equal variances:\n")
        print(levene_result)
      }
    }

    if (plots && length(values) >= 3) {
      qq_plot_with_ci(values, main = paste("Q-Q Plot:", variable))
    }

    results[[variable]] <- list(
      n = length(values),
      overall_shapiro = overall_shapiro,
      group_shapiro = group_shapiro,
      levene = levene_result
    )
  }

  invisible(results)
}


# ---------------------------------------------------------------
# 6. TRANSFORMATION WORKFLOW
# ---------------------------------------------------------------

transformation_workflow <- function(data = NULL,
                                    variable = NULL,
                                    alpha = 0.05,
                                    plots = TRUE,
                                    save_to = NULL) {
  .toolkit_require_inputs("transformation_workflow()",
                          list(data = data, variable = variable),
                          c(data = "the imported worksheet", variable = "the numeric column to compare across transformations"))
  .toolkit_check_alpha(alpha, "transformation_workflow()")
  transformed_data <- add_transformations(data, variable, save_to = save_to)

  variables_to_check <- c(
    variable,
    paste0(variable, "_sqrt"),
    paste0(variable, "_square"),
    paste0(variable, "_log"),
    paste0(variable, "_log10")
  )

  test_results <- lapply(variables_to_check, function(current_variable) {
    values <- transformed_data[[current_variable]]
    values <- values[is.finite(values)]

    if (length(values) < 3 || length(values) > 5000) {
      return(data.frame(
        variable = current_variable,
        n = length(values),
        shapiro_p_value = NA_real_,
        interpretation = "Shapiro-Wilk unavailable",
        stringsAsFactors = FALSE
      ))
    }

    shapiro_result <- shapiro.test(values)
    interpretation <- if (shapiro_result$p.value < alpha) {
      "Evidence against normality"
    } else {
      "No significant evidence against normality"
    }

    if (plots) {
      qq_plot_with_ci(values, main = paste("Q-Q Plot:", current_variable))
    }

    data.frame(
      variable = current_variable,
      n = length(values),
      shapiro_p_value = shapiro_result$p.value,
      interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  })

  test_results <- do.call(rbind, test_results)
  cat("\nTransformation normality comparison:\n")
  print(test_results, row.names = FALSE)

  invisible(list(
    data = transformed_data,
    normality_results = test_results
  ))
}

compare_transformations <- function(data = NULL,
                                    variable = NULL,
                                    alpha = 0.05,
                                    plots = TRUE,
                                    save_to = NULL) {
  transformation_workflow(
    data = data,
    variable = variable,
    alpha = alpha,
    plots = plots,
    save_to = save_to
  )
}


# ---------------------------------------------------------------
# 7. SHARED MODEL REPORTING
# ---------------------------------------------------------------

enhanced_model_report <- function(model,
                                  diagnostic_plots = TRUE,
                                  effect_sizes = TRUE,
                                  model_checks = TRUE) {
  if (!inherits(model, "lm") && !inherits(model, "aov")) {
    received <- if (is.null(model)) "nothing (NULL)" else paste(class(model), collapse = "/")
    stop("enhanced_model_report(): I received ", received, ", not a fitted linear model.\n",
         "What to try: first fit a model with lm(...) or aov(...), then pass the result, for example enhanced_model_report(my_model).",
         call. = FALSE)
  }
  options <- list(diagnostic_plots = diagnostic_plots,
                  effect_sizes = effect_sizes,
                  model_checks = model_checks)
  for (option_name in names(options)) {
    option_value <- options[[option_name]]
    if (!is.logical(option_value) || length(option_value) != 1 || is.na(option_value)) {
      stop("enhanced_model_report(): `", option_name,
           "` must be TRUE or FALSE. Use TRUE to include that output or FALSE to skip it.", call. = FALSE)
    }
  }

  output <- list()

  if (effect_sizes && requireNamespace("effectsize", quietly = TRUE)) {
    output$eta_squared <- effectsize::eta_squared(model, ci = 0.95)
    output$omega_squared <- effectsize::omega_squared(model, ci = 0.95)

    cat("\nEnhanced effect sizes:\n")
    cat("Eta-squared with 95% confidence intervals:\n")
    print(output$eta_squared)
    cat("\nOmega-squared with 95% confidence intervals:\n")
    print(output$omega_squared)
  }

  if (model_checks && requireNamespace("performance", quietly = TRUE)) {
    output$normality_check <- performance::check_normality(model)
    output$homogeneity_check <- performance::check_heteroscedasticity(model)

    cat("\nEnhanced model checks from performance:\n")
    cat("Normality check:\n")
    print(output$normality_check)
    cat("\nHeteroscedasticity check:\n")
    print(output$homogeneity_check)
  }

  if (diagnostic_plots && requireNamespace("ggResidpanel", quietly = TRUE)) {
    output$residual_panel <- ggResidpanel::resid_panel(
      model,
      plots = c("resid", "qq", "yvp"),
      qqbands = TRUE
    )
    print(output$residual_panel)
  }

  if (length(output) == 0) {
    cat("\nNo enhanced package results were available.\n")
    cat("Run install_biostats_packages() to install the optional packages.\n")
  }

  invisible(output)
}


# ---------------------------------------------------------------
# 8. ANALYSIS GUIDE, CORRELATION, AND REGRESSION
# ---------------------------------------------------------------

# Start here when a new worksheet is imported. This function reports what R
# can see, suggests a test family from selected columns, and lists the study
# details that still need a human answer.
analysis_start <- function(data = NULL,
                           outcome = NULL,
                           predictors = NULL,
                           goal = NULL,
                           design = NULL,
                           question = NULL) {
  data <- .toolkit_prepare_data(data, "analysis_start()")
  if (nrow(data) == 0 || ncol(data) == 0) {
    stop("analysis_start(): This table has no rows or no columns. Recheck that you imported the worksheet you intended.", call. = FALSE)
  }
  valid_goals <- c("difference", "association", "prediction", "goodness_of_fit")
  valid_designs <- c("independent", "paired", "repeated", "one_sample", "unknown")
  if (!is.null(goal) && (!is.character(goal) || length(goal) != 1 ||
                         is.na(goal) || !goal %in% valid_goals)) {
    stop("analysis_start(): `goal` should be one of \"difference\", \"association\", \"prediction\", or \"goodness_of_fit\". Pick the one that matches the wording of the research question.", call. = FALSE)
  }
  if (!is.null(design) && (!is.character(design) || length(design) != 1 ||
                           is.na(design) || !design %in% valid_designs)) {
    stop("analysis_start(): `design` should be \"independent\", \"paired\", \"repeated\", \"one_sample\", or \"unknown\". This comes from how observations were collected.", call. = FALSE)
  }
  if (!is.null(question) && (!is.character(question) || length(question) != 1 ||
                             is.na(question) || !nzchar(trimws(question)))) {
    stop("analysis_start(): `question` is optional. If supplied, write one sentence in quotation marks describing what the problem asks, for example question = \"Does attraction differ between the two light conditions?\".", call. = FALSE)
  }
  looks_like_id <- function(name, x) {
    nonmissing <- x[!is.na(x)]
    clean_name <- tolower(gsub("[ .-]+", "_", name))
    name_hint <- grepl(
      "(^id$|_id$|^(subject|individual|sample|specimen|tree|plant|plot|transect|quadrat)(_?(id|number|no|code))?$)",
      clean_name
    )
    unique_values <- unique(nonmissing)
    sequence_id <- FALSE
    if (is.numeric(x) && length(nonmissing) >= 3 &&
        length(unique_values) == length(nonmissing)) {
      ordered_values <- sort(as.numeric(unique_values))
      sequence_id <- identical(ordered_values, as.numeric(seq_along(ordered_values))) ||
        identical(ordered_values, as.numeric(seq_along(ordered_values) - 1L))
    }
    many_unique_labels <- (is.character(x) || is.factor(x)) &&
      length(nonmissing) >= 5 &&
      length(unique_values) >= 0.9 * length(nonmissing)
    name_hint || sequence_id || many_unique_labels
  }
  describe_column <- function(name, x) {
    n_nonmissing <- sum(!is.na(x))
    unique_n <- length(unique(x[!is.na(x)]))
    examples <- unique(x[!is.na(x)])
    examples <- paste(utils::head(as.character(examples), 3), collapse = ", ")
    id_hint <- looks_like_id(name, x)
    if (inherits(x, "Date") || inherits(x, "POSIXt")) {
      kind <- "date/time"
      interpretation <- "May mark when something happened; decide whether time is part of the question."
    } else if (is.ordered(x)) {
      kind <- "ordered category (rank)"
      interpretation <- "Values have an order; do not assume the gaps between ranks are equal."
    } else if (is.factor(x) || is.character(x) || is.logical(x)) {
      kind <- if (is.logical(x) || unique_n == 2) "category (2 levels)" else "category/label"
      interpretation <- "Use as a group or outcome only if its labels match the research question."
      if (is.character(x) && n_nonmissing > 0) {
        numeric_text <- grepl("^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$",
                              trimws(as.character(x[!is.na(x)])))
        if (all(numeric_text)) {
          kind <- "number-like text"
          interpretation <- "Looks like numbers stored as text. Check decimal marks, units, and symbols before converting."
        }
      }
    } else if (is.numeric(x)) {
      if (id_hint) {
        kind <- "possibly an ID (stored as a number)"
        interpretation <- "Likely identifies rows/subjects; usually do not analyze it as a measurement. Confirm from the study notes."
      } else if (unique_n <= 10 && all(x[!is.na(x)] == floor(x[!is.na(x)]))) {
        category_name_hint <- grepl("(group|treat|condition|sex|species|habitat|category|class|level)",
                                   tolower(name))
        if (category_name_hint) {
          kind <- "possible numeric category code"
          interpretation <- "The name and values could indicate numeric labels for groups. Confirm from the study notes; do not convert automatically."
        } else {
          kind <- "whole-number count or measurement"
          interpretation <- "May be a count or a measurement recorded in whole numbers. Use the study's definitions; do not turn it into a category just because it has few values."
        }
      } else {
        kind <- "numeric measurement"
        interpretation <- "Could be an outcome or predictor; check units and whether the values make biological sense."
      }
    } else {
      kind <- paste("other:", paste(class(x), collapse = "/"))
      interpretation <- "Check this column's meaning before choosing a test."
    }
    if (id_hint && !grepl("ID", kind, fixed = TRUE)) {
      interpretation <- paste(interpretation, "Its name or values also look ID-like; verify before analyzing it.")
    }
    data.frame(column = name, R_type = paste(class(x), collapse = "/"),
               suggested_kind = kind, distinct_values = unique_n,
               missing = sum(is.na(x)), examples = examples,
               plain_English_note = interpretation,
               stringsAsFactors = FALSE)
  }
  variable_map <- do.call(rbind, lapply(names(data), function(name) {
    describe_column(name, data[[name]])
  }))
  rownames(variable_map) <- NULL

  cat("\nNEW WORKSHEET: FIRST LOOK\n")
  cat("Rows (observations):", nrow(data), " | Columns (variables):", ncol(data), "\n")
  cat("\nColumn guide (these are clues, not final decisions):\n")
  print(variable_map, row.names = FALSE, right = FALSE)
  cat("\nHOW TO IDENTIFY THE RESPONSE VARIABLE\n")
  cat("The response (also called outcome or dependent variable) is the result the question asks you to measure, explain, compare, or predict. It is the Y / left side in a formula such as Y ~ X. Column type alone cannot tell us this.\n")
  cat("Example: ‘Does plant height differ by soil type?’ -> response = plant height; predictor/group = soil type.\n")
  cat("Example: ‘Are species and habitat associated?’ -> neither is naturally the response; treat both as categories in an association test.\n")
  cat("\nPRE-TEST CALCULATION CHECK\n")
  cat("Look in the problem wording for a defined quantity you must calculate first: a paired difference, product (such as area = length × width), rate/ratio, or total. The worksheet alone cannot tell which formula the biology question intends; do not create a product or ratio just because columns are available.\n")
  cat("A calculated product and a regression interaction are different: if the question asks whether one predictor's effect changes with another, use outcome ~ predictor1 * predictor2 in a regression model; R includes both predictors and their interaction.\n")
  cat("\nQUICK TEST PATTERN CLUES (confirm the question and design)\n")
  cat("- Numeric measurements in two columns for the same subjects -> paired t-test; check the within-subject differences.\n")
  cat("- Numeric measurement across two unrelated groups -> Welch two-sample t-test; across 3+ unrelated groups -> ANOVA.\n")
  cat("- Two numeric measurements -> correlation for association, regression for a slope/equation/prediction.\n")
  cat("- Two categorical variables -> chi-square/G test of association; one set of category counts vs a stated ratio -> goodness-of-fit.\n")
  cat("- Yes/no outcome with predictors -> logistic regression; one yes/no count vs a stated chance -> binomial test.\n")
  cat("- Same subjects across 3+ conditions -> repeated-measures method, not ordinary independent-groups ANOVA.\n")
  cat("- A quantity defined as length × width -> calculate that product first; do not substitute a regression interaction.\n")
  cat("- ‘per’, ‘per unit’, concentration, or a defined index -> check whether the problem defines a ratio, sum, or special formula before testing.\n")
  cat("- Species/categories repeated within each transect/sample -> calculate one diversity index per independent sample, then compare those sample-level values.\n")
  cat("- Several individuals measured within each gall/nest/plot -> individuals may be grouped within a shared unit; ordinary independent-group tests can be wrong. Check for a nested/random-effects method.\n")
  cat("- Counts in categories vs an expected proportion/ratio -> goodness-of-fit; counts for two categorical variables -> association test. First build the observed-count table.\n")
  cat("- Distance and an angle in a tree-height lab -> use the stated trigonometric formula and confirm whether the angle is degrees; R's tan() needs radians.\n")
  if (!is.null(question)) {
    question_lower <- tolower(question)
    cat("Question you supplied: ", question, "\n", sep = "")
    if (grepl("interaction|depends on|depending on|varies with|moderates", question_lower)) {
      cat("QUESTION CLUE: This wording may describe an interaction. If that matches the intended biology, put both predictors and * in the regression formula, for example linear_regression_workflow(data, Y ~ X1 * X2). Do not manually multiply columns unless the assignment specifically requests a calculated product.\n")
    }
    if (grepl("product|multiply|multiplied|area|volume", question_lower)) {
      cat("QUESTION CLUE: The wording may require a calculated product. Confirm which two columns and the exact meaning/units, then use add_product_variable(data, first = \"COLUMN1\", second = \"COLUMN2\", name = \"NewProduct\"). Replace the placeholders with the specified columns and output name.\n")
    }
    if (grepl("difference|change|before|after|polarized|depolarized|same (male|female|animal|subject|individual)", question_lower)) {
      cat("QUESTION CLUE: This may be a within-subject change. If each row links the same individual across the two conditions, create after-minus-before with add_paired_difference(data, before = \"FIRST_CONDITION\", after = \"SECOND_CONDITION\") and test the paired differences. Confirm the row matching first.\n")
    }
    if (grepl("rate|ratio|per (unit|individual|day|hour|area)|proportion|percentage|percent", question_lower)) {
      cat("QUESTION CLUE: The wording may require a rate or proportion. Identify the numerator and denominator from the definition, then use add_ratio_variable(data, numerator = \"NUMERATOR\", denominator = \"DENOMINATOR\", name = \"Rate\", scale = 1). Use scale = 100 only for a percentage. Check the units; zero denominators become NA.\n")
    }
    if (grepl("total|sum of|combined amount", question_lower)) {
      cat("QUESTION CLUE: The wording may require a total. Confirm that the listed components should be added and use the exact columns named by the problem; do not sum unrelated measurements.\n")
    }
    if (grepl("shannon|diversity|species richness|relative abundance", question_lower)) {
      cat("QUESTION CLUE: This may require one Shannon diversity value per independent sample. For long-format data (one row per sample and category), use shannon_by_sample(data, sample_unit = \"SAMPLE_ID\", category = \"SPECIES\", abundance = \"COUNT_OR_ABUNDANCE\", group = \"GROUP\"). Replace names; omit group = ... if there is no comparison group. Never treat the category rows within a sample as independent replicates.\n")
    }
    if (grepl("phenotypic gender|gender index|male.*female.*cone|female.*male.*cone", question_lower)) {
      cat("QUESTION CLUE: If the assignment defines the index as (individual male count / total male count) divided by [(individual female count / total female count) + (individual male count / total male count)], use add_gender_contribution_index(data, male_count = \"MALE_COUNT\", female_count = \"FEMALE_COUNT\"). Check that these columns and this exact formula match the assignment.\n")
    }
    if (grepl("tree height|height.*distance|angle.*height|tangent", question_lower)) {
      cat("QUESTION CLUE: Your lab may define height as distance × tan(angle) + observer height. If that exact formula applies, use add_tree_height(data, distance = \"DISTANCE\", angle = \"ANGLE\", observer_height = \"OBSERVER_HEIGHT\", angle_unit = \"degrees\"). Check the lab formula, angle units, and matching distance units first.\n")
    }
    if (grepl("coefficient of variation|\\bcv\\b|variance among|spread by|variability", question_lower)) {
      cat("QUESTION CLUE: The problem may ask for spread rather than a difference in means. Use dispersion_summary(data, response = \"MEASUREMENT\", group = \"GROUP\") to show variance and CV by group; then follow the assigned equal-variance test (Bartlett assumes normality; Levene is more robust).\n")
    }
    if (grepl("gall|nest|litter|within.*(group|gall|nest)|model ii|heritability", question_lower)) {
      cat("DESIGN CLUE: Multiple measurements from the same gall/nest/litter are clustered, not automatically independent. Identify the shared-unit ID and ask whether the course expects a nested/random-effects variance-components model. Do not call an ICC ‘heritability’ unless the assignment's design justifies that interpretation.\n")
    }
    if (grepl("expected ratio|expected proportion|goodness.of.fit|fit.*(1:1|9:3:3:1|3:1)|poisson", question_lower)) {
      cat("QUESTION CLUE: First total the observed counts into clearly labeled categories, then enter that vector and the matching expected probabilities in goodness_of_fit_workflow(). If a Poisson parameter is estimated from these data, account for that in the degrees of freedom as the course specifies.\n")
    }
  } else {
    cat("If you add question = \"...\" to analysis_start(data), I can highlight wording clues such as product, rate, total, change, or interaction. These are prompts to verify against the problem, not automatic decisions.\n")
  }
  code_columns <- variable_map$column[
    grepl("possible numeric category code", variable_map$suggested_kind, fixed = TRUE)
  ]
  if (length(code_columns) > 0) {
    cat("\nPOSSIBLE NUMERIC CATEGORY CODES\n")
    cat("These columns' names and values could indicate numeric group labels: ",
        paste(code_columns, collapse = ", "), ".\n", sep = "")
    cat("Only if the numbers are labels (for example, 1 = control and 2 = treatment), confirm from your notes, then convert to a category before comparing groups:\n")
    for (code_column in code_columns) {
      cat("  data[[\"", code_column, "\"]] <- factor(data[[\"",
          code_column, "\"]])\n", sep = "")
    }
    cat("If the values are actual counts or measurements, leave them numeric.\n")
  }
  whole_number_columns <- variable_map$column[
    grepl("whole-number count or measurement", variable_map$suggested_kind, fixed = TRUE)
  ]
  if (length(whole_number_columns) > 0) {
    cat("\nWHOLE-NUMBER COLUMNS ARE NOT AUTOMATICALLY CATEGORIES\n")
    cat("These remain numeric for now: ", paste(whole_number_columns, collapse = ", "),
        ". They could be counts or measurements. The values alone do not show that they are group codes.\n", sep = "")
  }
  paired_hint <- .toolkit_paired_wide_hint(data)
  if (!is.null(paired_hint)) {
    before_column <- paired_hint$measurements[1]
    after_column <- paired_hint$measurements[2]
    cat("\nPOSSIBLE PAIRED TWO-CONDITION DATA\n")
    if (length(paired_hint$id) > 0) {
      cat("Possible individual ID column(s): ", paste(paired_hint$id, collapse = ", "),
          ". A sequence of labels may identify subjects; it is usually not a measurement. Confirm this from the study description.\n", sep = "")
    }
    cat("I found two numeric columns, ", before_column, " and ", after_column,
        ". If each row is the same individual measured in both conditions, this is a paired design. The worksheet alone cannot confirm the row matching.\n", sep = "")
    cat("For a within-individual difference column and the paired test, use:\n")
    cat("  data <- add_paired_difference(data, before = \"", before_column,
        "\", after = \"", after_column, "\")\n", sep = "")
    cat("  View(data)  # the extra difference column is optional for the test\n")
    cat("  paired_t_test(data, before = \"", before_column,
        "\", after = \"", after_column, "\")\n", sep = "")
  }
  cat("\nBefore a test: check what one row represents, whether measurements are independent or matched, and whether columns are IDs, counts, ranks, or measurements. Those facts come from the study description, not just the worksheet.\n")

  if (is.null(outcome)) {
    cat("\nWHAT I STILL NEED TO SUGGEST A TEST\n")
    if (is.null(paired_hint)) {
      cat("RESPONSE / OUTCOME: NOT SELECTED YET. Choose the column that represents the result named in your question; do not pick a column just because it is numeric.\n")
      cat("1. Which column is that response? For an association question between two categories, say so—there may be no single response.\n")
    } else {
      cat("POSSIBLE PAIRED COMPARISON: If each row is the same individual in both condition columns, the response for the comparison is the within-individual difference. You do not need to choose one condition column as the sole outcome.\n")
      cat("1. Does each row link the same individual across both columns? If yes, the paired test compares the average difference; if no, tell me how the groups were sampled.\n")
    }
    cat("2. Which column(s) contain the group or possible explanation? Leave out ID columns unless the question is about subject-level differences.\n")
    cat("3. What are you asking: a difference, an association, a prediction, or whether counts fit an expected ratio?\n")
    cat("4. Were observations from different independent units, or were the same/matched units measured more than once?\n")
    if (!is.null(paired_hint)) {
      cat("For the possible paired layout above, answer this specific check: does each row's ID refer to the same individual in both condition columns? If yes, compare the within-individual differences with paired_t_test(); you do not need to invent a single outcome column first.\n")
    }
    if (is.null(paired_hint)) {
      cat("In the example below, replace YOUR_OUTCOME with the result column and YOUR_GROUP with the column used to explain or compare it:\n")
      cat('analysis_start(data, outcome = "YOUR_OUTCOME", predictors = "YOUR_GROUP", goal = "difference", design = "independent")\n')
    } else {
      cat("If you confirm the rows are matched, use the exact column names in this paired test (it calculates after minus before):\n")
      cat("  paired_t_test(data, before = \"", before_column,
          "\", after = \"", after_column, "\")\n", sep = "")
      cat("You can also add the difference column with add_paired_difference() and open View(data); that extra column is optional for the test.\n")
    }
    cat("Use the exact column spelling shown in the table. Do not guess the study design from the data values.\n")
    return(invisible(list(data = data, variable_map = variable_map,
                          outcome = NULL,
                          question = question,
                          suggestions = if (is.null(paired_hint)) character(0) else
                            "If rows are matched by individual, use paired_t_test() on the two condition columns.",
                          missing_inputs = if (is.null(paired_hint)) {
                            c("outcome", "research question/goal", "predictor(s), if any",
                              "independent or paired/repeated design")
                          } else {
                            c("confirm that each row is the same individual in both condition columns",
                              "confirm the comparison goal and hypothesis direction")
                          })))
  }

  .toolkit_check_column(data, outcome, "outcome", "analysis_start()")
  cat("\nRESPONSE VARIABLE (OUTCOME) SELECTED FOR THIS PLAN: ", outcome, "\n", sep = "")
  cat("I will treat this as the measured result (Y). You selected it based on your research question; the worksheet alone cannot establish that role.\n")
  if (is.null(predictors)) {
    predictors <- setdiff(names(data), outcome)
    likely_ids <- variable_map$column[
      grepl("ID", variable_map$suggested_kind, fixed = TRUE) |
        grepl("ID-like", variable_map$plain_English_note, fixed = TRUE)
    ]
    predictors <- setdiff(predictors, likely_ids)
    cat("\nNo predictor list was given, so I will show possibilities from the other non-ID-looking columns. Please remove columns that do not belong to your research question.\n")
  }
  if (length(predictors) > 0) {
    if (!is.character(predictors) || anyNA(predictors) || any(!nzchar(predictors))) {
      stop("analysis_start(): `predictors` should be one or more column names in quotation marks, for example c(\"Treatment\", \"Temperature\").", call. = FALSE)
    }
    for (predictor in predictors) {
      .toolkit_check_column(data, predictor, "predictors", "analysis_start()")
    }
    if (outcome %in% predictors) {
      stop("analysis_start(): The outcome also appears in `predictors`. Remove it from the predictor list; the outcome and predictor need to be different columns.", call. = FALSE)
    }
  }

  outcome_kind <- variable_map$suggested_kind[match(outcome, variable_map$column)]
  y <- data[[outcome]]
  y_is_numeric <- is.numeric(y) && !grepl("possibly an ID", outcome_kind, fixed = TRUE)
  y_is_binary <- length(unique(y[!is.na(y)]) ) == 2 &&
    (is.factor(y) || is.character(y) || is.logical(y) ||
       (is.numeric(y) && all(unique(y[!is.na(y)]) %in% c(0, 1))))
  predictor_kinds <- if (length(predictors) == 0) character(0) else
    variable_map$suggested_kind[match(predictors, variable_map$column)]
  predictor_is_numeric <- if (length(predictors) == 0) logical(0) else
    vapply(predictors, function(name) is.numeric(data[[name]]) &&
             !grepl("possibly an ID", variable_map$suggested_kind[
               match(name, variable_map$column)], fixed = TRUE), logical(1))
  predictor_is_category <- if (length(predictors) == 0) logical(0) else
    vapply(predictors, function(name) is.factor(data[[name]]) ||
             is.character(data[[name]]) || is.logical(data[[name]]) ||
             grepl("possible numeric category code", variable_map$suggested_kind[
               match(name, variable_map$column)], fixed = TRUE), logical(1))

  suggestions <- character(0)
  hypothesis <- NULL
  next_step <- NULL
  unresolved <- character(0)
  if (length(predictors) == 0) {
    if (is.numeric(y) && identical(goal, "goodness_of_fit")) {
      suggestions <- "For counts compared with a proposed probability, use binomial_test_workflow() for two outcomes or goodness_of_fit_workflow() for several categories."
      unresolved <- c(unresolved, "expected probability/ratio", "whether any model parameters were estimated from these same data")
    } else {
      suggestions <- "This is a one-column question only if you have a reference value or expected pattern to compare against."
      unresolved <- c(unresolved, "reference mean or expected count pattern", "research question/goal")
    }
  } else if (y_is_binary) {
    if (all(predictor_is_category)) {
      suggestions <- c("For a question about association between categories, start with contingency_test_workflow().",
                       "For a question about how predictors relate to the probability of one chosen outcome, use logistic_regression_workflow().")
      next_step <- paste0('contingency_test_workflow(data, row = "', outcome,
                          '", column = "', predictors[1], '")')
    } else {
      suggestions <- "A two-outcome result is usually modeled with logistic regression when the question asks whether predictors relate to the chance of an event."
      next_step <- paste0('logistic_regression_workflow(data, ', outcome, ' ~ ',
                          paste(predictors, collapse = " + "), ', event = "CHOOSE_EVENT")')
      unresolved <- c(unresolved, "which outcome label counts as the event")
    }
    hypothesis <- paste0("H0: the outcome's event probability is not related to ",
                         paste(predictors, collapse = " or "),
                         ". HA: it is related. Check whether your question predicts a direction.")
  } else if (y_is_numeric && all(predictor_is_numeric)) {
    if (identical(goal, "prediction")) {
      suggestions <- "Use linear regression when you want an equation to predict the numeric outcome from the predictor(s)."
      next_step <- paste0("linear_regression_workflow(data, ", outcome, " ~ ",
                          paste(predictors, collapse = " + "), ")")
      hypothesis <- paste0("H0: the predictor(s) do not improve prediction of ", outcome,
                           ". HA: at least one predictor helps predict it.")
    } else {
      suggestions <- c("Use correlation to describe the direction and strength of a relationship between two numeric variables.",
                       "Use linear regression if your question asks how the outcome changes with a predictor, or asks for prediction.")
      if (length(predictors) == 1) {
        next_step <- paste0('correlation_workflow(data, x = "', outcome,
                            '", y = "', predictors[1], '")')
      }
      hypothesis <- paste0("H0: there is no relationship between ", outcome, " and ",
                           paste(predictors, collapse = ", "),
                           ". HA: there is a relationship. Correlation by itself does not show cause.")
    }
    if (length(predictors) > 1 && identical(goal, "association")) {
      unresolved <- c(unresolved, "whether predictors should be analyzed separately (as some course assignments require) or together")
    }
  } else if (y_is_numeric && any(predictor_is_numeric) &&
             any(predictor_is_category)) {
    suggestions <- "Your numeric outcome is paired with both measurement and category predictors. A linear model can include both, but confirm the assignment asks you to adjust for them together and make sure category columns are treated as factors."
    next_step <- paste0("linear_regression_workflow(data, ", outcome, " ~ ",
                        paste(predictors, collapse = " + "), ")")
    hypothesis <- paste0("H0: the predictor(s), considered in this model, do not explain variation in ",
                         outcome, ". HA: at least one predictor is related to it.")
  } else if (y_is_numeric && all(predictor_is_category)) {
    group_counts <- vapply(predictors, function(name) {
      length(unique(data[[name]][!is.na(data[[name]])]))
    }, integer(1))
    if (length(predictors) == 1 && group_counts[1] == 2) {
      if (!is.null(design) && design %in% c("paired", "repeated")) {
        suggestions <- "For two matched measurements stored in separate columns, use paired_t_test(). If your data are in long format (one outcome column plus a condition column), ask before reshaping or use an appropriate repeated-measures model."
        next_step <- "paired_t_test(data, before = \"FIRST_CONDITION_COLUMN\", after = \"SECOND_CONDITION_COLUMN\")"
      } else if (identical(design, "independent")) {
        suggestions <- "For two separate groups, use independent_t_test(). It defaults to Welch's t-test, which does not assume equal group variances."
        next_step <- paste0('independent_t_test(data, outcome = "', outcome,
                            '", group = "', predictors[1], '")')
      } else {
        suggestions <- c("If the groups contain different, unrelated study units, use independent_t_test().",
                         "If the same or deliberately matched units appear in both groups, use paired_t_test().")
        unresolved <- c(unresolved, "independent groups or matched/paired measurements")
      }
      hypothesis <- paste0("H0: the mean ", outcome, " is the same across the two conditions. HA: the means differ (unless your prior prediction named a direction).")
    } else if (length(predictors) == 1 && group_counts[1] >= 3) {
      if (!is.null(design) && design %in% c("paired", "repeated")) {
        suggestions <- "For the same subjects measured at three or more conditions, use repeated_measures_anova_workflow() with the subject ID column."
        unresolved <- c(unresolved, "subject ID column")
        next_step <- paste0('repeated_measures_anova_workflow(data, outcome = "',
                            outcome, '", subject = "SUBJECT_ID", within = "',
                            predictors[1], '")')
      } else {
        suggestions <- c("For three or more independent groups, use anova_workflow() to compare means.",
                         "If the data remain strongly non-normal or are ranks, ask whether Kruskal-Wallis is appropriate; it compares rank patterns, not simply means.")
        if (identical(design, "independent")) {
          next_step <- paste0("anova_workflow(data, ", outcome, " ~ ", predictors[1], ")")
        } else {
          unresolved <- c(unresolved, "independent groups or repeated measurements")
        }
      }
      hypothesis <- paste0("H0: the mean ", outcome, " is equal across all groups. HA: at least one group mean differs.")
    } else {
      suggestions <- "For a numeric outcome and categorical predictors, ANOVA/regression may fit, but this data alone cannot tell whether factors are fixed, random, or repeated. Check the sampling/design description first."
      unresolved <- c(unresolved, "which factor levels were selected as treatments versus sampled from a larger population", "whether observations are independent or repeated")
    }
  } else if (!y_is_numeric && all(predictor_is_category)) {
    suggestions <- "For two categorical variables, use contingency_test_workflow() to test for an association. Its output checks expected counts and suggests Fisher's exact test or simulation when counts are sparse."
    next_step <- paste0('contingency_test_workflow(data, row = "', outcome,
                        '", column = "', predictors[1], '")')
    hypothesis <- paste0("H0: ", outcome, " and ", predictors[1],
                         " are not associated. HA: they are associated.")
  } else {
    suggestions <- "This combination does not map cleanly to the toolkit's introductory tests. Check if a numeric column was imported as text, whether a numeric group code should be a factor, and whether the outcome has more than two categories."
    unresolved <- c(unresolved, "meaning of each selected variable", "whether number-like text or numeric codes are measurements or labels")
  }

  cat("\nQUESTION-BASED STARTING POINT\n")
  cat("Outcome:", outcome, "(", outcome_kind, ")\n", sep = " ")
  cat("Predictors considered:", if (length(predictors)) paste(predictors, collapse = ", ") else "none supplied", "\n")
  if (is.null(goal)) {
    cat("Goal not supplied: decide whether you want a difference, association, prediction, or fit to an expected count pattern. The same columns can answer different questions.\n")
    unresolved <- c(unresolved, "research question/goal")
  }
  if (length(predictors) > 0 && (is.null(design) || identical(design, "unknown"))) {
    unresolved <- c(unresolved, "whether rows are independent or measurements are paired/repeated (check the study description)")
  }
  if (length(suggestions)) {
    cat("\nPossible route(s):\n")
    for (suggestion in suggestions) cat("- ", suggestion, "\n", sep = "")
  }
  if (!is.null(hypothesis)) {
    cat("\nPlain-English hypothesis template (check it against the study question):\n")
    cat("  ", hypothesis, "\n", sep = "")
    cat("This template is two-sided (different/associated). Use a one-direction hypothesis only if your research prediction named that direction before you examined the data.\n")
  }
  if (!is.null(next_step)) {
    cat("\nLikely next call:\n  ", next_step, "\n", sep = "")
  }
  cat("\nDesign reminder: a test cannot fix dependent observations, pseudoreplication, confounding, or biased sampling. Check that one row represents the correct independent unit. If the factor levels are a sample from a broader population rather than the specific treatments of interest, ask whether a random-effects model is intended.\n")
  if (length(predictors) > 1) {
    cat("Your course materials sometimes analyze predictors in separate models and sometimes together. Follow the assignment wording: for separate models, fit one predictor at a time; use ` + ` in one formula only when the question asks for a combined model.\n")
  }
  unresolved <- unique(unresolved)
  if (length(unresolved) > 0) {
    cat("\nImportant detail(s) still to confirm:\n")
    for (item in unresolved) cat("- ", item, "\n", sep = "")
    cat("If the study description does not answer one of these, ask your instructor before choosing the test.\n")
  } else {
    cat("\nThe column types suggest a starting route. Before running it, confirm the outcome, the scientific question, and the study design in the assignment or field/lab notes.\n")
  }
  invisible(list(data = data, variable_map = variable_map, outcome = outcome,
                 question = question,
                 predictors = predictors, goal = goal, design = design,
                 suggestions = suggestions, hypothesis_template = hypothesis,
                 next_step = next_step,
                 missing_inputs = unresolved))
}

# Suggest a starting analysis based on the outcome and predictor column types.
# This is a guide, not a substitute for the research question or study design.
analysis_guide <- function(data = NULL, outcome = NULL, predictors = NULL) {
  .toolkit_require_inputs("analysis_guide()", list(data = data, outcome = outcome),
                          c(data = "the imported worksheet", outcome = "the column you want to explain or compare"))
  data <- .toolkit_prepare_data(data, "analysis_guide()")
  .toolkit_check_column(data, outcome, "outcome", "analysis_guide()")

  if (is.null(predictors)) {
    predictors <- setdiff(names(data), outcome)
    message("No predictors were listed, so I considered every other column. If the dataset has IDs or unrelated columns, pass only the predictors that match your question.")
  }
  if (!is.character(predictors) || length(predictors) == 0 || anyNA(predictors)) {
    stop("analysis_guide(): `predictors` should be one or more column names in quotation marks. If left blank, all other columns are considered.", call. = FALSE)
  }
  for (predictor in predictors) {
    .toolkit_check_column(data, predictor, "predictors", "analysis_guide()")
  }

  y <- data[[outcome]]
  y_values <- unique(na.omit(y))
  categorical_outcome <- is.factor(y) || is.character(y) || is.logical(y)
  is_binary <- (is.factor(y) || is.character(y) || is.logical(y)) &&
    length(y_values) == 2
  if (is.numeric(y) && length(y_values) == 2 &&
      all(y_values %in% c(0, 1))) is_binary <- TRUE

  predictor_types <- vapply(predictors, function(name) {
    value <- data[[name]]
    if (is.numeric(value)) {
      "numeric measurement"
    } else if (is.factor(value) || is.character(value) || is.logical(value)) {
      "category/label"
    } else {
      paste("other:", paste(class(value), collapse = "/"))
    }
  }, character(1))

  cat("\n====================================\n")
  cat("WHICH ANALYSIS SHOULD I START WITH?\n")
  cat("====================================\n")
  cat("Outcome (the result you want to understand):", outcome, "\n")
  cat("Predictor(s) (possible explanations):", paste(predictors, collapse = ", "), "\n")
  cat("Predictor types:", paste(names(predictor_types), predictor_types,
                                  sep = " = ", collapse = "; "), "\n\n")

  if (is_binary) {
    if (all(predictor_types == "category/label")) {
      cat("Both the outcome and predictors are categories. Use a contingency test for a question about whether two categories are associated; use logistic regression when you want to model the probability of a chosen outcome category.\n")
      cat("Next function for association: contingency_test_workflow(data, row = \"", outcome,
          "\", column = \"", predictors[1], "\").\n", sep = "")
      cat("Next function for modeling event probability: logistic_regression_workflow(data, ", outcome,
          " ~ ", paste(predictors, collapse = " + "), ", event = \"choose the event label\").\n", sep = "")
    } else {
      cat("The outcome has two possible values. Use logistic regression if you want to model how numeric or categorical predictors relate to the chance of the event.\n")
      cat("Next function: logistic_regression_workflow(data, outcome ~ predictor). Specify which outcome value counts as the event when you run it.\n")
    }
  } else if (is.numeric(y) && all(predictor_types == "numeric measurement")) {
    if (length(predictors) == 1) {
      cat("Both columns are numeric. Use correlation for the strength/direction of a two-variable association; it does not fit a regression model, so there are no regression residuals to check. Use linear regression when you want a fitted equation, prediction, or the outcome's change per unit of the predictor; that function fits the model and automatically shows residual checks. These answer related but different questions.\n")
      cat("Next functions: correlation_workflow(data, \"", outcome, "\", \"", predictors[1], "\") or linear_regression_workflow(data, ", outcome, " ~ ", predictors[1], ").\n", sep = "")
    } else {
      cat("The outcome and predictors are numeric. Use multiple linear regression to examine the outcome in relation to several predictors together. The regression function fits the model and automatically shows residual checks. Correlation only describes two variables at a time and does not adjust for the others.\n")
      cat("Next function: linear_regression_workflow(data, ", outcome, " ~ ",
          paste(predictors, collapse = " + "), ").\n", sep = "")
    }
  } else if (is.numeric(y) && all(predictor_types == "category/label")) {
    if (length(predictors) == 1 && length(unique(na.omit(data[[predictors[1]]]))) == 2) {
      cat("The outcome is numeric and the predictor has two groups. If the groups are different, unrelated subjects, use independent_t_test(). If the same subjects or matched pairs are measured twice, use paired_t_test() with one column per measurement.\n")
      cat("For independent groups: independent_t_test(data, outcome = \"", outcome,
          "\", group = \"", predictors[1], "\").\n", sep = "")
      cat("If the predictor has 3+ independent groups, use anova_workflow().\n")
    } else {
      cat("The outcome is numeric and the predictors are categories. Use ANOVA to compare the outcome across groups.\n")
      cat("Next function: anova_workflow(data, ", outcome, " ~ ",
          paste(predictors, collapse = " + "), ").\n", sep = "")
    }
  } else if (is.numeric(y) &&
             all(predictor_types %in% c("numeric measurement", "category/label"))) {
    cat("The outcome is numeric and your predictors mix measurements and categories. Use linear regression to include them together; the model estimates numeric slopes and category differences in one analysis, then automatically shows residual checks.\n")
    cat("Next function: linear_regression_workflow(data, ", outcome, " ~ ",
        paste(predictors, collapse = " + "), ").\n", sep = "")
  } else if (categorical_outcome && all(predictor_types == "category/label")) {
    cat("The outcome and predictors are categories. Use a chi-square or G-test of independence to test whether they are associated. If expected counts are sparse, consider Fisher's exact test or simulation.\n")
    cat("Next function: contingency_test_workflow(data, row = \"", outcome,
        "\", column = \"", predictors[1], "\").\n", sep = "")
  } else {
    cat("I cannot match this outcome to one of the toolkit's basic workflows. Numeric outcomes can use ANOVA or linear regression; two-category outcomes can use logistic regression. For an outcome with 3+ categories, a different model may be needed.\n")
  }

  cat("\nBefore choosing: identify what counts as the outcome, whether the question is about association or prediction, and whether observations are independent or paired/repeated. A data-type match alone cannot determine the study design.\n")
  cat("If your assignment asks you to analyze predictors separately, use one predictor per model formula. Combine predictors only when the question and course method call for considering them together.\n")
  cat("A numeric column containing group codes (such as 1 and 2) may need to be changed to a factor before it is treated as a category. Use factor(data[[\"column_name\"]]) only when those numbers are labels, not measurements.\n")
  invisible(list(outcome = outcome, predictors = predictors,
                 outcome_class = class(y), predictor_types = predictor_types,
                 binary_outcome = is_binary))
}


correlation_workflow <- function(data = NULL, x = NULL, y = NULL,
                                 method = "pearson",
                                 alpha = 0.05,
                                 plots = TRUE) {
  .toolkit_require_inputs("correlation_workflow()",
                          list(data = data, x = x, y = y),
                          c(data = "the imported worksheet", x = "the first numeric column", y = "the second numeric column"))
  data <- .toolkit_prepare_data(data, "correlation_workflow()")
  .toolkit_check_column(data, x, "x", "correlation_workflow()")
  .toolkit_check_column(data, y, "y", "correlation_workflow()")
  .toolkit_check_alpha(alpha, "correlation_workflow()")
  if (!is.character(method) || length(method) != 1 || is.na(method) ||
      !method %in% c("pearson", "spearman", "kendall")) {
    stop("correlation_workflow(): Choose `method = \"pearson\"`, `\"spearman\"`, or `\"kendall\"`. Pearson checks a straight-line pattern; Spearman and Kendall check rank-based increasing/decreasing patterns.", call. = FALSE)
  }
  if (!is.logical(plots) || length(plots) != 1 || is.na(plots)) {
    stop("correlation_workflow(): `plots` should be TRUE to draw the scatterplot or FALSE to skip it.", call. = FALSE)
  }
  if (!is.numeric(data[[x]]) || !is.numeric(data[[y]])) {
    stop("correlation_workflow(): Both selected columns must be numeric. I received ", x, " as ",
         paste(class(data[[x]]), collapse = "/"), " and ", y, " as ",
         paste(class(data[[y]]), collapse = "/"), ". Choose numeric measurements. If values are numbers stored as text, inspect them before converting.", call. = FALSE)
  }

  complete <- is.finite(data[[x]]) & is.finite(data[[y]])
  x_values <- data[[x]][complete]
  y_values <- data[[y]][complete]
  n <- length(x_values)
  if (n < 3) {
    stop("correlation_workflow(): I found fewer than 3 rows with usable values in both columns. Check missing values and make sure you selected the right columns.", call. = FALSE)
  }
  if (sd(x_values) == 0 || sd(y_values) == 0) {
    stop("correlation_workflow(): One selected column has the same value in every usable row. Correlation needs variation in both columns.", call. = FALSE)
  }

  test <- tryCatch(
    cor.test(x_values, y_values, method = method, exact = FALSE),
    error = function(e) stop("correlation_workflow(): R could not calculate this correlation. Check that both columns contain enough varied numeric values. R's message was: ", conditionMessage(e), call. = FALSE)
  )
  estimate <- unname(test$estimate)

  cat("\n====================================\n")
  cat("CORRELATION WORKFLOW\n")
  cat("====================================\n")
  cat("Variables:", x, "and", y, "| Method:", method, "| Complete pairs:", n, "\n")
  cat("This is a correlation analysis, not a fitted regression model. Correlation does not require residual checks. The optional red line on the Pearson plot is only a visual straight-line guide.\n")
  print(test)
  cat("\nHow to read this: the correlation estimate ranges from -1 to +1. The sign gives direction; values closer to either end indicate a stronger pattern of association. The p-value tests for evidence of association under this method.\n")
  cat("Compare the p-value with alpha =", alpha,
      ". A smaller p-value is evidence of association under the test assumptions; it does not measure how important the association is.\n")
  if (test$p.value < alpha) {
    cat("At this alpha level, the test gives evidence of an association. Check its direction, size, and scatterplot before describing what it means.\n")
  } else {
    cat("At this alpha level, the test does not give strong evidence of an association. This does not prove the variables are unrelated.\n")
  }
  if (method == "pearson") {
    cat("Pearson measures a straight-line relationship and can be strongly affected by outliers. Check the scatterplot before interpreting it.\n")
  } else {
    cat("This rank-based method looks for a consistently increasing or decreasing pattern; it does not require the pattern to be a straight line. Check the scatterplot for the shape and outliers.\n")
  }
  cat("Correlation does not show that one variable caused the other. If your question is prediction, or you need to account for other predictors, consider linear_regression_workflow().\n")

  if (plots) {
    plot(x_values, y_values, xlab = x, ylab = y,
         main = paste("", method, "correlation:", x, "and", y), pch = 19)
    if (method == "pearson") {
      abline(lm(y_values ~ x_values), col = "red", lwd = 2)
      legend("topleft", legend = "straight-line guide only", col = "red",
             lty = 1, bty = "n", cex = 0.8)
    }
  }

  invisible(list(test = test, estimate = estimate, n = n,
                 method = method, x = x, y = y))
}


linear_regression_workflow <- function(data = NULL, formula = NULL,
                                       alpha = 0.05,
                                       plots = TRUE) {
  .toolkit_require_inputs("linear_regression_workflow()",
                          list(data = data, formula = formula),
                          c(data = "the imported worksheet", formula = "a formula with the numeric outcome on the left, such as Height ~ Treatment"))
  data <- .toolkit_prepare_data(data, "linear_regression_workflow()")
  .toolkit_check_alpha(alpha, "linear_regression_workflow()")
  if (is.character(formula) && length(formula) == 1 && !is.na(formula)) {
    formula <- tryCatch(as.formula(formula), error = function(e) NULL)
    if (!is.null(formula)) message("I converted the text you supplied into an R formula.")
  }
  if (!inherits(formula, "formula")) {
    stop("linear_regression_workflow(): Give a formula that names the numeric outcome first and predictor(s) after ~, such as RatDnsty ~ ShrbCvr + SeedPrd.", call. = FALSE)
  }
  if (!is.logical(plots) || length(plots) != 1 || is.na(plots)) {
    stop("linear_regression_workflow(): `plots` should be TRUE to show diagnostic plots or FALSE to skip them.", call. = FALSE)
  }
  variables <- all.vars(formula)
  missing_variables <- setdiff(variables, names(data))
  if (length(missing_variables) > 0) {
    stop("linear_regression_workflow(): I cannot find these formula column(s): ",
         paste(missing_variables, collapse = ", "), ". Available columns: ",
         paste(names(data), collapse = ", "), ". Check spelling and capitalization.", call. = FALSE)
  }
  if (length(variables) < 2) {
    stop("linear_regression_workflow(): Include at least one predictor after ~, such as outcome ~ predictor.", call. = FALSE)
  }
  outcome <- variables[1]
  if (!is.numeric(data[[outcome]])) {
    stop("linear_regression_workflow(): The outcome column \"", outcome,
         "\" must be numeric. If the outcome has only two categories, use logistic_regression_workflow() instead.", call. = FALSE)
  }

  model <- tryCatch(
    lm(formula, data = data, na.action = na.omit),
    error = function(e) stop("linear_regression_workflow(): R could not fit this model. Check that the formula columns are correct and that the data have enough usable observations. R's message was: ", conditionMessage(e), call. = FALSE)
  )
  model_summary <- summary(model)
  model_frame <- model.frame(model)
  residual_values <- residuals(model)
  fitted_values <- fitted(model)

  cat("\n====================================\n")
  cat("LINEAR REGRESSION WORKFLOW\n")
  cat("====================================\n")
  cat("I fitted a linear model automatically with lm(). The left side of the formula is the response; the right side contains predictor(s).\n")
  cat("Model:", deparse(formula), "\n")
  cat("Complete observations used:", nrow(model_frame), "\n\n")
  print(model_summary)
  cat("\nHow to read this: each numeric coefficient is the predicted change in the outcome for a one-unit increase in that predictor, holding the other predictors in this model constant. A factor coefficient compares that category with R's reference category.\n")
  cat("R-squared is the fraction of outcome variation described by this model in these data; it is not prediction accuracy on new data. The overall F-test asks whether the model explains more variation than an intercept-only model.\n")
  cat("For a coefficient's p-value, compare it with alpha =", alpha,
      ". A small p-value is evidence of a relationship in this model; also consider the coefficient size, confidence interval, and study design.\n")
  cat("I am checking this model's residuals below because regression assumptions apply to the fitted model, not to the raw columns by themselves. Residuals are observed values minus model-predicted values. Check that the residuals have no curve or funnel, the Q-Q plot is roughly straight, and there are no extreme influential points. Independence comes from how the study was collected; a plot cannot prove it.\n")

  shapiro_result <- NULL
  if (length(residual_values) >= 3 && length(residual_values) <= 5000) {
    shapiro_result <- tryCatch(shapiro.test(residual_values), error = function(e) NULL)
    if (is.null(shapiro_result)) {
      cat("\nShapiro-Wilk residual check could not be calculated (for example, residuals may all be identical). Use the Q-Q plot and residual spread plot instead.\n")
    } else {
      cat("\nShapiro-Wilk check of model residuals (a clue, not a pass/fail verdict):\n")
      print(shapiro_result)
      if (shapiro_result$p.value < alpha) {
        cat("This gives evidence that residuals depart from a normal shape. Look at the Q-Q plot and sample size; this alone does not tell you which fix is right.\n")
      } else {
        cat("This does not find strong evidence against a normal residual shape; it does not prove normality. Also inspect the Q-Q plot.\n")
      }
    }
  } else {
    cat("\nThe Shapiro-Wilk residual check was skipped: it requires 3 to 5,000 residuals. Use the Q-Q plot as a visual check.\n")
  }

  variance_result <- NULL
  if (requireNamespace("car", quietly = TRUE)) {
    variance_result <- tryCatch(car::ncvTest(model), error = function(e) NULL)
    if (!is.null(variance_result)) {
      cat("\nNon-constant variance check (car package):\n")
      print(variance_result)
      if (!is.null(variance_result$p) && variance_result$p < alpha) {
        cat("This gives evidence that residual spread changes with fitted values. Check the residual plot for a funnel and discuss a suitable model or variance approach with your instructor.\n")
      } else {
        cat("This check did not find strong evidence that residual spread changes. Still inspect the residual plot; this is not proof of equal spread.\n")
      }
    }
  } else {
    cat("The optional car package is unavailable, so the formal changing-spread check was skipped; the residual plot is still shown.\n")
  }

  if (plots) {
    old_par <- par(mfrow = c(2, 2))
    on.exit(par(old_par), add = TRUE)
    formula_predictors <- all.vars(delete.response(terms(model)))
    single_numeric_predictor <- length(formula_predictors) == 1 &&
      is.numeric(model_frame[[formula_predictors[1]]])
    if (single_numeric_predictor) {
      predictor_values <- model_frame[[formula_predictors[1]]]
      plot(predictor_values, model.response(model_frame),
           xlab = formula_predictors[1], ylab = outcome,
           main = "Observed data + fitted line", pch = 19)
      abline(model, col = "red", lwd = 2)
      if (length(unique(predictor_values)) >= 4) {
        lines(lowess(predictor_values, model.response(model_frame)),
              col = "blue", lwd = 2)
        legend("topleft", legend = c("linear model", "curved LOWESS guide"),
               col = c("red", "blue"), lty = 1, bty = "n", cex = 0.75)
        cat("Plot guide: the first panel shows the observed points, fitted straight line (red), and a flexible curved LOWESS guide (blue). A clear curved mismatch can suggest the relationship is not straight-line; LOWESS is a visual clue, not a replacement test.\n")
      } else {
        cat("Plot guide: the first panel shows the observed points and fitted straight line. There are too few distinct predictor values for a useful curved guide.\n")
      }
    } else {
      plot(fitted_values, model.response(model_frame),
           xlab = "Model-predicted values", ylab = "Observed outcome",
           main = "Observed vs. predicted", pch = 19)
      abline(0, 1, col = "red", lwd = 2)
      cat("Plot guide: with multiple or nonnumeric predictors there is no single X-versus-Y view; this panel compares observed outcomes with model predictions. The red line marks perfect prediction.\n")
    }
    plot(fitted_values, residual_values,
         xlab = "Model-predicted values", ylab = "Residuals",
         main = "Residuals vs. predicted", pch = 19)
    abline(h = 0, col = "red", lwd = 2)
    if (length(unique(fitted_values)) >= 4) {
      lines(lowess(fitted_values, residual_values), col = "blue", lwd = 2)
    }
    cat("Residual plot guide: look for a random band around zero. A curve suggests the straight-line model may miss a pattern; a funnel suggests changing residual spread.\n")
    qqnorm(residual_values, main = "Q-Q plot of residuals", pch = 19)
    qqline(residual_values, col = "red", lwd = 2)
    plot(fitted_values, sqrt(abs(rstandard(model))),
         xlab = "Model-predicted values", ylab = "Sqrt |standardized residual|",
         main = "Residual spread", pch = 19)
    if (length(unique(fitted_values)) >= 4) {
      lines(lowess(fitted_values, sqrt(abs(rstandard(model)))), col = "blue", lwd = 2)
    }
  }

  invisible(list(model = model, summary = model_summary,
                 coefficients = model_summary$coefficients,
                 observations = nrow(model_frame), residuals = residual_values,
                 fitted = fitted_values, shapiro = shapiro_result,
                 nonconstant_variance = variance_result))
}


logistic_regression_workflow <- function(data = NULL, formula = NULL,
                                         event = NULL,
                                         alpha = 0.05) {
  .toolkit_require_inputs("logistic_regression_workflow()",
                          list(data = data, formula = formula),
                          c(data = "the imported worksheet", formula = "a formula with a two-category outcome on the left"))
  data <- .toolkit_prepare_data(data, "logistic_regression_workflow()")
  .toolkit_check_alpha(alpha, "logistic_regression_workflow()")
  if (is.character(formula) && length(formula) == 1 && !is.na(formula)) {
    formula <- tryCatch(as.formula(formula), error = function(e) NULL)
    if (!is.null(formula)) message("I converted the text you supplied into an R formula.")
  }
  if (!inherits(formula, "formula") || !is.name(formula[[2]])) {
    stop("logistic_regression_workflow(): Give a formula with a simple two-category outcome column on the left, such as Survival ~ BodySize. Do not put a calculation on the left side.", call. = FALSE)
  }
  variables <- all.vars(formula)
  missing_variables <- setdiff(variables, names(data))
  if (length(missing_variables) > 0) {
    stop("logistic_regression_workflow(): I cannot find these formula column(s): ",
         paste(missing_variables, collapse = ", "), ". Available columns: ",
         paste(names(data), collapse = ", "), ". Check spelling and capitalization.", call. = FALSE)
  }
  if (length(variables) < 2) {
    stop("logistic_regression_workflow(): Include at least one predictor after ~, such as event ~ predictor.", call. = FALSE)
  }

  outcome <- as.character(formula[[2]])
  original_outcome <- data[[outcome]]
  observed_values <- unique(na.omit(original_outcome))
  is_01 <- is.numeric(original_outcome) && length(observed_values) == 2 &&
    all(observed_values %in% c(0, 1))

  if (is_01) {
    event_value <- if (is.null(event)) "1" else as.character(event)
    if (length(event_value) != 1 || is.na(event_value) ||
        !event_value %in% as.character(observed_values)) {
      stop("logistic_regression_workflow(): For a numeric 0/1 outcome, `event` must be either 0 or 1. By default, 1 is treated as the event.", call. = FALSE)
    }
    binary_outcome <- as.integer(as.character(original_outcome) == event_value)
  } else if (is.factor(original_outcome) || is.character(original_outcome) ||
             is.logical(original_outcome)) {
    outcome_levels <- if (is.factor(original_outcome)) {
      levels(droplevels(original_outcome))
    } else {
      sort(unique(as.character(na.omit(original_outcome))))
    }
    if (length(outcome_levels) != 2) {
      stop("logistic_regression_workflow(): The outcome must have exactly two observed categories. I found ",
           length(outcome_levels), ". For 3+ categories, use a multinomial model; this toolkit's logistic function handles two-category outcomes only.", call. = FALSE)
    }
    event_value <- if (is.null(event)) outcome_levels[2] else as.character(event)
    if (length(event_value) != 1 || is.na(event_value) ||
        !event_value %in% outcome_levels) {
      stop("logistic_regression_workflow(): `event` must exactly match one outcome label: ",
           paste(outcome_levels, collapse = " or "), ". Check spelling and capitalization.", call. = FALSE)
    }
    if (is.null(event)) {
      message("No `event` was specified. I treated \"", event_value,
              "\" (the second outcome level) as the event. Set event = \"",
              event_value, "\" explicitly to make this choice clear.")
    }
    binary_outcome <- as.integer(as.character(original_outcome) == event_value)
  } else {
    stop("logistic_regression_workflow(): The outcome must be a factor, text, or TRUE/FALSE column with exactly two categories, or numeric values coded 0 and 1. I received class ",
         paste(class(original_outcome), collapse = "/"), ". Check the outcome column or recode it intentionally.", call. = FALSE)
  }

  analysis_data <- data
  analysis_data[[outcome]] <- binary_outcome
  model <- tryCatch(
    glm(formula, data = analysis_data, family = binomial(), na.action = na.omit),
    error = function(e) stop("logistic_regression_workflow(): R could not fit this model. Check that each predictor has usable values and that the data contain both outcome categories. R's message was: ", conditionMessage(e), call. = FALSE)
  )
  model_summary <- summary(model)
  coefficient_table <- model_summary$coefficients
  odds_ratio <- exp(coefficient_table[, "Estimate"])
  lower_95 <- exp(coefficient_table[, "Estimate"] - 1.96 * coefficient_table[, "Std. Error"])
  upper_95 <- exp(coefficient_table[, "Estimate"] + 1.96 * coefficient_table[, "Std. Error"])
  odds_table <- data.frame(
    term = rownames(coefficient_table),
    log_odds_coefficient = coefficient_table[, "Estimate"],
    odds_ratio = odds_ratio,
    lower_95_Wald = lower_95,
    upper_95_Wald = upper_95,
    p_value = coefficient_table[, "Pr(>|z|)"],
    row.names = NULL,
    check.names = FALSE
  )
  pseudo_r_squared <- if (model$null.deviance > 0) {
    1 - model$deviance / model$null.deviance
  } else {
    NA_real_
  }
  likelihood_ratio <- model$null.deviance - model$deviance
  likelihood_df <- model$df.null - model$df.residual
  overall_p <- if (likelihood_df > 0) {
    pchisq(likelihood_ratio, df = likelihood_df, lower.tail = FALSE)
  } else {
    NA_real_
  }

  cat("\n====================================\n")
  cat("LOGISTIC REGRESSION WORKFLOW\n")
  cat("====================================\n")
  cat("Model:", deparse(formula), "\n")
  cat("Event treated as 1:", event_value, "| Complete observations:", nobs(model), "\n")
  cat("\nOdds ratios (with approximate 95% Wald intervals):\n")
  print(odds_table, row.names = FALSE)
  cat("\nOverall likelihood-ratio test p-value:", format.pval(overall_p), "\n")
  cat("McFadden pseudo-R-squared:", round(pseudo_r_squared, 4), "\n")
  cat("\nHow to read this: an odds ratio above 1 means higher predictor values/categories are associated with higher odds of the selected event; below 1 means lower odds. For numeric predictors, this is per one-unit increase. For factors, it compares with the reference category.\n")
  cat("Compare a term's p-value with alpha =", alpha,
      "and check whether its 95% interval includes 1. These are related checks; the odds-ratio size and uncertainty matter too.\n")
  cat("The intervals and p-values describe uncertainty in this fitted model. They do not show causation. Check that observations are independent, that numeric predictors relate roughly linearly to log-odds, and that the model converged without separation or extreme estimates.\n")
  cat("For a prediction goal, evaluate predictions on new or held-out data and choose a decision threshold for the actual purpose; the default 0.5 cutoff is not automatically best.\n")
  if (!model$converged) {
    cat("\nWARNING: The model did not converge. Treat these estimates as unreliable; check sparse groups, too few events, and predictors that perfectly separate the outcome.\n")
  }
  if (any(!is.finite(coefficient_table[, "Estimate"]))) {
    cat("\nWARNING: At least one coefficient is not finite. This often happens when a predictor perfectly separates the two outcomes or some categories have too few observations.\n")
  }

  invisible(list(model = model, summary = model_summary,
                 odds_ratios = odds_table,
                 event = event_value,
                 pseudo_r_squared = pseudo_r_squared,
                 overall_p_value = overall_p,
                 observations = nobs(model)))
}

# ---------------------------------------------------------------
# 9. ANOVA WORKFLOW
# ---------------------------------------------------------------

anova_workflow <- function(data = NULL,
                           formula = NULL,
                           alpha = 0.05,
                           plots = TRUE,
                           run_tukey = TRUE,
                           enhanced = TRUE) {
  .toolkit_require_inputs("anova_workflow()",
                          list(data = data, formula = formula),
                          c(data = "the imported worksheet", formula = "a formula with the numeric outcome on the left and group(s) on the right"))
  data <- .toolkit_prepare_data(data, "anova_workflow()")
  .toolkit_check_alpha(alpha, "anova_workflow()")
  options <- list(plots = plots, run_tukey = run_tukey, enhanced = enhanced)
  for (option_name in names(options)) {
    option_value <- options[[option_name]]
    if (!is.logical(option_value) || length(option_value) != 1 || is.na(option_value)) {
      stop("anova_workflow(): `", option_name,
           "` must be TRUE or FALSE. Use TRUE to include that output or FALSE to skip it.", call. = FALSE)
    }
  }
  if (is.character(formula) && length(formula) == 1 && !is.na(formula)) {
    formula <- tryCatch(as.formula(formula), error = function(e) NULL)
    if (!is.null(formula)) message("I converted the text you supplied into an R formula.")
  }
  if (!inherits(formula, "formula")) {
    stop("anova_workflow(): `formula` describes which columns to compare. Give it as y ~ group, for example RatDnsty ~ Treatment.\n",
         "What I received was not an R formula. You can also supply the formula as text in quotation marks, such as \"RatDnsty ~ Treatment\".",
         call. = FALSE)
  }

  formula_variables <- all.vars(formula)
  missing_variables <- setdiff(formula_variables, names(data))
  if (length(missing_variables) > 0) {
    stop("anova_workflow(): The formula refers to column(s) I cannot find: ",
         paste(missing_variables, collapse = ", "), ".\n",
         "Columns in your dataset: ", paste(names(data), collapse = ", "), ".\n",
         "Check spelling and capitalization in the formula, or use the correct dataset.", call. = FALSE)
  }
  if (length(formula_variables) < 2) {
    stop("anova_workflow(): The formula has no predictor/group to compare. Include a response and at least one predictor, like response ~ group.", call. = FALSE)
  }
  response_name <- formula_variables[1]
  if (!is.numeric(data[[response_name]])) {
    stop("anova_workflow(): The response column \"", response_name, "\" is stored as ",
         paste(class(data[[response_name]]), collapse = "/"), ", but this workflow expects a numeric measurement.\n",
         "Choose a numeric outcome. If its values are numbers stored as text, inspect them before converting.", call. = FALSE)
  }

  for (predictor_name in setdiff(formula_variables, response_name)) {
    if (is.character(data[[predictor_name]])) {
      data[[predictor_name]] <- factor(data[[predictor_name]])
      message("anova_workflow(): I treated the text column \"", predictor_name,
              "\" as a group/category variable by converting it to a factor.")
    }
  }

  model <- tryCatch(
    lm(formula, data = data, na.action = na.omit),
    error = function(e) {
      stop("anova_workflow(): R could not fit this model.\n",
           "R's message was: ", conditionMessage(e), "\n\n",
           "Common things to check: every formula name matches a dataset column; the response is numeric; and each group you compare has enough usable observations and variation.",
           call. = FALSE)
    }
  )
  model_frame <- model.frame(model)
  predictor_names <- all.vars(delete.response(terms(model)))
  response_name <- all.vars(formula)[1]
  residual_values <- residuals(model)
  fitted_values <- fitted(model)

  numeric_predictors <- predictor_names[
    vapply(model_frame[predictor_names], is.numeric, logical(1))
  ]
  if (length(numeric_predictors) > 0) {
    message(
      "I found numeric predictor(s): ", paste(numeric_predictors, collapse = ", "),
      ". Numeric predictors are measurements, not group labels. If your goal is to estimate their relationship with the outcome, use linear_regression_workflow() instead."
    )
  }

  cat("\n====================================\n")
  cat("ANOVA WORKFLOW\n")
  cat("====================================\n")
  cat("I fitted the ANOVA linear model automatically with lm(). I am checking its residuals and group spread below because these assumptions belong to the fitted model.\n")
  cat("Model: ", deparse(formula), "\n", sep = "")
  cat("Complete observations used:", nrow(model_frame), "\n")

  anova_table <- anova(model)
  cat("\nANOVA table:\n")
  print(anova_table)

  cat("\nANOVA interpretation:\n")
  anova_p_values <- anova_table[["Pr(>F)"]]
  for (term_number in seq_along(anova_p_values)) {
    term_name <- rownames(anova_table)[term_number]
    p_value <- anova_p_values[term_number]
    is_numeric_term <- term_name %in% numeric_predictors

    if (is.na(p_value)) {
      next
    }

    if (p_value < alpha) {
      if (is_numeric_term) {
        cat("-", term_name, ": p =", format.pval(p_value),
            "-> evidence this numeric predictor is associated with the response.\n")
      } else {
        cat("-", term_name, ": p =", format.pval(p_value),
            "-> evidence that means differ across this term's groups.\n")
      }
    } else {
      if (is_numeric_term) {
        cat("-", term_name, ": p =", format.pval(p_value),
            "-> no strong evidence this numeric predictor is associated with the response.\n")
      } else {
        cat("-", term_name, ": p =", format.pval(p_value),
            "-> no strong evidence that means differ across this term's groups.\n")
      }
    }
  }

  residual_diagnostics <- data.frame(
    observation = seq_along(residual_values),
    fitted_value = fitted_values,
    residual = residual_values,
    standardized_residual = rstandard(model)
  )

  cat("\nResidual diagnostics:\n")
  cat("Residuals are observed values minus model-predicted values.\n")
  cat("A residual near 0 means the model predicted that observation well.\n")
  cat("Large positive or negative residuals are observations the model predicted poorly.\n")
  cat("\nFirst residuals:\n")
  print(head(residual_diagnostics, 10), row.names = FALSE)
  cat("\nResidual summary:\n")
  print(summary(residual_values))

  shapiro_result <- NULL
  if (length(residual_values) >= 3 && length(residual_values) <= 5000) {
    shapiro_result <- shapiro.test(residual_values)
    cat("\nShapiro-Wilk test on model residuals:\n")
    print(shapiro_result)

    if (shapiro_result$p.value < alpha) {
      cat("Interpretation: evidence against normally distributed residuals.\n")
      cat("This suggests inspecting the Q-Q plot and considering transformation or a different model.\n")
    } else {
      cat("Interpretation: no significant evidence against normal residuals.\n")
      cat("This assumption appears reasonable based on this test.\n")
    }
  }

  grouping_data <- model_frame[predictor_names]
  group_label <- do.call(
    interaction,
    c(grouping_data, list(drop = TRUE, sep = " × "))
  )

  response_values <- model_frame[[response_name]]
  levene_result <- NULL
  levene_p_value <- NA_real_
  if (requireNamespace("car", quietly = TRUE)) {
    levene_data <- data.frame(response = response_values, group = group_label)
    levene_result <- car::leveneTest(response ~ group, data = levene_data)
    cat("\nLevene's test for equal variances:\n")
    print(levene_result)
    levene_p_value <- levene_result[["Pr(>F)"]][1]
    if (levene_p_value < alpha) {
      cat("Interpretation: evidence that group variances are unequal.\n")
      if (length(predictor_names) == 1 &&
          (is.factor(model_frame[[predictor_names]]) ||
           is.character(model_frame[[predictor_names]]))) {
        cat("Next step for this one-factor design: consider Welch's one-way ANOVA with oneway.test(formula, data = data, var.equal = FALSE). The usual Tukey HSD assumes equal variances; use a suitable unequal-variance pairwise method, such as Games-Howell, if follow-up comparisons are needed.\n")
      } else {
        cat("The usual equal-variance ANOVA follow-up may not fit this design. Discuss an unequal-variance approach with your instructor or use a model suited to this study design.\n")
      }
    } else {
      cat("Interpretation: no significant evidence that group variances differ.\n")
    }
  } else {
    cat("\nLevene's test skipped. To enable it, run:\n")
    cat('install.packages("car")\n')
  }

  split_values <- split(response_values, group_label)
  group_summary <- do.call(rbind, lapply(names(split_values), function(group_name) {
    x <- split_values[[group_name]]
    n <- length(x)
    mean_value <- mean(x)
    sd_value <- if (n > 1) sd(x) else NA_real_
    se_value <- if (n > 1) sd_value / sqrt(n) else NA_real_
    margin <- if (n > 1) qt(.975, df = n - 1) * se_value else NA_real_

    data.frame(
      group = group_name,
      n = n,
      mean = mean_value,
      sd = sd_value,
      se = se_value,
      lower_95_CI = mean_value - margin,
      upper_95_CI = mean_value + margin
    )
  }))
  rownames(group_summary) <- NULL

  cat("\nGroup summary:\n")
  print(group_summary, row.names = FALSE)

  anova_ss <- anova_table[["Sum Sq"]]
  eta_squared <- sum(anova_ss[-length(anova_ss)]) / sum(anova_ss)
  cat("\nOverall eta-squared:", round(eta_squared, 4), "\n")
  eta_interpretation <- if (eta_squared < 0.01) {
    "very small effect"
  } else if (eta_squared < 0.06) {
    "small effect"
  } else if (eta_squared < 0.14) {
    "medium effect"
  } else {
    "large effect"
  }
  cat("Approximate effect-size interpretation:", eta_interpretation, "\n")

  tukey_result <- NULL
  emmeans_result <- NULL
  if (run_tukey) {
    if (length(predictor_names) != 1 ||
        !(is.factor(model_frame[[predictor_names[1]]]) ||
          is.character(model_frame[[predictor_names[1]]]))) {
      cat("\nTukey follow-up was skipped. This workflow only runs it automatically for a one-factor ANOVA with a categorical group. For multiple factors or interactions, use planned comparisons or estimated marginal means for the specific question.\n")
    } else {
      factor_name <- predictor_names[1]
      number_of_groups <- nlevels(droplevels(as.factor(model_frame[[factor_name]])))
      term_row <- match(factor_name, rownames(anova_table))
      omnibus_p <- if (length(term_row) == 1 && !is.na(term_row)) {
        anova_p_values[term_row]
      } else {
        NA_real_
      }

      if (number_of_groups < 3) {
        cat("\nTukey follow-up was skipped: this factor has only two groups, so there is only one pair to compare. The ANOVA result already tests that comparison.\n")
      } else if (!is.na(levene_p_value) && levene_p_value < alpha) {
        cat("\nTukey follow-up was skipped because Levene's test found evidence of unequal group variances. For a one-factor comparison, consider Welch's ANOVA and a pairwise method designed for unequal variances, such as Games-Howell.\n")
      } else if (is.na(omnibus_p) || omnibus_p >= alpha) {
        cat("\nTukey follow-up was skipped because the overall ANOVA was not statistically significant at alpha = ", alpha, ". First interpret the overall test; it did not provide evidence of a difference among group means. If you specified a particular comparison before looking at the results, that planned comparison is a separate question from an all-pairs Tukey test.\n", sep = "")
      } else if (requireNamespace("emmeans", quietly = TRUE)) {
        if (is.na(levene_p_value)) {
          cat("\nI could not run Levene's variance check because the optional `car` package is unavailable. Review the residual plots and equal-variance assumption before interpreting Tukey comparisons.\n")
        }
        cat("\nThe overall ANOVA was significant. Tukey follow-up checks which group pairs differ while adjusting for the fact that all pairs are being checked. This assumes the ANOVA design and residual checks are reasonably suitable.\n")
        emmeans_result <- list(
          estimated_means = emmeans::emmeans(model, specs = factor_name),
          pairwise_comparisons = emmeans::pairs(
            emmeans::emmeans(model, specs = factor_name),
            adjust = "tukey"
          )
        )
        cat("\nEstimated marginal means:\n")
        print(emmeans_result$estimated_means)
        cat("\nTukey-adjusted pairwise comparisons:\n")
        print(emmeans_result$pairwise_comparisons)
        cat("Read the adjusted p-value column: pairs below alpha =", alpha,
            "show evidence that those group means differ. Tukey adjustment accounts for checking all pairs in this factor.\n")
      } else {
        if (is.na(levene_p_value)) {
          cat("\nI could not run Levene's variance check because the optional `car` package is unavailable. Review the residual plots and equal-variance assumption before interpreting Tukey comparisons.\n")
        }
        cat("\nThe overall ANOVA was significant. Tukey follow-up checks which group pairs differ while adjusting for the fact that all pairs are being checked. This assumes the ANOVA design and residual checks are reasonably suitable.\n")
        tukey_result <- TukeyHSD(aov(model))
        cat("\nTukey HSD post-hoc comparisons:\n")
        print(tukey_result)
        cat("Read the adjusted p-value column: pairs below alpha =", alpha,
            "show evidence that those group means differ. Tukey adjustment accounts for checking all pairs in this factor.\n")
      }
    }
  }

  enhanced_results <- NULL
  if (enhanced) {
    enhanced_results <- enhanced_model_report(
      model,
      diagnostic_plots = plots,
      effect_sizes = TRUE,
      model_checks = TRUE
    )

    # rstatix provides an additional transparent one-way decision check.
    if (requireNamespace("rstatix", quietly = TRUE) &&
        length(predictor_names) == 1 &&
        (is.factor(model_frame[[predictor_names]]) ||
         is.character(model_frame[[predictor_names]]))) {
      rstatix_result <- tryCatch(
        rstatix::check_test_assumptions(
          data = model_frame,
          formula = formula
        ),
        error = function(e) {
          message("rstatix assumption recommendation was unavailable: ", e$message)
          NULL
        }
      )

      if (!is.null(rstatix_result)) {
        enhanced_results$rstatix_assumptions <- rstatix_result
        cat("\nrstatix assumption-based recommendation:\n")
        print(rstatix_result)
      }
    }
  }

  if (plots) {
    old_par <- par(mfrow = c(1, 2))
    on.exit(par(old_par), add = TRUE)

    qq_plot_with_ci(
      residual_values,
      main = "Q-Q Plot of ANOVA Residuals"
    )
    plot(fitted_values, residual_values,
         main = "Residuals vs. Fitted Values",
         xlab = "Fitted values",
         ylab = "Residuals")
    abline(h = 0, col = "red", lwd = 2)

    cat("\nPlot interpretation:\n")
    cat("- Q-Q plot: residuals should generally follow the red line and remain near the blue envelope.\n")
    cat("- Residuals-versus-fitted plot: residuals should be randomly scattered around zero without a funnel or curve.\n")
  }

  invisible(list(
    model = model,
    anova = anova_table,
    residual_diagnostics = residual_diagnostics,
    shapiro = shapiro_result,
    levene = levene_result,
    group_summary = group_summary,
    eta_squared = eta_squared,
    tukey = tukey_result,
    emmeans = emmeans_result,
    enhanced = enhanced_results
  ))
}


anova_assumption_report <- function(data = NULL,
                                    formula = NULL,
                                    alpha = 0.05,
                                    plots = TRUE) {
  anova_workflow(
    data = data,
    formula = formula,
    alpha = alpha,
    plots = plots,
    run_tukey = FALSE
  )
}


# ---------------------------------------------------------------
# 10. T TESTS
# ---------------------------------------------------------------

t_test_guide <- function() {
  cat("\nT TEST PICKER\n")
  cat("1. One-sample t-test: one numeric sample compared with a reference value (for example, mean body size versus 10 cm).\n")
  cat("2. Independent two-sample t-test: compare the means of two separate, unrelated groups (for example, treatment A versus treatment B). This toolkit uses Welch's version by default because it does not assume equal group variances. The equal-variance Student version is available with var.equal = TRUE when that assumption is justified.\n")
  cat("3. Paired t-test: compare two measurements from the same individuals or deliberately matched pairs (for example, before versus after). Each row must keep the two measurements for the same individual together.\n")
  cat("If you have three or more independent groups, start with ANOVA instead of running many t-tests. If the same individuals are measured at three or more times, use a repeated-measures method; ordinary one-way ANOVA is not a substitute.\n")
  cat("Choose a one-sided alternative only when its direction was part of the hypothesis before looking at the results. Otherwise use the two-sided default.\n")
  invisible(c("one_sample_t_test", "independent_t_test", "paired_t_test"))
}


one_sample_t_test <- function(x = NULL, mu = 0,
                              alternative = "two.sided",
                              alpha = 0.05,
                              plots = TRUE) {
  .toolkit_require_inputs("one_sample_t_test()", list(x = x),
                          c(x = "one numeric measurement column/vector"))
  .toolkit_check_alpha(alpha, "one_sample_t_test()")
  x <- .toolkit_clean_numeric_vector(x, "one_sample_t_test()")
  if (!is.numeric(mu) || length(mu) != 1 || is.na(mu) || !is.finite(mu)) {
    stop("one_sample_t_test(): `mu` is the reference mean you are comparing against. Give one finite number, such as mu = 10.", call. = FALSE)
  }
  if (!is.character(alternative) || length(alternative) != 1 || is.na(alternative) ||
      !alternative %in% c("two.sided", "less", "greater")) {
    stop("one_sample_t_test(): `alternative` must be \"two.sided\", \"less\", or \"greater\", based on your hypothesis.", call. = FALSE)
  }
  if (!is.logical(plots) || length(plots) != 1 || is.na(plots)) {
    stop("one_sample_t_test(): `plots` should be TRUE to show plots or FALSE to skip them.", call. = FALSE)
  }
  if (length(x) < 2) {
    stop("one_sample_t_test(): At least two usable numeric observations are needed.", call. = FALSE)
  }
  if (sd(x) == 0) {
    stop("one_sample_t_test(): All usable values are identical, so the sample has no spread for estimating a t-test.", call. = FALSE)
  }

  result <- t.test(x, mu = mu, alternative = alternative,
                   conf.level = 1 - alpha)
  cat("\nONE-SAMPLE T-TEST\n")
  cat("Question: does the sample mean differ from the reference value", mu, "?\n")
  print(result)
  cat("The estimate and confidence interval describe the population mean; compare the interval with your reference value. The p-value tests whether the mean differs from that reference.\n")
  cat("The t-test assumes independent observations and, especially with a small sample, a roughly normal population without extreme outliers.\n")
  cat("Hypothesis direction:", alternative,
      if (alternative == "greater") " means the sample is predicted to be above the reference.\n" else if (alternative == "less") " means the sample is predicted to be below the reference.\n" else " checks for a difference in either direction.\n")
  if (result$p.value < alpha) {
    cat("At alpha =", alpha, "there is evidence that the population mean differs from the reference value in the direction of the alternative hypothesis.\n")
  } else {
    cat("At alpha =", alpha, "there is not strong evidence that the population mean differs from the reference value. This does not prove the means are equal.\n")
  }
  cat("Next check whether observations are independent and whether the data have severe outliers or extreme skew. The t-test is about a population mean, not whether every observation is near the reference value.\n")
  if (plots) {
    old_par <- par(mfrow = c(1, 2))
    on.exit(par(old_par), add = TRUE)
    hist(x, main = "Sample values", xlab = "Value")
    qqnorm(x, main = "Q-Q plot of sample values")
    qqline(x, col = "red", lwd = 2)
  }
  invisible(list(test = result, n = length(x), mean = mean(x), reference = mu))
}


independent_t_test <- function(data = NULL, outcome = NULL, group = NULL,
                               var.equal = FALSE,
                               alternative = "two.sided",
                               alpha = 0.05,
                               plots = TRUE) {
  .toolkit_require_inputs("independent_t_test()",
                          list(data = data, outcome = outcome, group = group),
                          c(data = "the imported worksheet", outcome = "the numeric measurement column", group = "the two-group column"))
  data <- .toolkit_prepare_data(data, "independent_t_test()")
  .toolkit_check_column(data, outcome, "outcome", "independent_t_test()")
  .toolkit_check_column(data, group, "group", "independent_t_test()")
  .toolkit_check_alpha(alpha, "independent_t_test()")
  if (!is.numeric(data[[outcome]])) {
    stop("independent_t_test(): The outcome column must contain numeric measurements. I received class ",
         paste(class(data[[outcome]]), collapse = "/"), ". If it contains numbers stored as text, inspect it before converting.", call. = FALSE)
  }
  if (!is.logical(var.equal) || length(var.equal) != 1 || is.na(var.equal)) {
    stop("independent_t_test(): `var.equal` must be TRUE or FALSE. FALSE (the default) runs Welch's test, which does not assume the two groups have equal variances.", call. = FALSE)
  }
  if (!is.character(alternative) || length(alternative) != 1 || is.na(alternative) ||
      !alternative %in% c("two.sided", "less", "greater")) {
    stop("independent_t_test(): `alternative` must be \"two.sided\", \"less\", or \"greater\", based on your hypothesis.", call. = FALSE)
  }
  if (!is.logical(plots) || length(plots) != 1 || is.na(plots)) {
    stop("independent_t_test(): `plots` should be TRUE to show a boxplot or FALSE to skip it.", call. = FALSE)
  }

  response <- data[[outcome]]
  group_values <- data[[group]]
  complete <- is.finite(response) & !is.na(group_values)
  if (sum(!complete) > 0) {
    message("independent_t_test(): I skipped ", sum(!complete),
            " row(s) with a missing group or non-finite outcome.")
  }
  response <- response[complete]
  group_values <- droplevels(as.factor(group_values[complete]))
  if (nlevels(group_values) != 2) {
    stop("independent_t_test(): This test needs exactly two groups, but I found ",
         nlevels(group_values), ". For 3+ independent groups, use anova_workflow(). For repeated measurements from the same subjects, use a paired/repeated-measures method instead; paired_t_test() handles two measurements only.", call. = FALSE)
  }
  group_names <- levels(group_values)
  first_group <- response[group_values == group_names[1]]
  second_group <- response[group_values == group_names[2]]
  if (length(first_group) < 2 || length(second_group) < 2) {
    stop("independent_t_test(): Each group needs at least two usable observations. Check missing values and the group column.", call. = FALSE)
  }

  result <- t.test(response ~ group_values,
                   var.equal = var.equal,
                   alternative = alternative,
                   conf.level = 1 - alpha)
  method_name <- if (var.equal) "Student two-sample t-test (equal-variance version)" else "Welch two-sample t-test (default; does not assume equal variances)"
  summary_table <- data.frame(
    group = group_names,
    n = c(length(first_group), length(second_group)),
    mean = c(mean(first_group), mean(second_group)),
    sd = c(sd(first_group), sd(second_group)),
    stringsAsFactors = FALSE
  )
  cat("\nINDEPENDENT TWO-SAMPLE T-TEST\n")
  cat("Method:", method_name, "\n")
  cat("The groups must contain different, unrelated observational units.\n")
  print(summary_table, row.names = FALSE)
  print(result)
  cat("The estimated difference is mean of ", group_names[1], " minus mean of ",
      group_names[2], ". The confidence interval is for that difference.\n", sep = "")
  cat("Hypothesis direction:", alternative,
      if (alternative == "greater") paste0(" means ", group_names[1], " is predicted to have the larger mean.\n") else if (alternative == "less") paste0(" means ", group_names[1], " is predicted to have the smaller mean.\n") else " checks for a difference in either direction.\n")
  if (result$p.value < alpha) {
    cat("At alpha =", alpha, "there is evidence that the two population means differ in the direction tested.\n")
  } else {
    cat("At alpha =", alpha, "there is not strong evidence that the two population means differ. This does not prove the means are equal.\n")
  }
  cat("Next check the group distributions for extreme outliers and severe skew, especially when group sizes are small, and confirm each observation belongs to only one group. For this toolkit's Shapiro and Levene checks, run assumption_report(data, response = \"", outcome, "\", group = \"", group, "\"). If the same subjects were measured twice, use paired_t_test() instead.\n", sep = "")
  if (plots) {
    boxplot(response ~ group_values, xlab = group, ylab = outcome,
            main = "Outcome by independent group")
  }
  invisible(list(test = result, group_summary = summary_table,
                 method = method_name, observations = length(response)))
}


paired_t_test <- function(data = NULL, before = NULL, after = NULL,
                          mu = 0,
                          alternative = "two.sided",
                          alpha = 0.05,
                          plots = TRUE) {
  .toolkit_require_inputs("paired_t_test()",
                          list(data = data, before = before, after = after),
                          c(data = "the worksheet with matched measurements", before = "the first measurement column", after = "the second measurement column"))
  data <- .toolkit_prepare_data(data, "paired_t_test()")
  .toolkit_check_column(data, before, "before", "paired_t_test()")
  .toolkit_check_column(data, after, "after", "paired_t_test()")
  if (identical(before, after)) {
    stop("paired_t_test(): `before` and `after` must name two different measurement columns.", call. = FALSE)
  }
  .toolkit_check_alpha(alpha, "paired_t_test()")
  if (!is.numeric(data[[before]]) || !is.numeric(data[[after]])) {
    stop("paired_t_test(): Both before/after columns must contain numeric measurements. Check their classes with str(data).", call. = FALSE)
  }
  if (!is.numeric(mu) || length(mu) != 1 || is.na(mu) || !is.finite(mu)) {
    stop("paired_t_test(): `mu` is the reference mean change (usually 0). Give one finite number.", call. = FALSE)
  }
  if (!is.character(alternative) || length(alternative) != 1 || is.na(alternative) ||
      !alternative %in% c("two.sided", "less", "greater")) {
    stop("paired_t_test(): `alternative` must be \"two.sided\", \"less\", or \"greater\", based on your hypothesis about after minus before.", call. = FALSE)
  }
  if (!is.logical(plots) || length(plots) != 1 || is.na(plots)) {
    stop("paired_t_test(): `plots` should be TRUE to show plots or FALSE to skip them.", call. = FALSE)
  }

  before_values <- data[[before]]
  after_values <- data[[after]]
  complete <- is.finite(before_values) & is.finite(after_values)
  if (sum(!complete) > 0) {
    message("paired_t_test(): I skipped ", sum(!complete),
            " row(s) without a usable value in both columns. Each row must represent the same subject/pair in both columns.")
  }
  before_values <- before_values[complete]
  after_values <- after_values[complete]
  differences <- after_values - before_values
  if (length(differences) < 2) {
    stop("paired_t_test(): At least two complete pairs are needed. Check that each row has both measurements.", call. = FALSE)
  }
  if (sd(differences) == 0) {
    stop("paired_t_test(): All paired differences are identical, so there is no variation for estimating a t-test.", call. = FALSE)
  }

  result <- t.test(after_values, before_values, paired = TRUE, mu = mu,
                   alternative = alternative, conf.level = 1 - alpha)
  cat("\nPAIRED T-TEST\n")
  cat("Complete pairs:", length(differences), "| Difference calculated as", after,
      "minus", before, "\n")
  print(result)
  cat("This test is performed on the within-pair differences, not as if the two columns came from unrelated groups.\n")
  cat("Hypothesis direction:", alternative,
      if (alternative == "greater") " means after minus before is predicted to exceed the reference change.\n" else if (alternative == "less") " means after minus before is predicted to be below the reference change.\n" else " checks for a difference in either direction.\n")
  if (result$p.value < alpha) {
    cat("At alpha =", alpha, "there is evidence of an average change/difference in the direction tested.\n")
  } else {
    cat("At alpha =", alpha, "there is not strong evidence of an average change/difference. This does not prove the measurements are identical.\n")
  }
  cat("Next check the distribution of the paired differences for severe outliers or skew; with small samples, it should be roughly normal. The pairing must be correct: values in the same row must belong to the same individual or matched pair.\n")
  if (plots) {
    old_par <- par(mfrow = c(1, 2))
    on.exit(par(old_par), add = TRUE)
    boxplot(differences, ylab = paste(after, "minus", before),
            main = "Paired differences")
    abline(h = mu, col = "red", lwd = 2)
    qqnorm(differences, main = "Q-Q plot of paired differences")
    qqline(differences, col = "red", lwd = 2)
  }
  invisible(list(test = result, differences = differences,
                 n_pairs = length(differences)))
}


# ---------------------------------------------------------------
# 11. CHI-SQUARE AND G-TESTS
# ---------------------------------------------------------------

# Exact test for one binary outcome compared with a stated probability.
binomial_test_workflow <- function(successes = NULL,
                                   trials = NULL,
                                   probability = 0.5,
                                   alternative = "two.sided",
                                   alpha = 0.05,
                                   success_label = "success") {
  .toolkit_require_inputs("binomial_test_workflow()",
                          list(successes = successes, trials = trials),
                          c(successes = "how many trials had the outcome", trials = "the total number of trials"))
  .toolkit_check_alpha(alpha, "binomial_test_workflow()")
  whole_number <- function(x) is.numeric(x) && length(x) == 1 &&
    !is.na(x) && is.finite(x) && x >= 0 && x == floor(x)
  if (!whole_number(trials) || trials < 1) {
    stop("binomial_test_workflow(): `trials` is the total number of independent yes/no outcomes. Give one positive whole number, such as 20.", call. = FALSE)
  }
  if (!whole_number(successes) || successes > trials) {
    stop("binomial_test_workflow(): `successes` is the number of trials with the outcome you are counting. Give a whole number from 0 through `trials`.", call. = FALSE)
  }
  if (!is.numeric(probability) || length(probability) != 1 ||
      is.na(probability) || !is.finite(probability) ||
      probability < 0 || probability > 1) {
    stop("binomial_test_workflow(): `probability` is the expected success chance under the null hypothesis. Give one number from 0 through 1, such as 0.5 or 0.95.", call. = FALSE)
  }
  if (!is.character(alternative) || length(alternative) != 1 ||
      is.na(alternative) || !alternative %in% c("two.sided", "less", "greater")) {
    stop("binomial_test_workflow(): `alternative` should be \"two.sided\", \"less\", or \"greater\". Choose a one-direction test only if that direction was predicted before examining the results.", call. = FALSE)
  }
  if (!is.character(success_label) || length(success_label) != 1 ||
      is.na(success_label) || !nzchar(success_label)) {
    stop("binomial_test_workflow(): `success_label` should be a short name for the outcome being counted, such as \"trees with fruit\".", call. = FALSE)
  }
  result <- stats::binom.test(successes, trials, p = probability,
                              alternative = alternative,
                              conf.level = 1 - alpha)
  cat("\nEXACT BINOMIAL TEST\n")
  cat("Question: is the chance of", success_label,
      "different from the expected probability?\n")
  cat("Counted:", successes, "of", trials, "trials | observed proportion:",
      round(successes / trials, 4), "| expected proportion:", probability, "\n")
  cat("H0: the true chance of", success_label, "is", probability, ".\n")
  if (alternative == "two.sided") {
    cat("HA: the true chance is different from", probability, "in either direction.\n")
  } else if (alternative == "less") {
    cat("HA: the true chance is lower than", probability, ".\n")
  } else {
    cat("HA: the true chance is higher than", probability, ".\n")
  }
  print(result)
  if (result$p.value < alpha) {
    cat("Plain-English result: the data provide evidence that the true chance differs from the expected value in the direction tested.\n")
  } else {
    cat("Plain-English result: the data do not provide strong evidence that the true chance differs from the expected value. This does not prove the expected value is exact.\n")
  }
  cat("Check: each trial should be an independent yes/no outcome with the same expected chance. If you are counting several categories, use goodness_of_fit_workflow() instead.\n")
  invisible(list(result = result, successes = successes, trials = trials,
                 observed_proportion = successes / trials,
                 expected_probability = probability,
                 alternative = alternative))
}

categorical_test_guide <- function() {
  cat("\nCHI-SQUARE / G-TEST PICKER\n")
  cat("1. One categorical variable versus an expected pattern: use goodness_of_fit_workflow() with observed counts and expected counts or probabilities. Example: do four color categories fit a 1:1:1:1 ratio?\n")
  cat("2. Two categorical variables: use contingency_test_workflow() to test whether they are associated. Give raw observations with row/column names, or give an already-counted table.\n")
  cat("Pearson chi-square and the G-test are two ways to compare observed counts with counts expected under the null hypothesis. They usually give similar conclusions when expected counts are adequate; follow the method required by your class.\n")
  cat("Sparse expected counts can make the usual chi-square/G p-value approximation unreliable. For a small 2x2 table consider Fisher's exact test; for larger tables consider simulation. Combine categories only when the categories can reasonably be combined scientifically.\n")
  cat("The rows being counted should be independent. Do not treat repeated measurements from the same individual as independent counts.\n")
  invisible(c("goodness_of_fit_workflow", "contingency_test_workflow"))
}


# Repeated-measures ANOVA for the same subjects measured at 2+ conditions.
# This follows the professor's afex::aov_ez() / emmeans workflow.
repeated_measures_anova_workflow <- function(data = NULL,
                                             outcome = NULL,
                                             subject = NULL,
                                             within = NULL,
                                             alpha = 0.05,
                                             posthoc = TRUE) {
  .toolkit_require_inputs("repeated_measures_anova_workflow()",
                          list(data = data, outcome = outcome,
                               subject = subject, within = within),
                          c(data = "the imported worksheet", outcome = "the numeric measurement column", subject = "the ID that links repeated measurements", within = "the column naming each repeated condition"))
  data <- .toolkit_prepare_data(data, "repeated_measures_anova_workflow()")
  .toolkit_check_alpha(alpha, "repeated_measures_anova_workflow()")
  for (item in c("outcome", "subject", "within")) {
    .toolkit_check_column(data, get(item), item,
                          "repeated_measures_anova_workflow()")
  }
  if (length(unique(c(outcome, subject, within))) != 3) {
    stop("repeated_measures_anova_workflow(): `outcome`, `subject`, and `within` must be three different columns. The subject column identifies who/what was measured more than once.", call. = FALSE)
  }
  if (!is.numeric(data[[outcome]])) {
    stop("repeated_measures_anova_workflow(): The outcome column must contain numeric measurements. Check whether numbers were imported as text.", call. = FALSE)
  }
  if (!is.logical(posthoc) || length(posthoc) != 1 || is.na(posthoc)) {
    stop("repeated_measures_anova_workflow(): `posthoc` should be TRUE to check condition pairs after a significant overall test, or FALSE to skip pair comparisons.", call. = FALSE)
  }
  if (!requireNamespace("afex", quietly = TRUE)) {
    stop("repeated_measures_anova_workflow(): The `afex` package is needed for this analysis. Run install_biostats_packages(), then source the toolkit again.", call. = FALSE)
  }
  selected <- data[, c(outcome, subject, within), drop = FALSE]
  complete <- stats::complete.cases(selected)
  if (any(!complete)) {
    message("repeated_measures_anova_workflow(): I skipped ", sum(!complete),
            " row(s) with a missing outcome, subject ID, or condition.")
  }
  selected <- selected[complete, , drop = FALSE]
  selected$.toolkit_subject <- as.factor(selected[[subject]])
  selected$.toolkit_condition <- as.factor(selected[[within]])
  if (nlevels(droplevels(selected$.toolkit_condition)) < 2) {
    stop("repeated_measures_anova_workflow(): I found fewer than two measured conditions. Check the `within` column and missing values.", call. = FALSE)
  }
  counts <- table(selected$.toolkit_subject, selected$.toolkit_condition)
  if (any(counts != 1)) {
    stop("repeated_measures_anova_workflow(): Each subject needs exactly one usable measurement in every condition for this simple repeated-measures workflow. Check duplicate rows, missing conditions, and the subject ID column. If your design is more complex, ask your instructor about a mixed-effects model.", call. = FALSE)
  }
  model <- afex::aov_ez(data = selected, dv = outcome,
                        id = ".toolkit_subject",
                        within = ".toolkit_condition")
  cat("\nREPEATED-MEASURES ANOVA\n")
  cat("Same subjects measured across:", within, "| Outcome:", outcome,
      "| Subjects:", nlevels(selected$.toolkit_subject), "\n")
  cat("H0: the mean outcome is the same across all conditions.\n")
  cat("HA: at least one condition has a different mean outcome.\n")
  cat("This analysis accounts for repeated measurements from the same subject. Check the assignment's assumptions and the sphericity correction shown in the output.\n\n")
  print(model)
  anova_table <- model$anova_table
  posthoc_result <- NULL
  p_column <- grep("Pr\\(>F\\)", names(anova_table), value = TRUE)
  omnibus_p <- if (length(p_column) > 0 && nrow(anova_table) > 0) {
    as.numeric(anova_table[[p_column[1]]][1])
  } else {
    NA_real_
  }
  if (posthoc && nlevels(selected$.toolkit_condition) > 2 &&
      !is.na(omnibus_p) && omnibus_p < alpha) {
    if (requireNamespace("emmeans", quietly = TRUE)) {
      means <- emmeans::emmeans(model, specs = ".toolkit_condition")
      posthoc_result <- emmeans::pairs(means, adjust = "tukey")
      cat("\nThe overall test is below alpha =", alpha,
          ". Tukey-adjusted condition comparisons:\n")
      print(posthoc_result)
      cat("Use the adjusted p-values to see which pairs differ.\n")
    } else {
      cat("\nThe overall test is below alpha =", alpha,
          ", so follow-up comparisons may be useful. To use the professor's Tukey-style workflow, install the optional `emmeans` package and rerun with posthoc = TRUE.\n")
    }
  } else if (posthoc && nlevels(selected$.toolkit_condition) > 2) {
    cat("\nNo Tukey pair comparisons were run because the overall condition test was not shown to be significant.\n")
  }
  invisible(list(model = model, anova = anova_table,
                 posthoc = posthoc_result, data_used = selected,
                 omnibus_p_value = omnibus_p))
}


# Random-intercept ANOVA: estimates how much variation occurs among sampled
# groups (for example, genets), rather than comparing each group as a treatment.
random_effects_anova_workflow <- function(data = NULL,
                                          outcome = NULL,
                                          random_group = NULL,
                                          alpha = 0.05) {
  .toolkit_require_inputs("random_effects_anova_workflow()",
                          list(data = data, outcome = outcome,
                               random_group = random_group),
                          c(data = "the imported worksheet", outcome = "the numeric measurement column", random_group = "the grouping factor sampled from a broader population"))
  data <- .toolkit_prepare_data(data, "random_effects_anova_workflow()")
  .toolkit_check_alpha(alpha, "random_effects_anova_workflow()")
  .toolkit_check_column(data, outcome, "outcome",
                        "random_effects_anova_workflow()")
  .toolkit_check_column(data, random_group, "random_group",
                        "random_effects_anova_workflow()")
  if (outcome == random_group) {
    stop("random_effects_anova_workflow(): `outcome` and `random_group` must be different columns.", call. = FALSE)
  }
  if (!is.numeric(data[[outcome]])) {
    stop("random_effects_anova_workflow(): The outcome must be a numeric measurement. Check whether numbers were imported as text.", call. = FALSE)
  }
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("random_effects_anova_workflow(): The `lme4` package is needed for a random-effects model. Run install_biostats_packages(), then source the toolkit again.", call. = FALSE)
  }
  selected <- data.frame(y = data[[outcome]],
                         group = as.factor(data[[random_group]]))
  complete <- stats::complete.cases(selected) & is.finite(selected$y)
  if (any(!complete)) {
    message("random_effects_anova_workflow(): I skipped ", sum(!complete),
            " row(s) with a missing/non-finite outcome or missing group label.")
  }
  selected <- selected[complete, , drop = FALSE]
  selected$group <- droplevels(selected$group)
  if (nlevels(selected$group) < 2) {
    stop("random_effects_anova_workflow(): I need at least two group levels to estimate variation among groups.", call. = FALSE)
  }
  if (any(table(selected$group) < 2)) {
    stop("random_effects_anova_workflow(): Every sampled group needs at least two observations for this simple variance-components workflow. Check the group labels or ask about a model for an unbalanced design.", call. = FALSE)
  }
  model <- lme4::lmer(y ~ 1 + (1 | group), data = selected,
                      REML = FALSE)
  null_model <- stats::lm(y ~ 1, data = selected)
  lr_statistic <- max(0, 2 * (as.numeric(stats::logLik(model)) -
                                as.numeric(stats::logLik(null_model))))
  course_p <- stats::pchisq(lr_statistic, df = 1, lower.tail = FALSE)
  boundary_p <- if (lr_statistic > 0) course_p / 2 else 1
  variance_table <- as.data.frame(lme4::VarCorr(model))
  group_variance <- variance_table$vcov[variance_table$grp == "group"][1]
  residual_variance <- variance_table$vcov[variance_table$grp == "Residual"][1]
  icc <- if (is.finite(group_variance + residual_variance) &&
             group_variance + residual_variance > 0) {
    group_variance / (group_variance + residual_variance)
  } else {
    NA_real_
  }
  cat("\nRANDOM-EFFECTS ANOVA (RANDOM INTERCEPT)\n")
  cat("Outcome:", outcome, "| Random grouping variable:", random_group,
      "| Groups:", nlevels(selected$group), "\n")
  cat("Use this when the groups were sampled from a broader population and you want to estimate among-group variation. If these are specific treatments you want to compare, use ordinary ANOVA instead.\n")
  cat("H0: there is no variation among the sampled groups.\n")
  cat("HA: some variation in the outcome is associated with group membership.\n\n")
  print(summary(model))
  cat("\nVariance components:\n")
  print(variance_table, row.names = FALSE)
  cat("Intraclass correlation (share of estimated variation among groups):",
      if (is.na(icc)) "not available" else round(icc, 4), "\n")
  cat("Likelihood-ratio statistic:", round(lr_statistic, 4),
      "| course-style chi-square p-value (1 df):", signif(course_p, 4), "\n")
  cat("Approximate boundary-adjusted p-value:", signif(boundary_p, 4),
      "\n")
  cat("The random-effect variance is tested at a boundary (variance cannot be below zero), so p-values are approximate. Your course script uses the ordinary 1-df chi-square comparison; follow your professor's required reporting method.\n")
  if (course_p < alpha) {
    cat("At alpha =", alpha,
        "the course-style test gives evidence of among-group variation. The ICC describes its estimated share, not cause.\n")
  } else {
    cat("At alpha =", alpha,
        "the course-style test does not give strong evidence of among-group variation. This does not prove the group variance is exactly zero.\n")
  }
  invisible(list(model = model, null_model = null_model,
                 variance_components = variance_table, ICC = icc,
                 likelihood_ratio = lr_statistic,
                 p_value_course = course_p,
                 p_value_boundary_adjusted = boundary_p,
                 data_used = selected))
}


# PERMANOVA with the matching group-dispersion check shown in the course
# scripts. Use numeric measurement columns only and interpret dispersion too.
permanova_workflow <- function(data = NULL,
                               variables = NULL,
                               group = NULL,
                               permutations = 999,
                               seed = 123,
                               standardize = FALSE) {
  .toolkit_require_inputs("permanova_workflow()",
                          list(data = data, variables = variables, group = group),
                          c(data = "the imported worksheet", variables = "one or more numeric measurement columns", group = "the grouping column"))
  data <- .toolkit_prepare_data(data, "permanova_workflow()")
  .toolkit_check_column(data, group, "group", "permanova_workflow()")
  if (!is.character(variables) || length(variables) == 0 || anyNA(variables)) {
    stop("permanova_workflow(): `variables` should list one or more numeric measurement columns in quotation marks. Do not include the group or ID column.", call. = FALSE)
  }
  for (variable in variables) {
    .toolkit_check_column(data, variable, "variables", "permanova_workflow()")
  }
  if (group %in% variables) {
    stop("permanova_workflow(): The grouping column cannot also be an outcome measurement. Remove it from `variables`.", call. = FALSE)
  }
  if (any(!vapply(data[variables], is.numeric, logical(1)))) {
    bad <- variables[!vapply(data[variables], is.numeric, logical(1))]
    stop("permanova_workflow(): These selected outcome columns are not stored as numbers: ",
         paste(bad, collapse = ", "), ". Check their meaning and import type before continuing.", call. = FALSE)
  }
  if (!is.numeric(permutations) || length(permutations) != 1 ||
      is.na(permutations) || !is.finite(permutations) ||
      permutations < 99 || permutations != floor(permutations)) {
    stop("permanova_workflow(): `permutations` is how many shuffled versions of the data to compare. Use a whole number of at least 99; 999 is a common starting value.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1 || is.na(seed) ||
      !is.finite(seed) || seed != floor(seed)) {
    stop("permanova_workflow(): `seed` should be one whole number so the random shuffling can be reproduced. The default 123 is fine.", call. = FALSE)
  }
  if (!is.logical(standardize) || length(standardize) != 1 || is.na(standardize)) {
    stop("permanova_workflow(): `standardize` should be TRUE or FALSE. Use TRUE only when variables measured on different scales should contribute equally.", call. = FALSE)
  }
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("permanova_workflow(): The `vegan` package is needed. Run install_biostats_packages(), then source the toolkit again.", call. = FALSE)
  }
  selected <- data[, c(variables, group), drop = FALSE]
  complete <- stats::complete.cases(selected)
  if (any(!complete)) {
    message("permanova_workflow(): I skipped ", sum(!complete),
            " row(s) with missing values in a selected measurement or group.")
  }
  selected <- selected[complete, , drop = FALSE]
  selected$.toolkit_group <- droplevels(as.factor(selected[[group]]))
  if (nlevels(selected$.toolkit_group) < 2 ||
      any(table(selected$.toolkit_group) < 2)) {
    stop("permanova_workflow(): I need at least two groups and at least two observations in every group.", call. = FALSE)
  }
  x <- as.matrix(selected[variables])
  if (standardize) x <- scale(x)
  if (any(!is.finite(x))) {
    stop("permanova_workflow(): The selected columns contain non-finite values or a constant column that cannot be standardized. Review the data first.", call. = FALSE)
  }
  distance <- stats::dist(x)
  had_random_seed <- exists(".Random.seed", envir = .GlobalEnv,
                            inherits = FALSE)
  if (had_random_seed) {
    previous_random_seed <- get(".Random.seed", envir = .GlobalEnv,
                                inherits = FALSE)
  }
  on.exit({
    if (had_random_seed) {
      assign(".Random.seed", previous_random_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv,
                      inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  fit <- vegan::adonis2(distance ~ .toolkit_group, data = selected,
                       permutations = permutations)
  dispersion <- vegan::betadisper(distance, selected$.toolkit_group)
  dispersion_anova <- stats::anova(dispersion)
  dispersion_permutation <- vegan::permutest(dispersion,
                                              permutations = permutations)
  cat("\nPERMANOVA (PERMUTATION-BASED GROUP COMPARISON)\n")
  cat("Measurements:", paste(variables, collapse = ", "),
      "| Group:", group, "| Complete rows:", nrow(selected), "\n")
  cat("H0: the groups do not differ in their multivariable distance pattern.\n")
  cat("HA: at least one group differs in that pattern.\n\n")
  print(fit)
  cat("\nGroup dispersion check (a significant result suggests groups differ in spread):\n")
  print(dispersion_anova)
  print(dispersion_permutation)
  cat("\nInterpret both outputs together: PERMANOVA can respond to differences in group centers and to differences in spread. If the dispersion check is significant, do not describe the PERMANOVA result as a clean difference in group centers. Confirm the permutation scheme matches the study design, especially for paired, blocked, or spatially related observations.\n")
  invisible(list(permanova = fit, dispersion_anova = dispersion_anova,
                 dispersion_permutation = dispersion_permutation,
                 distance = distance, data_used = selected,
                 seed = seed, permutations = permutations))
}

goodness_of_fit_workflow <- function(observed = NULL,
                                     expected = NULL,
                                     probabilities = NULL,
                                     test = "chisq",
                                     estimated_parameters = 0,
                                     alpha = 0.05) {
  .toolkit_require_inputs("goodness_of_fit_workflow()", list(observed = observed),
                          c(observed = "one count for each category"))
  .toolkit_check_alpha(alpha, "goodness_of_fit_workflow()")
  if (!is.character(test) || length(test) != 1 || is.na(test) ||
      !test %in% c("chisq", "g")) {
    stop("goodness_of_fit_workflow(): Choose test = \"chisq\" for Pearson's chi-square test or test = \"g\" for the G-test.", call. = FALSE)
  }
  if (!is.null(expected) && !is.null(probabilities)) {
    stop("goodness_of_fit_workflow(): Provide expected counts OR expected probabilities, not both.", call. = FALSE)
  }
  if (!is.numeric(estimated_parameters) || length(estimated_parameters) != 1 ||
      is.na(estimated_parameters) || !is.finite(estimated_parameters) ||
      estimated_parameters < 0 || estimated_parameters != floor(estimated_parameters)) {
    stop("goodness_of_fit_workflow(): `estimated_parameters` is how many model values (such as a Poisson mean) were estimated from these same data. Use a whole number 0 or larger; use 0 when none were estimated.", call. = FALSE)
  }

  category_names <- names(observed)
  if (is.matrix(observed)) {
    stop("goodness_of_fit_workflow(): `observed` should be one vector of category counts, not a two-dimensional table. Use contingency_test_workflow() for a table with two categorical variables.", call. = FALSE)
  }
  if (is.table(observed)) observed <- as.vector(observed)
  if (!is.numeric(observed) || !is.null(dim(observed))) {
    stop("goodness_of_fit_workflow(): `observed` must be a numeric vector of category counts, such as c(12, 8, 10). Do not give raw individual labels; count them first.", call. = FALSE)
  }
  if (length(observed) < 2 || anyNA(observed) || any(!is.finite(observed)) ||
      any(observed < 0) || any(observed != floor(observed))) {
    stop("goodness_of_fit_workflow(): Observed counts must be at least two nonnegative whole numbers with no missing values.", call. = FALSE)
  }
  if (is.null(category_names) || length(category_names) != length(observed)) {
    category_names <- paste("Category", seq_along(observed))
  }

  total <- sum(observed)
  if (total <= 0) {
    stop("goodness_of_fit_workflow(): The observed counts add up to zero. At least one count is needed.", call. = FALSE)
  }
  if (!is.null(probabilities)) {
    if (!is.numeric(probabilities) || length(probabilities) != length(observed) ||
        anyNA(probabilities) || any(!is.finite(probabilities)) ||
        any(probabilities < 0) || abs(sum(probabilities) - 1) > 1e-8) {
      stop("goodness_of_fit_workflow(): `probabilities` must have one nonnegative value per observed category and add up to 1. Example: c(0.5, 0.3, 0.2).", call. = FALSE)
    }
    expected <- probabilities * total
  }
  if (is.null(expected)) {
    stop("goodness_of_fit_workflow(): Tell me what counts are expected, either with `expected = c(...)` or `probabilities = c(...)`.", call. = FALSE)
  }
  if (!is.numeric(expected) || length(expected) != length(observed) ||
      anyNA(expected) || any(!is.finite(expected)) || any(expected <= 0)) {
    stop("goodness_of_fit_workflow(): Expected counts must be positive numbers, one for each observed category. If you have probabilities, supply them with `probabilities =` instead.", call. = FALSE)
  }
  if (abs(sum(expected) - total) > 1e-7 * max(1, total)) {
    stop("goodness_of_fit_workflow(): Expected counts must add up to the same total as observed counts (", total, "). If you know proportions instead, use `probabilities =`.", call. = FALSE)
  }

  degrees_freedom <- length(observed) - 1 - estimated_parameters
  if (degrees_freedom <= 0) {
    stop("goodness_of_fit_workflow(): The degrees of freedom are zero or below. Check the number of categories and how many parameters were estimated from the data.", call. = FALSE)
  }
  if (test == "chisq") {
    contribution <- (observed - expected)^2 / expected
  } else {
    contribution <- numeric(length(observed))
    positive_observed <- observed > 0
    contribution[positive_observed] <- 2 * observed[positive_observed] *
      log(observed[positive_observed] / expected[positive_observed])
  }
  statistic <- sum(contribution)
  p_value <- pchisq(statistic, df = degrees_freedom, lower.tail = FALSE)
  comparison <- data.frame(
    category = category_names,
    observed = observed,
    expected = expected,
    contribution_to_test_statistic = contribution,
    stringsAsFactors = FALSE
  )

  cat("\nGOODNESS-OF-FIT ", if (test == "chisq") "CHI-SQUARE" else "G-TEST", "\n", sep = "")
  cat("Question: do the observed category counts fit the expected counts/proportions?\n")
  cat("Pearson chi-square and G-tests ask the same general fit question using different formulas. Use the method your assignment or instructor specifies; both rely on an adequate expected-count approximation.\n")
  cat("Total observations:", total, "| Categories:", length(observed),
      "| Parameters estimated from these data:", estimated_parameters, "\n")
  print(comparison, row.names = FALSE)
  cat("Test statistic:", statistic, "| df:", degrees_freedom,
      "| p-value:", format.pval(p_value), "\n")
  if (p_value < alpha) {
    cat("At alpha =", alpha, "there is evidence that the observed counts do not fit the expected pattern. The contribution column shows which categories contribute most to the difference; it does not by itself identify the cause.\n")
  } else {
    cat("At alpha =", alpha, "there is not strong evidence that the observed counts differ from the expected pattern. This does not prove the expected pattern is exactly true.\n")
  }
  small_expected <- expected < 5
  if (any(expected < 1) || mean(small_expected) > 0.20) {
    cat("Caution: several expected counts are small. The usual chi-square/G-test p-value approximation may be poor. Combine categories only when that makes biological sense, or ask about an exact/simulation approach.\n")
  }
  cat("If you estimated a parameter from these same observations (for example, a Poisson mean), count that in `estimated_parameters` so the degrees of freedom are reduced appropriately.\n")

  invisible(list(test = test, statistic = statistic, df = degrees_freedom,
                 p_value = p_value, comparison = comparison,
                 expected = expected, observed = observed))
}


contingency_test_workflow <- function(data = NULL,
                                      row = NULL,
                                      column = NULL,
                                      counts = NULL,
                                      test = "chisq",
                                      simulate = FALSE,
                                      B = 2000,
                                      alpha = 0.05) {
  .toolkit_check_alpha(alpha, "contingency_test_workflow()")
  if (!is.character(test) || length(test) != 1 || is.na(test) ||
      !test %in% c("chisq", "g", "fisher")) {
    stop("contingency_test_workflow(): Choose test = \"chisq\" for Pearson's chi-square, test = \"g\" for the G-test, or test = \"fisher\" for Fisher's exact test.", call. = FALSE)
  }
  if (!is.logical(simulate) || length(simulate) != 1 || is.na(simulate)) {
    stop("contingency_test_workflow(): `simulate` must be TRUE or FALSE. TRUE estimates a p-value using repeated random tables when the chosen test supports it.", call. = FALSE)
  }
  if (!is.numeric(B) || length(B) != 1 || is.na(B) || !is.finite(B) ||
      B < 100 || B != floor(B)) {
    stop("contingency_test_workflow(): `B` is the number of random tables for simulation. Use a whole number of at least 100, such as 2000.", call. = FALSE)
  }
  if (test == "g" && simulate) {
    stop("contingency_test_workflow(): The G-test in this toolkit uses its usual chi-square approximation and does not simulate a p-value. For sparse tables, choose test = \"chisq\", simulate = TRUE, or use Fisher's exact test.", call. = FALSE)
  }

  if (is.table(data) || is.matrix(data)) {
    if (!is.null(row) || !is.null(column) || !is.null(counts)) {
      stop("contingency_test_workflow(): You supplied a table/matrix of counts, so leave `row`, `column`, and `counts` blank. Or supply raw data with row and column names instead.", call. = FALSE)
    }
    if (length(dim(data)) != 2) {
      stop("contingency_test_workflow(): The count table must have exactly two dimensions: one for each categorical variable.", call. = FALSE)
    }
    observed_table <- as.matrix(data)
    if (!is.numeric(observed_table)) {
      stop("contingency_test_workflow(): The table must contain numeric counts, not text labels.", call. = FALSE)
    }
  } else {
    data <- .toolkit_prepare_data(data, "contingency_test_workflow()")
    if (is.null(row) || is.null(column)) {
      stop("contingency_test_workflow(): For raw observations, name both categorical columns, such as row = \"Sex\", column = \"Survival\". For an existing count table, pass the table alone.", call. = FALSE)
    }
    .toolkit_check_column(data, row, "row", "contingency_test_workflow()")
    .toolkit_check_column(data, column, "column", "contingency_test_workflow()")
    if (identical(row, column)) {
      stop("contingency_test_workflow(): `row` and `column` must be two different categorical variables.", call. = FALSE)
    }

    if (is.null(counts)) {
      complete <- !is.na(data[[row]]) & !is.na(data[[column]])
      if (sum(!complete) > 0) {
        message("contingency_test_workflow(): I skipped ", sum(!complete),
                " row(s) with a missing category.")
      }
      row_values <- droplevels(as.factor(data[[row]][complete]))
      column_values <- droplevels(as.factor(data[[column]][complete]))
      observed_table <- table(row_values, column_values)
    } else {
      .toolkit_check_column(data, counts, "counts", "contingency_test_workflow()")
      if (!is.numeric(data[[counts]])) {
        stop("contingency_test_workflow(): The `counts` column must contain numeric counts.", call. = FALSE)
      }
      complete <- !is.na(data[[row]]) & !is.na(data[[column]]) &
        is.finite(data[[counts]])
      if (sum(!complete) > 0) {
        message("contingency_test_workflow(): I skipped ", sum(!complete),
                " row(s) with a missing category or unusable count.")
      }
      count_values <- data[[counts]][complete]
      if (any(count_values < 0) || any(count_values != floor(count_values))) {
        stop("contingency_test_workflow(): The `counts` column must contain nonnegative whole-number frequencies.", call. = FALSE)
      }
      row_values <- droplevels(as.factor(data[[row]][complete]))
      column_values <- droplevels(as.factor(data[[column]][complete]))
      observed_table <- matrix(0, nrow = nlevels(row_values),
                               ncol = nlevels(column_values),
                               dimnames = list(levels(row_values), levels(column_values)))
      for (i in seq_along(count_values)) {
        row_name <- as.character(row_values[i])
        column_name <- as.character(column_values[i])
        observed_table[row_name, column_name] <-
          observed_table[row_name, column_name] + count_values[i]
      }
    }
  }

  if (!is.numeric(observed_table) || anyNA(observed_table) ||
      any(!is.finite(observed_table)) || any(observed_table < 0) ||
      any(observed_table != floor(observed_table))) {
    stop("contingency_test_workflow(): Every cell in the table must be a nonnegative whole-number count, with no missing cells.", call. = FALSE)
  }
  if (nrow(observed_table) < 2 || ncol(observed_table) < 2) {
    stop("contingency_test_workflow(): The table needs at least two categories for each variable. Check that you selected two categorical columns with more than one observed category each.", call. = FALSE)
  }
  if (sum(observed_table) == 0 || any(rowSums(observed_table) == 0) ||
      any(colSums(observed_table) == 0)) {
    stop("contingency_test_workflow(): The table has a category with zero total count. Remove unused categories or check how the table was built.", call. = FALSE)
  }

  expected_table <- outer(rowSums(observed_table), colSums(observed_table)) /
    sum(observed_table)
  dimnames(expected_table) <- dimnames(observed_table)
  df <- (nrow(observed_table) - 1) * (ncol(observed_table) - 1)

  if (test == "chisq") {
    result <- chisq.test(observed_table, correct = FALSE,
                         simulate.p.value = simulate, B = B)
    statistic <- unname(result$statistic)
    p_value <- result$p.value
    expected_table <- result$expected
  } else if (test == "fisher") {
    result <- tryCatch(
      fisher.test(observed_table, simulate.p.value = simulate, B = B),
      error = function(e) stop("contingency_test_workflow(): Fisher's exact calculation could not finish for this table. Try `test = \"chisq\", simulate = TRUE` for a simulated chi-square p-value, or ask your instructor which method to use. R's message was: ", conditionMessage(e), call. = FALSE)
    )
    statistic <- NA_real_
    p_value <- result$p.value
  } else {
    contribution <- matrix(0, nrow = nrow(observed_table),
                           ncol = ncol(observed_table),
                           dimnames = dimnames(observed_table))
    positive_observed <- observed_table > 0
    contribution[positive_observed] <-
      2 * observed_table[positive_observed] *
      log(observed_table[positive_observed] / expected_table[positive_observed])
    statistic <- sum(contribution)
    p_value <- pchisq(statistic, df = df, lower.tail = FALSE)
    result <- list(statistic = statistic, parameter = df,
                   p.value = p_value, expected = expected_table,
                   method = "G-test of independence")
  }

  residual_table <- (observed_table - expected_table) / sqrt(expected_table)
  cat("\nCATEGORICAL ASSOCIATION TEST\n")
  cat("Question: are the two categorical variables associated?\n")
  cat("Test:", test, "| Total observations/counts:", sum(observed_table), "\n")
  if (test != "fisher") cat("Degrees of freedom:", df, "\n")
  cat("Observed counts:\n")
  print(observed_table)
  cat("\nExpected counts if the variables were unrelated:\n")
  print(round(expected_table, 3))
  cat("\nTest result:\n")
  print(result)
  if (test != "fisher") {
    cat("Pearson residuals by cell (positive = more observed than expected; negative = fewer):\n")
    print(round(residual_table, 2))
  }
  if (p_value < alpha) {
    cat("At alpha =", alpha, "there is evidence that the variables are associated. The test does not show which variable caused the other. Review the observed-versus-expected counts to see where the pattern occurs.\n")
  } else {
    cat("At alpha =", alpha, "there is not strong evidence of an association. This does not prove the variables are independent.\n")
  }
  small_cells <- expected_table < 5
  if (any(expected_table < 1) || mean(small_cells) > 0.20) {
    if (test == "fisher") {
      cat("The expected counts are sparse; Fisher's exact test is designed for this situation.\n")
    } else {
      cat("Caution: expected cell counts are sparse. For a small 2x2 table, consider Fisher's exact test. For larger tables, consider a simulated chi-square or Fisher p-value, or combine categories only when scientifically appropriate.\n")
    }
  }
  if (simulate) {
    cat("The p-value was estimated by simulation using", B,
        "random tables. More simulations make the estimate more stable but take longer.\n")
  }
  cat("The observations should be independent. If the same individuals contribute to multiple cells or measurements, this ordinary association test is not the right design.\n")

  invisible(list(test = test, result = result,
                 observed = observed_table, expected = expected_table,
                 residuals = residual_table, p_value = p_value,
                 df = df, simulated = simulate, B = if (simulate) B else NULL))
}


# Load optional packages when this master file is sourced.
load_biostats_packages()
