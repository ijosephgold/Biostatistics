# ==============================================================================
# Gbiostats_toolkit.r (v6) - UNIVERSAL CSUN BIOSTATISTICS, QUESTION-NLP,
#                            DERIVED-MATH, AUTO-TRANSFORM & SELF-CHECK ENGINE
# ==============================================================================
# NEW IN v6 (FOOLPROOF DERIVED-VARIABLE & QUESTION-READING INTELLIGENCE):
#   * Pass the Question Text Directly (`question = "..."` or `biostat_ask(data, "...")`):
#       biostat_wizard(data, "calculate plant size as the product of height and width")
#       -> Reads English cues ("product", "area", "plant size", "ratio", "difference",
#          "total", "Model III", "paired", "log-transform"), computes the new column
#          (e.g., Height * Width), drops the raw component columns so they don't
#          trigger a false ANCOVA, and runs the exact right test!
#   * Inline Formula Arithmetic (`Height * Width ~ Treatment`):
#       biostat_wizard(data, Height * Width ~ Treatment)
#       -> Automatically evaluates math on the left-hand side of `~` (`*`, `/`, `+`, `-`)
#          as a derived response variable instead of misinterpreting `*` as interaction!
#   * Auto-Detects Columns You Already Created (`erodiumdata$productHW`):
#       erodiumdata$productHW <- erodiumdata$Width * erodiumdata$Height
#       biostat_wizard(erodiumdata)
#       -> Automatically detects that `productHW` is the exact product of `Width` &
#          `Height`, selects `productHW` as Y, and ignores raw `Width` and `Height`!
# Built for zero-input execution (`biostat_wizard(data)`) + Automated Data
# Transformation Rescue (`log10`, `ln`, `sqrt`, `arcsine-sqrt`) + Full Self-Checking
# Transparency across ALL course modules & Problem Set 7 datasets:
#
#   * MASTER SUPERPOWER:
#     - biostat_wizard(data)         : Auto-detects 1, 2, or 3+ variables, checks
#                                      assumptions on RAW data, tests if a basic
#                                      transformation (log10, ln, sqrt, arcsine)
#                                      rescues normality/homoscedasticity BEFORE
#                                      jumping to non-parametric tests, runs BOTH
#                                      raw & transformed models when relevant,
#                                      plots diagnostic + final graphs, AND prints
#                                      copy-pasteable R verification code!
#
#   * NEW TRANSFORMATION & ASSUMPTION RESCUE ENGINE:
#     - biostat_transform(data)      : Tests Raw vs. log10(x), log(x), sqrt(x), &
#                                      asin(sqrt(x)) on your model residuals,
#                                      groups, or X-Y covariates; prints a side-by-
#                                      side comparison table of Shapiro-Wilk W/p &
#                                      Fligner/Bartlett variance p-values; returns
#                                      the best transformed dataframe; and plots
#                                      Before-vs-After Q-Q & distribution panels!
#
#   * SELF-CHECKING & STEP-BY-STEP CLARIFICATION FUNCTIONS:
#     - biostat_code(data)           : Prints a clean, copy-pasteable Base-R script
#                                      (including any log/sqrt transformation steps)
#                                      with side-by-side comments to check work!
#     - biostat_explain(data)        : Narrative breakdown of WHY a test/transform
#                                      was chosen & what each step did!
#     - biostat_check_math(data)     : Shows the raw under-the-hood math formulas
#                                      plugged in with your dataset's numbers!
#
#   * DEDICATED MULTI-VARIABLE & PROBLEM-SET SOLVERS:
#     - biostat_logistic(data)       : [PS7 Q1] Logistic GLM, Odds Factor exp(b1),
#                                      McFadden R2, LRT G2, & S-curve plot
#     - biostat_factorial(data)      : [PS7 Q2] Two-Way Model I ANOVA (Checks A*B ->
#                                      Separate t-tests if sig OR Additive A+B if NS;
#                                      auto-checks sqrt/log transform if skewed!)
#     - biostat_model3(data)         : [PS7 Q3] Mixed-Model III ANOVA (lme4 LLRT
#                                      Chi-Sq for Random Interaction + Fixed Effect;
#                                      compares Raw vs. sqrt/log transformed!)
#     - biostat_ancova(data)         : [PS7 Q4] ANCOVA (Tests X*Group slope int ->
#                                      Separate regressions if sig OR X+Group if NS;
#                                      auto-checks log-transformed covariates!)
#     - biostat_multireg(data)       : [PS7 Q5] Multiple Regression, Collinearity
#                                      SPLOM, All-Subsets AIC Ranking, & Best Model
# ==============================================================================

.has_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

.biostat_last_audit <- new.env(parent = emptyenv())

.sym_q <- function(col_name) {
  if (grepl("[^a-zA-Z0-9._]", col_name)) {
    paste0("`", col_name, "`")
  } else {
    col_name
  }
}

# ------------------------------------------------------------------------------
# 0.1. SKEWNESS, KURTOSIS, R/R2 INTERPRETATION & PROFESSOR-STYLE QQ HELPERS
# ------------------------------------------------------------------------------
.calc_skew_kurt <- function(x) {
  xc <- na.omit(as.numeric(x))
  if (length(xc) < 3 || sd(xc) == 0) {
    return(list(skew = NA, kurt = NA, skew_note = "Insufficient variation", kurt_note = "Insufficient variation"))
  }
  if (.has_pkg("moments")) {
    sk <- moments::skewness(xc)
    ku <- moments::kurtosis(xc)
  } else {
    n <- length(xc)
    m <- mean(xc)
    s <- sqrt(mean((xc - m)^2))
    sk <- mean(((xc - m) / s)^3)
    ku <- mean(((xc - m) / s)^4)
  }

  # Plain-English Skewness interpretation (0 = perfectly symmetric)
  skew_note <- if (abs(sk) < 0.5) {
    "Approximately symmetric (close to 0 -> good for normality)"
  } else if (abs(sk) < 1.0) {
    ifelse(sk > 0, "Moderately RIGHT-skewed (tail stretches high; consider sqrt/log if Shapiro fails)",
                   "Moderately LEFT-skewed (tail stretches low)")
  } else {
    ifelse(sk > 0, "Highly RIGHT-skewed (|skew| >= 1 -> strong candidate for log10/sqrt transform)",
                   "Highly LEFT-skewed (|skew| >= 1 -> strong departure from normality)")
  }

  # Plain-English Pearson Kurtosis interpretation (Normal distribution = 3.0)
  kurt_note <- if (ku >= 2.5 && ku <= 3.5) {
    "Mesokurtic (near 3.0 -> normal-like tail weight & peak)"
  } else if (ku > 3.5) {
    "Leptokurtic (> 3.5 -> sharper peak / heavier tails & potential outliers)"
  } else {
    "Platykurtic (< 2.5 -> flatter peak / lighter tails than normal)"
  }

  list(skew = sk, kurt = ku, skew_note = skew_note, kurt_note = kurt_note)
}

.interpret_r <- function(r_val) {
  ar <- abs(r_val)
  dir_str <- ifelse(r_val > 0, "positive", ifelse(r_val < 0, "negative", "zero"))
  mag_str <- if (ar < 0.1) "negligible" else if (ar < 0.3) "weak" else if (ar < 0.5) "moderate" else if (ar < 0.7) "moderately strong" else "strong"
  sprintf("%s %s linear association", tools::toTitleCase(mag_str), dir_str)
}

.plot_prof_qq <- function(x, main_title = "Normal Q-Q Plot", ylab_str = "Residuals") {
  xc <- na.omit(as.numeric(x))
  if (.has_pkg("car")) {
    # Exact professor-style car::qqPlot with red center line + shaded 95% confidence envelope!
    tryCatch(
      car::qqPlot(xc, main = main_title, ylab = ylab_str, col.lines = "firebrick",
                  envelope = list(style = "filled", col = "gray65", alpha = 0.28), pch = 19),
      error = function(e) {
        car::qqPlot(xc, main = main_title, ylab = ylab_str, col.lines = "firebrick", pch = 19)
      }
    )
  } else {
    qqnorm(xc, main = main_title, ylab = ylab_str, pch = 19, col = "steelblue")
    qqline(xc, col = "firebrick", lwd = 2)
  }
}

.plot_hist_normal <- function(x, main_title = "Histogram vs. Normal Curve", xlab_str = "Values") {
  xc <- na.omit(as.numeric(x))
  h <- hist(xc, freq = FALSE, col = "lightblue", border = "white",
            main = main_title, xlab = xlab_str)
  if (length(xc) >= 3 && sd(xc) > 0) {
    curve(dnorm(x, mean = mean(xc), sd = sd(xc)), add = TRUE, col = "firebrick", lwd = 2.2)
    legend("topright", legend = "Normal Curve", col = "firebrick", lwd = 2.2, bty = "n", cex = 0.85)
  }
  invisible(h)
}



# ------------------------------------------------------------------------------
# 0A. NATURAL LANGUAGE QUESTION & DERIVED VARIABLE ENGINE (.parse_question_and_derived)
# ------------------------------------------------------------------------------
# Reads plain-English problem prompts (or inline formula math / newly created columns)
# to automatically construct derived variables (Area = Height * Width, Ratio = A / B,
# Difference = Final - Initial, Sum = A + B, BMI, Density, Rate) and infer study flags!
# ------------------------------------------------------------------------------
.parse_question_and_derived <- function(df, response = NULL, predictor = NULL,
                                        question = NULL, model = "auto",
                                        paired = NULL, transform = "auto", verbose = TRUE) {
  df <- as.data.frame(df)
  cols <- names(df)
  num_cols <- cols[sapply(cols, function(c) is.numeric(df[[c]]) && length(unique(na.omit(df[[c]]))) > 2)]
  cat_cols <- cols[sapply(cols, function(c) is.factor(df[[c]]) || is.character(df[[c]]) || length(unique(na.omit(df[[c]]))) <= 6)]
  dropped_components <- character(0)
  derived_note <- NULL

  # 1. If the user passed a long plain-English sentence as `response` or `predictor`
  #    (e.g., biostat_wizard(data, "calculate plant size as the product of height and width")),
  #    treat it as the `question` text automatically!
  if (is.character(response) && length(response) == 1 && !response %in% cols &&
      (nchar(response) > 25 || grepl("\\s+", trimws(response)))) {
    # Check if it's an inline math string like "Height * Width" vs a full English sentence
    parsed_expr <- tryCatch(parse(text = response)[[1]], error = function(e) NULL)
    expr_vars <- if (!is.null(parsed_expr)) all.vars(parsed_expr) else character(0)
    if (!is.null(parsed_expr) && length(expr_vars) >= 1 && all(expr_vars %in% cols) &&
        grepl("[\\*\\+\\-/\\^]", response)) {
      new_col <- make.names(gsub("\\s+", "", response))
      df[[new_col]] <- eval(parsed_expr, envir = df)
      dropped_components <- setdiff(expr_vars, new_col)
      derived_note <- sprintf("Computed derived response `%s <- %s` (excluded component columns: %s)",
                              new_col, response, paste(dropped_components, collapse = ", "))
      response <- new_col
    } else {
      question <- paste(c(question, response), collapse = " ")
      response <- NULL
    }
  }
  if (is.character(predictor) && length(predictor) == 1 && !predictor %in% cols &&
      (nchar(predictor) > 25 || grepl("\\s+", trimws(predictor)))) {
    question <- paste(c(question, predictor), collapse = " ")
    predictor <- NULL
  }

  # 2. Parse Natural Language `question` Prompt (if provided)
  if (!is.null(question) && is.character(question) && nzchar(trimws(question))) {
    q_low <- tolower(question)

    # 2A. Detect Model III / Random Factor cues from question text
    if (model == "auto" && grepl("(model iii|model 3|mixed[- ]model|random effect|random factor|llrt|log[- ]likelihood)", q_low)) {
      model <- "III"
      if (verbose) cat("[QUESTION NLP]: Detected 'Model III / Random Effect' cue -> setting `model = 'III'`.\n")
    }

    # 2B. Detect Paired / Repeated Measures cues from question text
    if (is.null(paired) && grepl("(paired|before and after|pre and post|same individual|matched pair|repeated measure)", q_low)) {
      paired <- TRUE
      if (verbose) cat("[QUESTION NLP]: Detected 'paired / matched' cue -> setting `paired = TRUE`.\n")
    }

    # 2C. Detect explicit Transformation cues from question text
    if (transform == "auto") {
      if (grepl("(log[- ]transform|log10|logarithm)", q_low)) {
        transform <- "log10"
        if (verbose) cat("[QUESTION NLP]: Detected 'log-transform' cue in question text.\n")
      } else if (grepl("(square[- ]root|sqrt)", q_low)) {
        transform <- "sqrt"
        if (verbose) cat("[QUESTION NLP]: Detected 'square-root transform' cue in question text.\n")
      } else if (grepl("(arcsine|angular transform)", q_low)) {
        transform <- "asin"
        if (verbose) cat("[QUESTION NLP]: Detected 'arcsine transform' cue in question text.\n")
      }
    }

    # 2D. Identify which dataset columns are mentioned in the question text
    mentioned_num <- num_cols[sapply(num_cols, function(cn) {
      # Match column name or base stem (e.g., "height", "width", "length", "mass")
      stem <- tolower(gsub("(_mm|_cm|_m|_g|_kg|_n|_s|1|2)$", "", cn, ignore.case = TRUE))
      grepl(paste0("\\b", tolower(cn), "\\b"), q_low) || (nchar(stem) >= 3 && grepl(paste0("\\b", stem), q_low))
    })]

    # If 2 numeric columns aren't explicitly named, check for standard dimension pairs (Height/Width, Length/Width, Mass/Volume)
    if (length(mentioned_num) < 2 && length(num_cols) >= 2) {
      dim_cands <- num_cols[grepl("(height|width|length|breadth|depth|diameter|radius|mass|weight|volume|area|time|initial|final|start|end|pre|post)",
                                  num_cols, ignore.case = TRUE)]
      if (length(dim_cands) >= 2) mentioned_num <- dim_cands[1:2]
    }

    # 2E. Detect Mathematical Operations in the Question Prompt!
    if (is.null(response) && length(mentioned_num) >= 2) {
      c1 <- mentioned_num[1]; c2 <- mentioned_num[2]

      # Operation 1: PRODUCT / AREA / SIZE / VOLUME (e.g., "product of height and width", "area", "height times width")
      if (grepl("(product|multipl|times|\\barea\\b|plant size|leaf size|basal area|surface area|volume|height\\s*\\*\\s*width|width\\s*\\*\\s*height|length\\s*\\*\\s*width)", q_low)) {
        if (length(mentioned_num) >= 3 && grepl("(volume|3d|three[- ]dimension|length.*width.*height)", q_low)) {
          c3 <- mentioned_num[3]
          new_col <- paste0("Product_", c1, "_", c2, "_", c3)
          df[[new_col]] <- df[[c1]] * df[[c2]] * df[[c3]]
          dropped_components <- c(c1, c2, c3)
          derived_note <- sprintf("[QUESTION NLP]: Question asks for Volume/3D Product -> Computed `%s = %s * %s * %s`",
                                  new_col, c1, c2, c3)
        } else {
          new_col <- paste0("Product_", c1, "_", c2)
          df[[new_col]] <- df[[c1]] * df[[c2]]
          dropped_components <- c(c1, c2)
          derived_note <- sprintf("[QUESTION NLP]: Question asks for Area/Size/Product -> Computed `%s = %s * %s` (excluded raw %s & %s)",
                                  new_col, c1, c2, c1, c2)
        }
        response <- new_col

      # Operation 2: RATIO / PROPORTION / DENSITY / RATE / PER (e.g., "ratio of height to width", "force per area", "divided by")
      } else if (grepl("(ratio|divided by|quotient|\\bper\\b|aspect ratio|density|rate of|proportion of)", q_low)) {
        new_col <- paste0("Ratio_", c1, "_to_", c2)
        df[[new_col]] <- df[[c1]] / df[[c2]]
        dropped_components <- c(c1, c2)
        derived_note <- sprintf("[QUESTION NLP]: Question asks for Ratio/Rate -> Computed `%s = %s / %s` (excluded raw %s & %s)",
                                new_col, c1, c2, c1, c2)
        response <- new_col

      # Operation 3: DIFFERENCE / CHANGE / GAIN / LOSS (e.g., "difference between male and female", "change in mass", "final minus initial")
      } else if (grepl("(difference between|change in|net change|gain in|loss in|minus|subtract)", q_low)) {
        new_col <- paste0("Diff_", c1, "_minus_", c2)
        df[[new_col]] <- df[[c1]] - df[[c2]]
        dropped_components <- c(c1, c2)
        derived_note <- sprintf("[QUESTION NLP]: Question asks for Difference/Change -> Computed `%s = %s - %s` (excluded raw %s & %s)",
                                new_col, c1, c2, c1, c2)
        response <- new_col

      # Operation 4: SUM / TOTAL / COMBINED (e.g., "total biomass", "sum of", "combined length")
      } else if (grepl("(sum of|total|combined|added to|plus)", q_low)) {
        new_col <- paste0("Total_", c1, "_plus_", c2)
        df[[new_col]] <- df[[c1]] + df[[c2]]
        dropped_components <- c(c1, c2)
        derived_note <- sprintf("[QUESTION NLP]: Question asks for Sum/Total -> Computed `%s = %s + %s` (excluded raw %s & %s)",
                                new_col, c1, c2, c1, c2)
        response <- new_col
      }
    }
  }

  # 3. AUTO-DETECT USER-CREATED DERIVED COLUMNS ALREADY IN THE DATAFRAME!
  #    (e.g., You already ran `erodiumdata$productHW <- erodiumdata$Width * erodiumdata$Height`
  #     and then called `biostat_wizard(erodiumdata)` without specifying column names!)
  if (is.null(response) && length(num_cols) >= 3) {
    for (target_col in rev(num_cols)) {
      other_nums <- setdiff(num_cols, target_col)
      if (length(other_nums) >= 2) {
        pairs <- combn(other_nums, 2, simplify = FALSE)
        for (pr in pairs) {
          v_t <- df[[target_col]]
          v1  <- df[[pr[1]]]
          v2  <- df[[pr[2]]]
          ok  <- complete.cases(v_t, v1, v2)
          if (sum(ok) >= 3) {
            tol <- 1e-5 * max(abs(v_t[ok]), 1)
            is_prod  <- all(abs(v_t[ok] - (v1[ok] * v2[ok])) <= tol)
            is_sum   <- all(abs(v_t[ok] - (v1[ok] + v2[ok])) <= tol)
            is_diff1 <- all(abs(v_t[ok] - (v1[ok] - v2[ok])) <= tol)
            is_diff2 <- all(abs(v_t[ok] - (v2[ok] - v1[ok])) <= tol)
            is_rat1  <- all(v2[ok] != 0) && all(abs(v_t[ok] - (v1[ok] / v2[ok])) <= tol)
            is_rat2  <- all(v1[ok] != 0) && all(abs(v_t[ok] - (v2[ok] / v1[ok])) <= tol)

            if (is_prod || is_sum || is_diff1 || is_diff2 || is_rat1 || is_rat2) {
              op_sym <- if (is_prod) "*" else if (is_sum) "+" else if (is_diff1 || is_diff2) "-" else "/"
              response <- target_col
              dropped_components <- c(pr[1], pr[2])
              derived_note <- sprintf("[AUTO-DETECTED DERIVED COLUMN]: Column '%s' is the exact mathematical combination (%s %s %s)!\n  -> Automatically selected '%s' as Response (Y) and excluded raw component columns ('%s', '%s')!",
                                      target_col, pr[1], op_sym, pr[2], target_col, pr[1], pr[2])
              break
            }
          }
        }
      }
      if (!is.null(response)) break
    }
  }

  # 4. SMART MORPHOMETRIC DIMENSION DETECTOR (When no `question` and no `productHW` column was created yet!)
  #    If a dataset has 1 Categorical Treatment/Group column + 2 Morphometric Dimension columns
  #    like (`Height` & `Width`) or (`Length` & `Width`), alert the user prominently that
  #    `Height * Width` (Area/Size) is available with 1 flag or `biostat_derive()`!
  morph_hint <- NULL
  if (is.null(response) && length(dropped_components) == 0 && length(num_cols) == 2 && length(cat_cols) == 1) {
    has_hw <- all(sapply(c("height|length", "width|breadth|diameter"), function(pat) any(grepl(pat, num_cols, ignore.case = TRUE))))
    has_trt <- grepl("(treat|trt|group|plot|site|cond|water|shade|wind|control)", cat_cols[1], ignore.case = TRUE)
    if (has_hw && has_trt) {
      morph_hint <- sprintf(
        paste0("\n======================================================================\n",
               " [SMART DIMENSION ALERT]: Detected '%s' & '%s' alongside '%s'!\n",
               " * If your problem asks for PLANT/LEAF AREA or SIZE (Product = %s * %s)\n",
               "   rather than ANCOVA, simply run ONE of these:\n",
               "     1) biostat_wizard(data, %s * %s ~ %s)\n",
               "     2) biostat_wizard(data, question = 'product of %s and %s')\n",
               "======================================================================\n"),
        num_cols[1], num_cols[2], cat_cols[1],
        num_cols[1], num_cols[2],
        num_cols[1], num_cols[2], cat_cols[1],
        tolower(num_cols[1]), tolower(num_cols[2])
      )
    }
  }

  if (verbose && !is.null(derived_note)) {
    cat("----------------------------------------------------------------------\n")
    cat(derived_note, "\n")
    cat("----------------------------------------------------------------------\n")
  }
  if (verbose && !is.null(morph_hint)) {
    cat(morph_hint)
  }

  # Remove raw component columns from the working dataframe when a derived response is active
  if (length(dropped_components) > 0 && is.null(predictor)) {
    keep_cols <- setdiff(names(df), dropped_components)
    df <- df[, keep_cols, drop = FALSE]
  }

  list(
    df = df, response = response, predictor = predictor,
    model = model, paired = paired, transform = transform,
    dropped_components = dropped_components, derived_note = derived_note
  )
}

# ------------------------------------------------------------------------------
# 0B. UNIVERSAL STRUCTURAL & ROLE INFERENCE ENGINE (.auto_parse_roles)
# ------------------------------------------------------------------------------
.auto_parse_roles <- function(data, response = NULL, predictor = NULL, block = NULL,
                              paired = NULL, model = "auto", random = NULL,
                              question = NULL, transform = "auto",
                              wide_mode = "auto", verbose = TRUE) {
  has_interaction_fml <- FALSE

  # A1. Formula passed as first argument: e.g., biostat_wizard(Height * Width ~ Treatment, data)
  if (inherits(data, "formula")) {
    fml <- data
    if (is.null(response) || !is.data.frame(response)) {
      stop("When passing a formula first, provide the dataframe second: e.g., biostat_wizard(y ~ x, data)")
    }
    data <- as.data.frame(response)
    response <- fml
  }

  # A2. Formula passed as second argument: e.g., biostat_wizard(data, Height * Width ~ Treatment)
  if (inherits(response, "formula")) {
    fml <- response
    df_tmp <- as.data.frame(data)
    lhs_expr <- paste(deparse(fml[[2]]), collapse = " ")
    rhs_expr <- paste(deparse(fml[[3]]), collapse = " ")
    lhs_vars <- all.vars(fml[[2]])
    rhs_vars <- all.vars(fml[[3]])

    # If LHS is a mathematical expression of existing columns (like `Height * Width` or `I(Height * Width)` or `log10(Height)`):
    if (!lhs_expr %in% names(df_tmp) && length(lhs_vars) >= 1 && all(lhs_vars %in% names(df_tmp))) {
      clean_lhs <- gsub("^I\\((.*)\\)$", "\\1", lhs_expr)
      new_col_name <- make.names(gsub("\\s+", "", clean_lhs))
      df_tmp[[new_col_name]] <- eval(parse(text = clean_lhs), envir = df_tmp)
      # Remove raw component columns if they are not on the RHS of the formula!
      drop_c <- setdiff(lhs_vars, c(rhs_vars, new_col_name))
      if (length(drop_c) > 0) {
        if (verbose) {
          cat("----------------------------------------------------------------------\n")
          cat(sprintf("[INLINE FORMULA MATH]: Computed `%s = %s` as Response (Y)\n  (Excluded raw component columns: %s)\n",
                      new_col_name, clean_lhs, paste(drop_c, collapse = ", ")))
          cat("----------------------------------------------------------------------\n")
        }
        df_tmp <- df_tmp[, setdiff(names(df_tmp), drop_c), drop = FALSE]
      }
      data <- df_tmp
      response <- new_col_name
    } else {
      response <- lhs_vars[1]
    }
    predictor <- if (length(rhs_vars) > 0) rhs_vars else NULL
    has_interaction_fml <- grepl("\\*|:", rhs_expr)
  }

  # A3. Run Natural Language Question & Auto-Derived Column Detector on Dataframes
  if (is.data.frame(data)) {
    q_res <- .parse_question_and_derived(
      df = data, response = response, predictor = predictor,
      question = question, model = model, paired = paired,
      transform = transform, verbose = verbose
    )
    data      <- q_res$df
    response  <- q_res$response
    predictor <- q_res$predictor
    model     <- q_res$model
    paired    <- q_res$paired
    transform <- q_res$transform
  }

  if (is.table(data) || (is.matrix(data) && is.numeric(data) && nrow(data) <= 15 && ncol(data) <= 15)) {
    if (verbose) cat("[AUTO-DETECT]: Input is a pre-tabulated Contingency Table / Matrix.\n")
    return(list(mode = "table", table = as.table(data)))
  }

  if (is.atomic(data) && is.null(dim(data))) {
    vec_name <- deparse(substitute(data))
    df <- data.frame(Value = data)
    names(df) <- make.names(vec_name)
    return(list(
      mode = "dataframe", df = df, orig_df = df,
      response = names(df)[1], predictor = NULL, extra_predictors = NULL,
      block = NULL, random = NULL, paired = FALSE, wide_paired = FALSE, wide_kgroup = FALSE,
      model = model
    ))
  }

  df <- as.data.frame(data)
  orig_df <- df
  cols <- names(df)
  n_rows <- nrow(df)

  col_meta <- lapply(cols, function(cn) {
    v <- df[[cn]]
    v_clean <- na.omit(v)
    n_u <- length(unique(v_clean))
    n_clean <- length(v_clean)

    is_id_name <- grepl("^(id|sub|subject|patient|sample|specimen|tank|individual|mouse|rat|male|female|bird|fish|wing|butterfly|animal|plant|tree|leaf|person|pair|dyad|family|litter|clutch|nest|trial|rep|replicate|obs|index|row|colony)s?(_id|_num|_no)?$",
                        cn, ignore.case = TRUE)
    is_seq_int <- (is.numeric(v) && n_clean >= 3 && n_u == n_clean &&
                   all(v_clean == floor(v_clean)) &&
                   all(sort(v_clean) == seq(min(v_clean), max(v_clean))))
    is_unique_char <- ((is.character(v) || is.factor(v)) && n_clean >= 3 && n_u == n_clean)
    is_rep_col <- grepl("^(rep|replicate|subsample|lawn|quadrat|transect|sub|trial|obs|index)s?(_id|_num|_no)?$", cn, ignore.case = TRUE) && (is.numeric(v) || n_u > 12)

    if (is_seq_int || is_unique_char || (is_id_name && n_u == n_clean)) {
      role_type <- "id_or_block"
    } else if (is_rep_col && ncol(df) >= 4) {
      role_type <- "replicate_index"
    } else if (is.factor(v) || is.character(v) || is.logical(v)) {
      if (n_u == 2) {
        role_type <- "binary_cat"
      } else if (n_u <= 15) {
        role_type <- "multi_cat"
      } else {
        role_type <- "high_card_cat"
      }
    } else if (is.numeric(v)) {
      if (n_u == 2 && all(v_clean %in% c(0, 1))) {
        role_type <- "binary_num"
      } else if (n_u <= 4 && all(v_clean == floor(v_clean)) &&
                 grepl("(temp|co2|trt|treat|dose|level|group|plot|disturb|range|site|block)", cn, ignore.case = TRUE)) {
        role_type <- if (n_u == 2) "binary_cat" else "multi_cat"
      } else if (all(v_clean == floor(v_clean)) && min(v_clean) >= 0 && n_u <= 6) {
        role_type <- "small_integer"
      } else {
        role_type <- "continuous"
      }
    } else {
      role_type <- "other"
    }

    is_resp_name <- grepl("(resp|outcome|yield|growth|area|basal|weight|mass|height|size|score|rate|conc|density|count|freq|number|fruits|flowers|noflrs|seed|offspring|fecundity|fitness|survival|mortality|fertil|mating|mated|success|breed|bred|force|strength|break|newton|chirp|pulse|diversity|shannon|^h$|cells|nocells|chloroplast|symbiont|photosyn|abundance|productivity|bp|expression|activity|preference|biomass|speed|value|^y$)",
                          cn, ignore.case = TRUE)
    is_pred_name <- grepl("(group|treat|trt|trmt|clipped|weeded|clip|weed|dose|dosage|species|sp|sex|gender|diet|cond|condition|temp|temperature|co2|acid|ph|day|age|site|habitat|location|loc|range|origin|native|invas|disturb|plot|zone|depth|drug|fertilizer|genotype|strain|status|exposure|thickness|thallus|stipe|length|claw|diff|red|green|brown|blue|light|wave|reflect|elev|elevation|tpi|twi|slope|cover|pctcover|^x$)",
                          cn, ignore.case = TRUE)
    is_strong_resp <- grepl("(growth|area|fruits|flowers|noflrs|offspring|fitness|fertil|mating|mated|success|breed|bred|force|strength|break|chirp|pulse|diversity|shannon|cells|nocells|chloroplast|abundance|productivity|response|outcome|^y$)",
                            cn, ignore.case = TRUE)
    is_random_name <- grepl("(site|location|loc|range|origin|pop|population|block|plot|lawn|transect|region|field|meadow|forest|island|reef|stream|lake)",
                            cn, ignore.case = TRUE)
    is_freq_name <- grepl("^(freq|frequency|count|counts|n|num|number|observed|obs_count)$", cn, ignore.case = TRUE)

    list(
      name = cn, type = role_type, n_unique = n_u,
      is_id_name = is_id_name, is_seq_int = is_seq_int, is_rep_col = is_rep_col,
      is_resp_name = is_resp_name, is_strong_resp = is_strong_resp,
      is_pred_name = is_pred_name, is_random_name = is_random_name,
      is_freq_name = is_freq_name,
      mean = if (is.numeric(v) && n_clean > 0) mean(v_clean) else NA,
      sd   = if (is.numeric(v) && n_clean > 1) sd(v_clean) else NA
    )
  })
  names(col_meta) <- cols

  if (is.null(block)) {
    id_cands <- cols[sapply(col_meta, function(m) m$type == "id_or_block")]
    if (length(id_cands) > 0) block <- id_cands[1]
  }
  rep_cands <- cols[sapply(col_meta, function(m) m$type == "replicate_index")]

  active_cols <- setdiff(cols, c(block, rep_cands))
  if (length(active_cols) == 0) active_cols <- cols

  if (is.null(response) && is.null(predictor) && length(active_cols) == 3) {
    cat_c <- active_cols[sapply(active_cols, function(c) col_meta[[c]]$type %in% c("binary_cat", "multi_cat"))]
    freq_c <- active_cols[sapply(active_cols, function(c) col_meta[[c]]$is_freq_name && is.numeric(df[[c]]))]
    if (length(cat_c) == 2 && length(freq_c) == 1) {
      xt <- xtabs(as.formula(paste(.sym_q(freq_c), "~", paste(sapply(cat_c, .sym_q), collapse = " + "))), data = df)
      if (verbose) {
        cat(sprintf("[AUTO-DETECT]: Pre-aggregated frequency dataframe detected ('%s' counts across '%s' x '%s').\n",
                    freq_c, cat_c[1], cat_c[2]))
      }
      return(list(mode = "table", table = as.table(xt)))
    }
  }

  wide_paired <- FALSE
  wide_kgroup <- FALSE
  if (is.null(response) && is.null(predictor) && wide_mode != "long") {
    all_active_num <- all(sapply(active_cols, function(c) is.numeric(df[[c]]) && col_meta[[c]]$type %in% c("continuous", "small_integer")))

    if (all_active_num && length(active_cols) == 2) {
      c1 <- active_cols[1]; c2 <- active_cols[2]
      has_subject_id <- !is.null(block)
      has_paired_names <- any(grepl("(pre|post|before|after|baseline|followup|t1|t2|time1|time2|left|right|twin|control|treated|polar|depolar|light|dark|wet|dry|high|low|warm|cold|sun|shade|morning|evening|day|night|dominant|subordinate|infected|uninfected|mated|virgin)",
                                    active_cols, ignore.case = TRUE))
      one_strong_resp <- xor(col_meta[[c1]]$is_strong_resp, col_meta[[c2]]$is_strong_resp)
      m1 <- col_meta[[c1]]$mean; m2 <- col_meta[[c2]]$mean
      comparable_scale <- (!is.na(m1) && !is.na(m2) && max(abs(c(m1, m2))) > 0 &&
                           (min(abs(c(m1, m2))) / max(abs(c(m1, m2))) > 0.15))

      if (isTRUE(paired) || has_paired_names || (has_subject_id && !one_strong_resp && comparable_scale)) {
        wide_paired <- TRUE
        complete_rows <- complete.cases(df[, c(c1, c2)])
        n_comp <- sum(complete_rows)
        subj_ids <- if (!is.null(block)) df[[block]][complete_rows] else seq_len(n_comp)
        id_label <- if (!is.null(block)) block else "Pair_ID"

        df_long <- data.frame(
          Subject_ID = factor(rep(subj_ids, 2)),
          Condition  = factor(rep(c(c1, c2), each = n_comp), levels = c(c1, c2)),
          Value      = c(df[[c1]][complete_rows], df[[c2]][complete_rows])
        )
        names(df_long)[1] <- id_label

        if (verbose) {
          cat("----------------------------------------------------------------------\n")
          cat("AUTO-INFERRED STUDY STRUCTURE: WIDE PAIRED / MATCHED DESIGN\n")
          cat(sprintf("  * Subject / Unit ID   : %s (n = %d matched pairs)\n", id_label, n_comp))
          cat(sprintf("  * Paired Conditions   : '%s' vs. '%s'\n", c1, c2))
          cat(sprintf("  * Paired Difference   : d = %s - %s\n", c1, c2))
          cat("  * Note: To force Regression/Correlation instead, pass `paired = FALSE`\n")
          cat("----------------------------------------------------------------------\n")
        }
        return(list(
          mode = "dataframe", df = df_long, orig_df = orig_df,
          response = "Value", predictor = "Condition", extra_predictors = NULL,
          block = id_label, random = NULL, paired = TRUE, wide_paired = TRUE, wide_kgroup = FALSE,
          orig_cols = c(c1, c2), model = model
        ))
      }
    }
  }

  cont_cols <- active_cols[sapply(active_cols, function(c) col_meta[[c]]$type %in% c("continuous", "small_integer"))]
  cat_cols  <- active_cols[sapply(active_cols, function(c) col_meta[[c]]$type %in% c("binary_cat", "multi_cat", "binary_num"))]

  if (is.null(response)) {
    if (length(active_cols) == 1) {
      response <- active_cols[1]
      predictor <- NULL

    } else if (length(active_cols) == 2) {
      if (length(cont_cols) == 1 && length(cat_cols) == 1) {
        c_cont <- cont_cols[1]
        c_cat  <- cat_cols[1]
        cat_is_binary <- (col_meta[[c_cat]]$n_unique == 2)
        cat_looks_like_outcome <- col_meta[[c_cat]]$is_strong_resp || col_meta[[c_cat]]$type == "binary_num" ||
                                  grepl("(fertil|mating|mated|success|breed|bred|surviv|dead|alive|status|infected|disease|presence|outcome|resp|^y$)", c_cat, ignore.case = TRUE)
        cont_looks_like_pred <- col_meta[[c_cont]]$is_pred_name && !col_meta[[c_cont]]$is_strong_resp

        if (cat_is_binary && (cat_looks_like_outcome && (cont_looks_like_pred || col_meta[[c_cat]]$type == "binary_num"))) {
          response  <- c_cat
          predictor <- c_cont
        } else {
          response  <- c_cont
          predictor <- c_cat
        }
      } else if (length(cont_cols) == 2) {
        s_resp <- cont_cols[sapply(cont_cols, function(c) col_meta[[c]]$is_strong_resp)]
        resp_n <- cont_cols[sapply(cont_cols, function(c) col_meta[[c]]$is_resp_name && !col_meta[[c]]$is_pred_name)]
        pred_n <- cont_cols[sapply(cont_cols, function(c) col_meta[[c]]$is_pred_name && !col_meta[[c]]$is_strong_resp)]
        if (length(s_resp) == 1) {
          response  <- s_resp[1]
          predictor <- setdiff(cont_cols, response)[1]
        } else if (length(resp_n) >= 1) {
          response  <- resp_n[1]
          predictor <- setdiff(cont_cols, response)[1]
        } else if (length(pred_n) >= 1) {
          predictor <- pred_n[1]
          response  <- setdiff(cont_cols, predictor)[1]
        } else {
          response  <- cont_cols[2]
          predictor <- cont_cols[1]
        }
      } else {
        response  <- active_cols[2]
        predictor <- active_cols[1]
      }

    } else if (length(active_cols) >= 3) {
      if (length(cont_cols) == 1 && length(cat_cols) >= 2) {
        response  <- cont_cols[1]
        predictor <- cat_cols
      } else if (length(cont_cols) == 2 && length(cat_cols) == 1) {
        c_sresp <- cont_cols[sapply(cont_cols, function(c) col_meta[[c]]$is_strong_resp)]
        c_resp  <- cont_cols[sapply(cont_cols, function(c) col_meta[[c]]$is_resp_name && !col_meta[[c]]$is_pred_name)]
        c_pred  <- cont_cols[sapply(cont_cols, function(c) col_meta[[c]]$is_pred_name && !col_meta[[c]]$is_strong_resp)]

        if (length(c_sresp) == 1) {
          response <- c_sresp[1]
        } else if (length(c_resp) >= 1) {
          response <- c_resp[1]
        } else if (length(c_pred) == 1) {
          response <- setdiff(cont_cols, c_pred[1])[1]
        } else {
          response <- cont_cols[length(cont_cols)]
        }
        covar_col <- setdiff(cont_cols, response)[1]
        predictor <- c(covar_col, cat_cols[1])

      } else if (length(cont_cols) >= 3 && length(cat_cols) == 0) {
        c_sresp <- cont_cols[sapply(cont_cols, function(c) col_meta[[c]]$is_strong_resp)]
        c_resp  <- cont_cols[sapply(cont_cols, function(c) col_meta[[c]]$is_resp_name && !col_meta[[c]]$is_pred_name)]
        if (length(c_sresp) == 1) {
          response <- c_sresp[1]
        } else if (length(c_resp) == 1) {
          response <- c_resp[1]
        } else {
          if (col_meta[[cont_cols[1]]]$is_resp_name && !col_meta[[cont_cols[length(cont_cols)]]]$is_resp_name) {
            response <- cont_cols[1]
          } else {
            response <- cont_cols[length(cont_cols)]
          }
        }
        predictor <- setdiff(cont_cols, response)
      } else if (length(cont_cols) >= 2 && length(cat_cols) == 1) {
        # E.g. Byssal threads data.csv: Site (categorical), Strength (strong response), Length (covariate), Thread.count (extra continuous)
        c_sresp <- cont_cols[sapply(cont_cols, function(c) col_meta[[c]]$is_strong_resp)]
        response <- if (length(c_sresp) >= 1) c_sresp[1] else cont_cols[1]
        rem_cont <- setdiff(cont_cols, response)
        c_pred <- rem_cont[sapply(rem_cont, function(c) col_meta[[c]]$is_pred_name && !col_meta[[c]]$is_strong_resp)]
        covar_col <- if (length(c_pred) >= 1) c_pred[1] else rem_cont[1]
        predictor <- c(covar_col, cat_cols[1])
      } else {
        c_sresp <- cont_cols[sapply(cont_cols, function(c) col_meta[[c]]$is_strong_resp)]
        response <- if (length(c_sresp) >= 1) c_sresp[1] else if (length(cont_cols) >= 1) cont_cols[1] else active_cols[1]
        predictor <- setdiff(active_cols, response)
      }
    }
  } else if (is.null(predictor) && length(active_cols) >= 2) {
    predictor <- setdiff(active_cols, response)
  }

  extra_predictors <- setdiff(cols, c(response, predictor, block))

  if (length(predictor) == 2 && all(predictor %in% cat_cols)) {
    # Auto-detect Model III if user set model='III' OR if one factor is Location/Region/Lawn and the other is Genotype/Allele
    has_genotype_loc <- any(grepl("(genotype|allele|locus|clone)", predictor, ignore.case = TRUE)) &&
                        any(grepl("(location|region|lawn|meadow)", predictor, ignore.case = TRUE))
    if (is.null(random) && (has_genotype_loc || identical(model, "III") || identical(model, 3) || identical(tolower(model), "model3"))) {
      r_cands <- predictor[sapply(predictor, function(p) col_meta[[p]]$is_random_name)]
      random <- if (length(r_cands) >= 1) r_cands[1] else predictor[2]
      model <- "III"
    }
  }

  if (is.null(paired)) {
    paired <- FALSE
    if (!is.null(block) && !is.null(predictor) && length(predictor) == 1) {
      b_vec <- df[[block]]
      x_vec <- df[[predictor]]
      if (length(unique(na.omit(x_vec))) >= 2 && length(unique(na.omit(x_vec))) <= 10) {
        ct <- table(x_vec, b_vec)
        if (ncol(ct) >= 3 && all(ct == 1)) {
          paired <- TRUE
        }
      }
    }
  }

  if (verbose) {
    cat("----------------------------------------------------------------------\n")
    cat("AUTO-INFERRED VARIABLE ROLES (Override anytime via formula Y ~ X):\n")
    cat(sprintf("  * Dependent Response (Y) : %s\n", ifelse(is.null(response), "None", response)))
    cat(sprintf("  * Independent Factor(s) X: %s\n", ifelse(is.null(predictor), "None (One-Sample)", paste(predictor, collapse = ", "))))
    if (length(rep_cands) > 0) {
      cat(sprintf("  * Replicate/Subsample Col: %s (Ignored as predictor; used for design check)\n", paste(rep_cands, collapse = ", ")))
    }
    if (!is.null(random)) {
      fixed_f <- setdiff(predictor, random)
      cat(sprintf("  * Model III Structure    : Fixed = '%s' | Random = '%s'\n", paste(fixed_f, collapse = ", "), random))
    }
    if (!is.null(block)) {
      cat(sprintf("  * Subject / ID Column    : %s (Paired/Blocked = %s)\n", block, ifelse(paired, "TRUE", "FALSE")))
    }
    cat("----------------------------------------------------------------------\n")
  }

  list(
    mode = "dataframe", df = df, orig_df = orig_df,
    response = response, predictor = predictor,
    extra_predictors = extra_predictors, rep_cols = rep_cands,
    block = block, random = random, paired = paired,
    wide_paired = FALSE, wide_kgroup = FALSE,
    model = model, transform = transform,
    has_interaction_fml = has_interaction_fml,
    col_meta = col_meta
  )
}

# ------------------------------------------------------------------------------
# 1. INTERACTIVE MENU & CSUN COURSE WORKFLOW GUIDE
# ------------------------------------------------------------------------------
biostat_guide <- function() {
  cat("======================================================================\n")
  cat("   UNIVERSAL BIOSTATISTICS TOOLKIT v5 (AUTO-TRANSFORM & AUDIT EDITION)\n")
  cat("======================================================================\n")
  cat("1. ONE-LINE MASTER SUPERPOWER (Pass just `data`!):\n")
  cat("   - biostat_wizard(data)               : Auto-detects study design, checks\n")
  cat("                                          assumptions, tests log/sqrt/arcsine\n")
  cat("                                          transformations BEFORE jumping to\n")
  cat("                                          non-parametric tests, runs & plots!\n")
  cat("   - biostat_wizard(data, model = 'III'): Runs Model III Mixed-Effects ANOVA\n")
  cat("   - biostat_wizard(data, transform = 'log10') : Force a specific transform\n")
  cat("                                          ('auto', 'none', 'log10', 'ln', 'sqrt', 'asin')\n\n")
  cat("2. DATA TRANSFORMATION & ASSUMPTION RESCUE ENGINE:\n")
  cat("   - biostat_transform(data)            : Evaluates Raw vs. log10, ln, sqrt, &\n")
  cat("                                          arcsine-sqrt on model residuals & groups;\n")
  cat("                                          prints comparison table & Q-Q plots;\n")
  cat("                                          and returns the transformed dataframe!\n\n")
  cat("3. SELF-CHECKING, TRANSPARENCY & MANUAL VERIFICATION FUNCTIONS:\n")
  cat("   - biostat_code(data)                 : Prints copy-pasteable Base-R code\n")
  cat("                                          (with transformation steps!) to check work\n")
  cat("   - biostat_explain(data)              : Step-by-step plain-English walkthrough\n")
  cat("   - biostat_check_math(data)           : Plugs your dataset's numbers into formulas\n\n")
  cat("4. DEDICATED PROBLEM-SET 7 & ADVANCED GLM SOLVERS:\n")
  cat("   - biostat_logistic(data)  : [PS7 Q1] Binary Logistic GLM, exp(b1), McFadden R2\n")
  cat("   - biostat_factorial(data) : [PS7 Q2] Two-Way Model I ANOVA (+ transform check)\n")
  cat("   - biostat_model3(data)    : [PS7 Q3] Model III Mixed ANOVA (+ transform check)\n")
  cat("   - biostat_ancova(data)    : [PS7 Q4] ANCOVA (Raw + Log-Log covariate check!)\n")
  cat("   - biostat_multireg(data)  : [PS7 Q5] Multiple Regression, Collinearity & AIC\n")
  cat("======================================================================\n")
}

# ------------------------------------------------------------------------------
# 2. DATA TRANSFORMATION & ASSUMPTION RESCUE ENGINE (`biostat_transform`)
# ------------------------------------------------------------------------------
# Tests standard biological transformations when raw data/residuals violate
# normality (Shapiro-Wilk p < alpha) or homoscedasticity (Fligner/Bartlett p < alpha):
#   - log10(y)   [or log10(y + 1) if zeros exist] : Right-skewed sizes, weights, forces, lengths
#   - log(y)     [natural log ln(y) or ln(y + 1)] : Allometric / exponential growth data
#   - sqrt(y)    [square-root]                    : Count data (fruits, flowers, cells, colonies)
#   - asin(sqrt(y)) [arcsine-square-root]         : Proportions / percentages (0 to 1 or 0 to 100%)
#
# Also checks if a continuous predictor X in Regression/ANCOVA is strongly skewed
# (e.g., Thickness in Stipe data.csv!) and tests log10(X) alongside log10(Y)!
# ------------------------------------------------------------------------------
biostat_transform <- function(df, response = NULL, predictor = NULL, paired = NULL,
                              alpha = 0.05, plot = TRUE, verbose = TRUE) {
  parsed <- .auto_parse_roles(df, response = response, predictor = predictor, paired = paired, verbose = FALSE)
  if (parsed$mode == "table") {
    if (verbose) cat("Contingency tables use count frequencies directly; transformations do not apply.\n")
    return(invisible(NULL))
  }

  df_use    <- parsed$df
  response  <- parsed$response
  predictors <- parsed$predictor
  predictor <- predictors[1]
  paired    <- parsed$paired

  y_raw <- df_use[[response]]
  if (!is.numeric(y_raw) || (length(unique(na.omit(y_raw))) == 2 && all(na.omit(y_raw) %in% c(0, 1)))) {
    if (verbose) cat(sprintf("Response '%s' is binary/categorical; transformations do not apply.\n", response))
    return(invisible(NULL))
  }

  y_clean <- na.omit(y_raw)
  min_y   <- min(y_clean)
  max_y   <- max(y_clean)
  has_zero <- (min_y == 0)
  is_nonneg <- (min_y >= 0)
  is_prop   <- (min_y >= 0 && max_y <= 1)
  is_pct    <- (min_y >= 0 && max_y > 1 && max_y <= 100 && grepl("(pct|percent|prop|freq|rate)", response, ignore.case = TRUE))

  # Build candidate transformation functions for Y
  cand_list <- list(
    list(name = "Raw (Untransformed)", code = "none", label = response,
         r_expr = .sym_q(response),
         fn = function(v) v)
  )

  if (min_y > 0) {
    cand_list[[length(cand_list) + 1]] <- list(
      name = "Log10: log10(Y)", code = "log10", label = paste0("log10(", response, ")"),
      r_expr = sprintf("log10(%s)", .sym_q(response)),
      fn = function(v) log10(v)
    )
    cand_list[[length(cand_list) + 1]] <- list(
      name = "Natural Log: log(Y)", code = "ln", label = paste0("log(", response, ")"),
      r_expr = sprintf("log(%s)", .sym_q(response)),
      fn = function(v) log(v)
    )
  } else if (has_zero) {
    cand_list[[length(cand_list) + 1]] <- list(
      name = "Log10(Y + 1)", code = "log10", label = paste0("log10(", response, " + 1)"),
      r_expr = sprintf("log10(%s + 1)", .sym_q(response)),
      fn = function(v) log10(v + 1)
    )
    cand_list[[length(cand_list) + 1]] <- list(
      name = "Natural Log(Y + 1)", code = "ln", label = paste0("log(", response, " + 1)"),
      r_expr = sprintf("log(%s + 1)", .sym_q(response)),
      fn = function(v) log(v + 1)
    )
  }

  if (is_nonneg) {
    cand_list[[length(cand_list) + 1]] <- list(
      name = "Square-Root: sqrt(Y)", code = "sqrt", label = paste0("sqrt(", response, ")"),
      r_expr = sprintf("sqrt(%s)", .sym_q(response)),
      fn = function(v) sqrt(v)
    )
  }

  if (is_prop) {
    cand_list[[length(cand_list) + 1]] <- list(
      name = "Arcsine-Sqrt: asin(sqrt(Y))", code = "asin", label = paste0("asin(sqrt(", response, "))"),
      r_expr = sprintf("asin(sqrt(%s))", .sym_q(response)),
      fn = function(v) asin(sqrt(v))
    )
  } else if (is_pct) {
    cand_list[[length(cand_list) + 1]] <- list(
      name = "Arcsine-Sqrt: asin(sqrt(Y/100))", code = "asin", label = paste0("asin(sqrt(", response, "/100))"),
      r_expr = sprintf("asin(sqrt(%s / 100))", .sym_q(response)),
      fn = function(v) asin(sqrt(v / 100))
    )
  }

  # Helper to evaluate normality and variance homogeneity for a transformed Y vector
  eval_candidate <- function(y_tr, x_tr_override = NULL) {
    if (any(is.infinite(y_tr) | is.nan(y_tr))) {
      return(list(sw_w = NA, sw_p = NA, grp_min_p = NA, flig_p = NA, bart_p = NA,
                  norm_ok = FALSE, var_ok = FALSE, both_ok = FALSE, resids = y_tr))
    }

    if (is.null(predictor)) {
      yc <- na.omit(y_tr)
      sw <- if (length(yc) >= 3 && length(yc) <= 5000) shapiro.test(yc) else list(statistic = NA, p.value = NA)
      n_ok <- (!is.na(sw$p.value) && sw$p.value >= alpha)
      return(list(sw_w = as.numeric(sw$statistic), sw_p = sw$p.value, grp_min_p = sw$p.value,
                  flig_p = NA, bart_p = NA, norm_ok = n_ok, var_ok = TRUE, both_ok = n_ok, resids = yc))
    }

    # Multi-predictor or single-predictor evaluation
    tmp_df <- df_use
    tmp_df$.Y_TR <- y_tr
    if (!is.null(x_tr_override)) {
      tmp_df[[predictors[1]]] <- x_tr_override
    }
    tmp_df <- tmp_df[complete.cases(tmp_df[, c(".Y_TR", predictors)]), ]

    # Paired 2-group case: evaluate paired differences!
    x1 <- tmp_df[[predictor]]
    is_cat_x1 <- is.factor(x1) || is.character(x1) || (is.numeric(x1) && length(unique(x1)) <= 6)

    if (length(predictors) == 1 && is_cat_x1 && length(unique(x1)) == 2 && isTRUE(paired)) {
      xf <- as.factor(x1)
      g1 <- tmp_df$.Y_TR[xf == levels(xf)[1]]
      g2 <- tmp_df$.Y_TR[xf == levels(xf)[2]]
      diffs <- g1 - g2
      sw <- if (length(diffs) >= 3 && length(diffs) <= 5000 && length(unique(diffs)) > 1) shapiro.test(diffs) else list(statistic = NA, p.value = NA)
      n_ok <- (!is.na(sw$p.value) && sw$p.value >= alpha)
      return(list(sw_w = as.numeric(sw$statistic), sw_p = sw$p.value, grp_min_p = sw$p.value,
                  flig_p = NA, bart_p = NA, norm_ok = n_ok, var_ok = TRUE, both_ok = n_ok, resids = diffs))
    }

    # Build linear model to extract residuals
    if (length(predictors) == 1) {
      fml_e <- as.formula(paste(".Y_TR ~", .sym_q(predictor)))
      fit_e <- lm(fml_e, data = tmp_df)
      res_e <- residuals(fit_e)
      sw <- if (length(res_e) >= 3 && length(res_e) <= 5000 && length(unique(res_e)) > 1) shapiro.test(res_e) else list(statistic = NA, p.value = NA)

      if (is_cat_x1) {
        xf <- as.factor(tmp_df[[predictor]])
        grp_p <- tapply(tmp_df$.Y_TR, xf, function(v) {
          if (length(v) >= 3 && length(v) <= 5000 && length(unique(v)) > 1) shapiro.test(v)$p.value else NA
        })
        grp_min <- if (all(is.na(grp_p))) NA else min(grp_p, na.rm = TRUE)
        flig_p <- tryCatch(fligner.test(tmp_df$.Y_TR ~ xf)$p.value, error = function(e) NA)
        bart_p <- tryCatch(bartlett.test(tmp_df$.Y_TR ~ xf)$p.value, error = function(e) NA)
        n_ok <- (!is.na(sw$p.value) && sw$p.value >= alpha) && (is.na(grp_min) || grp_min >= alpha)
        v_ok <- (is.na(flig_p) || flig_p >= alpha)
      } else {
        # Continuous X: also check Shapiro on Y and X
        sw_y <- if (nrow(tmp_df) >= 3 && nrow(tmp_df) <= 5000) shapiro.test(tmp_df$.Y_TR)$p.value else NA
        grp_min <- sw_y
        flig_p <- NA; bart_p <- NA
        n_ok <- (!is.na(sw$p.value) && sw$p.value >= alpha)
        v_ok <- TRUE
      }
      return(list(sw_w = as.numeric(sw$statistic), sw_p = sw$p.value, grp_min_p = grp_min,
                  flig_p = flig_p, bart_p = bart_p, norm_ok = n_ok, var_ok = v_ok, both_ok = n_ok && v_ok, resids = res_e))

    } else {
      # Multi-predictor (2-Way ANOVA, Model III, ANCOVA, Multiple Regression)
      is_cats <- sapply(predictors, function(p) {
        v <- tmp_df[[p]]
        is.factor(v) || is.character(v) || length(unique(na.omit(v))) <= 6
      })
      if (length(predictors) == 2 && any(is_cats)) {
        fml_e <- as.formula(paste(".Y_TR ~", .sym_q(predictors[1]), "*", .sym_q(predictors[2])))
      } else {
        fml_e <- as.formula(paste(".Y_TR ~", paste(sapply(predictors, .sym_q), collapse = " + ")))
      }
      fit_e <- lm(fml_e, data = tmp_df)
      res_e <- residuals(fit_e)
      sw <- if (length(res_e) >= 3 && length(res_e) <= 5000 && length(unique(res_e)) > 1) shapiro.test(res_e) else list(statistic = NA, p.value = NA)

      if (all(is_cats)) {
        int_fac <- interaction(tmp_df[[predictors[1]]], tmp_df[[predictors[2]]])
        flig_p <- tryCatch(fligner.test(tmp_df$.Y_TR ~ int_fac)$p.value, error = function(e) NA)
        bart_p <- tryCatch(bartlett.test(tmp_df$.Y_TR ~ int_fac)$p.value, error = function(e) NA)
      } else {
        flig_p <- NA; bart_p <- NA
      }
      sw_y <- if (nrow(tmp_df) >= 3 && nrow(tmp_df) <= 5000) shapiro.test(tmp_df$.Y_TR)$p.value else NA
      n_ok <- (!is.na(sw$p.value) && sw$p.value >= alpha)
      v_ok <- (is.na(flig_p) || flig_p >= alpha)
      return(list(sw_w = as.numeric(sw$statistic), sw_p = sw$p.value, grp_min_p = sw_y,
                  flig_p = flig_p, bart_p = bart_p, norm_ok = n_ok, var_ok = v_ok, both_ok = n_ok && v_ok, resids = res_e))
    }
  }

  # Special check for ANCOVA or Simple Regression: what if BOTH Y and continuous X are log10-transformed?
  # (Exactly like LogLength in the Mussels lecture example or Thickness in Stipe data.csv!)
  has_cont_x <- (!is.null(predictor) && is.numeric(df_use[[predictors[1]]]) &&
                 length(unique(na.omit(df_use[[predictors[1]]]))) > 6 &&
                 min(df_use[[predictors[1]]], na.rm = TRUE) > 0)
  if (has_cont_x && min_y > 0) {
    x_name <- predictors[1]
    cand_list[[length(cand_list) + 1]] <- list(
      name = sprintf("Log-Log: log10(Y) & log10(%s)", x_name),
      code = "log_xy",
      label = paste0("log10(", response, ")"),
      x_label = paste0("log10(", x_name, ")"),
      r_expr = sprintf("log10(%s) [with log10(%s)]", .sym_q(response), .sym_q(x_name)),
      fn = function(v) log10(v),
      x_fn = function(vx) log10(vx)
    )
  }

  # Evaluate all candidates
  eval_results <- lapply(cand_list, function(cand) {
    y_tr <- suppressWarnings(cand$fn(df_use[[response]]))
    x_tr <- if (!is.null(cand$x_fn)) suppressWarnings(cand$x_fn(df_use[[predictors[1]]])) else NULL
    ev <- eval_candidate(y_tr, x_tr_override = x_tr)
    list(cand = cand, ev = ev, y_tr = y_tr, x_tr = x_tr)
  })

  # Build summary comparison table
  comp_rows <- lapply(eval_results, function(item) {
    c_info <- item$cand
    ev     <- item$ev
    status_str <- if (isTRUE(ev$both_ok)) {
      "PASS (All Met!)"
    } else if (isTRUE(ev$norm_ok) && !isTRUE(ev$var_ok)) {
      "Norm PASS / Var FAIL"
    } else if (!isTRUE(ev$norm_ok) && isTRUE(ev$var_ok) && !is.na(ev$flig_p)) {
      "Var PASS / Norm FAIL"
    } else {
      "FAIL"
    }

    data.frame(
      Transformation = c_info$name,
      Code           = c_info$code,
      Resid_Shapiro_W = ifelse(is.na(ev$sw_w), NA, round(ev$sw_w, 4)),
      Resid_Shapiro_p = ifelse(is.na(ev$sw_p), NA, signif(ev$sw_p, 4)),
      Y_or_Grp_Min_p  = ifelse(is.na(ev$grp_min_p), NA, signif(ev$grp_min_p, 4)),
      Var_Fligner_p   = ifelse(is.na(ev$flig_p), NA, signif(ev$flig_p, 4)),
      Status          = status_str,
      stringsAsFactors = FALSE
    )
  })
  comp_df <- do.call(rbind, comp_rows)

  # Select best transformation if Raw fails (or if user inspects transformations)
  raw_ev <- eval_results[[1]]$ev
  best_idx <- 1

  # Check if Raw failed residual normality, group/Y normality, or variance homogeneity
  raw_needs_help <- (!isTRUE(raw_ev$both_ok)) || (!is.na(raw_ev$grp_min_p) && raw_ev$grp_min_p < alpha)
  if (has_cont_x) {
    x_raw_sw <- shapiro.test(na.omit(df_use[[predictors[1]]]))$p.value
    if (x_raw_sw < alpha) raw_needs_help <- TRUE
  }

  if (raw_needs_help && length(eval_results) > 1) {
    # Score each candidate:
    # Priority 1: Both Residual Normality & Variance Homogeneity PASS (+100 pts)
    # Priority 2: Y/Group Normality PASS (+40 pts)
    # Priority 3: Variance Homogeneity PASS (+30 pts)
    # Priority 4: Highest Shapiro-Wilk W statistic on residuals (+ 20 * W)
    scores <- sapply(seq_along(eval_results), function(i) {
      ev <- eval_results[[i]]$ev
      if (is.na(ev$sw_w)) return(-999)
      sc <- 20 * ev$sw_w
      if (isTRUE(ev$both_ok)) sc <- sc + 100
      if (!is.na(ev$grp_min_p) && ev$grp_min_p >= alpha) sc <- sc + 40
      if (isTRUE(ev$var_ok)) sc <- sc + 30
      # Slight preference for keeping Raw if Raw already passes everything
      if (i == 1 && isTRUE(ev$both_ok) && (is.na(ev$grp_min_p) || ev$grp_min_p >= alpha)) sc <- sc + 200
      sc
    })
    best_idx <- which.max(scores)
  }

  best_item <- eval_results[[best_idx]]
  best_code <- best_item$cand$code
  rescued   <- (best_idx > 1) && (isTRUE(best_item$ev$both_ok) ||
               (!isTRUE(raw_ev$var_ok) && isTRUE(best_item$ev$var_ok)) ||
               (best_item$ev$sw_w - raw_ev$sw_w > 0.025))

  if (verbose) {
    cat("\n======================================================================\n")
    cat("       DATA TRANSFORMATION & ASSUMPTION RESCUE DIAGNOSTICS            \n")
    cat("======================================================================\n")
    print(comp_df, row.names = FALSE)

    if (best_idx == 1 && !raw_needs_help) {
      cat("\n[VERDICT]: Raw (untransformed) data already meets parametric assumptions! No transformation needed.\n")
    } else if (rescued) {
      cat(sprintf("\n[TRANSFORMATION RESCUE SUCCESS]: '%s' (`transform = '%s'`)\n",
                  best_item$cand$name, best_code))
      cat(sprintf("  * Residual Normality improved from W = %.4f (p = %.4g) -> W = %.4f (p = %.4g)\n",
                  raw_ev$sw_w, raw_ev$sw_p, best_item$ev$sw_w, best_item$ev$sw_p))
      if (!is.na(raw_ev$flig_p)) {
        cat(sprintf("  * Variance Homogeneity (Fligner p): %.4g -> %.4g\n",
                    raw_ev$flig_p, best_item$ev$flig_p))
      }
      cat(sprintf("  * R Code to create transformed column: data$%s <- %s\n",
                  make.names(best_item$cand$label), best_item$cand$r_expr))
    } else {
      cat("\n[VERDICT]: Standard transformations (log10, ln, sqrt) improved skewness slightly\n")
      cat(sprintf("           (Best = '%s', W = %.4f, p = %.4g), but residuals remain non-normal (p < %.2f).\n",
                  best_item$cand$name, best_item$ev$sw_w, best_item$ev$sw_p, alpha))
      cat("           -> Jumping to Non-Parametric test (or reporting both Transformed & Non-Parametric) is justified!\n")
    }
    cat("======================================================================\n")
  }

  # Build transformed dataframe for immediate downstream use
  df_trans <- df_use
  trans_y_col <- response
  trans_x_col <- if (!is.null(predictor)) predictors[1] else NULL

  if (best_idx > 1) {
    df_trans[[response]] <- best_item$y_tr
    if (!is.null(best_item$x_tr)) {
      df_trans[[predictors[1]]] <- best_item$x_tr
    }
  }

  # Plot Before vs After Transformation Diagnostics (Q-Q plots & Histograms)
  if (plot && best_idx > 1) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

    # 1. Raw Histogram
    hist(na.omit(raw_ev$resids), freq = FALSE, col = "mistyrose", border = "white",
         main = sprintf("1. Raw Residuals / Data\n(Shapiro W = %.3f, p = %.3g)", raw_ev$sw_w, raw_ev$sw_p),
         xlab = paste("Raw", response, "Residuals"))
    lines(density(na.omit(raw_ev$resids)), col = "firebrick", lwd = 2)

    # 2. Raw Normal Q-Q Plot
    qqnorm(na.omit(raw_ev$resids), pch = 19, col = "firebrick",
           main = "2. Raw Normal Q-Q Plot")
    qqline(na.omit(raw_ev$resids), col = "black", lwd = 2, lty = 2)

    # 3. Transformed Histogram
    hist(na.omit(best_item$ev$resids), freq = FALSE, col = "lightblue", border = "white",
         main = sprintf("3. Transformed: %s\n(Shapiro W = %.3f, p = %.3g)",
                        best_item$cand$name, best_item$ev$sw_w, best_item$ev$sw_p),
         xlab = paste(best_item$cand$label, "Residuals"))
    lines(density(na.omit(best_item$ev$resids)), col = "darkblue", lwd = 2)

    # 4. Transformed Normal Q-Q Plot
    qqnorm(na.omit(best_item$ev$resids), pch = 19, col = "steelblue4",
           main = paste("4. Transformed Q-Q Plot:", best_item$cand$code))
    qqline(na.omit(best_item$ev$resids), col = "black", lwd = 2, lty = 2)
  }

  invisible(list(
    table = comp_df,
    best_code = if (rescued || best_idx > 1) best_code else "none",
    best_name = best_item$cand$name,
    best_label = best_item$cand$label,
    best_x_label = best_item$cand$x_label,
    r_expr = best_item$cand$r_expr,
    rescued = rescued,
    raw_ok = !raw_needs_help,
    best_ev = best_item$ev,
    raw_ev = raw_ev,
    df_transformed = df_trans
  ))
}

# Helper to apply an explicit transformation code ("log10", "ln", "sqrt", "asin", "log_xy") to a dataframe
.apply_transform_df <- function(df, response, predictor = NULL, transform_code = "none") {
  if (is.null(transform_code) || transform_code %in% c("none", "raw", "FALSE")) {
    return(list(df = df, y_label = response, x_label = predictor, r_expr = NULL))
  }
  df_out <- df
  y <- df_out[[response]]
  min_y <- min(y, na.rm = TRUE)

  if (transform_code == "log10") {
    df_out[[response]] <- if (min_y > 0) log10(y) else log10(y + 1)
    y_lbl <- if (min_y > 0) paste0("log10(", response, ")") else paste0("log10(", response, "+1)")
    r_ex  <- if (min_y > 0) sprintf("log10(data$%s)", .sym_q(response)) else sprintf("log10(data$%s + 1)", .sym_q(response))
    return(list(df = df_out, y_label = y_lbl, x_label = predictor, r_expr = r_ex))
  } else if (transform_code %in% c("ln", "log")) {
    df_out[[response]] <- if (min_y > 0) log(y) else log(y + 1)
    y_lbl <- if (min_y > 0) paste0("ln(", response, ")") else paste0("ln(", response, "+1)")
    r_ex  <- if (min_y > 0) sprintf("log(data$%s)", .sym_q(response)) else sprintf("log(data$%s + 1)", .sym_q(response))
    return(list(df = df_out, y_label = y_lbl, x_label = predictor, r_expr = r_ex))
  } else if (transform_code == "sqrt") {
    df_out[[response]] <- sqrt(pmax(0, y))
    y_lbl <- paste0("sqrt(", response, ")")
    r_ex  <- sprintf("sqrt(data$%s)", .sym_q(response))
    return(list(df = df_out, y_label = y_lbl, x_label = predictor, r_expr = r_ex))
  } else if (transform_code == "asin") {
    max_y <- max(y, na.rm = TRUE)
    df_out[[response]] <- if (max_y <= 1) asin(sqrt(pmax(0, y))) else asin(sqrt(pmax(0, y / 100)))
    y_lbl <- paste0("asin(sqrt(", response, "))")
    r_ex  <- sprintf("asin(sqrt(data$%s))", .sym_q(response))
    return(list(df = df_out, y_label = y_lbl, x_label = predictor, r_expr = r_ex))
  } else if (transform_code %in% c("log_xy", "ln_xy") && !is.null(predictor)) {
    use_ln <- (transform_code == "ln_xy")
    fn_use <- if (use_ln) log else log10
    fn_nm  <- if (use_ln) "log" else "log10"
    df_out[[response]] <- if (min_y > 0) fn_use(y) else fn_use(y + 1)
    x <- df_out[[predictor[1]]]
    min_x <- min(x, na.rm = TRUE)
    df_out[[predictor[1]]] <- if (min_x > 0) fn_use(x) else fn_use(x + 1)
    y_lbl <- paste0(fn_nm, "(", response, ")")
    x_lbl <- paste0(fn_nm, "(", predictor[1], ")")
    r_ex  <- sprintf("data$%s <- %s(data$%s); data$%s <- %s(data$%s)",
                     .sym_q(response), fn_nm, .sym_q(response), .sym_q(predictor[1]), fn_nm, .sym_q(predictor[1]))
    return(list(df = df_out, y_label = y_lbl, x_label = x_lbl, r_expr = r_ex))
  } else if (FALSE) {
    df_out[[response]] <- if (min_y > 0) log10(y) else log10(y + 1)
    x <- df_out[[predictor[1]]]
    min_x <- min(x, na.rm = TRUE)
    df_out[[predictor[1]]] <- if (min_x > 0) log10(x) else log10(x + 1)
    y_lbl <- paste0("log10(", response, ")")
    x_lbl <- paste0("log10(", predictor[1], ")")
    r_ex  <- sprintf("data$%s <- log10(data$%s); data$%s <- log10(data$%s)",
                     .sym_q(response), .sym_q(response), .sym_q(predictor[1]), .sym_q(predictor[1]))
    return(list(df = df_out, y_label = y_lbl, x_label = x_lbl, r_expr = r_ex))
  }
  list(df = df_out, y_label = response, x_label = predictor, r_expr = NULL)
}

# ------------------------------------------------------------------------------
# 3. SELF-CHECKING ENGINE 1: COPY-PASTEABLE EQUIVALENT R CODE (`biostat_code`)
# ------------------------------------------------------------------------------
biostat_code <- function(df = NULL, response = NULL, predictor = NULL, paired = NULL,
                         model = "auto", random = NULL, transform = "auto", df_name = NULL, silent = FALSE) {
  if (is.null(df)) {
    if (exists("last_code", envir = .biostat_last_audit)) {
      cat(.biostat_last_audit$last_code)
      return(invisible(.biostat_last_audit$last_code))
    } else {
      stop("Pass your dataset `biostat_code(my_data)` (or run `biostat_wizard(my_data)` first).")
    }
  }

  # Auto-capture the user's REAL dataframe name (e.g., `pen_data` instead of generic `data`)!
  if (is.null(df_name)) {
    df_name <- paste(deparse(substitute(df)), collapse = "")
  }
  if (!nzchar(df_name) || df_name %in% c("df", "df_use", "parsed$df", "sub_df") || grepl("\\(|\\$|\\[", df_name)) {
    df_name <- "data"
  }

  inf <- biostat_infer(df, response = response, predictors = predictor,
                       paired = paired, model = model, random = random, verbose_roles = FALSE)
  prob_type <- inf$prob_type
  parsed    <- inf$parsed
  y_col     <- parsed$response
  x_cols    <- parsed$predictor

  tr_info <- tryCatch(biostat_transform(df, response = y_col, predictor = x_cols,
                                        paired = parsed$paired, plot = FALSE, verbose = FALSE),
                      error = function(e) NULL)

  lines <- c(
    "======================================================================",
    sprintf("   EQUIVALENT STANDALONE R CODE FOR `%s` (COPY & PASTE TO VERIFY!)   ", df_name),
    "======================================================================"
  )

  # Optional Data Transformation Rescue Block
  if (!is.null(tr_info) && !isTRUE(tr_info$raw_ok) && tr_info$best_code != "none") {
    lines <- c(lines,
      sprintf("# 0. Optional Data Transformation Rescue (Fixes Skewness / Unequal Variance!)"),
      sprintf("#    Best detected transform: %s (Residual Shapiro W = %.4f, p = %.4g)",
              tr_info$best_name, tr_info$best_ev$sw_w, tr_info$best_ev$sw_p),
      sprintf("#    To run the analysis below on transformed data, uncomment the next line:"),
      if (tr_info$best_code == "log_xy") {
        sprintf("# %s$%s <- log10(%s$%s); %s$%s <- log10(%s$%s)",
                df_name, .sym_q(y_col), df_name, .sym_q(y_col),
                df_name, .sym_q(x_cols[1]), df_name, .sym_q(x_cols[1]))
      } else {
        sprintf("# %s$%s <- %s", df_name, .sym_q(y_col),
                gsub("^`?([^`]+)`?", paste0(df_name, "$\\1"), tr_info$r_expr))
      },
      ""
    )
  }

  # ----------------------------------------------------------------------------
  # CASE 1: BINARY LOGISTIC REGRESSION
  # ----------------------------------------------------------------------------
  if (prob_type == "logistic") {
    x_col <- x_cols[1]
    yq <- .sym_q(y_col); xq <- .sym_q(x_col)
    lines <- c(lines,
      sprintf("# 1. Fit Binary Logistic Regression GLM (Binomial family, logit link)"),
      sprintf("mod <- glm(%s ~ %s, data = %s, family = binomial)", yq, xq, df_name),
      sprintf("summary(mod)                                  # Wald z-value and p-value (Pr(>|z|))"),
      sprintf(""),
      sprintf("# 2. Numeric Factor by which +1 unit of X multiplies Odds of Success"),
      sprintf("exp(coef(mod))                                # exp(b1) = Odds Ratio multiplier"),
      sprintf("exp(confint(mod))                             # 95%% Confidence Interval for Odds Ratio"),
      sprintf(""),
      sprintf("# 3. Strength of Association (McFadden's Pseudo-R2 & LRT Chi-Sq)"),
      sprintf("1 - (mod$deviance / mod$null.deviance)        # McFadden's Pseudo-R2 (prop. deviance explained)"),
      sprintf("anova(mod, test = 'Chisq')                    # Likelihood Ratio Test (Drop-in-Deviance G2)"),
      sprintf(""),
      sprintf("# 4. PRIMARY FINAL GRAPH: Sigmoidal Logistic Probability Curve (With Jitter)"),
      sprintf("#    WHY THIS GRAPH: Shows how the predicted probability smoothly transitions from 0 to 1."),
      sprintf("#    NOTE ON JITTER: `height = 0.03` nudges dots slightly above 1 / below 0 ONLY so"),
      sprintf("#                    overlapping 0s and 1s don't hide; actual data & GLM are strictly 0 and 1!"),
      sprintf("library(ggplot2)"),
      sprintf("ggplot(%s, aes(x = %s, y = %s)) +", df_name, xq, yq),
      sprintf("  geom_jitter(height = 0.03, width = 0, alpha = 0.65, size = 2.5, color = 'steelblue4') +"),
      sprintf("  stat_smooth(method = 'glm', method.args = list(family = 'binomial'),"),
      sprintf("              se = TRUE, color = 'firebrick', fill = 'pink', alpha = 0.3, linewidth = 1.1) +"),
      sprintf("  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(-0.05, 1.05)) +"),
      sprintf("  labs(title = 'Logistic Regression: %s vs. %s',", y_col, x_col),
      sprintf("       x = '%s', y = 'Probability of %s (1 = Success)') +", x_col, y_col),
      sprintf("  theme_minimal(base_size = 13)"),
      sprintf(""),
      sprintf("# 5. ALTERNATIVE CLEAN GRAPH: Zero Vertical Jitter (Dots Sit Strictly Flat on 0 and 1)"),
      sprintf("#    WHY: Uses semi-transparency (`alpha = 0.45`) instead of jitter so stacked points look darker."),
      sprintf("ggplot(%s, aes(x = %s, y = %s)) +", df_name, xq, yq),
      sprintf("  geom_point(alpha = 0.45, size = 2.8, color = 'steelblue4') +"),
      sprintf("  stat_smooth(method = 'glm', method.args = list(family = 'binomial'),"),
      sprintf("              se = TRUE, color = 'firebrick', fill = 'pink', alpha = 0.3) +"),
      sprintf("  labs(title = 'Logistic Regression (No Jitter): %s vs. %s',", y_col, x_col),
      sprintf("       x = '%s', y = 'Probability of %s') +", x_col, y_col),
      sprintf("  theme_classic(base_size = 13)")
    )

  # ----------------------------------------------------------------------------
  # CASE 2: TWO-WAY MODEL I FACTORIAL ANOVA
  # ----------------------------------------------------------------------------
  } else if (prob_type == "factorial_anova") {
    f1 <- .sym_q(x_cols[1]); f2 <- .sym_q(x_cols[2]); yq <- .sym_q(y_col)
    lines <- c(lines,
      sprintf("# 1. Fit Full Two-Way Factorial Model I ANOVA (WITH Interaction Term First!)"),
      sprintf("mod_int <- lm(%s ~ %s * %s, data = %s)", yq, f1, f2, df_name),
      sprintf("shapiro.test(residuals(mod_int))              # Formal normality test on residuals"),
      sprintf("anova(mod_int)                                # Read '%s:%s' interaction row FIRST!", x_cols[1], x_cols[2]),
      sprintf(""),
      sprintf("# 2. ASSUMPTION GRAPHS: Normal Q-Q Plot & Residuals vs. Fitted (Wedge-Shape Check)"),
      sprintf("#    WHY: Confirms residuals are normally distributed (dots on line) & equal variance (no funnel)."),
      sprintf("library(car)"),
      sprintf("par(mfrow = c(1, 3))"),
      sprintf("res <- residuals(mod_int); hist(res, freq = FALSE, col = 'lightblue', border = 'white', main = 'Histogram vs Normal'); curve(dnorm(x, mean(res), sd(res)), add = TRUE, col = 'firebrick', lwd = 2)"),
      sprintf("qqPlot(res, col.lines = 'firebrick', main = 'QQ Plot w/ 95%% Envelope') # Dots inside shaded gray/blue band = Normal!"),
      sprintf("plot(res ~ fitted(mod_int), pch = 19, xlab = 'Fitted Values', ylab = 'Residuals', main = 'Residuals vs Fitted'); abline(h = 0, lty = 2, col = 'red')"),
      sprintf("par(mfrow = c(1, 1))"),
      sprintf(""),
      sprintf("# 3A. IF INTERACTION IS NOT SIGNIFICANT (p >= 0.05): Drop interaction & test Main Effects"),
      sprintf("mod_add <- lm(%s ~ %s + %s, data = %s)        # Additive model (holds other factor constant)", yq, f1, f2, df_name),
      sprintf("anova(mod_add)                                # Report main effect F and p-values from here!"),
      sprintf(""),
      sprintf("# 3B. IF INTERACTION IS SIGNIFICANT (p < 0.05): Do NOT read main effects above; run t-tests"),
      sprintf("bartlett.test(%s ~ %s, data = %s)             # Check equal variance for %s", yq, f1, df_name, x_cols[1]),
      sprintf("t.test(%s ~ %s, data = %s)                    # Separate two-sample t-test for %s", yq, f1, df_name, x_cols[1]),
      sprintf("bartlett.test(%s ~ %s, data = %s)             # Check equal variance for %s", yq, f2, df_name, x_cols[2]),
      sprintf("t.test(%s ~ %s, data = %s)                    # Separate two-sample t-test for %s", yq, f2, df_name, x_cols[2]),
      sprintf(""),
      sprintf("# 4. PRIMARY FINAL GRAPH: Grouped Barplot of Means +/- SEM (CSUN Course Standard)"),
      sprintf("#    WHY: Standard biology format to compare group averages and precision (SEM) across 2 factors."),
      sprintf("library(ggplot2)"),
      sprintf("ggplot(%s, aes(x = factor(%s), y = %s, fill = factor(%s))) +", df_name, f1, yq, f2),
      sprintf("  stat_summary(fun = mean, geom = 'bar', position = position_dodge(0.75), width = 0.68, color = 'black', alpha = 0.85) +"),
      sprintf("  stat_summary(fun.data = mean_se, geom = 'errorbar', position = position_dodge(0.75), width = 0.22, linewidth = 0.75) +"),
      sprintf("  labs(title = '%s by %s and %s (Mean +/- SEM)', x = '%s', y = '%s (Mean +/- SEM)', fill = '%s') +",
              y_col, x_cols[1], x_cols[2], x_cols[1], y_col, x_cols[2]),
      sprintf("  theme_classic(base_size = 13)"),
      sprintf(""),
      sprintf("# 5. ALTERNATIVE GRAPH: Interaction Profile Plot (Best for Visually Spotting Interactions!)"),
      sprintf("#    WHY: Parallel lines = NO interaction; Crossing or diverging lines = SIGNIFICANT interaction!"),
      sprintf("ggplot(%s, aes(x = factor(%s), y = %s, color = factor(%s), group = factor(%s))) +", df_name, f1, yq, f2, f2),
      sprintf("  stat_summary(fun = mean, geom = 'point', size = 3.5) +"),
      sprintf("  stat_summary(fun = mean, geom = 'line', linewidth = 1.1) +"),
      sprintf("  stat_summary(fun.data = mean_se, geom = 'errorbar', width = 0.12) +"),
      sprintf("  labs(title = 'Interaction Plot: %s by %s & %s', x = '%s', y = '%s (Mean +/- SEM)', color = '%s') +",
              y_col, x_cols[1], x_cols[2], x_cols[1], y_col, x_cols[2]),
      sprintf("  theme_classic(base_size = 13)")
    )

  # ----------------------------------------------------------------------------
  # CASE 3: MODEL III MIXED-EFFECTS FACTORIAL ANOVA
  # ----------------------------------------------------------------------------
  } else if (prob_type == "model3_anova") {
    r_col <- if (!is.null(parsed$random)) parsed$random else x_cols[2]
    f_col <- setdiff(x_cols, r_col)[1]
    fq <- .sym_q(f_col); rq <- .sym_q(r_col); yq <- .sym_q(y_col)
    lines <- c(lines,
      sprintf("# 1. Load lme4 and Fit Full Mixed Model (WITH Random Interaction) vs Reduced Model"),
      sprintf("library(lme4)"),
      sprintf("mod_full   <- lmer(%s ~ %s + (1|%s) + (1|%s:%s), data = %s, REML = FALSE)", yq, fq, rq, fq, rq, df_name),
      sprintf("mod_no_int <- lmer(%s ~ %s + (1|%s), data = %s, REML = FALSE)", yq, fq, rq, df_name),
      sprintf(""),
      sprintf("# 2. Log-Likelihood Ratio Test (LLRT Chi-Square) for the Random Interaction Term"),
      sprintf("anova(mod_no_int, mod_full)                   # Check Chisq and Pr(>Chisq) for interaction!"),
      sprintf(""),
      sprintf("# 3. If Interaction is Significant: Test Fixed Effect in Model w/ Random Effect & NO Interaction"),
      sprintf("mod_null   <- lmer(%s ~ 1 + (1|%s), data = %s, REML = FALSE)", yq, rq, df_name),
      sprintf("anova(mod_null, mod_no_int)                   # LLRT for Fixed Effect ('%s')", f_col),
      sprintf("anova(lmer(%s ~ %s + (1|%s), data = %s))      # F-table for Fixed Effect ('%s')", yq, fq, rq, df_name, f_col),
      sprintf(""),
      sprintf("# 4. PRIMARY FINAL GRAPH: Grouped Barplot of Means +/- SEM (Fixed Factor on X, Random as Fill)"),
      sprintf("#    WHY: Illustrates how the Fixed effect varies across levels of the Random factor."),
      sprintf("library(ggplot2)"),
      sprintf("ggplot(%s, aes(x = factor(%s), y = %s, fill = factor(%s))) +", df_name, fq, yq, rq),
      sprintf("  stat_summary(fun = mean, geom = 'bar', position = position_dodge(0.75), width = 0.68, color = 'black', alpha = 0.85) +"),
      sprintf("  stat_summary(fun.data = mean_se, geom = 'errorbar', position = position_dodge(0.75), width = 0.22, linewidth = 0.75) +"),
      sprintf("  labs(title = 'Model III ANOVA: %s by %s & %s (Mean +/- SEM)', x = '%s', y = '%s (Mean +/- SEM)', fill = '%s') +",
              y_col, f_col, r_col, f_col, y_col, r_col),
      sprintf("  theme_classic(base_size = 13)")
    )

  # ----------------------------------------------------------------------------
  # CASE 4: ANCOVA
  # ----------------------------------------------------------------------------
  } else if (prob_type == "ancova") {
    cov_col <- x_cols[1]; grp_col <- x_cols[2]
    cq <- .sym_q(cov_col); gq <- .sym_q(grp_col); yq <- .sym_q(y_col)
    lvls <- levels(as.factor(parsed$df[[grp_col]]))
    lines <- c(lines,
      sprintf("# 1. Fit Full ANCOVA Model with Interaction (Tests Homogeneity of Slopes)"),
      sprintf("mod_int <- lm(%s ~ %s * %s, data = %s)", yq, cq, gq, df_name),
      sprintf("shapiro.test(residuals(mod_int))              # Check normality of residuals"),
      sprintf("anova(mod_int)                                # Read '%s:%s' interaction row FIRST!", cov_col, grp_col),
      sprintf(""),
      sprintf("# 2. ASSUMPTION GRAPHS: Normal Q-Q Plot & Residuals vs. Fitted (Check for Wedge Shape)"),
      sprintf("library(car)"),
      sprintf("par(mfrow = c(1, 3))"),
      sprintf("res <- residuals(mod_int); hist(res, freq = FALSE, col = 'lightblue', border = 'white', main = 'Histogram vs Normal'); curve(dnorm(x, mean(res), sd(res)), add = TRUE, col = 'firebrick', lwd = 2)"),
      sprintf("qqPlot(res, col.lines = 'firebrick', main = 'QQ Plot w/ 95%% Envelope') # Dots inside shaded gray/blue band = Normal!"),
      sprintf("plot(res ~ fitted(mod_int), pch = 19, xlab = 'Fitted Values', ylab = 'Residuals', main = 'Residuals vs Fitted'); abline(h = 0, lty = 2, col = 'red')"),
      sprintf("par(mfrow = c(1, 1))"),
      sprintf(""),
      sprintf("# 3A. IF INTERACTION IS SIGNIFICANT (p < 0.05): Stop & run separate regressions per group!"),
      sprintf("anova(lm(%s ~ %s, data = subset(%s, %s == '%s')))  # Regression for %s", yq, cq, df_name, gq, lvls[1], lvls[1]),
      sprintf("anova(lm(%s ~ %s, data = subset(%s, %s == '%s')))  # Regression for %s", yq, cq, df_name, gq, lvls[min(2, length(lvls))], lvls[min(2, length(lvls))]),
      sprintf(""),
      sprintf("# 3B. IF INTERACTION IS NOT SIGNIFICANT (p >= 0.05): Drop interaction & compare line heights!"),
      sprintf("mod_add <- lm(%s ~ %s + %s, data = %s)        # Parallel lines ANCOVA", yq, cq, gq, df_name),
      sprintf("anova(mod_add)                                # Row 1 = Slope (%s), Row 2 = Line Height (%s)", cov_col, grp_col),
      sprintf(""),
      sprintf("# 4. PRIMARY FINAL GRAPH: Scatterplot with Separate Trendlines by Group"),
      sprintf("#    WHY: Directly visualizes whether the two groups have different slopes (interaction) or"),
      sprintf("#         parallel slopes with different heights (main effect of grouping variable)."),
      sprintf("library(ggplot2)"),
      sprintf("ggplot(%s, aes(x = %s, y = %s, color = factor(%s))) +", df_name, cq, yq, gq),
      sprintf("  geom_point(size = 2.6, alpha = 0.8) +       # Raw observations colored by group (no jitter needed)"),
      sprintf("  geom_smooth(method = 'lm', se = FALSE, linewidth = 1.1) + # Separate OLS regression line per group"),
      sprintf("  labs(title = 'ANCOVA: %s vs. %s by %s', x = '%s', y = '%s', color = '%s') +",
              y_col, cov_col, grp_col, cov_col, y_col, grp_col),
      sprintf("  theme_classic(base_size = 13)")
    )

  # ----------------------------------------------------------------------------
  # CASE 5: MULTIPLE LINEAR REGRESSION & AIC
  # ----------------------------------------------------------------------------
  } else if (prob_type == "multiple_regression") {
    yq <- .sym_q(y_col)
    x_str <- paste(sapply(x_cols, .sym_q), collapse = " + ")
    lines <- c(lines,
      sprintf("# 1. Check Predictor Collinearity (Correlation Matrix & Scatterplot Matrix / SPLOM)"),
      sprintf("#    WHY THIS GRAPH: A SPLOM (`pairs`) reveals if predictors are strongly correlated (|r| > 0.7)"),
      sprintf("#                    with each other, which inflates variance and masks true effects."),
      sprintf("cor(%s[, c(%s)], use = 'complete.obs')", df_name, paste0("'", x_cols, "'", collapse = ", ")),
      sprintf("pairs(%s[, c('%s', %s)], main = 'Scatterplot Matrix (SPLOM)')",
              df_name, y_col, paste0("'", x_cols, "'", collapse = ", ")),
      sprintf(""),
      sprintf("# 2. Fit Full Additive Multiple Regression Model (No Interaction Terms)"),
      sprintf("mod_full <- lm(%s ~ %s, data = %s)", yq, x_str, df_name),
      sprintf("shapiro.test(residuals(mod_full))             # Check residual normality"),
      sprintf("anova(mod_full)                               # Sequential ANOVA table"),
      sprintf("summary(mod_full)                             # Partial slopes, t-tests, R2, and Adjusted R2"),
      sprintf(""),
      sprintf("# 3. ASSUMPTION GRAPHS: Normal Q-Q Plot & Residuals vs. Fitted Values"),
      sprintf("library(car)"),
      sprintf("par(mfrow = c(1, 3))"),
      sprintf("res <- residuals(mod_full); hist(res, freq = FALSE, col = 'lightblue', border = 'white', main = 'Histogram vs Normal'); curve(dnorm(x, mean(res), sd(res)), add = TRUE, col = 'firebrick', lwd = 2)"),
      sprintf("qqPlot(res, col.lines = 'firebrick', main = 'QQ Plot w/ 95%% Envelope') # Dots inside shaded gray/blue band = Normal!"),
      sprintf("plot(res ~ fitted(mod_full), pch = 19, xlab = 'Fitted Values', ylab = 'Residuals', main = 'Residuals vs Fitted'); abline(h = 0, lty = 2, col = 'red')"),
      sprintf("par(mfrow = c(1, 1))"),
      sprintf(""),
      sprintf("# 4. Compare Candidate Reduced Models Using Akaike's Information Criterion (AIC)"),
      sprintf("AIC(mod_full)                                 # Smaller (lower) AIC = BETTER parsimonious fit!")
    )
    for (xc in x_cols) {
      lines <- c(lines, sprintf("AIC(lm(%s ~ %s, data = %s))                   # Single-predictor model: %s only", yq, .sym_q(xc), df_name, xc))
    }

  # ----------------------------------------------------------------------------
  # CASE 6: PAIRED TWO-SAMPLE COMPARISON
  # ----------------------------------------------------------------------------
  } else if (prob_type == "paired_2group") {
    if (parsed$wide_paired) {
      c1 <- .sym_q(parsed$orig_cols[1]); c2 <- .sym_q(parsed$orig_cols[2])
      lines <- c(lines,
        sprintf("# 1. Compute Paired Differences & Check Normality (Required for Paired Design!)"),
        sprintf("diffs <- %s$%s - %s$%s", df_name, c1, df_name, c2),
        sprintf("shapiro.test(diffs)                           # p >= 0.05 -> Paired t-test; p < 0.05 -> Try transform or Wilcoxon"),
        sprintf(""),
        sprintf("# 2. ASSUMPTION GRAPH: Histogram & Normal Q-Q Plot of Paired Differences (d)"),
        sprintf("#    WHY: Paired t-tests assume the DIFFERENCES (not raw groups) are normally distributed!"),
        sprintf("par(mfrow = c(1, 2))"),
        sprintf("hist(diffs, col = 'lightblue', border = 'white', main = 'Paired Differences (%s - %s)', xlab = 'Difference'); abline(v = 0, col = 'red', lty = 2, lwd = 2)", parsed$orig_cols[1], parsed$orig_cols[2]),
        sprintf("qqnorm(diffs, pch = 19, col = 'steelblue'); qqline(diffs, col = 'firebrick', lwd = 2)"),
        sprintf("par(mfrow = c(1, 1))"),
        sprintf(""),
        sprintf("# 3. Execute Parametric Paired t-test (if normal) or Wilcoxon Signed-Rank (if non-normal)"),
        sprintf("t.test(%s$%s, %s$%s, paired = TRUE)           # Parametric Paired t-test", df_name, c1, df_name, c2),
        sprintf("wilcox.test(%s$%s, %s$%s, paired = TRUE)      # Non-Parametric Wilcoxon Signed-Rank test", df_name, c1, df_name, c2),
        sprintf(""),
        sprintf("# 4. PRIMARY FINAL GRAPH: Connected-Subject Paired Slopegraph + Boxplot"),
        sprintf("#    WHY: Gray lines connect each individual across both conditions so you can visually"),
        sprintf("#         see whether almost every subject shifts in the same direction!"),
        sprintf("library(ggplot2)"),
        sprintf("df_long <- data.frame(ID = factor(rep(seq_len(nrow(%s)), 2)),", df_name),
        sprintf("                      Condition = factor(rep(c('%s', '%s'), each = nrow(%s)), levels = c('%s', '%s')),",
                parsed$orig_cols[1], parsed$orig_cols[2], df_name, parsed$orig_cols[1], parsed$orig_cols[2]),
        sprintf("                      Value = c(%s$%s, %s$%s))", df_name, c1, df_name, c2),
        sprintf("ggplot(df_long, aes(x = Condition, y = Value)) +"),
        sprintf("  geom_boxplot(aes(fill = Condition), alpha = 0.35, width = 0.4, outlier.shape = NA) +"),
        sprintf("  geom_line(aes(group = ID), color = 'gray45', alpha = 0.6, linewidth = 0.7) +"),
        sprintf("  geom_point(aes(color = Condition), size = 2.8, alpha = 0.9) +"),
        sprintf("  stat_summary(fun = mean, geom = 'point', shape = 23, size = 4, fill = 'white', color = 'black') + # White Diamond = Mean"),
        sprintf("  labs(title = 'Paired Comparison: %s vs. %s', y = 'Measurement') +", parsed$orig_cols[1], parsed$orig_cols[2]),
        sprintf("  theme_minimal(base_size = 13) + theme(legend.position = 'none')")
      )
    } else {
      yq <- .sym_q(y_col); xq <- .sym_q(x_cols[1])
      lines <- c(lines,
        sprintf("# 1. Execute Paired t-test / Wilcoxon Signed-Rank Test"),
        sprintf("t.test(%s ~ %s, data = %s, paired = TRUE)", yq, xq, df_name),
        sprintf("wilcox.test(%s ~ %s, data = %s, paired = TRUE)", yq, xq, df_name)
      )
    }

  # ----------------------------------------------------------------------------
  # CASE 7: INDEPENDENT TWO-SAMPLE COMPARISON
  # ----------------------------------------------------------------------------
  } else if (prob_type == "indep_2group") {
    yq <- .sym_q(y_col); xq <- .sym_q(x_cols[1])
    lines <- c(lines,
      sprintf("# 1. Check Normality & Homogeneity of Variance (Raw vs Transformed)"),
      sprintf("mod <- lm(%s ~ %s, data = %s)", yq, xq, df_name),
      sprintf("shapiro.test(residuals(mod))                  # Normality check on model residuals"),
      sprintf("bartlett.test(%s ~ %s, data = %s)             # Bartlett's test for equal variance (Professor's primary test; p >= 0.05 = equal var)", yq, xq, df_name),
      sprintf("fligner.test(%s ~ %s, data = %s)              # Fligner-Killeen test (robust backup if non-normal)", yq, xq, df_name),
      sprintf("library(moments); skewness(residuals(mod)); kurtosis(residuals(mod)) # Skew near 0 = symmetric; Kurt near 3 = normal tails"),
      sprintf("biostat_transform(%s)                         # Check if log10(Y) or sqrt(Y) fixes assumptions!", df_name),
      sprintf(""),
      sprintf("# 2. ASSUMPTION GRAPHS: Normal Q-Q Plot & Residuals vs. Fitted Values"),
      sprintf("#    WHY: Verifies normality (dots follow diagonal line) and equal spread across both groups."),
      sprintf("library(car)"),
      sprintf("par(mfrow = c(1, 3))"),
      sprintf("res <- residuals(mod); hist(res, freq = FALSE, col = 'lightblue', border = 'white', main = 'Histogram vs Normal'); curve(dnorm(x, mean(res), sd(res)), add = TRUE, col = 'firebrick', lwd = 2)"),
      sprintf("qqPlot(res, col.lines = 'firebrick', main = 'QQ Plot w/ 95%% Envelope') # Dots inside shaded gray/blue band = Normal!"),
      sprintf("plot(res ~ fitted(mod), pch = 19, xlab = 'Fitted Values', ylab = 'Residuals', main = 'Residuals vs Fitted'); abline(h = 0, lty = 2, col = 'red')"),
      sprintf("par(mfrow = c(1, 1))"),
      sprintf(""),
      sprintf("# 3. Run Two-Sample Test Based on Assumptions"),
      sprintf("t.test(%s ~ %s, data = %s, var.equal = FALSE) # Welch's t-test (unequal variances / safe default)", yq, xq, df_name),
      sprintf("t.test(%s ~ %s, data = %s, var.equal = TRUE)  # Student's t-test (equal variances)", yq, xq, df_name),
      sprintf("wilcox.test(%s ~ %s, data = %s)               # Mann-Whitney U test (if transformations fail)", yq, xq, df_name),
      sprintf(""),
      sprintf("# 4. PRIMARY FINAL GRAPH: Boxplot + Horizontal Jitter + Mean Diamond (Exact Wizard Plot!)"),
      sprintf("#    WHY THIS GRAPH: Shows full raw distribution, median/IQR, skewness, AND group mean all at once."),
      sprintf("#    NOTE ON JITTER: `width = 0.15` nudges dots left/right ONLY so identical values don't overlap;"),
      sprintf("#                    vertical y-position is exact! `outlier.shape = NA` avoids drawing outliers twice."),
      sprintf("library(ggplot2)"),
      sprintf("ggplot(%s, aes(x = factor(%s), y = %s, fill = factor(%s))) +", df_name, xq, yq, xq),
      sprintf("  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.5) +"),
      sprintf("  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +          # Horizontal jitter only"),
      sprintf("  stat_summary(fun = mean, geom = 'point', shape = 23, size = 3.5, fill = 'white') + # White Diamond = Mean"),
      sprintf("  labs(title = '%s by %s',", y_col, x_cols[1]),
      sprintf("       subtitle = 'Boxes = Median/IQR | White Diamond = Group Mean | Points = Raw Data',"),
      sprintf("       x = '%s', y = '%s') +", x_cols[1], y_col),
      sprintf("  theme_minimal(base_size = 13) + theme(legend.position = 'none')"),
      sprintf(""),
      sprintf("# 5. ALTERNATIVE GRAPH: Barplot of Group Means +/- SEM (Classic Exam / Publication Format)"),
      sprintf("ggplot(%s, aes(x = factor(%s), y = %s, fill = factor(%s))) +", df_name, xq, yq, xq),
      sprintf("  stat_summary(fun = mean, geom = 'bar', width = 0.55, color = 'black', alpha = 0.8) +"),
      sprintf("  stat_summary(fun.data = mean_se, geom = 'errorbar', width = 0.18, linewidth = 0.8) +"),
      sprintf("  labs(title = '%s by %s (Mean +/- SEM)', x = '%s', y = '%s (Mean +/- SEM)') +", y_col, x_cols[1], x_cols[1], y_col),
      sprintf("  theme_classic(base_size = 13) + theme(legend.position = 'none')")
    )

  # ----------------------------------------------------------------------------
  # CASE 8: INDEPENDENT >2 GROUPS (ONE-WAY ANOVA / WELCH / KRUSKAL-WALLIS)
  # ----------------------------------------------------------------------------
  } else if (prob_type == "indep_kgroup") {
    yq <- .sym_q(y_col); xq <- .sym_q(x_cols[1])
    lines <- c(lines,
      sprintf("# 1. Fit One-Way ANOVA & Check Residual Normality / Variance Homogeneity"),
      sprintf("fit <- aov(%s ~ %s, data = %s)", yq, xq, df_name),
      sprintf("shapiro.test(residuals(fit))                  # Normality of residuals"),
      sprintf("bartlett.test(%s ~ %s, data = %s)             # Bartlett's test for equal variance (Professor's primary test; p >= 0.05 = equal var)", yq, xq, df_name),
      sprintf("fligner.test(%s ~ %s, data = %s)              # Fligner-Killeen test (robust backup if non-normal)", yq, xq, df_name),
      sprintf("library(moments); skewness(residuals(fit)); kurtosis(residuals(fit)) # Skew near 0 = symmetric; Kurt near 3 = normal tails"),
      sprintf("biostat_transform(%s)                         # Check if log10(Y) or sqrt(Y) fixes assumptions!", df_name),
      sprintf(""),
      sprintf("# 2. ASSUMPTION GRAPHS: Normal Q-Q Plot & Residuals vs. Fitted Values"),
      sprintf("library(car)"),
      sprintf("par(mfrow = c(1, 3))"),
      sprintf("res <- residuals(fit); hist(res, freq = FALSE, col = 'lightblue', border = 'white', main = 'Histogram vs Normal'); curve(dnorm(x, mean(res), sd(res)), add = TRUE, col = 'firebrick', lwd = 2)"),
      sprintf("qqPlot(res, col.lines = 'firebrick', main = 'QQ Plot w/ 95%% Envelope') # Dots inside shaded gray/blue band = Normal!"),
      sprintf("plot(res ~ fitted(fit), pch = 19, xlab = 'Fitted Values', ylab = 'Residuals', main = 'Residuals vs Fitted'); abline(h = 0, lty = 2, col = 'red')"),
      sprintf("par(mfrow = c(1, 1))"),
      sprintf(""),
      sprintf("# 3. Run Omnibus Test + Post-Hoc Comparisons"),
      sprintf("summary(fit)                                  # Standard One-Way ANOVA (if normal & equal var)"),
      sprintf("TukeyHSD(fit)                                 # Tukey's HSD post-hoc pairwise comparisons"),
      sprintf("oneway.test(%s ~ %s, data = %s, var.equal = FALSE) # Welch's ANOVA (if unequal var)", yq, xq, df_name),
      sprintf("kruskal.test(%s ~ %s, data = %s)              # Non-parametric Kruskal-Wallis (if transformations fail)", yq, xq, df_name),
      sprintf(""),
      sprintf("# 4. PRIMARY FINAL GRAPH: Boxplot + Horizontal Jitter + Mean Diamond (Exact Wizard Plot!)"),
      sprintf("#    NOTE ON JITTER: `width = 0.15` spreads points horizontally so duplicates don't overlap; y-values are exact."),
      sprintf("library(ggplot2)"),
      sprintf("ggplot(%s, aes(x = factor(%s), y = %s, fill = factor(%s))) +", df_name, xq, yq, xq),
      sprintf("  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.5) +"),
      sprintf("  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +"),
      sprintf("  stat_summary(fun = mean, geom = 'point', shape = 23, size = 3.5, fill = 'white') + # White Diamond = Mean"),
      sprintf("  labs(title = '%s by %s', subtitle = 'Boxes = Median/IQR | White Diamond = Mean | Points = Raw Data',", y_col, x_cols[1]),
      sprintf("       x = '%s', y = '%s') +", x_cols[1], y_col),
      sprintf("  theme_minimal(base_size = 13) + theme(legend.position = 'none')"),
      sprintf(""),
      sprintf("# 5. ALTERNATIVE GRAPH: Barplot of Group Means +/- SEM"),
      sprintf("ggplot(%s, aes(x = factor(%s), y = %s, fill = factor(%s))) +", df_name, xq, yq, xq),
      sprintf("  stat_summary(fun = mean, geom = 'bar', width = 0.6, color = 'black', alpha = 0.8) +"),
      sprintf("  stat_summary(fun.data = mean_se, geom = 'errorbar', width = 0.2, linewidth = 0.8) +"),
      sprintf("  labs(title = '%s by %s (Mean +/- SEM)', x = '%s', y = '%s (Mean +/- SEM)') +", y_col, x_cols[1], x_cols[1], y_col),
      sprintf("  theme_classic(base_size = 13) + theme(legend.position = 'none')")
    )

  # ----------------------------------------------------------------------------
  # CASE 9: SIMPLE LINEAR REGRESSION & CORRELATION
  # ----------------------------------------------------------------------------
  } else if (prob_type == "continuous_xy") {
    yq <- .sym_q(y_col); xq <- .sym_q(x_cols[1])
    lines <- c(lines,
      sprintf("# 1. Fit Ordinary Least Squares (OLS) Linear Regression"),
      sprintf("mod <- lm(%s ~ %s, data = %s)", yq, xq, df_name),
      sprintf("shapiro.test(residuals(mod))                  # Check normality of regression residuals"),
      sprintf("biostat_transform(%s)                         # Check if log10(Y) or log10(X) linearizes & normalizes!", df_name),
      sprintf("anova(mod)                                    # Regression ANOVA table (F-ratio & Mean Sq)"),
      sprintf("summary(mod)                                  # Intercept, Slope, R-squared, and p-value"),
      sprintf(""),
      sprintf("# 2. ASSUMPTION GRAPHS: Normal Q-Q Plot & Residuals vs. Fitted (Check Linearity & Homoscedasticity)"),
      sprintf("library(car)"),
      sprintf("par(mfrow = c(1, 3))"),
      sprintf("res <- residuals(mod); hist(res, freq = FALSE, col = 'lightblue', border = 'white', main = 'Histogram vs Normal'); curve(dnorm(x, mean(res), sd(res)), add = TRUE, col = 'firebrick', lwd = 2)"),
      sprintf("qqPlot(res, col.lines = 'firebrick', main = 'QQ Plot w/ 95%% Envelope') # Dots inside shaded gray/blue band = Normal!"),
      sprintf("plot(res ~ fitted(mod), pch = 19, xlab = 'Fitted Values', ylab = 'Residuals', main = 'Residuals vs Fitted'); abline(h = 0, lty = 2, col = 'red')"),
      sprintf("par(mfrow = c(1, 1))"),
      sprintf(""),
      sprintf("# 3. Correlation Analysis (Parametric Pearson r vs Non-Parametric Spearman rho)"),
      sprintf("cor.test(%s$%s, %s$%s, method = 'pearson')", df_name, xq, df_name, yq),
      sprintf("cor.test(%s$%s, %s$%s, method = 'spearman')", df_name, xq, df_name, yq),
      sprintf(""),
      sprintf("# 4. PRIMARY FINAL GRAPH: Scatterplot + OLS Regression Line + 95%% Confidence Band"),
      sprintf("#    WHY THIS GRAPH: Visualizes direction, strength, linearity, and uncertainty of the X -> Y relationship."),
      sprintf("library(ggplot2)"),
      sprintf("ggplot(%s, aes(x = %s, y = %s)) +", df_name, xq, yq),
      sprintf("  geom_point(alpha = 0.75, size = 2.5, color = 'steelblue') +"),
      sprintf("  geom_smooth(method = 'lm', se = TRUE, color = 'firebrick', fill = 'pink', alpha = 0.3) +"),
      sprintf("  labs(title = 'Linear Regression: %s vs. %s',", y_col, x_cols[1]),
      sprintf("       subtitle = 'Red line = OLS fit with 95%% confidence band', x = '%s', y = '%s') +", x_cols[1], y_col),
      sprintf("  theme_minimal(base_size = 13)")
    )

  # ----------------------------------------------------------------------------
  # CASE 10: CONTINGENCY TABLE (CHI-SQUARE / FISHER / G-TEST)
  # ----------------------------------------------------------------------------
  } else if (prob_type %in% c("contingency", "contingency_table")) {
    lines <- c(lines,
      sprintf("# 1. Build Contingency Table & Check Expected Counts"),
      if (prob_type == "contingency_table") {
        sprintf("tbl <- as.table(%s)", df_name)
      } else {
        sprintf("tbl <- table(%s$%s, %s$%s)", df_name, .sym_q(y_col), df_name, .sym_q(x_cols[1]))
      },
      sprintf("chisq.test(tbl, correct = FALSE)$expected     # All expected counts >= 5?"),
      sprintf("chisq.test(tbl, correct = FALSE)              # Pearson's Chi-Square Test (if all E >= 5)"),
      sprintf("fisher.test(tbl)                              # Fisher's Exact Test (if any E < 5)"),
      sprintf(""),
      sprintf("# 2. PRIMARY FINAL GRAPH: Mosaic Plot (Visualizes Categorical Proportions & Independence)"),
      sprintf("#    WHY THIS GRAPH: Tile areas are proportional to counts; misaligned splits reveal association."),
      sprintf("mosaicplot(tbl, main = 'Contingency Mosaic Plot', color = TRUE)")
    )

  # ----------------------------------------------------------------------------
  # CASE 11: ONE-SAMPLE CONTINUOUS
  # ----------------------------------------------------------------------------
  } else if (prob_type == "one_sample_continuous") {
    yq <- .sym_q(y_col)
    lines <- c(lines,
      sprintf("# 1. Check Normality & Run One-Sample Test"),
      sprintf("shapiro.test(%s$%s)                           # p >= 0.05 -> t.test; p < 0.05 -> wilcox.test", df_name, yq),
      sprintf("t.test(%s$%s, mu = 0)                         # Replace mu = 0 with your null hypothesis value", df_name, yq),
      sprintf("wilcox.test(%s$%s, mu = 0)                    # Non-parametric alternative", df_name, yq),
      sprintf(""),
      sprintf("# 2. ASSUMPTION & DISTRIBUTION GRAPHS: Histogram with Density Curve + Normal Q-Q Plot"),
      sprintf("par(mfrow = c(1, 2))"),
      sprintf("hist(%s$%s, freq = FALSE, col = 'lightblue', border = 'white', main = 'Distribution of %s', xlab = '%s')", df_name, yq, y_col, y_col),
      sprintf("lines(density(na.omit(%s$%s)), col = 'darkblue', lwd = 2)", df_name, yq),
      sprintf("qqnorm(%s$%s, pch = 19, col = 'steelblue'); qqline(%s$%s, col = 'firebrick', lwd = 2)", df_name, yq, df_name, yq),
      sprintf("par(mfrow = c(1, 1))")
    )
  }

  lines <- c(lines, "======================================================================\n")
  out_str <- paste(lines, collapse = "\n")
  if (!isTRUE(silent)) cat("\n", out_str, sep = "")
  invisible(out_str)
}
# ------------------------------------------------------------------------------
# 4. SELF-CHECKING ENGINE 2 & 3: `biostat_explain` & `biostat_check_math`
# ------------------------------------------------------------------------------
biostat_explain <- function(df = NULL, response = NULL, predictor = NULL, paired = NULL,
                            model = "auto", random = NULL, transform = "auto") {
  if (is.null(df)) {
    if (exists("last_explain", envir = .biostat_last_audit)) {
      cat(.biostat_last_audit$last_explain)
      return(invisible(.biostat_last_audit$last_explain))
    } else {
      stop("Pass your dataset `biostat_explain(data)` (or run `biostat_wizard(data)` first).")
    }
  }
  biostat_wizard(df, response = response, predictor = predictor, paired = paired,
                 model = model, random = random, transform = transform, plot = FALSE, show_code = FALSE)
}

biostat_check_math <- function(df, response = NULL, predictor = NULL, paired = NULL) {
  inf <- biostat_infer(df, response = response, predictors = predictor, paired = paired, verbose_roles = FALSE)
  prob_type <- inf$prob_type
  parsed <- inf$parsed
  df_use <- parsed$df
  y_col <- parsed$response
  x_cols <- parsed$predictor

  cat("\n======================================================================\n")
  cat("       UNDER-THE-HOOD MATH VERIFICATION (PLUGGED-IN FORMULAS)         \n")
  cat("======================================================================\n")

  if (prob_type == "logistic") {
    sub_df <- df_use[complete.cases(df_use[, c(y_col, x_cols[1])]), ]
    y_raw <- sub_df[[y_col]]
    y_num <- if (is.numeric(y_raw) && all(y_raw %in% c(0, 1))) y_raw else as.numeric(as.factor(y_raw)) - 1
    sub_df$.y_num <- y_num
    fit <- glm(as.formula(paste(".y_num ~", .sym_q(x_cols[1]))), data = sub_df, family = binomial)
    b0 <- coef(fit)[1]; b1 <- coef(fit)[2]
    se_b1 <- summary(fit)$coefficients[2, "Std. Error"]
    d0 <- fit$null.deviance; d1 <- fit$deviance
    cat("1. Logistic Regression Log-Odds Equation:\n")
    cat(sprintf("   ln(p / (1 - p)) = b0 + b1 * (%s) = %.4f + %.4f * (%s)\n\n", x_cols[1], b0, b1, x_cols[1]))
    cat("2. Wald Z-Statistic Formula:\n")
    cat(sprintf("   z = b1 / SE(b1) = %.4f / %.4f = %.4f (p = %.5g)\n\n", b1, se_b1, b1/se_b1, summary(fit)$coefficients[2, 4]))
    cat("3. Numeric Odds Multiplier Formula:\n")
    cat(sprintf("   Odds Ratio = e^(b1) = exp(%.4f) = %.4f\n", b1, exp(b1)))
    cat(sprintf("   (Each +1 unit of '%s' multiplies the odds of '%s' by %.4f)\n\n", x_cols[1], y_col, exp(b1)))
    cat("4. McFadden's Pseudo-R2 & Likelihood Ratio Test (G2) Formula:\n")
    cat(sprintf("   LRT G2 = Null_Deviance - Residual_Deviance = %.4f - %.4f = %.4f (df = 1, p = %.5g)\n",
                d0, d1, d0 - d1, pchisq(d0 - d1, 1, lower.tail = FALSE)))
    cat(sprintf("   McFadden R2 = 1 - (Residual_Deviance / Null_Deviance) = 1 - (%.4f / %.4f) = %.4f\n",
                d1, d0, 1 - d1/d0))

  } else if (prob_type %in% c("factorial_anova", "model3_anova", "ancova")) {
    f1 <- .sym_q(x_cols[1]); f2 <- .sym_q(x_cols[2]); yq <- .sym_q(y_col)
    fit_int <- lm(as.formula(paste(yq, "~", f1, "*", f2)), data = df_use)
    at <- anova(fit_int)
    cat(sprintf("1. Linear Model Partitioning of Variance (%s ~ %s * %s):\n", y_col, x_cols[1], x_cols[2]))
    for (r in seq_len(nrow(at))) {
      rn <- trimws(rownames(at)[r])
      ss <- at[r, "Sum Sq"]; df_r <- at[r, "Df"]; ms <- at[r, "Mean Sq"]
      cat(sprintf("   * %-22s: SS = %10.4f, df = %3d -> MS = SS/df = %10.4f\n", rn, ss, df_r, ms))
    }
    ms_err <- at[nrow(at), "Mean Sq"]
    cat("\n2. F-Ratio Formulas (MS_Effect / MS_Residual):\n")
    for (r in seq_len(nrow(at) - 1)) {
      rn <- trimws(rownames(at)[r])
      ms <- at[r, "Mean Sq"]; f_v <- at[r, "F value"]; p_v <- at[r, "Pr(>F)"]
      cat(sprintf("   * F(%-18s) = %.4f / %.4f = %.4f (p = %.5g)\n", rn, ms, ms_err, f_v, p_v))
    }

  } else if (prob_type == "paired_2group") {
    x_fac <- as.factor(df_use[[x_cols[1]]])
    y <- df_use[[y_col]]
    g1 <- y[x_fac == levels(x_fac)[1]]; g2 <- y[x_fac == levels(x_fac)[2]]
    d <- g1 - g2; n <- length(d); md <- mean(d); sdd <- sd(d); sed <- sdd / sqrt(n)
    t_v <- md / sed
    cat("1. Paired Differences Formula (d_i = Y_1i - Y_2i):\n")
    cat(sprintf("   Mean(d) = %.4f, SD(d) = %.4f, n = %d pairs\n", md, sdd, n))
    cat(sprintf("   SE(d)   = SD(d) / sqrt(n) = %.4f / sqrt(%d) = %.4f\n\n", sdd, n, sed))
    cat("2. Paired t-Statistic Formula:\n")
    cat(sprintf("   t = Mean(d) / SE(d) = %.4f / %.4f = %.4f (df = n - 1 = %d, p = %.5g)\n",
                md, sed, t_v, n - 1, 2 * pt(-abs(t_v), df = n - 1)))

  } else if (prob_type == "indep_2group") {
    x_fac <- as.factor(df_use[[x_cols[1]]])
    y <- df_use[[y_col]]
    g1 <- na.omit(y[x_fac == levels(x_fac)[1]]); g2 <- na.omit(y[x_fac == levels(x_fac)[2]])
    n1 <- length(g1); n2 <- length(g2); m1 <- mean(g1); m2 <- mean(g2); s1 <- sd(g1); s2 <- sd(g2)
    sp <- sqrt(((n1 - 1)*s1^2 + (n2 - 1)*s2^2) / (n1 + n2 - 2))
    se_s <- sp * sqrt(1/n1 + 1/n2)
    se_w <- sqrt(s1^2/n1 + s2^2/n2)
    cat(sprintf("1. Group Summary Stats:\n   %s: Mean1 = %.4f, SD1 = %.4f, n1 = %d\n   %s: Mean2 = %.4f, SD2 = %.4f, n2 = %d\n\n",
                levels(x_fac)[1], m1, s1, n1, levels(x_fac)[2], m2, s2, n2))
    cat("2. Student's Equal-Variance t-Formula:\n")
    cat(sprintf("   Pooled SD (s_p) = %.4f -> SE_pooled = %.4f -> t = (%.4f - %.4f) / %.4f = %.4f (df = %d)\n\n",
                sp, se_s, m1, m2, se_s, (m1 - m2)/se_s, n1 + n2 - 2))
    cat("3. Welch's Unequal-Variance t-Formula:\n")
    cat(sprintf("   SE_Welch = sqrt(s1^2/n1 + s2^2/n2) = %.4f -> t = (%.4f - %.4f) / %.4f = %.4f\n",
                se_w, m1, m2, se_w, (m1 - m2)/se_w))
  } else {
    cat("Run `biostat_code(data)` to view the exact R verification commands for this model.\n")
  }
  cat("======================================================================\n")
}

# ------------------------------------------------------------------------------
# 5. DATASET INSPECTOR & AUTO-ROLE CLASSIFIER
# ------------------------------------------------------------------------------
biostat_inspect <- function(df) {
  if (is.table(df) || is.matrix(df)) {
    cat("\n--- CONTINGENCY TABLE / MATRIX DETECTED ---\n")
    print(addmargins(as.table(df)))
    return(invisible(df))
  }
  if (is.atomic(df) && is.null(dim(df))) df <- data.frame(Value = df)
  df <- as.data.frame(df)
  parsed <- .auto_parse_roles(df, verbose = FALSE)

  cat("\n--- DATASET OVERVIEW & AUTO-INFERRED ROLES ---\n")
  cat(sprintf("Dimensions: %d observations (rows) x %d variables (columns)\n", nrow(df), ncol(df)))
  if (parsed$wide_paired) {
    cat(sprintf("Detected Structure: WIDE PAIRED DESIGN ('%s' vs '%s', Paired by '%s')\n\n",
                parsed$orig_cols[1], parsed$orig_cols[2], parsed$block))
  } else {
    cat(sprintf("Inferred Response (Y): %s | Predictor(s) (X): %s\n\n",
                ifelse(is.null(parsed$response), "None", parsed$response),
                ifelse(is.null(parsed$predictor), "None", paste(parsed$predictor, collapse = ", "))))
  }

  summary_list <- lapply(names(df), function(col_name) {
    x <- df[[col_name]]
    n_miss <- sum(is.na(x))
    pct_miss <- round(100 * n_miss / length(x), 1)
    n_uniq <- length(unique(na.omit(x)))

    if (parsed$wide_paired) {
      role_tag <- if (!is.null(parsed$block) && col_name == parsed$block) "-> SUBJECT / PAIR ID"
                  else if (col_name %in% parsed$orig_cols) "-> PAIRED COND (Wide)"
                  else "-"
    } else {
      role_tag <- if (!is.null(parsed$response) && col_name == parsed$response) "-> RESPONSE (Y)"
                  else if (!is.null(parsed$predictor) && col_name %in% parsed$predictor) {
                    idx <- which(parsed$predictor == col_name)
                    sprintf("-> PREDICTOR (X%d)", idx)
                  } else if (col_name %in% parsed$rep_cols) "-> REPLICATE INDEX"
                  else if (!is.null(parsed$block) && col_name == parsed$block) "-> ID / BLOCK"
                  else "-"
    }

    if (is.numeric(x)) {
      v_type <- if (!is.null(parsed$block) && col_name == parsed$block) "Sequential ID"
                else if (col_name %in% parsed$rep_cols) "Replicate Index"
                else if (n_uniq == 2 && all(na.omit(x) %in% c(0, 1))) "Binary (0/1)"
                else if (n_uniq <= 4 && all(na.omit(x) == floor(na.omit(x)))) sprintf("Numeric Factor (%d lvl)", n_uniq)
                else "Continuous"
      x_clean <- na.omit(x)
      mean_sd <- if (length(x_clean) > 0) sprintf("%.2f (%.2f)", mean(x_clean), sd(x_clean)) else "NA"
      med_iqr <- if (length(x_clean) > 0) sprintf("%.2f (%.2f)", median(x_clean), IQR(x_clean)) else "NA"
    } else {
      v_type <- if (n_uniq == 2) "Categorical (2 lvl)" else if (n_uniq <= 15) sprintf("Categorical (%d lvl)", n_uniq) else "High-Card / ID"
      mean_sd <- "-"
      med_iqr <- "-"
    }

    data.frame(
      Variable = col_name, Auto_Role = role_tag, Scale_Type = v_type,
      Unique = n_uniq, Missing = sprintf("%d (%.0f%%)", n_miss, pct_miss),
      Mean_SD = mean_sd, Median_IQR = med_iqr, stringsAsFactors = FALSE
    )
  })

  res_df <- do.call(rbind, summary_list)
  print(res_df, row.names = FALSE)
  cat("\nTip: Run `biostat_wizard(data)` (or `biostat_transform(data)` to test log/sqrt transformations)!\n")
  invisible(res_df)
}

# ------------------------------------------------------------------------------
# 6. DEDICATED SOLVER 1: BINARY LOGISTIC REGRESSION (PS7 Q1: Crayfish Mating)
# ------------------------------------------------------------------------------
biostat_logistic <- function(df, response = NULL, predictor = NULL, plot = TRUE, show_code = FALSE) {
  parsed <- .auto_parse_roles(df, response = response, predictor = predictor, paired = FALSE, verbose = FALSE)
  df_use <- parsed$df
  if (is.null(response)) response <- parsed$response
  if (is.null(predictor)) predictor <- parsed$predictor[1]

  sub_df <- df_use[complete.cases(df_use[, c(response, predictor)]), ]
  y_raw <- sub_df[[response]]
  x_vec <- sub_df[[predictor]]

  if (length(unique(y_raw)) > 2 && length(unique(x_vec)) == 2) {
    tmp_name <- response; response <- predictor; predictor <- tmp_name
    y_raw <- sub_df[[response]]; x_vec <- sub_df[[predictor]]
  }

  if (is.numeric(y_raw) && all(y_raw %in% c(0, 1))) {
    y_num <- y_raw
    pos_label <- "1 (Success)"
  } else {
    y_fac <- as.factor(y_raw)
    if (nlevels(y_fac) != 2) stop("biostat_logistic() requires a binary (2-level) response variable.")
    y_num <- as.numeric(y_fac) - 1
    pos_label <- levels(y_fac)[2]
  }
  sub_df$.y_num <- y_num

  cat("\n======================================================================\n")
  cat("        BINARY LOGISTIC REGRESSION & ODDS RATIO ANALYSIS              \n")
  cat("======================================================================\n")
  cat(sprintf("Dependent Variable (Binary Y) : '%s' (Modeling probability of: %s)\n", response, pos_label))
  cat(sprintf("Independent Variable (Cont X) : '%s' (n = %d observations)\n\n", predictor, nrow(sub_df)))

  fml <- as.formula(paste(".y_num ~", .sym_q(predictor)))
  fit_glm <- glm(fml, data = sub_df, family = binomial(link = "logit"))
  sum_glm <- summary(fit_glm)
  print(sum_glm)

  null_dev <- fit_glm$null.deviance
  res_dev  <- fit_glm$deviance
  df_diff  <- fit_glm$df.null - fit_glm$df.residual
  lrt_chisq <- null_dev - res_dev
  lrt_p     <- pchisq(lrt_chisq, df = df_diff, lower.tail = FALSE)

  b0 <- coef(fit_glm)[1]
  b1 <- coef(fit_glm)[2]
  wald_z <- sum_glm$coefficients[2, "z value"]
  wald_p <- sum_glm$coefficients[2, "Pr(>|z|)"]
  or_factor <- exp(b1)
  ci_b <- suppressMessages(confint(fit_glm))
  or_ci <- exp(ci_b[2, ])

  mcfadden_r2 <- 1 - (res_dev / null_dev)
  n_obs <- nrow(sub_df)
  cox_snell_r2 <- 1 - exp((res_dev - null_dev) / n_obs)
  nagelkerke_r2 <- cox_snell_r2 / (1 - exp(-null_dev / n_obs))
  preds <- fitted(fit_glm)
  tjur_r2 <- mean(preds[y_num == 1]) - mean(preds[y_num == 0])
  pb_cor <- suppressWarnings(cor.test(x_vec, y_num))

  cat("--- KEY EXAM & PROBLEM SET ANSWERS ---\n")
  cat(sprintf("1. Hypothesis Test (Effect of '%s' on '%s'):\n", predictor, response))
  cat(sprintf("   * Wald Z-test                : z = %.4f, p = %.5g (%s)\n",
              wald_z, wald_p, ifelse(wald_p < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT")))
  cat(sprintf("   * Likelihood Ratio Test (G2) : Chi-Sq = %.4f, df = %d, p = %.5g\n\n",
              lrt_chisq, df_diff, lrt_p))

  cat(sprintf("2. Odds Multiplier (Numeric Factor Increasing Odds of '%s'):\n", response))
  cat(sprintf("   * Log-Odds Slope (b1)        : %.4f\n", b1))
  cat(sprintf("   * Odds Ratio Factor exp(b1)  : %.4f (95%% CI: [%.4f, %.4f])\n",
              or_factor, or_ci[1], or_ci[2]))
  cat(sprintf("   -> INTERPRETATION: For each 1-unit increase in '%s', the odds of '%s'\n", predictor, response))
  cat(sprintf("      increase by a numeric factor of %.3f (a %.1f%% increase in odds per unit).\n\n",
              or_factor, (or_factor - 1) * 100))

  cat("3. Best-Fit Logistic Equation & Strength of Association:\n")
  cat(sprintf("   * Log-Odds Equation          : ln(p / (1 - p)) = %.4f %s %.4f * (%s)\n",
              b0, ifelse(b1 >= 0, "+", "-"), abs(b1), predictor))
  cat(sprintf("   * Predicted Probability Eq.  : p = 1 / (1 + exp(-(%.4f %s %.4f * %s)))\n",
              b0, ifelse(b1 >= 0, "+", "-"), abs(b1), predictor))
  mcf_note <- if (mcfadden_r2 >= 0.4) "Excellent fit (> 0.40)" else if (mcfadden_r2 >= 0.2) "Good / strong fit (0.20 - 0.40)" else "Modest fit (< 0.20)"
  cat(sprintf("   * McFadden's Pseudo-R2       : %.4f (%.1f%% reduction in deviance -> %s)\n",
              mcfadden_r2, mcfadden_r2 * 100, mcf_note))
  cat(sprintf("   * Nagelkerke's Pseudo-R2     : %.4f\n", nagelkerke_r2))
  cat(sprintf("   * Cox & Snell Pseudo-R2      : %.4f\n", cox_snell_r2))
  cat(sprintf("   * Tjur's Coefficient of Disc.: %.4f\n", tjur_r2))
  cat(sprintf("   * Point-Biserial Correlation : r = %.4f (p = %.5g)\n", pb_cor$estimate, pb_cor$p.value))

  .biostat_last_audit$last_code <- biostat_code(df_use, response = response, predictor = predictor, silent = !isTRUE(show_code))

  if (!isTRUE(show_code)) cat("\n[TIP]: Run `biostat_code()` anytime to print the standalone R code & graph code.\n")
  if (plot) {
    tryCatch({
      biostat_resid_plot(fit_int)
      cat("[PLOTS]: Generated 2 plots! Click '<' in the Plots pane for the 3-Panel Assumption Grid.\n")
    }, error = function(e) NULL)
    if (.has_pkg("ggplot2")) {
      p <- ggplot2::ggplot(sub_df, ggplot2::aes(x = .data[[predictor]], y = .y_num)) +
        ggplot2::geom_jitter(height = 0.03, width = 0, alpha = 0.65, size = 2.5, color = "steelblue4") +
        ggplot2::stat_smooth(method = "glm", method.args = list(family = "binomial"),
                             se = TRUE, color = "firebrick", fill = "pink", alpha = 0.3, linewidth = 1.1) +
        ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(-0.05, 1.05)) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::labs(
          title = paste("Logistic Regression:", response, "vs.", predictor),
          subtitle = sprintf("Odds Factor exp(b) = %.3f | Wald p = %.4g | McFadden R2 = %.3f",
                             or_factor, wald_p, mcfadden_r2),
          x = predictor,
          y = paste("Probability of", response, "(1 =", pos_label, ")")
        )
      print(p)
    } else {
      plot(x_vec, jitter(y_num, amount = 0.03), pch = 19, col = "steelblue",
           xlab = predictor, ylab = paste("Probability of", response),
           main = paste("Logistic Regression:", response, "vs.", predictor))
      x_grid <- seq(min(x_vec), max(x_vec), length.out = 200)
      nd <- data.frame(x_grid); names(nd) <- predictor
      lines(x_grid, predict(fit_glm, newdata = nd, type = "response"), col = "firebrick", lwd = 2.5)
    }
  }

  invisible(list(model = fit_glm, odds_ratio = or_factor, or_ci = or_ci,
                 mcfadden_r2 = mcfadden_r2, nagelkerke_r2 = nagelkerke_r2, lrt_p = lrt_p))
}

# ------------------------------------------------------------------------------
# 7. HELPER: BARPLOT OF MEANS +/- SEM FOR 1-WAY OR 2-WAY FACTORIAL DESIGNS
# ------------------------------------------------------------------------------
biostat_sem_barplot <- function(df, response, factor1, factor2 = NULL, title_str = NULL, y_label = NULL) {
  sub_df <- df[complete.cases(df[, c(response, factor1, factor2)]), ]
  sub_df[[factor1]] <- as.factor(sub_df[[factor1]])
  if (!is.null(factor2)) sub_df[[factor2]] <- as.factor(sub_df[[factor2]])
  y_lab_use <- ifelse(is.null(y_label), paste(response, "(Mean +/- SEM)"), paste(y_label, "(Mean +/- SEM)"))

  if (.has_pkg("ggplot2")) {
    if (is.null(factor2)) {
      p <- ggplot2::ggplot(sub_df, ggplot2::aes(x = .data[[factor1]], y = .data[[response]], fill = .data[[factor1]])) +
        ggplot2::stat_summary(fun = mean, geom = "bar", width = 0.6, color = "black", alpha = 0.8) +
        ggplot2::stat_summary(fun.data = function(y) {
          m <- mean(y); se <- sd(y) / sqrt(length(y))
          data.frame(y = m, ymin = m - se, ymax = m + se)
        }, geom = "errorbar", width = 0.2, linewidth = 0.8) +
        ggplot2::theme_classic(base_size = 13) +
        ggplot2::labs(title = ifelse(is.null(title_str), paste(y_lab_use, "by", factor1), title_str),
                      x = factor1, y = y_lab_use) +
        ggplot2::theme(legend.position = "none")
      print(p)
      return(invisible(p))
    } else {
      p <- ggplot2::ggplot(sub_df, ggplot2::aes(x = .data[[factor1]], y = .data[[response]], fill = .data[[factor2]])) +
        ggplot2::stat_summary(fun = mean, geom = "bar", position = ggplot2::position_dodge(width = 0.75),
                              width = 0.68, color = "black", alpha = 0.85) +
        ggplot2::stat_summary(fun.data = function(y) {
          m <- mean(y); se <- sd(y) / sqrt(length(y))
          data.frame(y = m, ymin = m - se, ymax = m + se)
        }, geom = "errorbar", position = ggplot2::position_dodge(width = 0.75),
        width = 0.22, linewidth = 0.75) +
        ggplot2::theme_classic(base_size = 13) +
        ggplot2::labs(
          title = ifelse(is.null(title_str), paste(y_lab_use, "by", factor1, "and", factor2), title_str),
          x = factor1, y = y_lab_use, fill = factor2
        )
      print(p)
      return(invisible(p))
    }
  } else {
    if (is.null(factor2)) {
      means <- tapply(sub_df[[response]], sub_df[[factor1]], mean)
      ses   <- tapply(sub_df[[response]], sub_df[[factor1]], function(v) sd(v)/sqrt(length(v)))
      bp <- barplot(means, ylim = c(0, max(means + ses, na.rm = TRUE) * 1.15),
                    col = "steelblue3", ylab = y_lab_use, xlab = factor1,
                    main = ifelse(is.null(title_str), paste(y_lab_use, "by", factor1), title_str))
      arrows(bp, means - ses, bp, means + ses, angle = 90, code = 3, length = 0.08, lwd = 1.5)
    } else {
      means <- tapply(sub_df[[response]], list(sub_df[[factor2]], sub_df[[factor1]]), mean)
      ses   <- tapply(sub_df[[response]], list(sub_df[[factor2]], sub_df[[factor1]]), function(v) sd(v)/sqrt(length(v)))
      bp <- barplot(means, beside = TRUE, ylim = c(0, max(means + ses, na.rm = TRUE) * 1.2),
                    col = c("salmon", "darkturquoise", "goldenrod", "steelblue")[seq_len(nrow(means))],
                    legend.text = TRUE, args.legend = list(title = factor2, x = "topright"),
                    xlab = factor1, ylab = y_lab_use,
                    main = ifelse(is.null(title_str), paste(y_lab_use, "by", factor1, "&", factor2), title_str))
      arrows(bp, means - ses, bp, means + ses, angle = 90, code = 3, length = 0.06, lwd = 1.5)
    }
  }
}

# ------------------------------------------------------------------------------
# 8. DEDICATED SOLVER 2: TWO-WAY MODEL I FACTORIAL ANOVA (PS7 Q2: Coral Growth)
# ------------------------------------------------------------------------------
biostat_factorial <- function(df, response = NULL, factor1 = NULL, factor2 = NULL,
                              transform = "auto", alpha = 0.05, plot = TRUE, show_code = FALSE) {
  parsed <- .auto_parse_roles(df, response = response, predictor = c(factor1, factor2), verbose = FALSE)
  df_use <- parsed$df
  if (is.null(response)) response <- parsed$response
  if (is.null(factor1))  factor1  <- parsed$predictor[1]
  if (is.null(factor2))  factor2  <- parsed$predictor[2]

  sub_df <- df_use[complete.cases(df_use[, c(response, factor1, factor2)]), ]
  sub_df[[factor1]] <- as.factor(sub_df[[factor1]])
  sub_df[[factor2]] <- as.factor(sub_df[[factor2]])

  # Explicit Transformation Option
  y_disp_label <- response
  if (!transform %in% c("auto", "none", "raw")) {
    tr_app <- .apply_transform_df(sub_df, response, c(factor1, factor2), transform)
    sub_df <- tr_app$df
    y_disp_label <- tr_app$y_label
  }

  cat("\n======================================================================\n")
  cat("          FACTORIAL TWO-WAY MODEL I ANOVA (FIXED EFFECTS)             \n")
  cat("======================================================================\n")
  cat(sprintf("Response Variable (Y) : '%s'\n", y_disp_label))
  cat(sprintf("Fixed Factor A        : '%s' (%s)\n", factor1, paste(levels(sub_df[[factor1]]), collapse = ", ")))
  cat(sprintf("Fixed Factor B        : '%s' (%s)\n\n", factor2, paste(levels(sub_df[[factor2]]), collapse = ", ")))

  cat("--- CELL MEANS +/- SEM SUMMARY TABLE ---\n")
  cell_stats <- aggregate(sub_df[[response]], by = list(sub_df[[factor1]], sub_df[[factor2]]),
                          FUN = function(v) c(n = length(v), Mean = mean(v), SD = sd(v), SEM = sd(v)/sqrt(length(v))))
  cell_df <- data.frame(
    Factor1 = cell_stats[[1]], Factor2 = cell_stats[[2]],
    n    = cell_stats$x[, "n"],
    Mean = round(cell_stats$x[, "Mean"], 3),
    SD   = round(cell_stats$x[, "SD"], 3),
    SEM  = round(cell_stats$x[, "SEM"], 3)
  )
  names(cell_df)[1:2] <- c(factor1, factor2)
  print(cell_df, row.names = FALSE)

  fml_int <- as.formula(paste(.sym_q(response), "~", .sym_q(factor1), "*", .sym_q(factor2)))
  fit_int <- lm(fml_int, data = sub_df)
  aov_int <- anova(fit_int)

  cat("\n--- STEP 1: FULL MODEL WITH INTERACTION (", deparse(fml_int), ") ---\n", sep = "")
  print(aov_int)

  sw_p <- if (nrow(sub_df) >= 3 && nrow(sub_df) <= 5000) shapiro.test(residuals(fit_int))$p.value else NA
  flig_p <- tryCatch(fligner.test(sub_df[[response]] ~ interaction(sub_df[[factor1]], sub_df[[factor2]]))$p.value, error = function(e) NA)
  cat(sprintf("Assumption Check -> Residual Shapiro-Wilk p = %.4g | Cell Variance Fligner p = %.4g\n", sw_p, flig_p))

  # Auto-Check Transformations if Residuals or Cell Variances Fail (like sqrt(NoFlrs) in Bromus Grass lecture!)
  if (transform == "auto" && ((!is.na(sw_p) && sw_p < alpha) || (!is.na(flig_p) && flig_p < alpha))) {
    cat("\n[AUTO-TRANSFORMATION CHECK]: Raw residuals or cell variances violated assumptions!\n")
    tr_res <- biostat_transform(sub_df, response = response, predictor = c(factor1, factor2),
                                alpha = alpha, plot = FALSE, verbose = TRUE)
    if (!is.null(tr_res) && isTRUE(tr_res$rescued)) {
      cat(sprintf("\n-> SWITCHING PRIMARY ANALYSIS & GRAPH TO RESCUED SCALE: %s\n", tr_res$best_label))
      sub_df <- tr_res$df_transformed
      y_disp_label <- tr_res$best_label
      fit_int <- lm(fml_int, data = sub_df)
      aov_int <- anova(fit_int)
      cat(sprintf("\n--- TRANSFORMED FULL MODEL WITH INTERACTION (%s ~ %s * %s) ---\n",
                  y_disp_label, factor1, factor2))
      print(aov_int)
    }
  }

  int_p <- aov_int[3, "Pr(>F)"]
  int_f <- aov_int[3, "F value"]

  cat("\n--- STEP 2: COURSE DECISION RULE ON INTERACTION TERM ---\n")
  if (!is.na(int_p) && int_p < alpha) {
    cat(sprintf("[INTERACTION IS SIGNIFICANT]: %s:%s F = %.3f, p = %.5g (< %.2f)\n", factor1, factor2, int_f, int_p, alpha))
    cat(sprintf("-> CONCLUSION: The effect of '%s' on '%s' DEPENDS UPON '%s'!\n", factor1, y_disp_label, factor2))
    cat("-> COURSE CONVENTION: Do NOT interpret the main effects from the interaction ANOVA table.\n")
    cat("   Instead, analyze the factors separately using two-sample t-tests (below):\n\n")

    for (fac in c(factor1, factor2)) {
      if (nlevels(sub_df[[fac]]) == 2) {
        f_fml <- as.formula(paste(.sym_q(response), "~", .sym_q(fac)))
        bart_p <- tryCatch(bartlett.test(f_fml, data = sub_df)$p.value, error = function(e) NA)
        flig_p <- tryCatch(fligner.test(f_fml, data = sub_df)$p.value, error = function(e) 0)
        veq <- if (!is.na(bart_p)) (bart_p >= alpha) else (!is.na(flig_p) && flig_p >= alpha)
        tt <- t.test(f_fml, data = sub_df, var.equal = veq)
        cat(sprintf("   * Separate Two-Sample t-test for '%s' (Bartlett p = %.4g -> var.equal = %s [%s]):\n",
                    fac, bart_p, veq, ifelse(veq, "Student's t-test", "Welch's t-test")))
        cat(sprintf("     t = %.3f, df = %.3f, p = %.5g (Means: %s)\n",
                    tt$statistic, tt$parameter, tt$p.value,
                    paste(names(tt$estimate), round(tt$estimate, 3), sep = " = ", collapse = ", ")))
      }
    }

    cat(sprintf("\n   * Simple Effects of '%s' within each level of '%s':\n", factor1, factor2))
    for (lvl in levels(sub_df[[factor2]])) {
      sub_lvl <- sub_df[sub_df[[factor2]] == lvl, ]
      if (nlevels(sub_lvl[[factor1]]) == 2) {
        f_fml <- as.formula(paste(.sym_q(response), "~", .sym_q(factor1)))
        var_p <- tryCatch(fligner.test(f_fml, data = sub_lvl)$p.value, error = function(e) 0)
        veq <- (!is.na(var_p) && var_p >= alpha)
        tt <- t.test(f_fml, data = sub_lvl, var.equal = veq)
        cat(sprintf("     - When %s = '%s': t = %.3f, df = %.2f, p = %.5g\n",
                    factor2, lvl, tt$statistic, tt$parameter, tt$p.value))
      }
    }
    final_model <- fit_int

  } else {
    cat(sprintf("[INTERACTION IS NOT SIGNIFICANT]: %s:%s F = %.3f, p = %.5g (>= %.2f)\n", factor1, factor2, int_f, int_p, alpha))
    cat(sprintf("-> CONCLUSION: The effect of '%s' on '%s' does NOT depend upon '%s'.\n", factor1, y_disp_label, factor2))
    cat("-> COURSE CONVENTION: Drop the interaction term and re-run the Additive Model\n")
    cat("   to test each main effect while holding the other factor constant:\n\n")

    fml_add <- as.formula(paste(.sym_q(response), "~", .sym_q(factor1), "+", .sym_q(factor2)))
    fit_add <- lm(fml_add, data = sub_df)
    aov_add <- anova(fit_add)
    cat("--- ADDITIVE TWO-WAY ANOVA TABLE (", deparse(fml_add), ") ---\n", sep = "")
    print(aov_add)
    final_model <- fit_add
  }

  .biostat_last_audit$last_code <- biostat_code(df_use, response = response, predictor = c(factor1, factor2), silent = !isTRUE(show_code))

  if (plot) {
    tryCatch({
      biostat_resid_plot(fit_int)
      cat("\n[PLOTS]: Generated 2 plots! Click '<' in the Plots pane for the 3-Panel Assumption Grid.\n")
    }, error = function(e) NULL)
    biostat_sem_barplot(sub_df, response, factor1, factor2, y_label = y_disp_label)
  }
  if (!isTRUE(show_code)) cat("[TIP]: Run `biostat_code()` anytime to print the standalone R code & graph code.\n")

  invisible(final_model)
}

# ------------------------------------------------------------------------------
# 9. DEDICATED SOLVER 3: MODEL III MIXED-EFFECTS FACTORIAL ANOVA (PS7 Q3: Medicago)
# ------------------------------------------------------------------------------
biostat_model3 <- function(df, response = NULL, fixed = NULL, random = NULL,
                           transform = "auto", alpha = 0.05, plot = TRUE, show_code = FALSE) {
  parsed <- .auto_parse_roles(df, response = response, predictor = c(fixed, random),
                              model = "III", random = random, verbose = FALSE)
  df_use <- parsed$df
  if (is.null(response)) response <- parsed$response

  preds <- parsed$predictor
  if (is.null(fixed) && is.null(random)) {
    random <- parsed$random
    if (is.null(random)) random <- preds[2]
    fixed <- setdiff(preds, random)[1]
  } else if (!is.null(fixed) && is.null(random)) {
    random <- setdiff(preds, fixed)[1]
  } else if (is.null(fixed) && !is.null(random)) {
    fixed <- setdiff(preds, random)[1]
  }

  sub_df <- df_use[complete.cases(df_use[, c(response, fixed, random)]), ]
  sub_df[[fixed]]  <- as.factor(sub_df[[fixed]])
  sub_df[[random]] <- as.factor(sub_df[[random]])

  y_disp_label <- response
  if (!transform %in% c("auto", "none", "raw")) {
    tr_app <- .apply_transform_df(sub_df, response, c(fixed, random), transform)
    sub_df <- tr_app$df
    y_disp_label <- tr_app$y_label
  }

  cat("\n======================================================================\n")
  cat("       MIXED-MODEL FACTORIAL ANOVA (MODEL III ANOVA + LLRT)           \n")
  cat("======================================================================\n")
  cat(sprintf("Response Variable (Y) : '%s'\n", y_disp_label))
  cat(sprintf("Fixed Effect Factor   : '%s' (%s)\n", fixed, paste(levels(sub_df[[fixed]]), collapse = ", ")))
  cat(sprintf("Random Effect Factor  : '%s' (%s)\n", random, paste(levels(sub_df[[random]]), collapse = ", ")))
  cat("(Note: To swap Fixed and Random roles, call: biostat_model3(data, fixed='...', random='...'))\n\n")

  cat("--- CELL MEANS +/- SEM SUMMARY TABLE ---\n")
  cell_stats <- aggregate(sub_df[[response]], by = list(sub_df[[fixed]], sub_df[[random]]),
                          FUN = function(v) c(n = length(v), Mean = mean(v), SD = sd(v), SEM = sd(v)/sqrt(length(v))))
  cell_df <- data.frame(
    Fixed = cell_stats[[1]], Random = cell_stats[[2]],
    n    = cell_stats$x[, "n"],
    Mean = round(cell_stats$x[, "Mean"], 3),
    SD   = round(cell_stats$x[, "SD"], 3),
    SEM  = round(cell_stats$x[, "SEM"], 3)
  )
  names(cell_df)[1:2] <- c(fixed, random)
  print(cell_df, row.names = FALSE)

  fq <- .sym_q(fixed); rq <- .sym_q(random); yq <- .sym_q(response)

  cat("\n--- CLASSICAL FACTORIAL ANOVA PARTITIONING (FOR F-RATIO REFERENCE) ---\n")
  fit_ols <- lm(as.formula(paste(yq, "~", fq, "*", rq)), data = sub_df)
  print(anova(fit_ols))

  sw_p <- if (nrow(sub_df) >= 3 && nrow(sub_df) <= 5000) shapiro.test(residuals(fit_ols))$p.value else NA
  flig_p <- tryCatch(fligner.test(sub_df[[response]] ~ interaction(sub_df[[fixed]], sub_df[[random]]))$p.value, error = function(e) NA)
  cat(sprintf("Assumption Check -> Residual Shapiro-Wilk p = %.4g | Cell Variance Fligner p = %.4g\n", sw_p, flig_p))

  if (.has_pkg("lme4")) {
    fml_full   <- as.formula(sprintf("%s ~ %s + (1|%s) + (1|%s:%s)", yq, fq, rq, fq, rq))
    fml_no_int <- as.formula(sprintf("%s ~ %s + (1|%s)", yq, fq, rq))
    fml_null   <- as.formula(sprintf("%s ~ 1 + (1|%s)", yq, rq))

    m_full   <- suppressMessages(lme4::lmer(fml_full, data = sub_df, REML = FALSE))
    m_no_int <- suppressMessages(lme4::lmer(fml_no_int, data = sub_df, REML = FALSE))
    m_null   <- suppressMessages(lme4::lmer(fml_null, data = sub_df, REML = FALSE))

    llrt_int <- anova(m_no_int, m_full)
    chisq_int <- llrt_int$Chisq[2]
    df_int    <- llrt_int$Df[2]
    p_int     <- llrt_int$`Pr(>Chisq)`[2]

    cat("\n--- STEP 1: LOG-LIKELIHOOD RATIO TEST (LLRT) FOR INTERACTION ---\n")
    cat(sprintf("Comparing Full Model [%s] vs. No-Interaction Model [%s]:\n",
                deparse(fml_full), deparse(fml_no_int)))
    print(llrt_int)
    cat(sprintf("\n-> INTERACTION LLRT RESULT: Chi-Sq = %.4f, df = %d, p = %.5g (%s)\n",
                chisq_int, df_int, p_int, ifelse(!is.na(p_int) && p_int < alpha, "SIGNIFICANT", "NOT SIGNIFICANT")))

    cat("\n--- STEP 2: FIXED EFFECT TEST (IN MODEL WITH RANDOM EFFECT & NO INTERACTION) ---\n")
    cat(sprintf("Model Formula: %s\n", deparse(fml_no_int)))
    m_no_int_reml <- suppressMessages(lme4::lmer(fml_no_int, data = sub_df, REML = TRUE))
    aov_fixed <- anova(m_no_int_reml)
    f_val <- aov_fixed[1, "F value"]
    den_df <- nrow(sub_df) - nlevels(sub_df[[fixed]]) - nlevels(sub_df[[random]]) + 1
    num_df <- nlevels(sub_df[[fixed]]) - 1
    p_f_fixed <- pf(f_val, df1 = num_df, df2 = den_df, lower.tail = FALSE)
    llrt_fixed <- anova(m_null, m_no_int)

    cat(sprintf("  * Fixed Effect ('%s') F-test : F(%d, %d) = %.4f, p = %.5g\n",
                fixed, num_df, den_df, f_val, p_f_fixed))
    cat(sprintf("  * Fixed Effect ('%s') LLRT   : Chi-Sq = %.4f, df = %d, p = %.5g\n",
                fixed, llrt_fixed$Chisq[2], llrt_fixed$Df[2], llrt_fixed$`Pr(>Chisq)`[2]))
    cat("  * Additive Fixed+Random OLS ANOVA Table (No Interaction):\n")
    print(anova(lm(as.formula(paste(yq, "~", fq, "+", rq)), data = sub_df)))
  }

  # Auto-Check Transformations if Count/Skewed Residuals or Unequal Variances are Detected!
  if (transform == "auto" && ((!is.na(sw_p) && sw_p < alpha) || (!is.na(flig_p) && flig_p < alpha))) {
    cat("\n[AUTO-TRANSFORMATION CHECK]: Count/skewed response violated normality or equal variance!\n")
    tr_res <- biostat_transform(sub_df, response = response, predictor = c(fixed, random),
                                alpha = alpha, plot = FALSE, verbose = TRUE)
    if (!is.null(tr_res) && tr_res$best_code != "none" && .has_pkg("lme4")) {
      cat(sprintf("\n--- TRANSFORMED MODEL III CHECK (%s: Fixes Homoscedasticity!) ---\n", tr_res$best_name))
      df_tr <- tr_res$df_transformed
      mt_full <- suppressMessages(lme4::lmer(fml_full, data = df_tr, REML = FALSE))
      mt_noint <- suppressMessages(lme4::lmer(fml_no_int, data = df_tr, REML = FALSE))
      mt_null  <- suppressMessages(lme4::lmer(fml_null, data = df_tr, REML = FALSE))
      llrt_tr  <- anova(mt_noint, mt_full)
      llrt_f_tr <- anova(mt_null, mt_noint)
      cat(sprintf("  * Transformed Interaction LLRT : Chi-Sq = %.4f, df = %d, p = %.5g\n",
                  llrt_tr$Chisq[2], llrt_tr$Df[2], llrt_tr$`Pr(>Chisq)`[2]))
      cat(sprintf("  * Transformed Fixed Effect LLRT: Chi-Sq = %.4f, df = %d, p = %.5g\n",
                  llrt_f_tr$Chisq[2], llrt_f_tr$Df[2], llrt_f_tr$`Pr(>Chisq)`[2]))
    }
  }

  .biostat_last_audit$last_code <- biostat_code(df_use, response = response, predictor = c(fixed, random),
                                                model = "III", random = random, silent = !isTRUE(show_code))

  if (plot) {
    biostat_sem_barplot(sub_df, response, fixed, random, y_label = y_disp_label,
                        title_str = paste("Model III ANOVA:", y_disp_label, "by", fixed, "&", random, "(Mean +/- SEM)"))
  }

  invisible(sub_df)
}

# ------------------------------------------------------------------------------
# 10. DEDICATED SOLVER 4: ANCOVA (PS7 Q4: Red Algae Stipes / Mussels / Crickets)
# ------------------------------------------------------------------------------
# Automatically runs the Raw ANCOVA decision tree AND checks if the continuous
# covariate (X) or response (Y) is skewed (like Thickness in Stipe data.csv or
# LogLength in the Mussels lecture example!), showing the Log-Transformed ANCOVA
# side by side when skewness is detected!
# ------------------------------------------------------------------------------
biostat_ancova <- function(df, response = NULL, covariate = NULL, group = NULL,
                           transform = "auto", alpha = 0.05, plot = TRUE, show_code = FALSE) {
  parsed <- .auto_parse_roles(df, response = response, predictor = c(covariate, group), verbose = FALSE)
  df_use <- parsed$df
  if (is.null(response)) response <- parsed$response

  preds <- parsed$predictor
  if (is.null(covariate) || is.null(group)) {
    is_cat <- sapply(preds, function(p) {
      v <- df_use[[p]]
      is.factor(v) || is.character(v) || length(unique(na.omit(v))) <= 6
    })
    if (any(is_cat) && any(!is_cat)) {
      group     <- preds[which(is_cat)[1]]
      covariate <- preds[which(!is_cat)[1]]
    } else {
      covariate <- preds[1]
      group     <- preds[2]
    }
  }

  sub_df <- df_use[complete.cases(df_use[, c(response, covariate, group)]), ]
  sub_df[[group]] <- as.factor(sub_df[[group]])

  y_disp <- response; x_disp <- covariate
  if (!transform %in% c("auto", "none", "raw")) {
    tr_app <- .apply_transform_df(sub_df, response, covariate, transform)
    sub_df <- tr_app$df
    y_disp <- tr_app$y_label
    if (!is.null(tr_app$x_label)) x_disp <- tr_app$x_label
  }

  yq <- .sym_q(response); cq <- .sym_q(covariate); gq <- .sym_q(group)

  cat("\n======================================================================\n")
  cat("        ANALYSIS OF COVARIANCE (ANCOVA: SLOPES & LINE HEIGHTS)        \n")
  cat("======================================================================\n")
  cat(sprintf("Dependent Variable (Cont Y)   : '%s'\n", y_disp))
  cat(sprintf("Covariate Predictor (Cont X)  : '%s'\n", x_disp))
  cat(sprintf("Grouping Variable (Factor G)  : '%s' (%s)\n\n", group, paste(levels(sub_df[[group]]), collapse = ", ")))

  fml_int <- as.formula(paste(yq, "~", cq, "*", gq))
  fit_int <- lm(fml_int, data = sub_df)
  aov_int <- anova(fit_int)

  cat("--- STEP 1: TEST HOMOGENEITY OF SLOPES (INTERACTION MODEL: ", deparse(fml_int), ") ---\n", sep = "")
  print(aov_int)

  sw_res <- if (nrow(sub_df) >= 3 && nrow(sub_df) <= 5000) shapiro.test(residuals(fit_int))$p.value else NA
  sw_y   <- if (nrow(sub_df) >= 3 && nrow(sub_df) <= 5000) shapiro.test(sub_df[[response]])$p.value else NA
  sw_x   <- if (nrow(sub_df) >= 3 && nrow(sub_df) <= 5000) shapiro.test(sub_df[[covariate]])$p.value else NA
  cat(sprintf("Normality Checks (Shapiro-Wilk) -> Residuals p = %.4g | Y ('%s') p = %.4g | X ('%s') p = %.4g\n",
              sw_res, y_disp, sw_y, x_disp, sw_x))

  int_p <- aov_int[3, "Pr(>F)"]
  int_f <- aov_int[3, "F value"]

  cat("\n--- STEP 2: COURSE DECISION RULE ON ANCOVA INTERACTION ---\n")
  if (!is.na(int_p) && int_p < alpha) {
    cat(sprintf("[INTERACTION IS SIGNIFICANT]: %s:%s F = %.3f, p = %.5g (< %.2f)\n", x_disp, group, int_f, int_p, alpha))
    cat(sprintf("-> CONCLUSION: The effect (slope) of '%s' on '%s' DEPENDS UPON '%s'!\n", x_disp, y_disp, group))
    cat("-> COURSE CONVENTION (PS7 Q4 / Crickets Example): Stop interpreting the combined ANOVA table\n")
    cat(sprintf("   and run SEPARATE linear regressions of '%s ~ %s' for each '%s' below:\n\n", y_disp, x_disp, group))

    for (lvl in levels(sub_df[[group]])) {
      sub_g <- sub_df[sub_df[[group]] == lvl, ]
      fml_sep <- as.formula(paste(yq, "~", cq))
      fit_sep <- lm(fml_sep, data = sub_g)
      aov_sep <- anova(fit_sep)
      sum_sep <- summary(fit_sep)

      cat(sprintf("   ================================================================\n"))
      cat(sprintf("   * Separate Linear Regression for %s = '%s' (n = %d):\n", group, lvl, nrow(sub_g)))
      cat(sprintf("     Equation : %s = %.4f + (%.4f * %s)\n",
                  y_disp, coef(fit_sep)[1], coef(fit_sep)[2], x_disp))
      cat(sprintf("     ANOVA    : F(1, %d) = %.3f, Mean Sq = %.3f, p = %.5g, R2 = %.4f\n",
                  fit_sep$df.residual, aov_sep[1, "F value"], aov_sep[1, "Mean Sq"],
                  aov_sep[1, "Pr(>F)"], sum_sep$r.squared))
      print(aov_sep)
      cat("\n")
    }
    final_model <- fit_int

  } else {
    cat(sprintf("[INTERACTION IS NOT SIGNIFICANT]: %s:%s F = %.3f, p = %.5g (>= %.2f)\n", x_disp, group, int_f, int_p, alpha))
    cat(sprintf("-> CONCLUSION: The slopes are parallel! The effect of '%s' on '%s' does NOT depend upon '%s'.\n",
                x_disp, y_disp, group))
    cat("-> COURSE CONVENTION (PS7 Q4 / Mussels Example): Drop the interaction term and run\n")
    cat("   the Additive ANCOVA model to test the common slope and compare the HEIGHTS of the lines:\n\n")

    fml_add <- as.formula(paste(yq, "~", cq, "+", gq))
    fit_add <- lm(fml_add, data = sub_df)
    aov_add <- anova(fit_add)
    cat("--- ADDITIVE ANCOVA TABLE (", deparse(fml_add), " — Comparing Line Heights) ---\n", sep = "")
    print(aov_add)
    cat("\nModel Coefficients (Common Slope & Intercept/Height Difference):\n")
    print(summary(fit_add)$coefficients)
    final_model <- fit_add
  }

  # Auto-Check Log-Transformation when Raw Y or Raw X is Skewed (like Mussels LogLength or Stipe data!)
  if (transform == "auto" && ((!is.na(sw_res) && sw_res < alpha) || (!is.na(sw_y) && sw_y < alpha) || (!is.na(sw_x) && sw_x < alpha)) &&
      min(sub_df[[response]], na.rm = TRUE) > 0 && min(sub_df[[covariate]], na.rm = TRUE) > 0) {
    cat("\n======================================================================\n")
    cat("   [AUTO-TRANSFORMATION CHECK]: LOG-TRANSFORMED ANCOVA COMPARISON     \n")
    cat("======================================================================\n")
    cat(sprintf("Because raw '%s' (Shapiro p = %.4g) or '%s' (Shapiro p = %.4g) is skewed,\n",
                response, sw_y, covariate, sw_x))
    cat("here is the Log10-Transformed ANCOVA (as in the Mytilus Mussels lecture example):\n")
    cat("(Tip: Run `biostat_ancova(data, transform = 'log_xy')` to plot the log-log trendlines!)\n\n")

    tr_res <- biostat_transform(sub_df, response = response, predictor = c(covariate, group),
                                alpha = alpha, plot = FALSE, verbose = TRUE)

    df_ln <- sub_df
    df_ln[[response]]  <- log(sub_df[[response]])
    df_ln[[covariate]] <- log(sub_df[[covariate]])
    fit_ln_int <- lm(fml_int, data = df_ln)
    aov_ln_int <- anova(fit_ln_int)
    cat(sprintf("\n1. Natural Log-Log Interaction ANCOVA (log(%s) ~ log(%s) * %s) [Shapiro p = %.4g]:\n",
                response, covariate, group, shapiro.test(residuals(fit_ln_int))$p.value))
    print(aov_ln_int)

    if (aov_ln_int[3, "Pr(>F)"] >= alpha) {
      cat(sprintf("\n   -> On the log-log scale, the interaction is NOT significant (F = %.3f, p = %.4f — Slopes are equal)!\n",
                  aov_ln_int[3, "F value"], aov_ln_int[3, "Pr(>F)"]))
      cat(sprintf("   -> Dropping interaction to compare line heights (log(%s) ~ log(%s) + %s):\n",
                  response, covariate, group))
      print(anova(lm(as.formula(paste(yq, "~", cq, "+", gq)), data = df_ln)))
    }
  }

  .biostat_last_audit$last_code <- biostat_code(df_use, response = response, predictor = c(covariate, group), silent = !isTRUE(show_code))

  if (plot) {
    if (.has_pkg("ggplot2")) {
      p <- ggplot2::ggplot(sub_df, ggplot2::aes(x = .data[[covariate]], y = .data[[response]], color = .data[[group]])) +
        ggplot2::geom_point(size = 2.6, alpha = 0.8) +
        ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 1.1) +
        ggplot2::theme_classic(base_size = 13) +
        ggplot2::labs(
          title = paste("ANCOVA:", y_disp, "vs.", x_disp, "by", group),
          subtitle = sprintf("Interaction %s:%s p = %.4g (%s)",
                             x_disp, group, int_p,
                             ifelse(!is.na(int_p) && int_p < alpha, "Separate Slopes", "Parallel Slopes")),
          x = x_disp, y = y_disp, color = group
        )
      print(p)
    } else {
      cols_g <- c("firebrick", "darkturquoise", "goldenrod", "navy")[seq_len(nlevels(sub_df[[group]]))]
      plot(sub_df[[covariate]], sub_df[[response]], col = cols_g[sub_df[[group]]], pch = 19,
           xlab = x_disp, ylab = y_disp, main = paste("ANCOVA:", y_disp, "vs.", x_disp, "by", group))
      for (i in seq_along(levels(sub_df[[group]]))) {
        lvl <- levels(sub_df[[group]])[i]
        sub_g <- sub_df[sub_df[[group]] == lvl, ]
        if (nrow(sub_g) >= 2) abline(lm(as.formula(paste(yq, "~", cq)), data = sub_g), col = cols_g[i], lwd = 2)
      }
      legend("topright", legend = levels(sub_df[[group]]), col = cols_g, pch = 19, lwd = 2, title = group)
    }
  }

  invisible(final_model)
}

# ------------------------------------------------------------------------------
# 11. DEDICATED SOLVER 5: MULTIPLE REGRESSION & AIC MODEL SELECTION (PS7 Q5: Symbionts)
# ------------------------------------------------------------------------------
biostat_multireg <- function(df, response = NULL, predictors = NULL,
                             transform = "auto", alpha = 0.05, plot = TRUE, show_code = FALSE) {
  parsed <- .auto_parse_roles(df, response = response, predictor = predictors, verbose = FALSE)
  df_use <- parsed$df
  if (is.null(response))   response   <- parsed$response
  if (is.null(predictors)) predictors <- parsed$predictor

  sub_df <- df_use[complete.cases(df_use[, c(response, predictors)]), c(response, predictors)]
  y_disp <- response
  if (!transform %in% c("auto", "none", "raw")) {
    tr_app <- .apply_transform_df(sub_df, response, predictors, transform)
    sub_df <- tr_app$df
    y_disp <- tr_app$y_label
  }
  yq <- .sym_q(response)

  cat("\n======================================================================\n")
  cat("     MULTIPLE LINEAR REGRESSION, COLLINEARITY & AIC MODEL SELECTION   \n")
  cat("======================================================================\n")
  cat(sprintf("Dependent Response (Y)     : '%s' (n = %d observations)\n", y_disp, nrow(sub_df)))
  cat(sprintf("Candidate Predictors (X_i) : %s\n\n", paste(predictors, collapse = ", ")))

  cat("--- STEP 1: COLLINEARITY CHECK & INDIVIDUAL PREDICTOR CORRELATIONS ---\n")
  cor_full <- round(cor(sub_df, use = "complete.obs"), 4)
  print(cor_full)

  cat("\n1A. Individual Linear Relationships with Response ('", y_disp, "'):\n", sep = "")
  for (xc in predictors) {
    ct_y <- suppressWarnings(cor.test(sub_df[[xc]], sub_df[[response]]))
    m_single <- lm(as.formula(paste(yq, "~", .sym_q(xc))), data = sub_df)
    cat(sprintf("  * %s ~ %-8s : r = %6.3f, R2 = %.4f, F = %6.2f, p = %.4g, AIC = %.2f\n",
                y_disp, xc, ct_y$estimate, summary(m_single)$r.squared,
                summary(m_single)$fstatistic[1], ct_y$p.value, AIC(m_single)))
  }

  if (length(predictors) >= 2) {
    pairs_idx <- combn(predictors, 2, simplify = FALSE)
    collinear_flags <- character(0)
    cat("\n1B. Pairwise Correlations AMONG Predictors (Collinearity Check):\n")
    for (pr in pairs_idx) {
      ct <- suppressWarnings(cor.test(sub_df[[pr[1]]], sub_df[[pr[2]]]))
      sig_str <- ifelse(ct$p.value < 0.001 && abs(ct$estimate) > 0.7, "*** (SEVERE COLLINEARITY!)",
                        ifelse(ct$p.value < 0.05, "* (Significant Correlation)", "NS (Independent)"))
      cat(sprintf("  * %-8s vs. %-8s : r = %.4f, p = %.4g -> %s\n",
                  pr[1], pr[2], ct$estimate, ct$p.value, sig_str))
      if (ct$p.value < 0.05 && abs(ct$estimate) > 0.5) {
        collinear_flags <- c(collinear_flags, sprintf("%s & %s (r=%.3f)", pr[1], pr[2], ct$estimate))
      }
    }
    if (length(collinear_flags) > 0) {
      cat("\n  [COLLINEARITY ALERT]: Predictors are strongly collinear: ",
          paste(collinear_flags, collapse = "; "), "\n", sep = "")
      cat("  -> Including collinear predictors together in one multiple regression inflates variance\n")
      cat("     and masks individual effects! Use the All-Subsets AIC Table below to select the\n")
      cat("     most parsimonious non-collinear model.\n")
    }
  }

  fml_full <- as.formula(paste(yq, "~", paste(sapply(predictors, .sym_q), collapse = " + ")))
  fit_full <- lm(fml_full, data = sub_df)
  sw_p <- if (nrow(sub_df) >= 3 && nrow(sub_df) <= 5000) shapiro.test(residuals(fit_full))$p.value else NA

  cat(sprintf("\n--- STEP 2: FULL MULTIPLE REGRESSION MODEL (%s) ---\n", deparse(fml_full)))
  cat(sprintf("Residual Normality (Shapiro-Wilk): p = %.4g\n\n", sw_p))
  cat("ANOVA Table (Full Model):\n")
  print(anova(fit_full))
  cat("\nRegression Coefficients & Overall R-Squared (Full Model):\n")
  print(summary(fit_full))

  # Check transformations if Full Model residuals are non-normal
  if (transform == "auto" && !is.na(sw_p) && sw_p < alpha) {
    cat("\n[AUTO-TRANSFORMATION CHECK]: Full model residuals are right-skewed (Shapiro p < 0.05)!\n")
    biostat_transform(sub_df, response = response, predictor = predictors, alpha = alpha, plot = FALSE, verbose = TRUE)
    cat("(Note: In AIC model selection, only compare AIC values among models sharing the EXACT same response scale!)\n\n")
  }

  cat("--- STEP 3: ALL-SUBSETS MODEL SELECTION (RANKED BY LOWEST AIC) ---\n")
  all_models_list <- list()
  for (k in seq_along(predictors)) {
    combos <- combn(predictors, k, simplify = FALSE)
    for (cmb in combos) {
      f_str <- paste(yq, "~", paste(sapply(cmb, .sym_q), collapse = " + "))
      m_k   <- lm(as.formula(f_str), data = sub_df)
      s_k   <- summary(m_k)
      fstat <- s_k$fstatistic
      m_p   <- if (!is.null(fstat)) pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE) else NA
      all_models_list[[length(all_models_list) + 1]] <- data.frame(
        Model_Formula = f_str,
        K_Preds       = k,
        AIC           = round(AIC(m_k), 2),
        R2            = round(s_k$r.squared, 4),
        Adj_R2        = round(s_k$adj.r.squared, 4),
        F_stat        = round(fstat[1], 2),
        Model_P       = signif(m_p, 4),
        stringsAsFactors = FALSE
      )
    }
  }
  aic_df <- do.call(rbind, all_models_list)
  aic_df <- aic_df[order(aic_df$AIC), ]
  aic_df$Delta_AIC <- round(aic_df$AIC - min(aic_df$AIC), 2)
  rownames(aic_df) <- NULL
  print(aic_df, row.names = FALSE)

  best_fml_str <- aic_df$Model_Formula[1]
  best_fit <- lm(as.formula(best_fml_str), data = sub_df)
  cf_b <- coef(best_fit)
  eq_terms <- paste(sapply(2:length(cf_b), function(i) {
    sprintf("%s %.4f * (%s)", ifelse(cf_b[i] >= 0, "+", "-"), abs(cf_b[i]), names(cf_b)[i])
  }), collapse = " ")
  cat(sprintf("\n--- STEP 4: BEST SELECTED MODEL BY AIC (%s) ---\n", best_fml_str))
  cat(sprintf("  * Best-Fit Equation : %s = %.4f %s\n", y_disp, cf_b[1], eq_terms))
  cat(sprintf("  * Lowest AIC = %.2f | Multiple R2 = %.4f (%.1f%% of variance in '%s' explained) | Adjusted R2 = %.4f\n\n",
              AIC(best_fit), summary(best_fit)$r.squared, summary(best_fit)$r.squared * 100, y_disp, summary(best_fit)$adj.r.squared))
  cat("Best Model ANOVA Table:\n")
  print(anova(best_fit))
  cat("\nBest Model Regression Summary:\n")
  print(summary(best_fit))

  .biostat_last_audit$last_code <- biostat_code(df_use, response = response, predictor = predictors, silent = !isTRUE(show_code))

  if (plot) {
    panel_cor <- function(x, y, digits = 3, ...) {
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r <- cor(x, y, use = "complete.obs")
      p <- suppressWarnings(cor.test(x, y)$p.value)
      stars <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
      txt <- sprintf("Corr:\n%.3f%s", r, stars)
      text(0.5, 0.5, txt, cex = 1.1, font = ifelse(p < 0.05, 2, 1),
           col = ifelse(p < 0.05, "firebrick", "gray30"))
    }
    panel_smooth_lm <- function(x, y, ...) {
      points(x, y, pch = 19, col = rgb(0.2, 0.4, 0.6, 0.7))
      abline(lm(y ~ x), col = "firebrick", lwd = 1.8)
    }
    pairs(sub_df, lower.panel = panel_smooth_lm, upper.panel = panel_cor,
          main = paste("Scatterplot Matrix (SPLOM):", y_disp, "& Predictors"))
  }

  invisible(list(full_model = fit_full, best_model = best_fit, aic_table = aic_df, cor_matrix = cor_full))
}

# ------------------------------------------------------------------------------
# 12. EXPERIMENTAL DESIGN & ASSUMPTION CHECKERS
# ------------------------------------------------------------------------------
biostat_design_check <- function(df, response = NULL, group = NULL, block = NULL) {
  parsed <- .auto_parse_roles(df, response = response, predictor = group, block = block, verbose = TRUE)
  if (parsed$mode == "table") {
    cat("Input is a contingency table. Use `biostat_contingency(data)`.\n")
    return(invisible(NULL))
  }
  df_use <- parsed$df
  response <- parsed$response
  group <- parsed$predictor[1]

  cat("\n--- EXPERIMENTAL DESIGN DIAGNOSTICS ---\n")
  if (is.null(group)) {
    cat("Single-variable dataset. N =", sum(!is.na(df_use[[response]])), "\n")
    return(invisible(NULL))
  }
  sub_df <- df_use[!is.na(df_use[[response]]) & !is.na(df_use[[group]]), ]
  g <- as.factor(sub_df[[group]])
  counts <- table(g)
  cat(sprintf("1. Replication across '%s' groups:\n", group))
  print(counts)
  is_balanced <- length(unique(counts)) == 1
  cat(sprintf("2. Design Balance: %s\n",
              ifelse(is_balanced, "BALANCED (Equal n across groups)", "UNBALANCED (Unequal n)")))
}

biostat_assumptions <- function(df, response = NULL, predictor = NULL, paired = NULL, alpha = 0.05, verbose_roles = TRUE) {
  parsed <- .auto_parse_roles(df, response = response, predictor = predictor, paired = paired, verbose = verbose_roles)
  if (parsed$mode == "table") {
    sq <- suppressWarnings(chisq.test(parsed$table, correct = FALSE))
    min_e <- min(sq$expected)
    return(invisible(list(parametric_ok = min_e >= 5, min_expected = min_e)))
  }

  df_use <- parsed$df
  response <- parsed$response
  predictor <- parsed$predictor[1]
  paired <- parsed$paired

  cat("\n--- STATISTICAL ASSUMPTION DIAGNOSTICS (RAW DATA) ---\n")
  y <- df_use[[response]]

  if (!is.numeric(y) || (length(unique(na.omit(y))) == 2 && all(na.omit(y) %in% c(0, 1)))) {
    cat(sprintf("Response '%s' is binary/categorical. Normality/homoscedasticity tests do not apply.\n", response))
    return(invisible(list(parametric_ok = FALSE, normality_ok = FALSE, variance_ok = FALSE)))
  }

  if (is.null(predictor)) {
    y_clean <- na.omit(y)
    sw <- if (length(y_clean) >= 3 && length(y_clean) <= 5000) shapiro.test(y_clean) else list(p.value = NA, statistic = NA)
    param_ok <- ifelse(is.na(sw$p.value), TRUE, sw$p.value >= alpha)
    cat(sprintf("1. Normality of '%s' (Shapiro-Wilk): W = %.4f, p = %.4g -> %s\n",
                sw$statistic, sw$p.value, ifelse(param_ok, "PASS (Normal)", "FAIL (Non-normal)")))
    return(invisible(list(parametric_ok = param_ok, normality_ok = param_ok, variance_ok = TRUE)))
  }

  complete_idx <- complete.cases(df_use[, c(response, predictor)])
  y_clean <- df_use[[response]][complete_idx]
  x_clean <- df_use[[predictor]][complete_idx]
  is_cat_x <- is.factor(x_clean) || is.character(x_clean) || (is.numeric(x_clean) && length(unique(x_clean)) <= 8)

  if (is_cat_x && length(unique(x_clean)) == 2 && isTRUE(paired)) {
    x_fac <- as.factor(x_clean)
    g1 <- y_clean[x_fac == levels(x_fac)[1]]
    g2 <- y_clean[x_fac == levels(x_fac)[2]]
    if (length(g1) == length(g2)) {
      diffs <- g1 - g2
      sw_d <- if (length(diffs) >= 3 && length(diffs) <= 5000) shapiro.test(diffs) else list(p.value = NA, statistic = NA)
      norm_ok <- ifelse(is.na(sw_d$p.value), TRUE, sw_d$p.value >= alpha)
      cat(sprintf("1. Paired Differences (d = %s - %s): Mean = %.3f, SD = %.3f, Median = %.3f\n",
                  levels(x_fac)[1], levels(x_fac)[2], mean(diffs), sd(diffs), median(diffs)))
      cat(sprintf("2. Normality of Paired Differences (Shapiro-Wilk): W = %.4f, p = %.4g -> %s\n",
                  sw_d$statistic, sw_d$p.value, ifelse(norm_ok, "PASS (Paired t-test valid)", "FAIL (Check Transformation / Wilcoxon)")))
      return(invisible(list(parametric_ok = norm_ok, normality_ok = norm_ok, variance_ok = TRUE, diffs = diffs)))
    }
  }

  if (is_cat_x) {
    x_fac <- as.factor(x_clean)
    fit <- lm(y_clean ~ x_fac)
    resids <- residuals(fit)
    sw_res <- if (length(resids) >= 3 && length(resids) <= 5000) shapiro.test(resids) else list(p.value = NA)
    grp_sw <- tapply(y_clean, x_fac, function(val) {
      if (length(val) >= 3 && length(val) <= 5000 && length(unique(val)) > 1) shapiro.test(val)$p.value else NA
    })
    bart_p <- tryCatch(bartlett.test(y_clean ~ x_fac)$p.value, error = function(e) NA)
    flig_p <- tryCatch(fligner.test(y_clean ~ x_fac)$p.value, error = function(e) NA)
    var_ok <- ifelse(!is.na(flig_p), flig_p >= alpha, TRUE)
    norm_ok <- ifelse(!is.na(sw_res$p.value), sw_res$p.value >= alpha, TRUE) && all(is.na(grp_sw) | grp_sw >= alpha)

    sk_ku <- .calc_skew_kurt(resids)
    cat(sprintf("1. Residual Normality (Shapiro-Wilk): p = %.4g -> %s\n",
                sw_res$p.value, ifelse(norm_ok, "PASS", "FAIL (Will test log/sqrt transformations!)")))
    cat(sprintf("   * Skewness = %.3f -> %s\n", sk_ku$skew, sk_ku$skew_note))
    cat(sprintf("   * Kurtosis = %.3f -> %s\n", sk_ku$kurt, sk_ku$kurt_note))
    cat("   Group-level Shapiro-Wilk p-values:\n")
    print(round(grp_sw, 4))
    cat(sprintf("2. Homogeneity of Variance (Bartlett & Fligner-Killeen):\n"))
    cat(sprintf("   * Bartlett's Test (Professor's primary test) : p = %.4g -> %s\n",
                bart_p, ifelse(!is.na(bart_p) && bart_p >= alpha, "PASS (Equal variances)", "FAIL (Unequal variances)")))
    cat(sprintf("   * Fligner-Killeen (Robust to non-normality)  : p = %.4g -> %s\n",
                flig_p, ifelse(var_ok, "PASS (Equal variances)", "FAIL (Unequal variances)")))
    return(invisible(list(parametric_ok = norm_ok && var_ok, normality_ok = norm_ok, variance_ok = var_ok)))
  } else {
    fit <- lm(y_clean ~ x_clean)
    resids <- residuals(fit)
    sw_res <- if (length(resids) >= 3 && length(resids) <= 5000) shapiro.test(resids) else list(p.value = NA)
    sw_x <- if (length(x_clean) >= 3 && length(x_clean) <= 5000) shapiro.test(x_clean) else list(p.value = NA)
    sw_y <- if (length(y_clean) >= 3 && length(y_clean) <= 5000) shapiro.test(y_clean) else list(p.value = NA)
    norm_ok <- (!is.na(sw_res$p.value) && sw_res$p.value >= alpha)
    bivar_ok <- (!is.na(sw_x$p.value) && sw_x$p.value >= alpha) && (!is.na(sw_y$p.value) && sw_y$p.value >= alpha)
    sk_ku <- .calc_skew_kurt(resids)
    cat(sprintf("1. Regression Residual Normality (Shapiro-Wilk): p = %.4g -> %s\n",
                sw_res$p.value, ifelse(norm_ok, "PASS", "FAIL (Will test log/sqrt transformations!)")))
    cat(sprintf("   * Residual Skewness = %.3f -> %s\n", sk_ku$skew, sk_ku$skew_note))
    cat(sprintf("   * Residual Kurtosis = %.3f -> %s\n", sk_ku$kurt, sk_ku$kurt_note))
    cat(sprintf("2. Bivariate Normality: Y p = %.4g | X p = %.4g -> %s\n",
                sw_y$p.value, sw_x$p.value, ifelse(bivar_ok, "PASS", "FAIL")))
    return(invisible(list(parametric_ok = norm_ok, normality_ok = norm_ok, bivariate_normal = bivar_ok)))
  }
}

# ------------------------------------------------------------------------------
# 13. PROBLEM TYPE INFERENCE ENGINE
# ------------------------------------------------------------------------------
biostat_infer <- function(df, response = NULL, predictors = NULL, paired = NULL,
                          model = "auto", random = NULL, question = NULL,
                          transform = "auto", mu = 0, verbose_roles = TRUE) {
  parsed <- .auto_parse_roles(df, response = response, predictor = predictors,
                              paired = paired, model = model, random = random,
                              question = question, transform = transform, verbose = verbose_roles)

  if (verbose_roles) {
    cat("\n======================================================================\n")
    cat("                 BIOSTATISTICS PROBLEM INFERENCE ENGINE               \n")
    cat("======================================================================\n")
  }

  if (parsed$mode == "table") {
    if (verbose_roles) cat("INFERRED PROBLEM: Contingency Table Analysis (Chi-Square / Fisher / G-Test)\n")
    return(invisible(list(prob_type = "contingency_table", parsed = parsed)))
  }

  df_use <- parsed$df
  response <- parsed$response
  predictors <- parsed$predictor
  paired <- parsed$paired

  classify_vec <- function(cn) {
    v <- df_use[[cn]]
    v_clean <- na.omit(v)
    n_u <- length(unique(v_clean))
    if (is.factor(v) || is.character(v) || is.logical(v)) {
      return(list(type = "categorical", levels = n_u))
    } else if (is.numeric(v)) {
      if (n_u == 2 && all(v_clean %in% c(0, 1))) return(list(type = "categorical", levels = 2))
      if (n_u <= 4 && all(v_clean == floor(v_clean)) &&
          grepl("(temp|co2|trt|treat|dose|level|group|plot|disturb|range|site|block)", cn, ignore.case = TRUE)) {
        return(list(type = "categorical", levels = n_u))
      }
      return(list(type = "continuous", levels = n_u))
    }
    return(list(type = "unknown", levels = n_u))
  }

  y_info <- classify_vec(response)
  if (verbose_roles) {
    cat(sprintf("Response Variable (Y) : '%s' -> [%s] (%d unique values)\n",
                response, toupper(y_info$type), y_info$levels))
  }

  if (is.null(predictors) || length(predictors) == 0) {
    if (verbose_roles) cat("Predictor Variable (X): None (One-Sample Problem)\n")
    prob <- if (y_info$type == "continuous") "one_sample_continuous" else "one_sample_categorical"
    return(invisible(list(prob_type = prob, parsed = parsed)))
  }

  if (length(predictors) == 1) {
    pred <- predictors[1]
    x_info <- classify_vec(pred)
    if (verbose_roles) {
      cat(sprintf("Predictor Variable (X): '%s' -> [%s] (%d unique values)\n",
                  pred, toupper(x_info$type), x_info$levels))
      cat(sprintf("Study Design          : %s\n\n", ifelse(paired, "PAIRED / MATCHED DESIGN", "INDEPENDENT OBSERVATIONS")))
    }

    if (y_info$type == "continuous" && x_info$type == "categorical" && x_info$levels == 2) {
      prob <- if (paired) "paired_2group" else "indep_2group"
      if (verbose_roles) {
        cat(sprintf("INFERRED PROBLEM: %s\n", ifelse(paired, "Paired Two-Sample Comparison (Paired t-test / Transform / Wilcoxon)",
                                                        "Independent Two-Sample Comparison (Student/Welch t-test / Transform / Mann-Whitney U)")))
      }
      return(invisible(list(prob_type = prob, parsed = parsed)))
    }
    if (y_info$type == "continuous" && x_info$type == "categorical" && x_info$levels > 2) {
      prob <- if (paired) "paired_kgroup" else "indep_kgroup"
      if (verbose_roles) {
        cat(sprintf("INFERRED PROBLEM: %s\n", ifelse(paired, "Blocked / Repeated-Measures ANOVA / Friedman Test",
                                                        "One-Way ANOVA / Welch ANOVA / Transform / Kruskal-Wallis")))
      }
      return(invisible(list(prob_type = prob, parsed = parsed)))
    }
    if (y_info$type == "continuous" && x_info$type == "continuous") {
      if (verbose_roles) cat("INFERRED PROBLEM: Simple Linear Regression & Correlation (OLS / Transform / Pearson r / Spearman rho)\n")
      return(invisible(list(prob_type = "continuous_xy", parsed = parsed)))
    }
    if (y_info$type == "categorical" && y_info$levels == 2 && x_info$type == "continuous") {
      if (verbose_roles) cat("INFERRED PROBLEM: Binary Logistic Regression (Binary Y ~ Continuous X, Odds Multiplier exp(b))\n")
      return(invisible(list(prob_type = "logistic", parsed = parsed)))
    }
    if (y_info$type == "categorical" && x_info$type == "categorical") {
      if (verbose_roles) cat("INFERRED PROBLEM: Contingency Table Association (Chi-Square / Fisher Exact / G-Test)\n")
      return(invisible(list(prob_type = "contingency", parsed = parsed)))
    }
  }

  if (length(predictors) >= 2) {
    pred_types <- sapply(predictors, function(p) classify_vec(p)$type)
    if (verbose_roles) {
      for (i in seq_along(predictors)) {
        cat(sprintf("Predictor X%d          : '%s' -> [%s]\n", i, predictors[i], toupper(pred_types[i])))
      }
      cat("\n")
    }
    if (y_info$type == "continuous") {
      if (length(predictors) == 2 && all(pred_types == "categorical")) {
        if (!is.null(parsed$random) || identical(model, "III") || identical(model, 3) || identical(tolower(model), "model3")) {
          if (verbose_roles) cat("INFERRED PROBLEM: Mixed-Model Factorial ANOVA (Model III ANOVA: 1 Fixed + 1 Random Factor w/ LLRT)\n")
          return(invisible(list(prob_type = "model3_anova", parsed = parsed)))
        } else {
          if (verbose_roles) {
            cat("INFERRED PROBLEM: Two-Way Factorial Model I ANOVA (Tests Interaction A*B -> Separate t-tests or A+B)\n")
            cat("  (Tip: If one factor is Random [Model III ANOVA], run: `biostat_wizard(data, model = 'III')`)\n")
          }
          return(invisible(list(prob_type = "factorial_anova", parsed = parsed)))
        }
      } else if (length(predictors) == 2 && sum(pred_types == "categorical") == 1 && sum(pred_types == "continuous") == 1) {
        if (verbose_roles) cat("INFERRED PROBLEM: Analysis of Covariance (ANCOVA: Tests Slope Interaction -> Separate Regressions or Line Heights)\n")
        return(invisible(list(prob_type = "ancova", parsed = parsed)))
      } else if (all(pred_types == "continuous")) {
        if (verbose_roles) cat("INFERRED PROBLEM: Multiple Linear Regression + Collinearity Check + All-Subsets AIC Selection\n")
        return(invisible(list(prob_type = "multiple_regression", parsed = parsed)))
      } else {
        if (verbose_roles) cat("INFERRED PROBLEM: General Linear Model (Multi-Factor)\n")
        return(invisible(list(prob_type = "multiple_regression", parsed = parsed)))
      }
    }
  }

  invisible(list(prob_type = "custom", parsed = parsed))
}

# ------------------------------------------------------------------------------
# 14. CONTINGENCY TABLE SOLVER (CHI-SQ, G-TEST, FISHER, MCNEMAR, OR, RR)
# ------------------------------------------------------------------------------
biostat_contingency <- function(df, var1 = NULL, var2 = NULL, paired = FALSE) {
  cat("\n--- CONTINGENCY TABLE & EPIDEMIOLOGICAL ANALYSIS ---\n")
  if (is.table(df) || is.matrix(df)) {
    tbl <- as.table(df)
  } else {
    parsed <- .auto_parse_roles(df, response = var1, predictor = var2, verbose = FALSE)
    if (parsed$mode == "table") {
      tbl <- parsed$table
    } else {
      v1 <- parsed$response; v2 <- parsed$predictor[1]
      cat(sprintf("Auto-Selected Variables: '%s' (Rows) by '%s' (Columns)\n\n", v1, v2))
      tbl <- table(parsed$df[[v1]], parsed$df[[v2]])
    }
  }

  cat("Observed Counts:\n")
  print(addmargins(tbl))

  sq <- suppressWarnings(chisq.test(tbl, correct = FALSE))
  cat("\nExpected Counts under Null Hypothesis (Independence):\n")
  print(round(sq$expected, 2))

  obs <- as.numeric(tbl); exp_c <- as.numeric(sq$expected)
  keep <- obs > 0 & exp_c > 0
  G_stat <- 2 * sum(obs[keep] * log(obs[keep] / exp_c[keep]))
  G_df   <- (nrow(tbl) - 1) * (ncol(tbl) - 1)
  G_p    <- pchisq(G_stat, df = G_df, lower.tail = FALSE)

  min_exp <- min(sq$expected)
  if (min_exp < 5) {
    cat(sprintf("\n[DECISION]: Minimum expected count is %.2f (< 5). Running Fisher's Exact Test:\n", min_exp))
    print(fisher.test(tbl, simulate.p.value = (nrow(tbl) > 2 || ncol(tbl) > 2)))
  } else {
    cat(sprintf("\n[DECISION]: All expected counts >= 5 (min = %.2f). Running Pearson's Chi-Square & G-Test:\n", min_exp))
    print(sq)
    cat(sprintf("Log-Likelihood Ratio G-Test: G = %.4f, df = %d, p-value = %.5g\n", G_stat, G_df, G_p))
  }

  if (all(dim(tbl) == c(2, 2))) {
    a <- tbl[1, 1]; b <- tbl[1, 2]; c <- tbl[2, 1]; d <- tbl[2, 2]
    if (any(tbl == 0)) { a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5 }
    or_val <- (a * d) / (b * c)
    or_ci <- exp(log(or_val) + c(-1, 1) * 1.96 * sqrt(1/a + 1/b + 1/c + 1/d))
    r1 <- a / (a + b); r2 <- c / (c + d)
    rr_val <- r1 / r2
    rr_ci <- exp(log(rr_val) + c(-1, 1) * 1.96 * sqrt((1 - r1)/a + (1 - r2)/c))
    cat("--- 2x2 Effect Sizes ---\n")
    cat(sprintf("Odds Ratio (OR)    : %.3f (95%% CI: %.3f to %.3f)\n", or_val, or_ci[1], or_ci[2]))
    cat(sprintf("Relative Risk (RR) : %.3f (95%% CI: %.3f to %.3f)\n", rr_val, rr_ci[1], rr_ci[2]))
  }
  invisible(tbl)
}

# ------------------------------------------------------------------------------
# 15. SMART PLOTTER (SUPPORTS RAW OR TRANSFORMED LABELS!)
# ------------------------------------------------------------------------------
biostat_plot <- function(df, response = NULL, predictor = NULL, paired = NULL,
                         y_label = NULL, x_label = NULL, verbose_roles = TRUE) {
  parsed <- .auto_parse_roles(df, response = response, predictor = predictor, paired = paired, verbose = verbose_roles)
  if (parsed$mode == "table") {
    mosaicplot(parsed$table, main = "Contingency Mosaic Plot", color = TRUE)
    return(invisible(NULL))
  }

  df_use <- parsed$df
  response <- parsed$response
  predictors <- parsed$predictor
  predictor <- predictors[1]
  paired <- parsed$paired
  block <- parsed$block

  y_lab <- ifelse(is.null(y_label), response, y_label)
  x_lab <- ifelse(is.null(x_label), ifelse(is.null(predictor), "", predictor), x_label)

  if (is.null(predictor)) {
    y <- na.omit(df_use[[response]])
    old_par <- par(no.readonly = TRUE); on.exit(par(old_par)); par(mfrow = c(1, 2))
    hist(y, freq = FALSE, main = paste("Distribution of", y_lab), xlab = y_lab, col = "lightblue", border = "white")
    lines(density(y), col = "darkblue", lwd = 2)
    qqnorm(y, main = paste("Normal Q-Q Plot:", y_lab), pch = 19, col = "steelblue"); qqline(y, col = "firebrick", lwd = 2)
    return(invisible(NULL))
  }

  sub_df <- df_use[complete.cases(df_use[, c(response, predictor)]), ]
  y <- sub_df[[response]]
  x <- sub_df[[predictor]]

  is_cat_x <- is.factor(x) || is.character(x) || (is.numeric(x) && length(unique(x)) <= 8)
  is_cat_y <- is.factor(y) || is.character(y) || (is.numeric(y) && length(unique(y)) <= 2)

  if (.has_pkg("ggplot2")) {
    if (!is_cat_y && is_cat_x && isTRUE(paired) && !is.null(block) && block %in% names(sub_df)) {
      sub_df[[predictor]] <- as.factor(sub_df[[predictor]])
      p <- ggplot2::ggplot(sub_df, ggplot2::aes(x = .data[[predictor]], y = .data[[response]])) +
        ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[predictor]]), alpha = 0.35, width = 0.4, outlier.shape = NA) +
        ggplot2::geom_line(ggplot2::aes(group = .data[[block]]), color = "gray45", alpha = 0.6, linewidth = 0.7) +
        ggplot2::geom_point(ggplot2::aes(color = .data[[predictor]]), size = 2.8, alpha = 0.9) +
        ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 4, fill = "white", color = "black") +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::labs(title = paste("Paired / Blocked Comparison:", y_lab, "by", x_lab),
                      subtitle = paste("Gray lines connect matched", block, "| White Diamond = Mean | Box = Median/IQR"),
                      x = x_lab, y = y_lab) +
        ggplot2::theme(legend.position = "none")
      print(p)
      return(invisible(p))
    } else if (!is_cat_y && is_cat_x) {
      sub_df[[predictor]] <- as.factor(sub_df[[predictor]])
      p <- ggplot2::ggplot(sub_df, ggplot2::aes(x = .data[[predictor]], y = .data[[response]], fill = .data[[predictor]])) +
        ggplot2::geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.5) +
        ggplot2::geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
        ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5, fill = "white") +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::labs(title = paste(y_lab, "by", x_lab),
                      subtitle = "Boxes = Median/IQR | White Diamond = Group Mean | Points = Raw Data",
                      x = x_lab, y = y_lab) +
        ggplot2::theme(legend.position = "none")
      print(p)
      return(invisible(p))
    } else if (!is_cat_y && !is_cat_x) {
      p <- ggplot2::ggplot(sub_df, ggplot2::aes(x = .data[[predictor]], y = .data[[response]])) +
        ggplot2::geom_point(alpha = 0.75, size = 2.5, color = "steelblue") +
        ggplot2::geom_smooth(method = "lm", se = TRUE, color = "firebrick", fill = "pink", alpha = 0.3) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::labs(title = paste("Linear Regression:", y_lab, "vs.", x_lab),
                      subtitle = "Red line = OLS fit with 95% confidence band",
                      x = x_lab, y = y_lab)
      print(p)
      return(invisible(p))
    }
  }

  if (!is_cat_y && is_cat_x) {
    x_fac <- as.factor(x)
    boxplot(y ~ x_fac, main = paste(y_lab, "by", x_lab), xlab = x_lab, ylab = y_lab, col = "lightblue", outline = FALSE)
    stripchart(y ~ x_fac, vertical = TRUE, method = "jitter", add = TRUE, pch = 19, col = "darkblue")
  } else if (!is_cat_y && !is_cat_x) {
    plot(x, y, main = paste(y_lab, "vs.", x_lab), xlab = x_lab, ylab = y_lab, pch = 19, col = "steelblue")
    abline(lm(y ~ x), col = "firebrick", lwd = 2)
  } else {
    mosaicplot(table(x, y), main = paste("Mosaic Plot:", x_lab, "vs", y_lab), xlab = x_lab, ylab = y_lab, color = TRUE)
  }
}

# ------------------------------------------------------------------------------
# 16. MASTER WIZARD: ONE-LINE EXECUTION + AUTO-TRANSFORMATION RESCUE + CODE
# ------------------------------------------------------------------------------
# Hierarchy for Parametric vs. Transformation vs. Non-Parametric:
#   Step 1: Check assumptions on RAW data.
#           - If Raw passes -> Run Parametric test on Raw data!
#   Step 2: If Raw fails normality or homoscedasticity (and transform = "auto"):
#           - Run `biostat_transform()` to test log10, ln, sqrt, & arcsine-sqrt!
#           - If a transformation RESCUES normality & variance -> Execute the
#             Parametric test on the TRANSFORMED data (and also print Raw +
#             Non-Parametric results below it for complete exam security!),
#             and plot the transformed data!
#           - If NO transformation rescues normality -> Execute the Non-Parametric
#             test (Mann-Whitney / Wilcoxon / Kruskal-Wallis / Spearman), while
#             still showing the best transformation attempt!
# ------------------------------------------------------------------------------
biostat_wizard <- function(df, response = NULL, predictor = NULL, paired = NULL,
                           model = "auto", random = NULL, transform = "auto",
                           question = NULL, mu = 0, alpha = 0.05, plot = TRUE, show_code = FALSE) {
  user_df_name <- paste(deparse(substitute(df)), collapse = "")
  inf <- biostat_infer(df, response = response, predictors = predictor,
                       paired = paired, model = model, random = random,
                       question = question, transform = transform, mu = mu, verbose_roles = TRUE)
  prob_type <- inf$prob_type
  parsed <- inf$parsed
  if (!is.null(parsed$transform) && transform == "auto" && parsed$transform != "auto") {
    transform <- parsed$transform
  }

  if (prob_type == "contingency_table") {
    if (plot) biostat_plot(parsed$table, verbose_roles = FALSE)
    return(biostat_contingency(parsed$table, paired = isTRUE(paired)))
  }
  if (prob_type == "logistic") {
    return(biostat_logistic(parsed$df, response = parsed$response, predictor = parsed$predictor[1],
                            plot = plot, show_code = show_code))
  }
  if (prob_type == "factorial_anova") {
    return(biostat_factorial(parsed$df, response = parsed$response,
                             factor1 = parsed$predictor[1], factor2 = parsed$predictor[2],
                             transform = transform, alpha = alpha, plot = plot, show_code = show_code))
  }
  if (prob_type == "model3_anova") {
    return(biostat_model3(parsed$df, response = parsed$response,
                          fixed = setdiff(parsed$predictor, parsed$random)[1],
                          random = parsed$random, transform = transform,
                          alpha = alpha, plot = plot, show_code = show_code))
  }
  if (prob_type == "ancova") {
    return(biostat_ancova(parsed$df, response = parsed$response,
                          covariate = parsed$predictor[1], group = parsed$predictor[2],
                          transform = transform, alpha = alpha, plot = plot, show_code = show_code))
  }
  if (prob_type == "multiple_regression") {
    return(biostat_multireg(parsed$df, response = parsed$response,
                            predictors = parsed$predictor, transform = transform,
                            alpha = alpha, plot = plot, show_code = show_code))
  }

  df_use <- parsed$df
  response <- parsed$response
  predictors <- parsed$predictor
  predictor <- predictors[1]
  paired <- parsed$paired
  block <- parsed$block

  # Step 1: Check assumptions on Raw Data
  assump <- biostat_assumptions(df_use, response = response, predictor = predictor,
                                paired = paired, alpha = alpha, verbose_roles = FALSE)

  # Step 2: Check Data Transformations BEFORE jumping to Non-Parametric!
  active_df <- df_use
  y_label_use <- response
  x_label_use <- predictor
  used_transform <- FALSE
  tr_res <- NULL

  if (!transform %in% c("auto", "none", "raw")) {
    # User explicitly requested a specific transformation (e.g., transform = "log10" or "sqrt")
    tr_app <- .apply_transform_df(df_use, response, predictor, transform)
    active_df <- tr_app$df
    y_label_use <- tr_app$y_label
    x_label_use <- tr_app$x_label
    used_transform <- TRUE
    cat(sprintf("\n[USER-SPECIFIED TRANSFORMATION APPLIED]: %s\n", y_label_use))
    assump <- biostat_assumptions(active_df, response = response, predictor = predictor,
                                  paired = paired, alpha = alpha, verbose_roles = FALSE)

  } else if (transform == "auto" && is.numeric(df_use[[response]]) &&
             (!isTRUE(assump$normality_ok) || !isTRUE(assump$variance_ok) ||
              (prob_type == "continuous_xy" && !isTRUE(assump$bivariate_normal)))) {
    # Raw data failed an assumption -> Run Transformation Rescue Engine!
    tr_res <- tryCatch(biostat_transform(df_use, response = response, predictor = predictor,
                                         paired = paired, alpha = alpha, plot = FALSE, verbose = TRUE),
                       error = function(e) NULL)
    if (!is.null(tr_res) && isTRUE(tr_res$rescued) && isTRUE(tr_res$best_ev$norm_ok)) {
      active_df <- tr_res$df_transformed
      y_label_use <- tr_res$best_label
      if (!is.null(tr_res$best_x_label)) x_label_use <- tr_res$best_x_label
      assump$normality_ok <- tr_res$best_ev$norm_ok
      assump$variance_ok  <- tr_res$best_ev$var_ok
      assump$parametric_ok <- tr_res$best_ev$both_ok
      assump$bivariate_normal <- tr_res$best_ev$norm_ok
      used_transform <- TRUE
    }
  }

  cat("\n======================================================================\n")
  cat("                 AUTOMATED STATISTICAL EXECUTION                      \n")
  cat("======================================================================\n")

  if (plot && length(predictors) <= 1) {
    # Plot 1: 3-Panel Assumption Diagnostic Grid (Histogram + Normal Curve, Professor's Shaded QQ Plot, Residuals vs Fitted)
    if (is.numeric(active_df[[response]]) && prob_type %in% c("indep_2group", "indep_kgroup", "continuous_xy", "paired_2group")) {
      tryCatch({
        biostat_resid_plot(active_df, response = response, predictor = predictor)
        cat("[PLOTS]: Generated 2 plots! Click the blue '<' back arrow in RStudio's Plots pane to view the\n")
        cat("         3-Panel Assumption Grid (Histogram+Normal, Shaded QQ Envelope, Residuals vs. Fitted).\n\n")
      }, error = function(e) NULL)
    }
    # Plot 2: Final Publication-Ready Results Graph (shown immediately on screen)
    tryCatch(biostat_plot(active_df, response = response, predictor = predictor,
                          paired = paired, y_label = y_label_use, x_label = x_label_use, verbose_roles = FALSE),
             error = function(e) NULL)
  }

  # 1. One-Sample Problem
  if (is.null(predictor)) {
    y <- na.omit(active_df[[response]])
    if (is.numeric(y)) {
      if (assump$parametric_ok) {
        cat(sprintf("-> Executing Parametric One-Sample t-test on %s (Normality PASS):\n", y_label_use))
        res <- t.test(y, mu = mu)
      } else {
        cat(sprintf("-> Transformations did not normalize data. Executing Non-Parametric One-Sample Wilcoxon Test:\n"))
        res <- wilcox.test(df_use[[response]], mu = mu, exact = FALSE)
      }
    } else {
      cat("-> Executing Chi-Square Goodness-of-Fit Test:\n")
      res <- chisq.test(table(y))
    }
    print(res)
    return(invisible(res))
  }

  sub_df <- active_df[complete.cases(active_df[, c(response, predictors)]), ]
  y <- sub_df[[response]]
  x <- sub_df[[predictor]]

  # 2. Independent 2-Group Problem
  if (prob_type == "indep_2group") {
    x_fac <- as.factor(x)
    g1 <- y[x_fac == levels(x_fac)[1]]
    g2 <- y[x_fac == levels(x_fac)[2]]
    s_pool <- sqrt(((length(g1) - 1)*var(g1) + (length(g2) - 1)*var(g2)) / (length(g1) + length(g2) - 2))
    cohens_d <- (mean(g1) - mean(g2)) / s_pool

    cat(sprintf("Scale Analyzed : %s\n", y_label_use))
    cat(sprintf("Group Summaries: %s (n=%d, Mean=%.3f, SD=%.3f, SEM=%.3f) vs. %s (n=%d, Mean=%.3f, SD=%.3f, SEM=%.3f)\n\n",
                levels(x_fac)[1], length(g1), mean(g1), sd(g1), sd(g1)/sqrt(length(g1)),
                levels(x_fac)[2], length(g2), mean(g2), sd(g2), sd(g2)/sqrt(length(g2))))

    if (assump$normality_ok) {
      var_eq <- assump$variance_ok
      cat(sprintf("-> Normality PASS (%s). Variance equality %s. Running %s t-test:\n",
                  ifelse(used_transform, paste("after", y_label_use, "transform"), "on raw data"),
                  ifelse(var_eq, "PASS", "FAIL"),
                  ifelse(var_eq, "Student's (var.equal = TRUE)", "Welch's (var.equal = FALSE)")))
      res <- t.test(y ~ x_fac, var.equal = var_eq)
      print(res)
      if (used_transform) {
        cat("(Non-parametric Mann-Whitney U on raw data also shown for reference):\n")
        print(wilcox.test(df_use[[response]] ~ as.factor(df_use[[predictor]]), exact = FALSE))
      }
    } else {
      cat("-> Neither raw data nor basic transformations met normality. Running Mann-Whitney U Test:\n")
      res <- wilcox.test(df_use[[response]] ~ as.factor(df_use[[predictor]]), exact = FALSE)
      print(res)
      cat("(Welch's t-test on raw data also shown for reference):\n")
      print(t.test(df_use[[response]] ~ as.factor(df_use[[predictor]]), var.equal = FALSE))
    }
    cat(sprintf("Effect Size -> Cohen's d = %.3f\n", cohens_d))
    .biostat_last_audit$last_code <- biostat_code(df, response = response, predictor = predictor, df_name = user_df_name, silent = !isTRUE(show_code))
    if (!isTRUE(show_code)) cat("\n[TIP]: Run `biostat_code()` anytime to print the standalone R code & graph code for this analysis.\n")
    return(invisible(res))

  # 3. Paired 2-Group Problem
  } else if (prob_type == "paired_2group") {
    x_fac <- as.factor(x)
    g1 <- y[x_fac == levels(x_fac)[1]]
    g2 <- y[x_fac == levels(x_fac)[2]]
    diffs <- g1 - g2
    dz <- mean(diffs) / sd(diffs)

    cat(sprintf("Scale Analyzed : %s\n", y_label_use))
    cat(sprintf("Condition 1 (%s): n = %d, Mean = %.3f, SD = %.3f, Median = %.3f\n",
                levels(x_fac)[1], length(g1), mean(g1), sd(g1), median(g1)))
    cat(sprintf("Condition 2 (%s): n = %d, Mean = %.3f, SD = %.3f, Median = %.3f\n",
                levels(x_fac)[2], length(g2), mean(g2), sd(g2), median(g2)))
    cat(sprintf("Paired Diff (%s - %s): Mean = %.3f, SD = %.3f, SEM = %.3f\n\n",
                levels(x_fac)[1], levels(x_fac)[2], mean(diffs), sd(diffs), sd(diffs)/sqrt(length(diffs))))

    if (assump$normality_ok) {
      cat(sprintf("-> Normality of Paired Differences PASS (%s). Running Parametric Paired t-test:\n",
                  ifelse(used_transform, paste("after", y_label_use, "transform"), "on raw data")))
      res <- t.test(g1, g2, paired = TRUE)
      print(res)
      cat("(Non-parametric Wilcoxon Signed-Rank also shown for reference):\n")
      print(suppressWarnings(wilcox.test(g1, g2, paired = TRUE, exact = FALSE)))
    } else {
      cat("-> Paired differences non-normal even after transformation check. Running Wilcoxon Signed-Rank Test:\n")
      res <- suppressWarnings(wilcox.test(g1, g2, paired = TRUE, exact = FALSE))
      print(res)
      cat("(Parametric Paired t-test also shown for reference):\n")
      print(t.test(g1, g2, paired = TRUE))
    }
    cat(sprintf("Paired Effect Size -> Cohen's dz = %.3f\n", dz))
    .biostat_last_audit$last_code <- biostat_code(df, response = response, predictor = predictor, paired = TRUE, df_name = user_df_name, silent = !isTRUE(show_code))
    if (!isTRUE(show_code)) cat("\n[TIP]: Run `biostat_code()` anytime to print the standalone R code & graph code for this analysis.\n")
    return(invisible(res))

  # 4. Paired / Blocked >2 Groups
  } else if (prob_type == "paired_kgroup") {
    x_fac <- as.factor(x); b_fac <- as.factor(sub_df[[block]])
    cat("-> Running Parametric Randomized Complete Block / Repeated-Measures ANOVA:\n")
    fit_rcbd <- aov(y ~ x_fac + b_fac)
    print(summary(fit_rcbd))
    cat("\n--- Post-Hoc: Tukey's HSD across Conditions ---\n")
    print(TukeyHSD(fit_rcbd, which = "x_fac"))
    cat("\n-> Running Non-Parametric Friedman Rank Sum Test (for reference):\n")
    print(friedman.test(y ~ x_fac | b_fac))
    return(invisible(fit_rcbd))

  # 5. Independent >2 Groups (One-Way ANOVA / Transform / Welch / Kruskal-Wallis)
  } else if (prob_type == "indep_kgroup") {
    x_fac <- as.factor(x)
    fit_aov <- aov(y ~ x_fac)
    ss <- summary(fit_aov)[[1]]
    eta_sq <- ss[1, "Sum Sq"] / sum(ss[, "Sum Sq"])

    cat(sprintf("Scale Analyzed : %s\n", y_label_use))
    if (assump$normality_ok && assump$variance_ok) {
      cat(sprintf("-> Assumptions PASS (%s). Running Parametric One-Way ANOVA + Tukey HSD Post-Hoc:\n",
                  ifelse(used_transform, paste("after", y_label_use, "transform"), "on raw data")))
      print(summary(fit_aov))
      cat(sprintf("\nEffect Size -> Eta-squared = %.3f\n", eta_sq))
      if (ss[1, "Pr(>F)"] < alpha) {
        cat("\n--- Post-Hoc: Tukey's HSD ---\n")
        print(TukeyHSD(fit_aov))
      }
      res <- fit_aov
    } else if (assump$normality_ok && !assump$variance_ok) {
      cat("-> Normality PASS, Equal Variance FAIL. Running Welch's One-Way ANOVA:\n")
      res <- oneway.test(y ~ x_fac, var.equal = FALSE)
      print(res)
      cat("\n--- Post-Hoc: Pairwise Welch t-tests (Holm-adjusted) ---\n")
      print(pairwise.t.test(y, x_fac, pool.sd = FALSE, p.adjust.method = "holm"))
    } else {
      cat("-> Normality FAIL (even after checking log/sqrt transformations). Running Kruskal-Wallis Test:\n")
      res <- kruskal.test(df_use[[response]] ~ as.factor(df_use[[predictor]]))
      print(res)
      if (res$p.value < alpha) {
        cat("\n--- Post-Hoc: Pairwise Wilcoxon Rank-Sum Tests (Holm-adjusted) ---\n")
        print(pairwise.wilcox.test(df_use[[response]], as.factor(df_use[[predictor]]), p.adjust.method = "holm", exact = FALSE))
      }
    }
    .biostat_last_audit$last_code <- biostat_code(df, response = response, predictor = predictor, df_name = user_df_name, silent = !isTRUE(show_code))
    if (!isTRUE(show_code)) cat("\n[TIP]: Run `biostat_code()` anytime to print the standalone R code & graph code for this analysis.\n")
    return(invisible(res))

  # 6. Continuous Y vs Continuous X (Simple Linear Regression & Correlation)
  } else if (prob_type == "continuous_xy") {
    cat(sprintf("-> Fitting Simple Linear Regression (OLS: %s ~ %s):\n", y_label_use, x_label_use))
    fit_lm <- lm(y ~ x)
    s_lm   <- summary(fit_lm)
    b0 <- coef(fit_lm)[1]; b1 <- coef(fit_lm)[2]
    r_val <- cor(x, y, use = "complete.obs")
    r2_val <- s_lm$r.squared; adj_r2 <- s_lm$adj.r.squared
    print(anova(fit_lm))
    print(s_lm)
    cat("\n--- BEST-FIT LINE EQUATION & STRENGTH OF ASSOCIATION ---\n")
    cat(sprintf("  * Best-Fit Line Equation : %s = %.4f %s %.4f * (%s)\n",
                y_label_use, b0, ifelse(b1 >= 0, "+", "-"), abs(b1), x_label_use))
    cat(sprintf("  * Pearson r              : r = %.4f -> %s\n", r_val, .interpret_r(r_val)))
    cat(sprintf("  * R-Squared (R2)         : R2 = %.4f (Adjusted R2 = %.4f)\n", r2_val, adj_r2))
    cat(sprintf("    -> INTERPRETATION: %.1f%% of the variance in '%s' is explained by '%s'.\n",
                r2_val * 100, y_label_use, x_label_use))

    if (assump$bivariate_normal) {
      cat(sprintf("\n-> Bivariate Normality PASS (%s). Running Pearson Correlation (r):\n",
                  ifelse(used_transform, "on transformed scale", "on raw scale")))
      print(cor.test(x, y, method = "pearson"))
    } else {
      cat("\n-> Bivariate Normality FAIL. Running Spearman Rank Correlation (rho):\n")
      print(suppressWarnings(cor.test(df_use[[predictor]], df_use[[response]], method = "spearman", exact = FALSE)))
      cat("(Pearson Correlation also shown for reference):\n")
      print(cor.test(x, y, method = "pearson"))
    }
    .biostat_last_audit$last_code <- biostat_code(df, response = response, predictor = predictor, df_name = user_df_name, silent = !isTRUE(show_code))
    if (!isTRUE(show_code)) cat("\n[TIP]: Run `biostat_code()` anytime to print the standalone R code & graph code for this analysis.\n")
    return(invisible(fit_lm))

  # 7. Contingency Table
  } else if (prob_type == "contingency") {
    return(biostat_contingency(sub_df, response, predictor, paired = isTRUE(paired)))
  }
}

# ------------------------------------------------------------------------------
# 17. SUMMARY STATISTICS & POWER SOLVERS
# ------------------------------------------------------------------------------
biostat_from_summary <- function(mean1, sd1, n1, mean2 = NULL, sd2 = NULL, n2 = NULL,
                                 mu = 0, var.equal = FALSE, conf.level = 0.95) {
  cat("\n--- SUMMARY-STATISTICS PARAMETRIC SOLVER ---\n")
  alpha <- 1 - conf.level
  if (is.null(mean2)) {
    se <- sd1 / sqrt(n1); df <- n1 - 1; t_val <- (mean1 - mu) / se
    p_val <- 2 * pt(-abs(t_val), df = df)
    ci <- mean1 + c(-1, 1) * qt(1 - alpha/2, df = df) * se
    cat(sprintf("One-Sample t-test (Null mu = %g): t = %.4f, df = %d, p = %.5g | 95%% CI: [%.4f, %.4f]\n",
                mu, t_val, df, p_val, ci[1], ci[2]))
    return(invisible(list(t = t_val, df = df, p.value = p_val, ci = ci)))
  } else {
    diff_m <- mean1 - mean2
    if (var.equal) {
      s_pool <- sqrt(((n1 - 1)*sd1^2 + (n2 - 1)*sd2^2) / (n1 + n2 - 2))
      se <- s_pool * sqrt(1/n1 + 1/n2); df <- n1 + n2 - 2
    } else {
      v1 <- sd1^2 / n1; v2 <- sd2^2 / n2; se <- sqrt(v1 + v2)
      df <- (v1 + v2)^2 / ((v1^2 / (n1 - 1)) + (v2^2 / (n2 - 1)))
    }
    t_val <- diff_m / se; p_val <- 2 * pt(-abs(t_val), df = df)
    ci <- diff_m + c(-1, 1) * qt(1 - alpha/2, df = df) * se
    cat(sprintf("Two-Sample t-test: t = %.4f, df = %.2f, p = %.5g | Diff = %.4f (95%% CI: [%.4f, %.4f])\n",
                t_val, df, p_val, diff_m, ci[1], ci[2]))
    return(invisible(list(t = t_val, df = df, p.value = p_val, ci = ci)))
  }
}

biostat_power_calc <- function(test_type = c("t.two", "t.paired", "t.one", "anova"),
                               effect_size = 0.5, n = NULL, power = 0.80, alpha = 0.05, k_groups = 3) {
  test_type <- match.arg(test_type)
  if (test_type %in% c("t.two", "t.paired", "t.one")) {
    type_map <- c("t.two" = "two.sample", "t.paired" = "paired", "t.one" = "one.sample")
    res <- power.t.test(n = n, delta = effect_size, sd = 1, sig.level = alpha, power = power, type = type_map[test_type])
    print(res); return(invisible(res))
  } else if (test_type == "anova") {
    between_var <- (effect_size^2) * (k_groups / (k_groups - 1))
    res <- power.anova.test(groups = k_groups, n = n, between.var = between_var, within.var = 1, sig.level = alpha, power = power)
    print(res); return(invisible(res))
  }
}


# ------------------------------------------------------------------------------
# 18. QUESTION-DRIVEN & DERIVED-COLUMN HELPERS (`biostat_ask` & `biostat_derive`)
# ------------------------------------------------------------------------------
# Paste the homework/exam question prompt directly in plain English!
biostat_ask <- function(df, question, plot = TRUE, show_code = FALSE) {
  biostat_wizard(df, question = question, plot = plot, show_code = show_code)
}

# Create a derived column (like Height * Width) AND automatically drop raw components
biostat_derive <- function(df, expr_or_question) {
  q_res <- .parse_question_and_derived(df, response = expr_or_question, verbose = TRUE)
  q_res$df
}


# ------------------------------------------------------------------------------
# 19. PROFESSOR-STYLE RESIDUAL DIAGNOSTIC PLOTTER (`biostat_resid_plot`)
# ------------------------------------------------------------------------------
# Reproduces the exact `qqPlot(resid(model))` + `plot(resid(model) ~ fitted(model))`
# wedge-shape / homoscedasticity check from your ANCOVA & Factorial ANOVA scripts!
biostat_resid_plot <- function(df_or_model, response = NULL, predictor = NULL, transform = "none") {
  if (inherits(df_or_model, c("lm", "aov", "glm"))) {
    fit <- df_or_model
    fml_str <- deparse(formula(fit))
  } else {
    parsed <- .auto_parse_roles(df_or_model, response = response, predictor = predictor, verbose = FALSE)
    df_use <- parsed$df
    resp   <- parsed$response
    preds  <- parsed$predictor
    if (!transform %in% c("none", "raw", "auto")) {
      tr_app <- .apply_transform_df(df_use, resp, preds, transform)
      df_use <- tr_app$df
    }
    if (length(preds) >= 2) {
      fml <- as.formula(paste(.sym_q(resp), "~", paste(sapply(preds[1:2], .sym_q), collapse = " * ")))
    } else if (length(preds) == 1) {
      fml <- as.formula(paste(.sym_q(resp), "~", .sym_q(preds[1])))
    } else {
      fml <- as.formula(paste(.sym_q(resp), "~ 1"))
    }
    fit <- lm(fml, data = df_use)
    fml_str <- deparse(fml)
  }

  res <- resid(fit)
  fit_vals <- fitted(fit)
  sw <- if (length(res) >= 3 && length(res) <= 5000) shapiro.test(res) else list(statistic = NA, p.value = NA)

  sk_ku <- .calc_skew_kurt(res)
  old_par <- par(no.readonly = TRUE); on.exit(par(old_par))
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.8, 1))

  # Panel 1: Histogram with Normal Bell Curve Overlay
  .plot_hist_normal(res,
                    main_title = sprintf("1. Histogram vs. Normal Curve\nSkew = %.2f | Kurt = %.2f", sk_ku$skew, sk_ku$kurt),
                    xlab_str = "Model Residuals")

  # Panel 2: Professor-style car::qqPlot (red center line + shaded 95% confidence containment!)
  .plot_prof_qq(res,
                main_title = sprintf("2. QQ Plot w/ 95%% Envelope\nShapiro W = %.3f, p = %.4g", sw$statistic, sw$p.value),
                ylab_str = "Residuals")

  # Panel 3: Residuals vs. Fitted Values (Professor's wedge-shape homoscedasticity check!)
  plot(res ~ fit_vals, pch = 19, col = "steelblue4",
       xlab = "Predicted / Fitted Y Values", ylab = "Residuals",
       main = "3. Homoscedasticity Check:\nResiduals vs. Fitted Y (Wedge Check)")
  abline(h = 0, col = "firebrick", lwd = 2, lty = 2)
  invisible(fit)
}

message("Universal Gbiostats_toolkit.r (v10 - Clean Wizard + Dual Plots + On-Demand Code) loaded!")