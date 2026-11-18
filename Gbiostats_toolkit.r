# ==============================================================================
# Gbiostats_toolkit.R - ZERO-FRICTION AUTO-INFERRING BIOSTATISTICS ENGINE
# ==============================================================================
# Designed so you can pass JUST a dataset (or vector/table) with ZERO column 
# names required:
#   biostat_wizard(df)       # Auto-detects Y, X, IDs, paired structure & runs all!
#   biostat_infer(df)        # Infers roles & prints exact parametric/non-parametric tests
#   biostat_assumptions(df)  # Runs Shapiro-Wilk, Levene/Fligner, & outlier checks
#   biostat_plot(df)         # Auto-plots the inferred relationship
#   biostat_inspect(df)      # Classifies all variables & shows guessed roles
# ==============================================================================

.has_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

# ------------------------------------------------------------------------------
# 0. CORE ENGINE: AUTOMATIC VARIABLE ROLE & STRUCTURE DETECTOR
# ------------------------------------------------------------------------------
# Accepts:
#   - Just a dataframe: biostat_wizard(df) -> guesses Response (Y), Predictor (X), & ID/Block
#   - A formula:        biostat_wizard(y ~ x, data = df) or biostat_wizard(df, y ~ x)
#   - A single vector:  biostat_wizard(df$y)
#   - A 2D table/matrix:biostat_wizard(my_table)
#   - Wide paired data: 2 numeric columns with similar names/scales (e.g., Before, After)
# ------------------------------------------------------------------------------
.auto_parse_roles <- function(data, response = NULL, predictor = NULL, block = NULL, paired = NULL, verbose = TRUE) {
  # 1. Handle Formula input as first argument: e.g., biostat_wizard(Sepal.Length ~ Species, iris)
  if (inherits(data, "formula")) {
    fml <- data
    if (is.null(response) || !is.data.frame(response)) {
      stop("When passing a formula first, provide the dataframe as the second argument: e.g., biostat_wizard(y ~ x, df)")
    }
    df <- as.data.frame(response)
    vars <- all.vars(fml)
    response <- vars[1]
    predictor <- if (length(vars) > 1) vars[-1] else NULL
    data <- df
  } else if (inherits(response, "formula")) {
    # Handle biostat_wizard(iris, Sepal.Length ~ Species)
    fml <- response
    vars <- all.vars(fml)
    response <- vars[1]
    predictor <- if (length(vars) > 1) vars[-1] else NULL
  }

  # 2. Handle Matrix / Contingency Table directly
  if (is.table(data) || (is.matrix(data) && is.numeric(data) && nrow(data) <= 15 && ncol(data) <= 15)) {
    if (verbose) cat("[AUTO-DETECT]: Input is a pre-tabulated Contingency Table / Matrix.\n")
    return(list(mode = "table", table = as.table(data)))
  }

  # 3. Handle Atomic Vector directly: e.g., biostat_wizard(iris$Sepal.Length)
  if (is.atomic(data) && is.null(dim(data))) {
    vec_name <- deparse(substitute(data))
    df <- data.frame(Value = data)
    names(df) <- make.names(vec_name)
    return(list(
      mode = "dataframe", df = df,
      response = names(df)[1], predictor = NULL, extra_predictors = NULL,
      block = NULL, paired = FALSE, wide_paired = FALSE
    ))
  }

  df <- as.data.frame(data)
  cols <- names(df)
  n_rows <- nrow(df)

  # Helper to classify each column
  col_meta <- lapply(cols, function(cn) {
    v <- df[[cn]]
    v_clean <- na.omit(v)
    n_u <- length(unique(v_clean))
    n_clean <- length(v_clean)

    # Detect ID / Subject / Block column by name or uniqueness
    is_id_name <- grepl("^(id|sub|subject|patient|sample|specimen|block|tank|plot|site|individual|mouse|rat|rep|replicate)s?(_id|_num|_no)?$",
                        cn, ignore.case = TRUE)
    is_unique_id <- (n_u == n_clean && n_clean > 2 && (is.character(v) || is.factor(v) || (is.numeric(v) && all(v_clean == floor(v_clean)) && all(diff(sort(v_clean)) == 1))))

    if (is_id_name || is_unique_id) {
      role_type <- "id_or_block"
    } else if (is.factor(v) || is.character(v) || is.logical(v)) {
      role_type <- if (n_u == 2) "binary_cat" else if (n_u <= 12) "multi_cat" else "high_card_cat"
    } else if (is.numeric(v)) {
      if (n_u == 2 && all(v_clean %in% c(0, 1))) {
        role_type <- "binary_num"
      } else if (all(v_clean == floor(v_clean)) && min(v_clean) >= 0 && n_u <= 8) {
        role_type <- "small_integer"
      } else {
        role_type <- "continuous"
      }
    } else {
      role_type <- "other"
    }

    # Semantic name cues for Response vs Predictor
    is_resp_name <- grepl("(resp|outcome|yield|growth|weight|mass|length|height|size|score|rate|conc|density|count|survival|mortality|bp|expression|value|^y$)",
                          cn, ignore.case = TRUE)
    is_pred_name <- grepl("(group|treat|trt|trmt|dose|dosage|species|sex|gender|diet|cond|condition|temp|time|day|age|site|habitat|drug|fertilizer|genotype|strain|^x$)",
                          cn, ignore.case = TRUE)

    list(
      name = cn, type = role_type, n_unique = n_u,
      is_id_name = is_id_name, is_resp_name = is_resp_name, is_pred_name = is_pred_name
    )
  })
  names(col_meta) <- cols

  # 4. Auto-detect Block / ID column if not provided
  if (is.null(block)) {
    id_candidates <- cols[sapply(col_meta, function(m) m$type == "id_or_block")]
    if (length(id_candidates) > 0) {
      block <- id_candidates[1]
    }
  }

  # Candidate pool excluding ID/Block columns
  active_cols <- setdiff(cols, block)
  if (length(active_cols) == 0) active_cols <- cols

  # 5. Special Check: 2-Column Wide Paired Data (e.g., df with "Before" and "After" or "Pre" and "Post")
  wide_paired <- FALSE
  if (is.null(response) && is.null(predictor) && length(active_cols) == 2) {
    t1 <- col_meta[[active_cols[1]]]$type
    t2 <- col_meta[[active_cols[2]]]$type
    paired_names <- any(grepl("(pre|post|before|after|baseline|followup|t1|t2|time1|time2|left|right|twin|control|treated)",
                              active_cols, ignore.case = TRUE))
    if (t1 == "continuous" && t2 == "continuous" && paired_names) {
      wide_paired <- TRUE
      # Reshape wide 2-column paired data into long format for seamless analysis
      c1 <- active_cols[1]; c2 <- active_cols[2]
      complete_rows <- complete.cases(df[, c(c1, c2)])
      df_long <- data.frame(
        Pair_ID = factor(rep(seq_len(sum(complete_rows)), 2)),
        Condition = factor(rep(c(c1, c2), each = sum(complete_rows)), levels = c(c1, c2)),
        Measurement = c(df[[c1]][complete_rows], df[[c2]][complete_rows])
      )
      if (verbose) {
        cat(sprintf("[AUTO-INFERRED STRUCTURE]: Wide Paired Dataset detected ('%s' vs '%s').\n", c1, c2))
      }
      return(list(
        mode = "dataframe", df = df_long,
        response = "Measurement", predictor = "Condition", extra_predictors = NULL,
        block = "Pair_ID", paired = TRUE, wide_paired = TRUE, orig_cols = c(c1, c2)
      ))
    }
  }

  # 6. Auto-infer Response (Y) and Predictor (X) if omitted!
  if (is.null(response)) {
    cont_cols <- active_cols[sapply(active_cols, function(c) col_meta[[c]]$type %in% c("continuous", "small_integer"))]
    cat_cols  <- active_cols[sapply(active_cols, function(c) col_meta[[c]]$type %in% c("binary_cat", "multi_cat", "binary_num"))]

    if (length(active_cols) == 1) {
      response <- active_cols[1]
      predictor <- NULL
    } else if (length(cont_cols) == 1 && length(cat_cols) >= 1) {
      # Classic Group Comparison (t-test / ANOVA): 1 Continuous Y, 1+ Categorical X
      response <- cont_cols[1]
      # Pick best categorical predictor (prefer explicit predictor names, otherwise first categorical)
      pred_named <- cat_cols[sapply(cat_cols, function(c) col_meta[[c]]$is_pred_name)]
      predictor <- if (length(pred_named) > 0) pred_named[1] else cat_cols[1]
    } else if (length(cont_cols) >= 2 && length(cat_cols) == 0) {
      # All continuous (Regression / Correlation):
      # Prefer column with response-like name as Y, or predictor-like name as X
      resp_named <- cont_cols[sapply(cont_cols, function(c) col_meta[[c]]$is_resp_name && !col_meta[[c]]$is_pred_name)]
      pred_named <- cont_cols[sapply(cont_cols, function(c) col_meta[[c]]$is_pred_name && !col_meta[[c]]$is_resp_name)]

      if (length(resp_named) >= 1) {
        response <- resp_named[1]
        rem <- setdiff(cont_cols, response)
        predictor <- if (length(pred_named) >= 1 && pred_named[1] != response) pred_named[1] else rem[1]
      } else if (length(pred_named) >= 1) {
        predictor <- pred_named[1]
        rem <- setdiff(cont_cols, predictor)
        response <- rem[1]
      } else {
        # Convention: in 2-column (X, Y) datasets, check if last column is Y or first is Y
        response <- cont_cols[length(cont_cols)]
        predictor <- cont_cols[1]
      }
    } else if (length(cont_cols) == 0 && length(cat_cols) >= 2) {
      # All categorical (Contingency Table / Chi-Square / Fisher)
      resp_named <- cat_cols[sapply(cat_cols, function(c) col_meta[[c]]$is_resp_name)]
      if (length(resp_named) >= 1) {
        response <- resp_named[1]
        predictor <- setdiff(cat_cols, response)[1]
      } else {
        response <- cat_cols[length(cat_cols)]
        predictor <- cat_cols[1]
      }
    } else {
      # Mixed multi-variable dataset (multiple continuous and categorical columns)
      resp_named <- cont_cols[sapply(cont_cols, function(c) col_meta[[c]]$is_resp_name)]
      response <- if (length(resp_named) >= 1) resp_named[1] else cont_cols[1]
      rem_cols <- setdiff(active_cols, response)
      pred_named <- rem_cols[sapply(rem_cols, function(c) col_meta[[c]]$is_pred_name)]
      # Prefer categorical grouping variable first if available, else first remaining
      rem_cat <- intersect(rem_cols, cat_cols)
      if (length(pred_named) >= 1) {
        predictor <- pred_named[1]
      } else if (length(rem_cat) >= 1) {
        predictor <- rem_cat[1]
      } else {
        predictor <- rem_cols[1]
      }
    }
  } else if (is.null(predictor) && length(active_cols) >= 2) {
    # User provided response, auto-infer predictor from remaining active columns
    rem_cols <- setdiff(active_cols, response)
    if (length(rem_cols) >= 1) {
      pred_named <- rem_cols[sapply(rem_cols, function(c) col_meta[[c]]$is_pred_name)]
      predictor <- if (length(pred_named) >= 1) pred_named[1] else rem_cols[1]
    }
  }

  # Identify any additional predictors in the dataset (for multi-factor / ANCOVA / Multiple Regression hints)
  extra_predictors <- setdiff(active_cols, c(response, predictor))

  # 7. Auto-detect Paired / Repeated Measures if block/subject ID has 1 obs per group
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
    cat("AUTO-INFERRED VARIABLE ROLES (Override anytime if needed):\n")
    cat(sprintf("  * Response (Y)       : %s\n", ifelse(is.null(response), "None", response)))
    cat(sprintf("  * Primary Factor (X) : %s\n", ifelse(is.null(predictor), "None (One-Sample)", paste(predictor, collapse = ", "))))
    if (length(extra_predictors) > 0) {
      cat(sprintf("  * Other Covariates   : %s (also in dataset)\n", paste(extra_predictors, collapse = ", ")))
    }
    if (!is.null(block)) {
      cat(sprintf("  * Subject/Block ID   : %s (Paired/Blocked = %s)\n", block, ifelse(paired, "TRUE", "FALSE")))
    }
    cat("----------------------------------------------------------------------\n")
  }

  list(
    mode = "dataframe", df = df,
    response = response, predictor = predictor,
    extra_predictors = extra_predictors,
    block = block, paired = paired, wide_paired = wide_paired,
    col_meta = col_meta
  )
}

# ------------------------------------------------------------------------------
# 1. INTERACTIVE MENU & ZERO-INPUT WORKFLOW GUIDE
# ------------------------------------------------------------------------------
biostat_guide <- function() {
  cat("======================================================================\n")
  cat("       ZERO-FRICTION BIOSTATISTICS TOOLKIT (Gbiostats_toolkit.R)      \n")
  cat("======================================================================\n")
  cat("PASS JUST YOUR DATASET — NO COLUMN NAMES REQUIRED!\n\n")
  cat("1. ONE-COMMAND MASTER SOLVER:\n")
  cat("   - biostat_wizard(df)      : Auto-detects Y & X, checks assumptions, runs\n")
  cat("                               parametric/non-parametric test, post-hoc & plot!\n\n")
  cat("2. STEP-BY-STEP AUTO-INFERENCE FUNCTIONS (All work with just `df`):\n")
  cat("   - biostat_inspect(df)     : Classifies variables, missingness, & guessed roles\n")
  cat("   - biostat_infer(df)       : Infers study design & prints exact R test commands\n")
  cat("   - biostat_assumptions(df) : Runs Shapiro-Wilk, Levene/Fligner, & outlier checks\n")
  cat("   - biostat_design_check(df): Checks sample balance & pseudoreplication\n")
  cat("   - biostat_plot(df)        : Auto-generates publication-ready plot for Y vs X\n")
  cat("   - biostat_contingency(df) : Chi-Square, Fisher's Exact, Odds Ratio & Rel. Risk\n\n")
  cat("3. OPTIONAL OVERRIDES (If your dataset has many columns):\n")
  cat("   - Formula syntax : biostat_wizard(df, Growth ~ Treatment)\n")
  cat("   - Named syntax   : biostat_wizard(df, 'Growth', 'Treatment', paired = TRUE)\n")
  cat("======================================================================\n")
}

# ------------------------------------------------------------------------------
# 2. DATASET INSPECTOR & AUTO-ROLE CLASSIFIER
# ------------------------------------------------------------------------------
biostat_inspect <- function(df) {
  if (is.table(df) || is.matrix(df)) {
    cat("\n--- CONTINGENCY TABLE / MATRIX DETECTED ---\n")
    print(addmargins(as.table(df)))
    return(invisible(df))
  }
  if (is.atomic(df) && is.null(dim(df))) {
    df <- data.frame(Value = df)
  }
  df <- as.data.frame(df)
  parsed <- .auto_parse_roles(df, verbose = FALSE)

  cat("\n--- DATASET OVERVIEW & AUTO-INFERRED ROLES ---\n")
  cat(sprintf("Dimensions: %d observations (rows) x %d variables (columns)\n\n", nrow(df), ncol(df)))

  summary_list <- lapply(names(df), function(col_name) {
    x <- df[[col_name]]
    n_miss <- sum(is.na(x))
    pct_miss <- round(100 * n_miss / length(x), 1)
    n_uniq <- length(unique(na.omit(x)))

    # Role tag
    role_tag <- if (!is.null(parsed$response) && col_name == parsed$response) {
      "-> RESPONSE (Y)"
    } else if (!is.null(parsed$predictor) && col_name %in% parsed$predictor) {
      "-> PREDICTOR (X1)"
    } else if (!is.null(parsed$block) && col_name == parsed$block) {
      "-> ID / BLOCK"
    } else if (col_name %in% parsed$extra_predictors) {
      "Covariate (X+)"
    } else {
      "-"
    }

    if (is.numeric(x)) {
      if (n_uniq == 2 && all(na.omit(x) %in% c(0, 1))) {
        v_type <- "Binary (0/1)"
      } else if (all(na.omit(x) == floor(na.omit(x))) && min(x, na.rm = TRUE) >= 0 && n_uniq <= 10) {
        v_type <- "Discrete / Count"
      } else {
        v_type <- "Continuous"
      }
      x_clean <- na.omit(x)
      mean_sd <- if (length(x_clean) > 0) sprintf("%.2f (%.2f)", mean(x_clean), sd(x_clean)) else "NA"
      med_iqr <- if (length(x_clean) > 0) sprintf("%.2f (%.2f)", median(x_clean), IQR(x_clean)) else "NA"
    } else {
      v_type <- if (n_uniq == 2) "Categorical (2 lvl)" else if (n_uniq <= 12) sprintf("Categorical (%d lvl)", n_uniq) else "High-Card / ID"
      mean_sd <- "-"
      med_iqr <- "-"
    }

    data.frame(
      Variable = col_name,
      Auto_Role = role_tag,
      Scale_Type = v_type,
      Unique = n_uniq,
      Missing = sprintf("%d (%.0f%%)", n_miss, pct_miss),
      Mean_SD = mean_sd,
      Median_IQR = med_iqr,
      stringsAsFactors = FALSE
    )
  })

  res_df <- do.call(rbind, summary_list)
  print(res_df, row.names = FALSE)
  cat("\nTip: Run `biostat_wizard(df)` to immediately analyze the inferred Y & X!\n")
  invisible(res_df)
}

# ------------------------------------------------------------------------------
# 3. EXPERIMENTAL DESIGN & BALANCE CHECKER (ZERO-INPUT READY)
# ------------------------------------------------------------------------------
biostat_design_check <- function(df, response = NULL, group = NULL, block = NULL) {
  parsed <- .auto_parse_roles(df, response = response, predictor = group, block = block, verbose = TRUE)
  if (parsed$mode == "table") {
    cat("Input is already a contingency table. Use `biostat_contingency(df)`.\n")
    return(invisible(NULL))
  }
  df <- parsed$df
  response <- parsed$response
  group <- parsed$predictor[1]
  block <- parsed$block

  cat("\n--- EXPERIMENTAL DESIGN DIAGNOSTICS ---\n")
  if (is.null(group)) {
    cat("Single-variable dataset (no grouping variable detected). N =", sum(!is.na(df[[response]])), "\n")
    return(invisible(NULL))
  }

  sub_df <- df[!is.na(df[[response]]) & !is.na(df[[group]]), ]
  g <- as.factor(sub_df[[group]])
  counts <- table(g)

  cat(sprintf("1. Replication across '%s' groups:\n", group))
  print(counts)

  if (any(counts < 3)) {
    cat("   [WARNING] Extremely low sample size (<3) in at least one group.\n")
  }

  is_balanced <- length(unique(counts)) == 1
  cat(sprintf("2. Design Balance: %s\n",
              ifelse(is_balanced, "BALANCED (Equal n across groups)",
                     "UNBALANCED (Unequal n; Welch's test or Type II/III SS recommended)")))

  if (!is.null(block) && block %in% names(sub_df)) {
    b <- as.factor(sub_df[[block]])
    cross_tab <- table(g, b)
    cat(sprintf("\n3. Blocking / Repeated Measures Check ('%s' x '%s'):\n", group, block))
    if (ncol(cross_tab) <= 12) print(cross_tab) else cat(sprintf("   [%d unique IDs/blocks detected]\n", ncol(cross_tab)))

    if (any(cross_tab > 1)) {
      cat("   [ALERT - PSEUDOREPLICATION / REPEATED MEASURES]: Multiple observations per ID within a group!\n")
      cat("   -> Average replicates per subject first, or fit a Mixed Model: lme4::lmer(y ~ group + (1|block))\n")
    } else if (all(cross_tab == 1)) {
      cat("   [INFO] Paired / Complete Randomized Block Design detected (1 obs per group/block).\n")
    }
  }
}

# ------------------------------------------------------------------------------
# 4. ASSUMPTION CHECKER (ZERO-INPUT READY)
# ------------------------------------------------------------------------------
biostat_assumptions <- function(df, response = NULL, predictor = NULL, alpha = 0.05, verbose_roles = TRUE) {
  parsed <- .auto_parse_roles(df, response = response, predictor = predictor, verbose = verbose_roles)
  if (parsed$mode == "table") {
    cat("\n--- CONTINGENCY TABLE ASSUMPTIONS ---\n")
    sq <- suppressWarnings(chisq.test(parsed$table, correct = FALSE))
    min_e <- min(sq$expected)
    cat(sprintf("Minimum Expected Cell Count = %.2f -> %s\n",
                min_e, ifelse(min_e >= 5, "PASS (Chi-Square Valid)", "FAIL (<5; Use Fisher's Exact Test)")))
    return(invisible(list(parametric_ok = min_e >= 5, min_expected = min_e)))
  }

  df <- parsed$df
  response <- parsed$response
  predictor <- parsed$predictor[1]

  cat("\n--- STATISTICAL ASSUMPTION DIAGNOSTICS ---\n")
  y <- df[[response]]

  if (!is.numeric(y)) {
    cat(sprintf("Response '%s' is categorical. Parametric normality/variance tests do not apply.\n", response))
    return(invisible(list(parametric_ok = FALSE, reason = "Categorical response")))
  }

  # Case A: Single numeric variable
  if (is.null(predictor)) {
    y_clean <- na.omit(y)
    sw <- if (length(y_clean) >= 3 && length(y_clean) <= 5000) shapiro.test(y_clean) else list(p.value = NA, statistic = NA)
    outliers <- boxplot.stats(y_clean)$out

    cat(sprintf("1. Normality of '%s' (Shapiro-Wilk): W = %.4f, p = %.4g -> %s\n",
                sw$statistic, sw$p.value,
                ifelse(is.na(sw$p.value), "Inspect Q-Q plot",
                       ifelse(sw$p.value >= alpha, "PASS (Normal)", "FAIL (Non-normal)"))))
    cat(sprintf("2. Outliers (1.5x IQR rule): %d detected\n", length(outliers)))

    param_ok <- ifelse(is.na(sw$p.value), TRUE, sw$p.value >= alpha)
    return(invisible(list(parametric_ok = param_ok, normality_ok = param_ok, shapiro_p = sw$p.value)))
  }

  # Case B: Response + Predictor
  complete_idx <- complete.cases(df[, c(response, predictor)])
  y_clean <- df[[response]][complete_idx]
  x_clean <- df[[predictor]][complete_idx]

  is_cat_x <- is.factor(x_clean) || is.character(x_clean) || (is.numeric(x_clean) && length(unique(x_clean)) <= 8)

  if (is_cat_x) {
    x_fac <- as.factor(x_clean)
    fit <- lm(y_clean ~ x_fac)
    resids <- residuals(fit)

    sw_res <- if (length(resids) >= 3 && length(resids) <= 5000) shapiro.test(resids) else list(p.value = NA)
    cat(sprintf("1. Normality of Residuals (Shapiro-Wilk): p = %.4g -> %s\n",
                sw_res$p.value,
                ifelse(is.na(sw_res$p.value), "Check Q-Q plot",
                       ifelse(sw_res$p.value >= alpha, "PASS (Normal residuals)", "FAIL (Non-normal residuals)"))))

    grp_sw <- tapply(y_clean, x_fac, function(val) {
      if (length(val) >= 3 && length(val) <= 5000 && length(unique(val)) > 1) shapiro.test(val)$p.value else NA
    })
    cat("   Group-level Shapiro-Wilk p-values:\n")
    print(round(grp_sw, 4))

    bart_p <- tryCatch(bartlett.test(y_clean ~ x_fac)$p.value, error = function(e) NA)
    flig_p <- tryCatch(fligner.test(y_clean ~ x_fac)$p.value, error = function(e) NA)

    cat("2. Homogeneity of Variance (Homoscedasticity):\n")
    cat(sprintf("   - Bartlett's Test (normal-sensitive) : p = %.4g -> %s\n",
                bart_p, ifelse(!is.na(bart_p) && bart_p >= alpha, "PASS (Equal var)", "FAIL (Unequal var)")))
    cat(sprintf("   - Fligner-Killeen (non-normal robust): p = %.4g -> %s\n",
                flig_p, ifelse(!is.na(flig_p) && flig_p >= alpha, "PASS (Equal var)", "FAIL (Unequal var)")))

    if (.has_pkg("car")) {
      lev_p <- car::leveneTest(y_clean ~ x_fac)$`Pr(>F)`[1]
      cat(sprintf("   - Levene's Test (car package)        : p = %.4g -> %s\n",
                  lev_p, ifelse(lev_p >= alpha, "PASS (Equal var)", "FAIL (Unequal var)")))
      var_ok <- lev_p >= alpha
    } else {
      var_ok <- ifelse(!is.na(flig_p), flig_p >= alpha, TRUE)
    }

    norm_ok <- ifelse(!is.na(sw_res$p.value), sw_res$p.value >= alpha, TRUE) &&
               all(is.na(grp_sw) | grp_sw >= alpha)

    if (!norm_ok) {
      cat("\n   [TIP]: Normality failed -> Wizard will automatically switch to Non-Parametric tests!\n")
    }

    return(invisible(list(
      parametric_ok = norm_ok && var_ok,
      normality_ok = norm_ok,
      variance_ok = var_ok,
      residual_shapiro_p = sw_res$p.value
    )))
  } else {
    fit <- lm(y_clean ~ x_clean)
    resids <- residuals(fit)
    sw_res <- if (length(resids) >= 3 && length(resids) <= 5000) shapiro.test(resids) else list(p.value = NA)
    sw_x <- if (length(x_clean) >= 3 && length(x_clean) <= 5000) shapiro.test(x_clean) else list(p.value = NA)
    sw_y <- if (length(y_clean) >= 3 && length(y_clean) <= 5000) shapiro.test(y_clean) else list(p.value = NA)

    cat(sprintf("1. Regression Residual Normality (Shapiro-Wilk): p = %.4g -> %s\n",
                sw_res$p.value, ifelse(sw_res$p.value >= alpha, "PASS", "FAIL")))
    cat(sprintf("2. Bivariate Normality (Pearson vs Spearman): Y p = %.4g | X p = %.4g\n",
                sw_y$p.value, sw_x$p.value))

    norm_ok <- (!is.na(sw_res$p.value) && sw_res$p.value >= alpha)
    bivar_ok <- (!is.na(sw_x$p.value) && sw_x$p.value >= alpha) && (!is.na(sw_y$p.value) && sw_y$p.value >= alpha)

    return(invisible(list(
      parametric_ok = norm_ok,
      normality_ok = norm_ok,
      bivariate_normal = bivar_ok,
      residual_shapiro_p = sw_res$p.value
    )))
  }
}

# ------------------------------------------------------------------------------
# 5. PROBLEM TYPE INFERENCE ENGINE (ZERO-INPUT READY)
# ------------------------------------------------------------------------------
biostat_infer <- function(df, response = NULL, predictors = NULL, paired = NULL, mu = 0, verbose_roles = TRUE) {
  parsed <- .auto_parse_roles(df, response = response, predictor = predictors, paired = paired, verbose = verbose_roles)

  cat("\n======================================================================\n")
  cat("                 BIOSTATISTICS PROBLEM INFERENCE ENGINE               \n")
  cat("======================================================================\n")

  if (parsed$mode == "table") {
    cat("INFERRED PROBLEM: Pre-Tabulated Contingency Table Analysis\n")
    cat("  * All Expected Counts >= 5 : chisq.test(tbl, correct = FALSE)\n")
    cat("  * Any Expected Count  <  5 : fisher.test(tbl)\n")
    return(invisible(list(prob_type = "contingency_table", parsed = parsed)))
  }

  df <- parsed$df
  response <- parsed$response
  predictors <- parsed$predictor
  paired <- parsed$paired

  classify_vec <- function(v) {
    v_clean <- na.omit(v)
    n_u <- length(unique(v_clean))
    if (is.factor(v) || is.character(v) || is.logical(v)) {
      return(list(type = "categorical", levels = n_u))
    } else if (is.numeric(v)) {
      if (n_u == 2) return(list(type = "categorical", levels = 2))
      if (all(v_clean == floor(v_clean)) && min(v_clean) >= 0 && n_u <= 8) {
        return(list(type = "count_or_ordinal", levels = n_u))
      }
      return(list(type = "continuous", levels = n_u))
    }
    return(list(type = "unknown", levels = n_u))
  }

  y_info <- classify_vec(df[[response]])
  cat(sprintf("Response Variable (Y) : '%s' -> [%s] (%d unique values)\n",
              response, toupper(y_info$type), y_info$levels))

  if (is.null(predictors) || length(predictors) == 0) {
    cat("Predictor Variable (X): None (One-Sample / Goodness-of-Fit problem)\n\n")
    if (y_info$type == "continuous") {
      cat("RECOMMENDED WORKFLOW:\n")
      cat(sprintf("  * Parametric     : One-Sample t-test -> t.test(df$%s, mu = %g)\n", response, mu))
      cat(sprintf("  * Non-Parametric : Wilcoxon Signed-Rank -> wilcox.test(df$%s, mu = %g)\n", response, mu))
      return(invisible(list(prob_type = "one_sample_continuous", parsed = parsed)))
    } else {
      cat("RECOMMENDED WORKFLOW:\n")
      cat(sprintf("  * Goodness-of-Fit: Chi-Square / Binomial -> chisq.test(table(df$%s))\n", response))
      return(invisible(list(prob_type = "one_sample_categorical", parsed = parsed)))
    }
  }

  if (length(predictors) == 1) {
    pred <- predictors[1]
    x_info <- classify_vec(df[[pred]])
    cat(sprintf("Predictor Variable (X): '%s' -> [%s] (%d unique values)\n",
                pred, toupper(x_info$type), x_info$levels))
    cat(sprintf("Study Design          : %s\n\n", ifelse(paired, "PAIRED / REPEATED MEASURES", "INDEPENDENT OBSERVATIONS")))

    if (y_info$type %in% c("continuous", "count_or_ordinal") && x_info$type == "categorical" && x_info$levels == 2) {
      if (paired) {
        cat("INFERRED PROBLEM: Comparing 2 Paired/Matched Groups (Continuous Outcome)\n")
        cat(sprintf("  * Parametric     : Paired t-test -> t.test(%s ~ %s, data = df, paired = TRUE)\n", response, pred))
        cat(sprintf("  * Non-Parametric : Wilcoxon Signed-Rank -> wilcox.test(%s ~ %s, data = df, paired = TRUE)\n", response, pred))
        return(invisible(list(prob_type = "paired_2group", parsed = parsed)))
      } else {
        cat("INFERRED PROBLEM: Comparing 2 Independent Groups (Continuous Outcome)\n")
        cat(sprintf("  * Parametric (Unequal Var) : Welch's t-test -> t.test(%s ~ %s, data = df, var.equal = FALSE)\n", response, pred))
        cat(sprintf("  * Parametric (Equal Var)   : Student's t-test -> t.test(%s ~ %s, data = df, var.equal = TRUE)\n", response, pred))
        cat(sprintf("  * Non-Parametric           : Mann-Whitney U -> wilcox.test(%s ~ %s, data = df)\n", response, pred))
        return(invisible(list(prob_type = "indep_2group", parsed = parsed)))
      }
    }

    if (y_info$type %in% c("continuous", "count_or_ordinal") && x_info$type == "categorical" && x_info$levels > 2) {
      if (paired) {
        cat("INFERRED PROBLEM: Comparing >2 Repeated Measures / Blocked Groups\n")
        cat("  * Parametric     : Repeated-Measures ANOVA -> aov(y ~ x + Error(block/x))\n")
        cat("  * Non-Parametric : Friedman Test -> friedman.test(y ~ x | block, data = df)\n")
        return(invisible(list(prob_type = "paired_kgroup", parsed = parsed)))
      } else {
        cat(sprintf("INFERRED PROBLEM: Comparing %d Independent Groups (One-Way ANOVA Design)\n", x_info$levels))
        cat(sprintf("  * Parametric (Equal Var)   : One-Way ANOVA -> aov(%s ~ %s, data = df) + TukeyHSD()\n", response, pred))
        cat(sprintf("  * Parametric (Unequal Var) : Welch's ANOVA -> oneway.test(%s ~ %s, data = df, var.equal = FALSE)\n", response, pred))
        cat(sprintf("  * Non-Parametric           : Kruskal-Wallis -> kruskal.test(%s ~ %s, data = df)\n", response, pred))
        return(invisible(list(prob_type = "indep_kgroup", parsed = parsed)))
      }
    }

    if (y_info$type == "continuous" && x_info$type %in% c("continuous", "count_or_ordinal")) {
      cat("INFERRED PROBLEM: Association / Prediction Between Two Quantitative Variables\n")
      cat(sprintf("  * Linear Regression (X -> Y) : summary(lm(%s ~ %s, data = df))\n", response, pred))
      cat(sprintf("  * Parametric Correlation     : cor.test(df$%s, df$%s, method = 'pearson')\n", pred, response))
      cat(sprintf("  * Non-Parametric Correlation : cor.test(df$%s, df$%s, method = 'spearman')\n", pred, response))
      return(invisible(list(prob_type = "continuous_xy", parsed = parsed)))
    }

    if (y_info$type == "categorical" && x_info$type == "categorical") {
      cat(sprintf("INFERRED PROBLEM: Contingency Table Analysis (%dx%d Categorical Association)\n", y_info$levels, x_info$levels))
      cat(sprintf("  * Expected Counts >= 5 : chisq.test(table(df$%s, df$%s), correct = FALSE)\n", response, pred))
      cat(sprintf("  * Expected Counts <  5 : fisher.test(table(df$%s, df$%s))\n", response, pred))
      return(invisible(list(prob_type = "contingency", parsed = parsed)))
    }

    if (y_info$type == "categorical" && y_info$levels == 2 && x_info$type == "continuous") {
      cat("INFERRED PROBLEM: Binary Outcome vs. Continuous Predictor\n")
      cat(sprintf("  * Logistic Regression : summary(glm(%s ~ %s, data = df, family = binomial))\n", response, pred))
      return(invisible(list(prob_type = "logistic", parsed = parsed)))
    }
  }

  invisible(list(prob_type = "custom", parsed = parsed))
}

# ------------------------------------------------------------------------------
# 6. CONTINGENCY TABLE SOLVER (ZERO-INPUT READY)
# ------------------------------------------------------------------------------
biostat_contingency <- function(df, var1 = NULL, var2 = NULL) {
  cat("\n--- CONTINGENCY TABLE & EPIDEMIOLOGICAL ANALYSIS ---\n")
  if (is.table(df) || is.matrix(df)) {
    tbl <- as.table(df)
  } else {
    parsed <- .auto_parse_roles(df, response = var1, predictor = var2, verbose = FALSE)
    v1 <- parsed$response
    v2 <- parsed$predictor[1]
    cat(sprintf("Auto-Selected Variables: '%s' (Rows) by '%s' (Columns)\n\n", v1, v2))
    tbl <- table(df[[v1]], df[[v2]])
  }

  cat("Observed Counts:\n")
  print(addmargins(tbl))

  cat("\nRow Proportions (%):\n")
  print(round(100 * prop.table(tbl, margin = 1), 1))

  sq <- suppressWarnings(chisq.test(tbl, correct = FALSE))
  cat("\nExpected Counts under Null Hypothesis (Independence):\n")
  print(round(sq$expected, 2))

  min_exp <- min(sq$expected)
  if (min_exp < 5) {
    cat(sprintf("\n[DECISION]: Minimum expected count is %.2f (< 5). Running Fisher's Exact Test:\n", min_exp))
    ft <- fisher.test(tbl, simulate.p.value = (nrow(tbl) > 2 || ncol(tbl) > 2))
    print(ft)
  } else {
    cat(sprintf("\n[DECISION]: All expected counts >= 5 (min = %.2f). Running Pearson's Chi-Square Test:\n", min_exp))
    print(sq)
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
# 7. SMART PLOTTER (ZERO-INPUT READY)
# ------------------------------------------------------------------------------
biostat_plot <- function(df, response = NULL, predictor = NULL, verbose_roles = TRUE) {
  parsed <- .auto_parse_roles(df, response = response, predictor = predictor, verbose = verbose_roles)
  if (parsed$mode == "table") {
    mosaicplot(parsed$table, main = "Contingency Mosaic Plot", color = TRUE)
    return(invisible(NULL))
  }

  df <- parsed$df
  response <- parsed$response
  predictor <- parsed$predictor[1]

  if (is.null(predictor)) {
    y <- na.omit(df[[response]])
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mfrow = c(1, 2))
    hist(y, freq = FALSE, main = paste("Distribution of", response),
         xlab = response, col = "lightblue", border = "white")
    lines(density(y), col = "darkblue", lwd = 2)
    qqnorm(y, main = paste("Normal Q-Q Plot:", response), pch = 19, col = "steelblue")
    qqline(y, col = "firebrick", lwd = 2)
    return(invisible(NULL))
  }

  sub_df <- df[complete.cases(df[, c(response, predictor)]), ]
  y <- sub_df[[response]]
  x <- sub_df[[predictor]]

  is_cat_x <- is.factor(x) || is.character(x) || (is.numeric(x) && length(unique(x)) <= 8)
  is_cat_y <- is.factor(y) || is.character(y) || (is.numeric(y) && length(unique(y)) <= 2)

  if (.has_pkg("ggplot2")) {
    if (!is_cat_y && is_cat_x) {
      sub_df[[predictor]] <- as.factor(sub_df[[predictor]])
      p <- ggplot2::ggplot(sub_df, ggplot2::aes(x = .data[[predictor]], y = .data[[response]], fill = .data[[predictor]])) +
        ggplot2::geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.5) +
        ggplot2::geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
        ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5, fill = "white") +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::labs(title = paste(response, "by", predictor),
                      subtitle = "Boxes = Median/IQR | White Diamond = Group Mean | Points = Raw Data") +
        ggplot2::theme(legend.position = "none")
      print(p)
      return(invisible(p))
    } else if (!is_cat_y && !is_cat_x) {
      p <- ggplot2::ggplot(sub_df, ggplot2::aes(x = .data[[predictor]], y = .data[[response]])) +
        ggplot2::geom_point(alpha = 0.75, size = 2.5, color = "steelblue") +
        ggplot2::geom_smooth(method = "lm", se = TRUE, color = "firebrick", fill = "pink", alpha = 0.3) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::labs(title = paste("Linear Regression:", response, "vs.", predictor),
                      subtitle = "Red line = OLS fit with 95% confidence band")
      print(p)
      return(invisible(p))
    }
  }

  if (!is_cat_y && is_cat_x) {
    x_fac <- as.factor(x)
    boxplot(y ~ x_fac, main = paste(response, "by", predictor),
            xlab = predictor, ylab = response, col = "lightblue", outline = FALSE)
    stripchart(y ~ x_fac, vertical = TRUE, method = "jitter", add = TRUE, pch = 19, col = "darkblue")
  } else if (!is_cat_y && !is_cat_x) {
    plot(x, y, main = paste(response, "vs.", predictor),
         xlab = predictor, ylab = response, pch = 19, col = "steelblue")
    abline(lm(y ~ x), col = "firebrick", lwd = 2)
  } else {
    mosaicplot(table(x, y), main = paste("Mosaic Plot:", predictor, "vs", response),
               xlab = predictor, ylab = response, color = TRUE)
  }
}

# ------------------------------------------------------------------------------
# 8. ZERO-INPUT MASTER WIZARD: CALL `biostat_wizard(df)` AND GET EVERYTHING!
# ------------------------------------------------------------------------------
biostat_wizard <- function(df, response = NULL, predictor = NULL, paired = NULL, alpha = 0.05, plot = TRUE) {
  # Step 1: Auto-infer roles and problem type
  inf <- biostat_infer(df, response = response, predictors = predictor, paired = paired, verbose_roles = TRUE)
  prob_type <- inf$prob_type
  parsed <- inf$parsed

  if (prob_type == "contingency_table") {
    if (plot) biostat_plot(parsed$table, verbose_roles = FALSE)
    return(biostat_contingency(parsed$table))
  }

  df_use <- parsed$df
  response <- parsed$response
  predictor <- parsed$predictor[1]
  paired <- parsed$paired

  # Step 2: Check assumptions automatically
  assump <- biostat_assumptions(df_use, response = response, predictor = predictor, alpha = alpha, verbose_roles = FALSE)

  cat("\n======================================================================\n")
  cat("                 AUTOMATED STATISTICAL EXECUTION                      \n")
  cat("======================================================================\n")

  # Plot automatically if requested
  if (plot) {
    tryCatch(biostat_plot(df_use, response = response, predictor = predictor, verbose_roles = FALSE),
             error = function(e) NULL)
  }

  if (is.null(predictor)) {
    y <- na.omit(df_use[[response]])
    if (is.numeric(y)) {
      if (assump$parametric_ok) {
        cat("-> Executing One-Sample t-test (Normality PASS):\n")
        res <- t.test(y)
      } else {
        cat("-> Executing One-Sample Wilcoxon Signed-Rank Test (Normality FAIL):\n")
        res <- wilcox.test(y, exact = FALSE)
      }
    } else {
      cat("-> Executing Chi-Square Goodness-of-Fit Test:\n")
      res <- chisq.test(table(y))
    }
    print(res)
    return(invisible(res))
  }

  sub_df <- df_use[complete.cases(df_use[, c(response, predictor)]), ]
  y <- sub_df[[response]]
  x <- sub_df[[predictor]]

  if (prob_type == "indep_2group") {
    x_fac <- as.factor(x)
    g1 <- y[x_fac == levels(x_fac)[1]]
    g2 <- y[x_fac == levels(x_fac)[2]]
    s_pool <- sqrt(((length(g1) - 1)*var(g1) + (length(g2) - 1)*var(g2)) / (length(g1) + length(g2) - 2))
    cohens_d <- (mean(g1) - mean(g2)) / s_pool

    if (assump$normality_ok) {
      var_eq <- assump$variance_ok
      cat(sprintf("-> Normality PASS. Variance equality %s. Running %s t-test:\n",
                  ifelse(var_eq, "PASS", "FAIL"),
                  ifelse(var_eq, "Student's (var.equal = TRUE)", "Welch's (var.equal = FALSE)")))
      res <- t.test(y ~ x_fac, var.equal = var_eq)
      print(res)
    } else {
      cat("-> Normality FAIL. Running Mann-Whitney U (Wilcoxon Rank-Sum) Test:\n")
      res <- wilcox.test(y ~ x_fac, exact = FALSE)
      print(res)
      cat("(Welch's t-test also shown for reference):\n")
      print(t.test(y ~ x_fac, var.equal = FALSE))
    }
    cat(sprintf("Effect Size -> Cohen's d = %.3f (%s)\n", cohens_d,
                ifelse(abs(cohens_d) < 0.2, "Negligible",
                       ifelse(abs(cohens_d) < 0.5, "Small",
                              ifelse(abs(cohens_d) < 0.8, "Medium", "Large")))))
    return(invisible(res))

  } else if (prob_type == "paired_2group") {
    x_fac <- as.factor(x)
    g1 <- y[x_fac == levels(x_fac)[1]]
    g2 <- y[x_fac == levels(x_fac)[2]]
    diffs <- g1 - g2
    sw_diff <- if (length(diffs) >= 3 && length(diffs) <= 5000) shapiro.test(diffs)$p.value else 1
    cat(sprintf("Normality of Paired Differences (Shapiro-Wilk): p = %.4g\n", sw_diff))
    if (sw_diff >= alpha) {
      cat("-> Running Paired t-test:\n")
      res <- t.test(g1, g2, paired = TRUE)
    } else {
      cat("-> Running Wilcoxon Signed-Rank Test on Paired Differences:\n")
      res <- wilcox.test(g1, g2, paired = TRUE, exact = FALSE)
    }
    print(res)
    return(invisible(res))

  } else if (prob_type == "indep_kgroup") {
    x_fac <- as.factor(x)
    fit_aov <- aov(y ~ x_fac)
    ss <- summary(fit_aov)[[1]]
    eta_sq <- ss[1, "Sum Sq"] / sum(ss[, "Sum Sq"])

    if (assump$normality_ok && assump$variance_ok) {
      cat("-> Assumptions PASS. Running Parametric One-Way ANOVA + Tukey HSD Post-Hoc:\n")
      print(summary(fit_aov))
      cat(sprintf("\nEffect Size -> Eta-squared = %.3f\n", eta_sq))
      if (ss[1, "Pr(>F)"] < alpha) {
        cat("\n--- Post-Hoc: Tukey's HSD ---\n")
        print(TukeyHSD(fit_aov))
      }
      return(invisible(fit_aov))
    } else if (assump$normality_ok && !assump$variance_ok) {
      cat("-> Normality PASS, Equal Variance FAIL. Running Welch's One-Way ANOVA:\n")
      res <- oneway.test(y ~ x_fac, var.equal = FALSE)
      print(res)
      cat("\n--- Post-Hoc: Pairwise Welch t-tests (Holm-adjusted) ---\n")
      print(pairwise.t.test(y, x_fac, pool.sd = FALSE, p.adjust.method = "holm"))
      return(invisible(res))
    } else {
      cat("-> Normality FAIL. Running Non-Parametric Kruskal-Wallis Test:\n")
      res <- kruskal.test(y ~ x_fac)
      print(res)
      if (res$p.value < alpha) {
        cat("\n--- Post-Hoc: Pairwise Wilcoxon Rank-Sum Tests (Holm-adjusted) ---\n")
        print(pairwise.wilcox.test(y, x_fac, p.adjust.method = "holm", exact = FALSE))
      }
      return(invisible(res))
    }

  } else if (prob_type == "continuous_xy") {
    cat("-> Fitting Linear Regression (OLS: y ~ x):\n")
    fit_lm <- lm(y ~ x)
    print(summary(fit_lm))

    if (assump$bivariate_normal) {
      cat("\n-> Bivariate Normality PASS. Running Pearson Correlation (r):\n")
      print(cor.test(x, y, method = "pearson"))
    } else {
      cat("\n-> Bivariate Normality FAIL. Running Spearman Rank Correlation (rho):\n")
      print(cor.test(x, y, method = "spearman", exact = FALSE))
    }
    return(invisible(fit_lm))

  } else if (prob_type == "contingency") {
    return(biostat_contingency(sub_df, response, predictor))

  } else if (prob_type == "logistic") {
    y_bin <- as.factor(y)
    cat("-> Fitting Binary Logistic Regression (GLM family = binomial):\n")
    fit_glm <- glm(y_bin ~ x, family = binomial)
    print(summary(fit_glm))
    cat("\nOdds Ratios (exp(beta)) and 95% CIs:\n")
    print(exp(cbind(OR = coef(fit_glm), suppressMessages(confint(fit_glm)))))
    return(invisible(fit_glm))
  }
}

# ------------------------------------------------------------------------------
# 9. SAMPLE SIZE & POWER CALCULATOR
# ------------------------------------------------------------------------------
biostat_power_calc <- function(test_type = c("t.two", "t.paired", "t.one", "anova"),
                               effect_size = 0.5, n = NULL, power = 0.80, alpha = 0.05, k_groups = 3) {
  test_type <- match.arg(test_type)
  cat("\n--- STATISTICAL POWER & SAMPLE SIZE CALCULATOR ---\n")
  if (test_type %in% c("t.two", "t.paired", "t.one")) {
    type_map <- c("t.two" = "two.sample", "t.paired" = "paired", "t.one" = "one.sample")
    res <- power.t.test(n = n, delta = effect_size, sd = 1, sig.level = alpha, power = power, type = type_map[test_type])
    print(res)
    return(invisible(res))
  } else if (test_type == "anova") {
    between_var <- (effect_size^2) * (k_groups / (k_groups - 1))
    res <- power.anova.test(groups = k_groups, n = n, between.var = between_var, within.var = 1, sig.level = alpha, power = power)
    print(res)
    return(invisible(res))
  }
}

message("Zero-Friction Gbiostats_toolkit.R loaded! Pass any dataset directly: e.g., biostat_wizard(df)")
