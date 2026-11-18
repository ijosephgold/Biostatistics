# ==============================================================================
# Gbiostats_toolkit.r (v4) - UNIVERSAL CSUN BIOSTATISTICS SUPERPOWER & AUDIT ENGINE
# ==============================================================================
# Built for zero-input execution (`biostat_wizard(data)`) + Full Self-Checking
# Transparency across ALL course modules & Problem Set 7 datasets:
#
#   * MASTER SUPERPOWER:
#     - biostat_wizard(data)         : Auto-detects 1, 2, or 3+ variables, runs
#                                      assumptions, executes decision tree, plots,
#                                      AND prints the Manual Verification Code!
#
#   * SELF-CHECKING & STEP-BY-STEP CLARIFICATION FUNCTIONS:
#     - biostat_explain(data)        : Narrative breakdown of WHY a test was chosen,
#                                      HOW assumptions were checked, & WHAT each
#                                      step of the decision tree did!
#     - biostat_code(data)           : Prints a clean, copy-pasteable Base-R script
#                                      with side-by-side comments so you can run
#                                      every single step manually to verify your work!
#     - biostat_check_math(data)     : Shows the raw under-the-hood math formulas
#                                      (Sums of Squares, Pooled SD, SE, t/F/Chi-Sq,
#                                      exp(b1), McFadden R2, AIC) plugged in with
#                                      your dataset's exact numbers!
#
#   * DEDICATED MULTI-VARIABLE & PROBLEM-SET SOLVERS:
#     - biostat_logistic(data)       : [PS7 Q1] Logistic GLM, Odds Factor exp(b1),
#                                      McFadden R2, LRT G2, & S-curve plot
#     - biostat_factorial(data)      : [PS7 Q2] Two-Way Model I ANOVA (Checks A*B ->
#                                      Separate t-tests if sig OR Additive A+B if NS)
#     - biostat_model3(data)         : [PS7 Q3] Mixed-Model III ANOVA (lme4 LLRT
#                                      Chi-Sq for Random Interaction + Fixed Effect)
#     - biostat_ancova(data)         : [PS7 Q4] ANCOVA (Checks X*Group slope int ->
#                                      Separate regressions if sig OR X+Group if NS)
#     - biostat_multireg(data)       : [PS7 Q5] Multiple Regression, Collinearity
#                                      SPLOM, All-Subsets AIC Ranking, & Best Model
# ==============================================================================

.has_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

# Store last run audit trail in toolkit environment so biostat_explain() / biostat_code()
# can be called even with zero arguments right after biostat_wizard(data)!
.biostat_last_audit <- new.env(parent = emptyenv())

# Safe column symbol wrapper for formulas (handles spaces like `Claw diff` if check.names=FALSE)
.sym_q <- function(col_name) {
  if (grepl("[^a-zA-Z0-9._]", col_name)) {
    paste0("`", col_name, "`")
  } else {
    col_name
  }
}

# ------------------------------------------------------------------------------
# 0. UNIVERSAL STRUCTURAL & ROLE INFERENCE ENGINE (.auto_parse_roles)
# ------------------------------------------------------------------------------
.auto_parse_roles <- function(data, response = NULL, predictor = NULL, block = NULL,
                              paired = NULL, model = "auto", random = NULL,
                              wide_mode = "auto", verbose = TRUE) {
  has_interaction_fml <- FALSE
  if (inherits(data, "formula")) {
    fml <- data
    if (is.null(response) || !is.data.frame(response)) {
      stop("When passing a formula first, provide the dataframe second: e.g., biostat_wizard(y ~ x, data)")
    }
    df <- as.data.frame(response)
    has_interaction_fml <- grepl("\\*|:", paste(deparse(fml), collapse = " "))
    vars <- all.vars(fml)
    response <- vars[1]
    predictor <- if (length(vars) > 1) vars[-1] else NULL
    data <- df
  } else if (inherits(response, "formula")) {
    fml <- response
    has_interaction_fml <- grepl("\\*|:", paste(deparse(fml), collapse = " "))
    vars <- all.vars(fml)
    response <- vars[1]
    predictor <- if (length(vars) > 1) vars[-1] else NULL
  }

  # Matrix / Contingency Table directly
  if (is.table(data) || (is.matrix(data) && is.numeric(data) && nrow(data) <= 15 && ncol(data) <= 15)) {
    if (verbose) cat("[AUTO-DETECT]: Input is a pre-tabulated Contingency Table / Matrix.\n")
    return(list(mode = "table", table = as.table(data)))
  }

  # Single Atomic Vector directly
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

  # Classify every column by structural pattern + semantic name cues
  col_meta <- lapply(cols, function(cn) {
    v <- df[[cn]]
    v_clean <- na.omit(v)
    n_u <- length(unique(v_clean))
    n_clean <- length(v_clean)

    # 1. Detect Sequential / Unique ID or Replicate/Subsample index columns
    # Note: In datasets like Medicago data.csv, `Replicate` (1..5) is a replicate index
    # when 2 other categorical treatment columns (Disturbance, Range) and 1 response (Fruits) exist!
    is_id_name <- grepl("^(id|sub|subject|patient|sample|specimen|tank|individual|mouse|rat|male|female|bird|fish|wing|butterfly|animal|plant|tree|leaf|person|pair|dyad|family|litter|clutch|nest|trial|rep|replicate|obs|index|row|colony)s?(_id|_num|_no)?$",
                        cn, ignore.case = TRUE)
    is_seq_int <- (is.numeric(v) && n_clean >= 3 && n_u == n_clean &&
                   all(v_clean == floor(v_clean)) &&
                   all(sort(v_clean) == seq(min(v_clean), max(v_clean))))
    is_unique_char <- ((is.character(v) || is.factor(v)) && n_clean >= 3 && n_u == n_clean)
    is_rep_col <- grepl("^(rep|replicate|subsample|trial|id|obs|index)s?(_id|_num|_no)?$", cn, ignore.case = TRUE) && is.numeric(v)

    if (is_seq_int || is_unique_char || (is_id_name && n_u == n_clean)) {
      role_type <- "id_or_block"
    } else if (is_rep_col && ncol(df) >= 4) {
      # E.g. Medicago data.csv has Disturbance, Range, Replicate (1..5), Fruits -> Replicate is a replicate index!
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

  # Auto-detect Pure Sequential/Unique Subject ID column if not explicitly supplied
  if (is.null(block)) {
    id_cands <- cols[sapply(col_meta, function(m) m$type == "id_or_block")]
    if (length(id_cands) > 0) block <- id_cands[1]
  }
  rep_cands <- cols[sapply(col_meta, function(m) m$type == "replicate_index")]

  active_cols <- setdiff(cols, c(block, rep_cands))
  if (length(active_cols) == 0) active_cols <- cols

  # Detect Pre-Aggregated Frequency Table in DataFrame Form (2 Categorical + 1 Count column)
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

  # Detect WIDE PAIRED Data (e.g., Butterfly Wings: Male | Polarized | Depolarized)
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

  # Intelligent Multi-Variable Role Inference (1, 2, or 3+ Active Columns!)
  cont_cols <- active_cols[sapply(active_cols, function(c) col_meta[[c]]$type %in% c("continuous", "small_integer"))]
  cat_cols  <- active_cols[sapply(active_cols, function(c) col_meta[[c]]$type %in% c("binary_cat", "multi_cat", "binary_num"))]

  if (is.null(response)) {
    if (length(active_cols) == 1) {
      response <- active_cols[1]
      predictor <- NULL

    } else if (length(active_cols) == 2) {
      # 2 ACTIVE COLUMNS
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
      # 3 OR MORE ACTIVE COLUMNS (Factorial ANOVA, Model III ANOVA, ANCOVA, or Multiple Regression!)
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

  # Auto-Identify Fixed vs. Random Factor if Model III ANOVA is requested or inferred
  if (length(predictor) == 2 && all(predictor %in% cat_cols)) {
    # If dataset has a Replicate column (like Medicago data.csv) or Range/Site/Location + model="III"
    if (is.null(random) && (identical(model, "III") || identical(model, 3) || identical(tolower(model), "model3"))) {
      r_cands <- predictor[sapply(predictor, function(p) col_meta[[p]]$is_random_name)]
      random <- if (length(r_cands) >= 1) r_cands[1] else predictor[2]
    }
  }

  # Auto-detect Paired / Blocked structure in 1-predictor Long Format
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
    model = model, has_interaction_fml = has_interaction_fml,
    col_meta = col_meta
  )
}

# ------------------------------------------------------------------------------
# 1. INTERACTIVE MENU & CSUN COURSE WORKFLOW GUIDE
# ------------------------------------------------------------------------------
biostat_guide <- function() {
  cat("======================================================================\n")
  cat("   UNIVERSAL BIOSTATISTICS TOOLKIT v4 (SUPERPOWER & SELF-CHECK EDITION)\n")
  cat("======================================================================\n")
  cat("1. ONE-LINE MASTER SUPERPOWER (Pass just `data`!):\n")
  cat("   - biostat_wizard(data)               : Auto-detects Paired, Logistic, Two-Way\n")
  cat("                                          ANOVA, ANCOVA, or Multiple Reg + AIC,\n")
  cat("                                          runs tests, plots, & prints R code!\n")
  cat("   - biostat_wizard(data, model = 'III'): Runs Model III Mixed-Effects ANOVA\n")
  cat("                                          (LLRT Chi-Sq for Random Interaction!)\n\n")
  cat("2. SELF-CHECKING, TRANSPARENCY & MANUAL VERIFICATION FUNCTIONS:\n")
  cat("   - biostat_code(data)                 : Prints copy-pasteable Base-R code with\n")
  cat("                                          side-by-side comments to check work!\n")
  cat("   - biostat_explain(data)              : Step-by-step plain-English walkthrough\n")
  cat("                                          of WHY & HOW the problem was solved!\n")
  cat("   - biostat_check_math(data)           : Shows the exact math formulas plugged\n")
  cat("                                          in with your dataset's numbers!\n\n")
  cat("3. DEDICATED PROBLEM-SET 7 & ADVANCED GLM SOLVERS:\n")
  cat("   - biostat_logistic(data)  : [PS7 Q1] Binary Logistic GLM, Odds Factor exp(b1),\n")
  cat("                               McFadden Pseudo-R2, & Sigmoidal Probability Curve\n")
  cat("   - biostat_factorial(data) : [PS7 Q2] Two-Way Model I ANOVA (Checks A*B first:\n")
  cat("                               if sig -> separate t-tests; if NS -> drop to A+B)\n")
  cat("   - biostat_model3(data)    : [PS7 Q3] Model III Mixed ANOVA (lme4 LLRT Chi-Sq\n")
  cat("                               for interaction; if sig -> fixed + random no-int)\n")
  cat("   - biostat_ancova(data)    : [PS7 Q4] ANCOVA (Tests X*Group slope interaction:\n")
  cat("                               if sig -> separate regressions; if NS -> X+Group)\n")
  cat("   - biostat_multireg(data)  : [PS7 Q5] Multiple Regression, Collinearity Matrix,\n")
  cat("                               All-Subsets AIC Table, & Best Reduced Model!\n")
  cat("======================================================================\n")
}

# ------------------------------------------------------------------------------
# 2. SELF-CHECKING ENGINE 1: COPY-PASTEABLE EQUIVALENT R CODE (`biostat_code`)
# ------------------------------------------------------------------------------
# Generates the exact standalone R code with side-by-side comments that a student
# would type by hand in RStudio to solve the problem and verify every number.
# ------------------------------------------------------------------------------
biostat_code <- function(df = NULL, response = NULL, predictor = NULL, paired = NULL,
                         model = "auto", random = NULL) {
  if (is.null(df)) {
    if (exists("last_code", envir = .biostat_last_audit)) {
      cat(.biostat_last_audit$last_code)
      return(invisible(.biostat_last_audit$last_code))
    } else {
      stop("Pass your dataset `biostat_code(data)` (or run `biostat_wizard(data)` first).")
    }
  }

  inf <- biostat_infer(df, response = response, predictors = predictor,
                       paired = paired, model = model, random = random, verbose_roles = FALSE)
  prob_type <- inf$prob_type
  parsed    <- inf$parsed
  y_col     <- parsed$response
  x_cols    <- parsed$predictor

  lines <- c(
    "======================================================================",
    "   EQUIVALENT STANDALONE R CODE (COPY & PASTE TO CHECK YOUR WORK!)    ",
    "======================================================================"
  )

  if (prob_type == "logistic") {
    x_col <- x_cols[1]
    yq <- .sym_q(y_col); xq <- .sym_q(x_col)
    lines <- c(lines,
      sprintf("# 1. Fit Binary Logistic Regression GLM (Binomial family, logit link)"),
      sprintf("mod <- glm(%s ~ %s, data = data, family = binomial)", yq, xq),
      sprintf("summary(mod)                                # Check Wald z-value and p-value (Pr(>|z|))"),
      sprintf(""),
      sprintf("# 2. Calculate Numeric Factor by which 1-unit increase in X increases Odds"),
      sprintf("exp(coef(mod))                              # exp(b1) = Odds Ratio multiplier"),
      sprintf("exp(confint(mod))                           # 95%% Confidence Interval for Odds Ratio"),
      sprintf(""),
      sprintf("# 3. Assess Strength of Association (McFadden's Pseudo-R2 & LRT Chi-Sq)"),
      sprintf("1 - (mod$deviance / mod$null.deviance)      # McFadden's Pseudo-R2 (proportion deviance explained)"),
      sprintf("anova(mod, test = 'Chisq')                  # Likelihood Ratio Test (Drop-in-Deviance G2)"),
      sprintf(""),
      sprintf("# 4. Plot Logistic Regression Curve with Raw Points"),
      sprintf("library(ggplot2)"),
      sprintf("ggplot(data, aes(x = %s, y = %s)) +", xq, yq),
      sprintf("  geom_jitter(height = 0.03, width = 0, alpha = 0.6) +"),
      sprintf("  stat_smooth(method = 'glm', method.args = list(family = 'binomial'), se = TRUE, color = 'red') +"),
      sprintf("  theme_classic()")
    )

  } else if (prob_type == "factorial_anova") {
    f1 <- .sym_q(x_cols[1]); f2 <- .sym_q(x_cols[2]); yq <- .sym_q(y_col)
    lines <- c(lines,
      sprintf("# 1. Fit Full Two-Way Factorial Model I ANOVA (WITH Interaction Term First!)"),
      sprintf("mod_int <- lm(%s ~ %s * %s, data = data)", yq, f1, f2),
      sprintf("shapiro.test(residuals(mod_int))            # Check normality of residuals"),
      sprintf("anova(mod_int)                              # Read the '%s:%s' interaction row FIRST!", x_cols[1], x_cols[2]),
      sprintf(""),
      sprintf("# 2A. IF INTERACTION IS NOT SIGNIFICANT (p >= 0.05): Drop interaction & test Main Effects"),
      sprintf("mod_add <- lm(%s ~ %s + %s, data = data)    # Additive model (holds other factor constant)", yq, f1, f2),
      sprintf("anova(mod_add)                              # Report main effect F and p-values from here!"),
      sprintf(""),
      sprintf("# 2B. IF INTERACTION IS SIGNIFICANT (p < 0.05): Do NOT read main effects above; run t-tests"),
      sprintf("t.test(%s ~ %s, data = data)                # Separate two-sample t-test for %s", yq, f1, x_cols[1]),
      sprintf("t.test(%s ~ %s, data = data)                # Separate two-sample t-test for %s", yq, f2, x_cols[2]),
      sprintf(""),
      sprintf("# 3. Plot Grouped Barplot of Means +/- SEM"),
      sprintf("library(ggplot2)"),
      sprintf("ggplot(data, aes(x = factor(%s), y = %s, fill = factor(%s))) +", f1, yq, f2),
      sprintf("  stat_summary(fun = mean, geom = 'bar', position = position_dodge(0.8), width = 0.7, color = 'black') +"),
      sprintf("  stat_summary(fun.data = mean_se, geom = 'errorbar', position = position_dodge(0.8), width = 0.2) +"),
      sprintf("  theme_classic()")
    )

  } else if (prob_type == "model3_anova") {
    r_col <- if (!is.null(parsed$random)) parsed$random else x_cols[2]
    f_col <- setdiff(x_cols, r_col)[1]
    fq <- .sym_q(f_col); rq <- .sym_q(r_col); yq <- .sym_q(y_col)
    lines <- c(lines,
      sprintf("# 1. Load lme4 and Fit Full Mixed Model (WITH Random Interaction) vs Reduced Model"),
      sprintf("library(lme4)"),
      sprintf("mod_full   <- lmer(%s ~ %s + (1|%s) + (1|%s:%s), data = data, REML = FALSE)", yq, fq, rq, fq, rq),
      sprintf("mod_no_int <- lmer(%s ~ %s + (1|%s), data = data, REML = FALSE)", yq, fq, rq),
      sprintf(""),
      sprintf("# 2. Run Log-Likelihood Ratio Test (LLRT Chi-Square) for the Interaction Term"),
      sprintf("anova(mod_no_int, mod_full)                 # Check Chisq and Pr(>Chisq) for interaction!"),
      sprintf(""),
      sprintf("# 3. If Interaction is Significant: Test Fixed Effect in Model w/ Random Effect & NO Interaction"),
      sprintf("mod_null   <- lmer(%s ~ 1 + (1|%s), data = data, REML = FALSE)", yq, rq),
      sprintf("anova(mod_null, mod_no_int)                 # LLRT for Fixed Effect ('%s')", f_col),
      sprintf("anova(lmer(%s ~ %s + (1|%s), data = data))  # F-table for Fixed Effect ('%s')", yq, fq, rq, f_col),
      sprintf(""),
      sprintf("# 4. Classical Two-Way ANOVA Table (for F-ratio reference)"),
      sprintf("anova(lm(%s ~ %s * %s, data = data))", yq, fq, rq)
    )

  } else if (prob_type == "ancova") {
    cov_col <- x_cols[1]; grp_col <- x_cols[2]
    cq <- .sym_q(cov_col); gq <- .sym_q(grp_col); yq <- .sym_q(y_col)
    lvls <- levels(as.factor(parsed$df[[grp_col]]))
    lines <- c(lines,
      sprintf("# 1. Fit Full ANCOVA Model with Interaction (Tests Homogeneity of Slopes)"),
      sprintf("mod_int <- lm(%s ~ %s * %s, data = data)", yq, cq, gq),
      sprintf("shapiro.test(residuals(mod_int))            # Check normality of residuals"),
      sprintf("anova(mod_int)                              # Read '%s:%s' interaction row FIRST!", cov_col, grp_col),
      sprintf(""),
      sprintf("# 2A. IF INTERACTION IS SIGNIFICANT (p < 0.05): Stop & run separate regressions per group!"),
      sprintf("anova(lm(%s ~ %s, data = subset(data, %s == '%s')))  # Regression for %s", yq, cq, gq, lvls[1], lvls[1]),
      sprintf("anova(lm(%s ~ %s, data = subset(data, %s == '%s')))  # Regression for %s", yq, cq, gq, lvls[min(2, length(lvls))], lvls[min(2, length(lvls))]),
      sprintf(""),
      sprintf("# 2B. IF INTERACTION IS NOT SIGNIFICANT (p >= 0.05): Drop interaction & compare line heights!"),
      sprintf("mod_add <- lm(%s ~ %s + %s, data = data)    # Parallel lines ANCOVA", yq, cq, gq),
      sprintf("anova(mod_add)                              # Row 1 = Slope (%s), Row 2 = Line Height (%s)", cov_col, grp_col),
      sprintf(""),
      sprintf("# 3. Scatterplot with Separate Trendlines by Group"),
      sprintf("library(ggplot2)"),
      sprintf("ggplot(data, aes(x = %s, y = %s, color = %s)) +", cq, yq, gq),
      sprintf("  geom_point(size = 2.5) + geom_smooth(method = 'lm', se = FALSE) + theme_classic()")
    )

  } else if (prob_type == "multiple_regression") {
    yq <- .sym_q(y_col)
    x_str <- paste(sapply(x_cols, .sym_q), collapse = " + ")
    lines <- c(lines,
      sprintf("# 1. Check Predictor Collinearity (Correlation Matrix & Scatterplot Matrix / SPLOM)"),
      sprintf("cor(data[, c(%s)], use = 'complete.obs')    # Check if predictors correlate strongly (|r| > 0.7)",
              paste0("'", x_cols, "'", collapse = ", ")),
      sprintf("pairs(data)                                 # Visual Scatterplot Matrix (SPLOM)"),
      sprintf(""),
      sprintf("# 2. Fit Full Additive Multiple Regression Model (No Interaction Terms)"),
      sprintf("mod_full <- lm(%s ~ %s, data = data)", yq, x_str),
      sprintf("shapiro.test(residuals(mod_full))           # Check residual normality"),
      sprintf("anova(mod_full)                             # Sequential ANOVA table"),
      sprintf("summary(mod_full)                           # Partial slopes, t-tests, R2, and Adjusted R2"),
      sprintf(""),
      sprintf("# 3. Compare Candidate Reduced Models Using Akaike's Information Criterion (AIC)"),
      sprintf("AIC(mod_full)                               # Smaller (more negative / lower) AIC = BETTER fit!")
    )
    for (xc in x_cols) {
      lines <- c(lines, sprintf("AIC(lm(%s ~ %s, data = data))               # Single-predictor model: %s only", yq, .sym_q(xc), xc))
    }

  } else if (prob_type == "paired_2group") {
    if (parsed$wide_paired) {
      c1 <- .sym_q(parsed$orig_cols[1]); c2 <- .sym_q(parsed$orig_cols[2])
      lines <- c(lines,
        sprintf("# 1. Compute Paired Differences & Check Normality (Required for Paired Design!)"),
        sprintf("diffs <- data$%s - data$%s", c1, c2),
        sprintf("shapiro.test(diffs)                         # p >= 0.05 -> Paired t-test; p < 0.05 -> Wilcoxon"),
        sprintf(""),
        sprintf("# 2. Execute Parametric Paired t-test (if normal) or Wilcoxon Signed-Rank (if non-normal)"),
        sprintf("t.test(data$%s, data$%s, paired = TRUE)     # Parametric Paired t-test", c1, c2),
        sprintf("wilcox.test(data$%s, data$%s, paired = TRUE)# Non-Parametric Wilcoxon Signed-Rank test", c1, c2)
      )
    } else {
      yq <- .sym_q(y_col); xq <- .sym_q(x_cols[1])
      lines <- c(lines,
        sprintf("# 1. Execute Paired t-test / Wilcoxon Signed-Rank Test"),
        sprintf("t.test(%s ~ %s, data = data, paired = TRUE)", yq, xq),
        sprintf("wilcox.test(%s ~ %s, data = data, paired = TRUE)", yq, xq)
      )
    }

  } else if (prob_type == "indep_2group") {
    yq <- .sym_q(y_col); xq <- .sym_q(x_cols[1])
    lines <- c(lines,
      sprintf("# 1. Check Normality & Homogeneity of Variance"),
      sprintf("shapiro.test(residuals(lm(%s ~ %s, data = data)))  # Normality check", yq, xq),
      sprintf("fligner.test(%s ~ %s, data = data)                 # Equal variance check (or bartlett.test)", yq, xq),
      sprintf(""),
      sprintf("# 2. Run Two-Sample Test Based on Assumptions"),
      sprintf("t.test(%s ~ %s, data = data, var.equal = FALSE)    # Welch's t-test (unequal variances / default)", yq, xq),
      sprintf("t.test(%s ~ %s, data = data, var.equal = TRUE)     # Student's t-test (equal variances)", yq, xq),
      sprintf("wilcox.test(%s ~ %s, data = data)                  # Mann-Whitney U test (non-normal data)", yq, xq)
    )

  } else if (prob_type == "indep_kgroup") {
    yq <- .sym_q(y_col); xq <- .sym_q(x_cols[1])
    lines <- c(lines,
      sprintf("# 1. Fit One-Way ANOVA & Check Residual Normality / Variance Homogeneity"),
      sprintf("fit <- aov(%s ~ %s, data = data)", yq, xq),
      sprintf("shapiro.test(residuals(fit))                # Normality of residuals"),
      sprintf("fligner.test(%s ~ %s, data = data)          # Homogeneity of variance", yq, xq),
      sprintf(""),
      sprintf("# 2. Run Omnibus Test + Post-Hoc Comparisons"),
      sprintf("summary(fit)                                # Standard One-Way ANOVA (if normal & equal var)"),
      sprintf("TukeyHSD(fit)                               # Tukey's HSD post-hoc pairwise comparisons"),
      sprintf("oneway.test(%s ~ %s, data = data, var.equal = FALSE) # Welch's ANOVA (if unequal var)", yq, xq),
      sprintf("kruskal.test(%s ~ %s, data = data)          # Non-parametric Kruskal-Wallis (if non-normal)", yq, xq)
    )

  } else if (prob_type == "continuous_xy") {
    yq <- .sym_q(y_col); xq <- .sym_q(x_cols[1])
    lines <- c(lines,
      sprintf("# 1. Fit Ordinary Least Squares (OLS) Linear Regression"),
      sprintf("mod <- lm(%s ~ %s, data = data)", yq, xq),
      sprintf("shapiro.test(residuals(mod))                # Check normality of regression residuals"),
      sprintf("anova(mod)                                  # Regression ANOVA table (F-ratio & Mean Sq)"),
      sprintf("summary(mod)                                # Intercept, Slope, R-squared, and p-value"),
      sprintf(""),
      sprintf("# 2. Correlation Analysis (Parametric Pearson r vs Non-Parametric Spearman rho)"),
      sprintf("cor.test(data$%s, data$%s, method = 'pearson')", xq, yq),
      sprintf("cor.test(data$%s, data$%s, method = 'spearman')", xq, yq)
    )
  }

  lines <- c(lines, "======================================================================\n")
  out_str <- paste(lines, collapse = "\n")
  cat("\n", out_str, sep = "")
  invisible(out_str)
}

# ------------------------------------------------------------------------------
# 3. SELF-CHECKING ENGINE 2: STEP-BY-STEP WALKTHROUGH (`biostat_explain`)
# ------------------------------------------------------------------------------
biostat_explain <- function(df = NULL, response = NULL, predictor = NULL, paired = NULL,
                            model = "auto", random = NULL) {
  if (is.null(df)) {
    if (exists("last_explain", envir = .biostat_last_audit)) {
      cat(.biostat_last_audit$last_explain)
      return(invisible(.biostat_last_audit$last_explain))
    } else {
      stop("Pass your dataset `biostat_explain(data)` (or run `biostat_wizard(data)` first).")
    }
  }
  # Running biostat_wizard with plot=FALSE populates and prints the full explanation!
  biostat_wizard(df, response = response, predictor = predictor, paired = paired,
                 model = model, random = random, plot = FALSE, show_code = TRUE)
}

# ------------------------------------------------------------------------------
# 4. SELF-CHECKING ENGINE 3: PLUGGED-IN MATH FORMULAS (`biostat_check_math`)
# ------------------------------------------------------------------------------
# Shows the exact hand-calculation equations plugged in with the dataset's numbers
# so you can verify how t, F, Odds Ratios, Pseudo-R2, or AIC were calculated!
# ------------------------------------------------------------------------------
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
  cat("\nTip: Run `biostat_wizard(data)` (or `biostat_wizard(data, model = 'III')` for mixed ANOVA)!\n")
  invisible(res_df)
}

# ------------------------------------------------------------------------------
# 6. DEDICATED SOLVER 1: BINARY LOGISTIC REGRESSION (PS7 Q1: Crayfish Mating)
# ------------------------------------------------------------------------------
biostat_logistic <- function(df, response = NULL, predictor = NULL, plot = TRUE, show_code = TRUE) {
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

  cat("3. Strength of Association Metrics:\n")
  cat(sprintf("   * McFadden's Pseudo-R2       : %.4f (%.1f%% reduction in deviance)\n", mcfadden_r2, mcfadden_r2 * 100))
  cat(sprintf("   * Nagelkerke's Pseudo-R2     : %.4f\n", nagelkerke_r2))
  cat(sprintf("   * Cox & Snell Pseudo-R2      : %.4f\n", cox_snell_r2))
  cat(sprintf("   * Tjur's Coefficient of Disc.: %.4f\n", tjur_r2))
  cat(sprintf("   * Point-Biserial Correlation : r = %.4f (p = %.5g)\n", pb_cor$estimate, pb_cor$p.value))

  if (show_code) {
    .biostat_last_audit$last_code <- biostat_code(df_use, response = response, predictor = predictor)
  }

  if (plot) {
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
biostat_sem_barplot <- function(df, response, factor1, factor2 = NULL, title_str = NULL) {
  sub_df <- df[complete.cases(df[, c(response, factor1, factor2)]), ]
  sub_df[[factor1]] <- as.factor(sub_df[[factor1]])
  if (!is.null(factor2)) sub_df[[factor2]] <- as.factor(sub_df[[factor2]])

  if (.has_pkg("ggplot2")) {
    if (is.null(factor2)) {
      p <- ggplot2::ggplot(sub_df, ggplot2::aes(x = .data[[factor1]], y = .data[[response]], fill = .data[[factor1]])) +
        ggplot2::stat_summary(fun = mean, geom = "bar", width = 0.6, color = "black", alpha = 0.8) +
        ggplot2::stat_summary(fun.data = function(y) {
          m <- mean(y); se <- sd(y) / sqrt(length(y))
          data.frame(y = m, ymin = m - se, ymax = m + se)
        }, geom = "errorbar", width = 0.2, linewidth = 0.8) +
        ggplot2::theme_classic(base_size = 13) +
        ggplot2::labs(title = ifelse(is.null(title_str), paste(response, "by", factor1, "(Mean +/- SEM)"), title_str),
                      x = factor1, y = paste(response, "(Mean +/- SEM)")) +
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
          title = ifelse(is.null(title_str), paste(response, "by", factor1, "and", factor2, "(Mean +/- SEM)"), title_str),
          x = factor1, y = paste(response, "(Mean +/- SEM)"), fill = factor2
        )
      print(p)
      return(invisible(p))
    }
  } else {
    if (is.null(factor2)) {
      means <- tapply(sub_df[[response]], sub_df[[factor1]], mean)
      ses   <- tapply(sub_df[[response]], sub_df[[factor1]], function(v) sd(v)/sqrt(length(v)))
      bp <- barplot(means, ylim = c(0, max(means + ses, na.rm = TRUE) * 1.15),
                    col = "steelblue3", ylab = paste(response, "(Mean +/- SEM)"), xlab = factor1,
                    main = ifelse(is.null(title_str), paste(response, "by", factor1), title_str))
      arrows(bp, means - ses, bp, means + ses, angle = 90, code = 3, length = 0.08, lwd = 1.5)
    } else {
      means <- tapply(sub_df[[response]], list(sub_df[[factor2]], sub_df[[factor1]]), mean)
      ses   <- tapply(sub_df[[response]], list(sub_df[[factor2]], sub_df[[factor1]]), function(v) sd(v)/sqrt(length(v)))
      bp <- barplot(means, beside = TRUE, ylim = c(0, max(means + ses, na.rm = TRUE) * 1.2),
                    col = c("salmon", "darkturquoise", "goldenrod", "steelblue")[seq_len(nrow(means))],
                    legend.text = TRUE, args.legend = list(title = factor2, x = "topright"),
                    xlab = factor1, ylab = paste(response, "(Mean +/- SEM)"),
                    main = ifelse(is.null(title_str), paste(response, "by", factor1, "&", factor2), title_str))
      arrows(bp, means - ses, bp, means + ses, angle = 90, code = 3, length = 0.06, lwd = 1.5)
    }
  }
}

# ------------------------------------------------------------------------------
# 8. DEDICATED SOLVER 2: TWO-WAY MODEL I FACTORIAL ANOVA (PS7 Q2: Coral Growth)
# ------------------------------------------------------------------------------
biostat_factorial <- function(df, response = NULL, factor1 = NULL, factor2 = NULL,
                              alpha = 0.05, plot = TRUE, show_code = TRUE) {
  parsed <- .auto_parse_roles(df, response = response, predictor = c(factor1, factor2), verbose = FALSE)
  df_use <- parsed$df
  if (is.null(response)) response <- parsed$response
  if (is.null(factor1))  factor1  <- parsed$predictor[1]
  if (is.null(factor2))  factor2  <- parsed$predictor[2]

  sub_df <- df_use[complete.cases(df_use[, c(response, factor1, factor2)]), ]
  sub_df[[factor1]] <- as.factor(sub_df[[factor1]])
  sub_df[[factor2]] <- as.factor(sub_df[[factor2]])

  cat("\n======================================================================\n")
  cat("          FACTORIAL TWO-WAY MODEL I ANOVA (FIXED EFFECTS)             \n")
  cat("======================================================================\n")
  cat(sprintf("Response Variable (Y) : '%s'\n", response))
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

  int_p <- aov_int[3, "Pr(>F)"]
  int_f <- aov_int[3, "F value"]

  cat("\n--- STEP 2: COURSE DECISION RULE ON INTERACTION TERM ---\n")
  if (!is.na(int_p) && int_p < alpha) {
    cat(sprintf("[INTERACTION IS SIGNIFICANT]: %s:%s F = %.3f, p = %.5g (< %.2f)\n", factor1, factor2, int_f, int_p, alpha))
    cat(sprintf("-> CONCLUSION: The effect of '%s' on '%s' DEPENDS UPON '%s'!\n", factor1, response, factor2))
    cat("-> COURSE CONVENTION: Do NOT interpret the main effects from the interaction ANOVA table.\n")
    cat("   Instead, analyze the factors separately using two-sample t-tests (below):\n\n")

    for (fac in c(factor1, factor2)) {
      if (nlevels(sub_df[[fac]]) == 2) {
        f_fml <- as.formula(paste(.sym_q(response), "~", .sym_q(fac)))
        var_p <- tryCatch(fligner.test(f_fml, data = sub_df)$p.value, error = function(e) 0)
        veq <- (!is.na(var_p) && var_p >= alpha)
        tt <- t.test(f_fml, data = sub_df, var.equal = veq)
        cat(sprintf("   * Separate Two-Sample t-test for '%s' (%s, var.equal = %s):\n",
                    fac, ifelse(veq, "Student's t-test", "Welch's t-test"), veq))
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
    cat(sprintf("-> CONCLUSION: The effect of '%s' on '%s' does NOT depend upon '%s'.\n", factor1, response, factor2))
    cat("-> COURSE CONVENTION: Drop the interaction term and re-run the Additive Model\n")
    cat("   to test each main effect while holding the other factor constant:\n\n")

    fml_add <- as.formula(paste(.sym_q(response), "~", .sym_q(factor1), "+", .sym_q(factor2)))
    fit_add <- lm(fml_add, data = sub_df)
    aov_add <- anova(fit_add)
    cat("--- ADDITIVE TWO-WAY ANOVA TABLE (", deparse(fml_add), ") ---\n", sep = "")
    print(aov_add)
    final_model <- fit_add
  }

  if (show_code) {
    .biostat_last_audit$last_code <- biostat_code(df_use, response = response, predictor = c(factor1, factor2))
  }

  if (plot) {
    biostat_sem_barplot(sub_df, response, factor1, factor2)
  }

  invisible(final_model)
}

# ------------------------------------------------------------------------------
# 9. DEDICATED SOLVER 3: MODEL III MIXED-EFFECTS FACTORIAL ANOVA (PS7 Q3: Medicago)
# ------------------------------------------------------------------------------
biostat_model3 <- function(df, response = NULL, fixed = NULL, random = NULL,
                           alpha = 0.05, plot = TRUE, show_code = TRUE) {
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

  cat("\n======================================================================\n")
  cat("       MIXED-MODEL FACTORIAL ANOVA (MODEL III ANOVA + LLRT)           \n")
  cat("======================================================================\n")
  cat(sprintf("Response Variable (Y) : '%s'\n", response))
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

  # Always show the classical ANOVA F-table for the interaction and main effects too
  cat("\n--- CLASSICAL FACTORIAL ANOVA PARTITIONING (FOR F-RATIO REFERENCE) ---\n")
  fit_ols <- lm(as.formula(paste(yq, "~", fq, "*", rq)), data = sub_df)
  print(anova(fit_ols))

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

  if (show_code) {
    .biostat_last_audit$last_code <- biostat_code(df_use, response = response, predictor = c(fixed, random),
                                                  model = "III", random = random)
  }

  if (plot) {
    biostat_sem_barplot(sub_df, response, fixed, random,
                        title_str = paste("Model III ANOVA:", response, "by", fixed, "&", random, "(Mean +/- SEM)"))
  }

  invisible(sub_df)
}

# ------------------------------------------------------------------------------
# 10. DEDICATED SOLVER 4: ANCOVA (PS7 Q4: Red Algae Stipes / Mussels / Crickets)
# ------------------------------------------------------------------------------
biostat_ancova <- function(df, response = NULL, covariate = NULL, group = NULL,
                           alpha = 0.05, plot = TRUE, show_code = TRUE) {
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
  yq <- .sym_q(response); cq <- .sym_q(covariate); gq <- .sym_q(group)

  cat("\n======================================================================\n")
  cat("        ANALYSIS OF COVARIANCE (ANCOVA: SLOPES & LINE HEIGHTS)        \n")
  cat("======================================================================\n")
  cat(sprintf("Dependent Variable (Cont Y)   : '%s'\n", response))
  cat(sprintf("Covariate Predictor (Cont X)  : '%s'\n", covariate))
  cat(sprintf("Grouping Variable (Factor G)  : '%s' (%s)\n\n", group, paste(levels(sub_df[[group]]), collapse = ", ")))

  fml_int <- as.formula(paste(yq, "~", cq, "*", gq))
  fit_int <- lm(fml_int, data = sub_df)
  aov_int <- anova(fit_int)

  cat("--- STEP 1: TEST HOMOGENEITY OF SLOPES (INTERACTION MODEL: ", deparse(fml_int), ") ---\n", sep = "")
  print(aov_int)

  sw_p <- if (nrow(sub_df) >= 3 && nrow(sub_df) <= 5000) shapiro.test(residuals(fit_int))$p.value else NA
  cat(sprintf("Residual Normality Check (Shapiro-Wilk): p = %.4g\n", sw_p))

  int_p <- aov_int[3, "Pr(>F)"]
  int_f <- aov_int[3, "F value"]

  cat("\n--- STEP 2: COURSE DECISION RULE ON ANCOVA INTERACTION ---\n")
  if (!is.na(int_p) && int_p < alpha) {
    cat(sprintf("[INTERACTION IS SIGNIFICANT]: %s:%s F = %.3f, p = %.5g (< %.2f)\n", covariate, group, int_f, int_p, alpha))
    cat(sprintf("-> CONCLUSION: The effect (slope) of '%s' on '%s' DEPENDS UPON '%s'!\n", covariate, response, group))
    cat("-> COURSE CONVENTION (PS7 Q4 / Crickets Example): Stop interpreting the combined ANOVA table\n")
    cat(sprintf("   and run SEPARATE linear regressions of '%s ~ %s' for each '%s' below:\n\n", response, covariate, group))

    for (lvl in levels(sub_df[[group]])) {
      sub_g <- sub_df[sub_df[[group]] == lvl, ]
      fml_sep <- as.formula(paste(yq, "~", cq))
      fit_sep <- lm(fml_sep, data = sub_g)
      aov_sep <- anova(fit_sep)
      sum_sep <- summary(fit_sep)

      cat(sprintf("   ================================================================\n"))
      cat(sprintf("   * Separate Linear Regression for %s = '%s' (n = %d):\n", group, lvl, nrow(sub_g)))
      cat(sprintf("     Equation : %s = %.4f + (%.4f * %s)\n",
                  response, coef(fit_sep)[1], coef(fit_sep)[2], covariate))
      cat(sprintf("     ANOVA    : F(1, %d) = %.3f, Mean Sq = %.3f, p = %.5g, R2 = %.4f\n",
                  fit_sep$df.residual, aov_sep[1, "F value"], aov_sep[1, "Mean Sq"],
                  aov_sep[1, "Pr(>F)"], sum_sep$r.squared))
      print(aov_sep)
      cat("\n")
    }
    final_model <- fit_int

  } else {
    cat(sprintf("[INTERACTION IS NOT SIGNIFICANT]: %s:%s F = %.3f, p = %.5g (>= %.2f)\n", covariate, group, int_f, int_p, alpha))
    cat(sprintf("-> CONCLUSION: The slopes are parallel! The effect of '%s' on '%s' does NOT depend upon '%s'.\n",
                covariate, response, group))
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

  if (show_code) {
    .biostat_last_audit$last_code <- biostat_code(df_use, response = response, predictor = c(covariate, group))
  }

  if (plot) {
    if (.has_pkg("ggplot2")) {
      p <- ggplot2::ggplot(sub_df, ggplot2::aes(x = .data[[covariate]], y = .data[[response]], color = .data[[group]])) +
        ggplot2::geom_point(size = 2.6, alpha = 0.8) +
        ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 1.1) +
        ggplot2::theme_classic(base_size = 13) +
        ggplot2::labs(
          title = paste("ANCOVA:", response, "vs.", covariate, "by", group),
          subtitle = sprintf("Interaction %s:%s p = %.4g (%s)",
                             covariate, group, int_p,
                             ifelse(!is.na(int_p) && int_p < alpha, "Separate Slopes", "Parallel Slopes")),
          x = covariate, y = response, color = group
        )
      print(p)
    } else {
      cols_g <- c("firebrick", "darkturquoise", "goldenrod", "navy")[seq_len(nlevels(sub_df[[group]]))]
      plot(sub_df[[covariate]], sub_df[[response]], col = cols_g[sub_df[[group]]], pch = 19,
           xlab = covariate, ylab = response, main = paste("ANCOVA:", response, "vs.", covariate, "by", group))
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
                             alpha = 0.05, plot = TRUE, show_code = TRUE) {
  parsed <- .auto_parse_roles(df, response = response, predictor = predictors, verbose = FALSE)
  df_use <- parsed$df
  if (is.null(response))   response   <- parsed$response
  if (is.null(predictors)) predictors <- parsed$predictor

  sub_df <- df_use[complete.cases(df_use[, c(response, predictors)]), c(response, predictors)]
  yq <- .sym_q(response)

  cat("\n======================================================================\n")
  cat("     MULTIPLE LINEAR REGRESSION, COLLINEARITY & AIC MODEL SELECTION   \n")
  cat("======================================================================\n")
  cat(sprintf("Dependent Response (Y)     : '%s' (n = %d observations)\n", response, nrow(sub_df)))
  cat(sprintf("Candidate Predictors (X_i) : %s\n\n", paste(predictors, collapse = ", ")))

  # Step 1: Collinearity Diagnostics (Pairwise Correlation Matrix & p-values)
  cat("--- STEP 1: COLLINEARITY CHECK & INDIVIDUAL PREDICTOR CORRELATIONS ---\n")
  cor_full <- round(cor(sub_df, use = "complete.obs"), 4)
  print(cor_full)

  cat("\n1A. Individual Linear Relationships with Response ('", response, "'):\n", sep = "")
  for (xc in predictors) {
    ct_y <- suppressWarnings(cor.test(sub_df[[xc]], sub_df[[response]]))
    m_single <- lm(as.formula(paste(yq, "~", .sym_q(xc))), data = sub_df)
    cat(sprintf("  * %s ~ %-8s : r = %6.3f, R2 = %.4f, F = %6.2f, p = %.4g, AIC = %.2f\n",
                response, xc, ct_y$estimate, summary(m_single)$r.squared,
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

  # Step 2: Full Additive Multiple Linear Regression
  fml_full <- as.formula(paste(yq, "~", paste(sapply(predictors, .sym_q), collapse = " + ")))
  fit_full <- lm(fml_full, data = sub_df)
  sw_p <- if (nrow(sub_df) >= 3 && nrow(sub_df) <= 5000) shapiro.test(residuals(fit_full))$p.value else NA

  cat(sprintf("\n--- STEP 2: FULL MULTIPLE REGRESSION MODEL (%s) ---\n", deparse(fml_full)))
  cat(sprintf("Residual Normality (Shapiro-Wilk): p = %.4g\n\n", sw_p))
  cat("ANOVA Table (Full Model):\n")
  print(anova(fit_full))
  cat("\nRegression Coefficients & Overall R-Squared (Full Model):\n")
  print(summary(fit_full))

  # Step 3: All-Subsets AIC Model Selection Table
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

  # Step 4: Fit & Report the Best (Lowest-AIC) Model
  best_fml_str <- aic_df$Model_Formula[1]
  best_fit <- lm(as.formula(best_fml_str), data = sub_df)
  cat(sprintf("\n--- STEP 4: BEST SELECTED MODEL BY AIC (%s) ---\n", best_fml_str))
  cat(sprintf("Lowest AIC = %.2f | Multiple R2 = %.4f (%.1f%% variance explained) | Adjusted R2 = %.4f\n\n",
              AIC(best_fit), summary(best_fit)$r.squared, summary(best_fit)$r.squared * 100, summary(best_fit)$adj.r.squared))
  cat("Best Model ANOVA Table:\n")
  print(anova(best_fit))
  cat("\nBest Model Regression Summary:\n")
  print(summary(best_fit))

  if (show_code) {
    .biostat_last_audit$last_code <- biostat_code(df_use, response = response, predictor = predictors)
  }

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
          main = paste("Scatterplot Matrix (SPLOM):", response, "& Predictors"))
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

  cat("\n--- STATISTICAL ASSUMPTION DIAGNOSTICS ---\n")
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
                  sw_d$statistic, sw_d$p.value, ifelse(norm_ok, "PASS (Paired t-test valid)", "FAIL (Use Wilcoxon Signed-Rank)")))
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

    cat(sprintf("1. Residual Normality (Shapiro-Wilk): p = %.4g -> %s\n",
                sw_res$p.value, ifelse(norm_ok, "PASS", "FAIL")))
    cat(sprintf("2. Homogeneity of Variance: Bartlett p = %.4g | Fligner-Killeen p = %.4g -> %s\n",
                bart_p, flig_p, ifelse(var_ok, "PASS (Equal var)", "FAIL (Unequal var -> Use Welch)")))
    return(invisible(list(parametric_ok = norm_ok && var_ok, normality_ok = norm_ok, variance_ok = var_ok)))
  } else {
    fit <- lm(y_clean ~ x_clean)
    resids <- residuals(fit)
    sw_res <- if (length(resids) >= 3 && length(resids) <= 5000) shapiro.test(resids) else list(p.value = NA)
    sw_x <- if (length(x_clean) >= 3 && length(x_clean) <= 5000) shapiro.test(x_clean) else list(p.value = NA)
    sw_y <- if (length(y_clean) >= 3 && length(y_clean) <= 5000) shapiro.test(y_clean) else list(p.value = NA)
    norm_ok <- (!is.na(sw_res$p.value) && sw_res$p.value >= alpha)
    bivar_ok <- (!is.na(sw_x$p.value) && sw_x$p.value >= alpha) && (!is.na(sw_y$p.value) && sw_y$p.value >= alpha)
    cat(sprintf("1. Regression Residual Normality (Shapiro-Wilk): p = %.4g -> %s\n",
                sw_res$p.value, ifelse(norm_ok, "PASS", "FAIL")))
    cat(sprintf("2. Bivariate Normality: Y p = %.4g | X p = %.4g\n", sw_y$p.value, sw_x$p.value))
    return(invisible(list(parametric_ok = norm_ok, normality_ok = norm_ok, bivariate_normal = bivar_ok)))
  }
}

# ------------------------------------------------------------------------------
# 13. PROBLEM TYPE INFERENCE ENGINE (HANDLES ALL SINGLE & MULTI-FACTOR DESIGNS)
# ------------------------------------------------------------------------------
biostat_infer <- function(df, response = NULL, predictors = NULL, paired = NULL,
                          model = "auto", random = NULL, mu = 0, verbose_roles = TRUE) {
  parsed <- .auto_parse_roles(df, response = response, predictor = predictors,
                              paired = paired, model = model, random = random, verbose = verbose_roles)

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
        cat(sprintf("INFERRED PROBLEM: %s\n", ifelse(paired, "Paired Two-Sample Comparison (Paired t-test / Wilcoxon Signed-Rank)",
                                                        "Independent Two-Sample Comparison (Student/Welch t-test / Mann-Whitney U)")))
      }
      return(invisible(list(prob_type = prob, parsed = parsed)))
    }
    if (y_info$type == "continuous" && x_info$type == "categorical" && x_info$levels > 2) {
      prob <- if (paired) "paired_kgroup" else "indep_kgroup"
      if (verbose_roles) {
        cat(sprintf("INFERRED PROBLEM: %s\n", ifelse(paired, "Blocked / Repeated-Measures ANOVA / Friedman Test",
                                                        "One-Way ANOVA / Welch ANOVA / Kruskal-Wallis")))
      }
      return(invisible(list(prob_type = prob, parsed = parsed)))
    }
    if (y_info$type == "continuous" && x_info$type == "continuous") {
      if (verbose_roles) cat("INFERRED PROBLEM: Simple Linear Regression & Correlation (OLS / Pearson r / Spearman rho)\n")
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
# 15. SMART PLOTTER (HANDLES 1-VAR, PAIRED, GROUP SEM/BOXPLOT, OLS, & MULTI-VAR)
# ------------------------------------------------------------------------------
biostat_plot <- function(df, response = NULL, predictor = NULL, paired = NULL, verbose_roles = TRUE) {
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

  if (is.null(predictor)) {
    y <- na.omit(df_use[[response]])
    old_par <- par(no.readonly = TRUE); on.exit(par(old_par)); par(mfrow = c(1, 2))
    hist(y, freq = FALSE, main = paste("Distribution of", response), xlab = response, col = "lightblue", border = "white")
    lines(density(y), col = "darkblue", lwd = 2)
    qqnorm(y, main = paste("Normal Q-Q Plot:", response), pch = 19, col = "steelblue"); qqline(y, col = "firebrick", lwd = 2)
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
        ggplot2::labs(title = paste("Paired / Blocked Comparison:", response, "by", predictor),
                      subtitle = paste("Gray lines connect matched", block, "| White Diamond = Mean | Box = Median/IQR")) +
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
    boxplot(y ~ x_fac, main = paste(response, "by", predictor), xlab = predictor, ylab = response, col = "lightblue", outline = FALSE)
    stripchart(y ~ x_fac, vertical = TRUE, method = "jitter", add = TRUE, pch = 19, col = "darkblue")
  } else if (!is_cat_y && !is_cat_x) {
    plot(x, y, main = paste(response, "vs.", predictor), xlab = predictor, ylab = response, pch = 19, col = "steelblue")
    abline(lm(y ~ x), col = "firebrick", lwd = 2)
  } else {
    mosaicplot(table(x, y), main = paste("Mosaic Plot:", predictor, "vs", response), xlab = predictor, ylab = response, color = TRUE)
  }
}

# ------------------------------------------------------------------------------
# 16. MASTER WIZARD: ONE-LINE EXECUTION + EQUIVALENT CODE FOR ANY DATASET
# ------------------------------------------------------------------------------
biostat_wizard <- function(df, response = NULL, predictor = NULL, paired = NULL,
                           model = "auto", random = NULL, mu = 0, alpha = 0.05,
                           plot = TRUE, show_code = TRUE) {
  inf <- biostat_infer(df, response = response, predictors = predictor,
                       paired = paired, model = model, random = random, mu = mu, verbose_roles = TRUE)
  prob_type <- inf$prob_type
  parsed <- inf$parsed

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
                             alpha = alpha, plot = plot, show_code = show_code))
  }
  if (prob_type == "model3_anova") {
    return(biostat_model3(parsed$df, response = parsed$response,
                          fixed = setdiff(parsed$predictor, parsed$random)[1],
                          random = parsed$random, alpha = alpha, plot = plot, show_code = show_code))
  }
  if (prob_type == "ancova") {
    return(biostat_ancova(parsed$df, response = parsed$response,
                          covariate = parsed$predictor[1], group = parsed$predictor[2],
                          alpha = alpha, plot = plot, show_code = show_code))
  }
  if (prob_type == "multiple_regression") {
    return(biostat_multireg(parsed$df, response = parsed$response,
                            predictors = parsed$predictor, alpha = alpha, plot = plot, show_code = show_code))
  }

  df_use <- parsed$df
  response <- parsed$response
  predictors <- parsed$predictor
  predictor <- predictors[1]
  paired <- parsed$paired
  block <- parsed$block

  assump <- biostat_assumptions(df_use, response = response, predictor = predictor,
                                paired = paired, alpha = alpha, verbose_roles = FALSE)

  cat("\n======================================================================\n")
  cat("                 AUTOMATED STATISTICAL EXECUTION                      \n")
  cat("======================================================================\n")

  if (plot && length(predictors) <= 1) {
    tryCatch(biostat_plot(df, response = response, predictor = predictor, paired = paired, verbose_roles = FALSE),
             error = function(e) NULL)
  }

  if (is.null(predictor)) {
    y <- na.omit(df_use[[response]])
    if (is.numeric(y)) {
      if (assump$parametric_ok) {
        cat(sprintf("-> Executing One-Sample t-test (mu = %g, Normality PASS):\n", mu))
        res <- t.test(y, mu = mu)
      } else {
        cat(sprintf("-> Executing One-Sample Wilcoxon Signed-Rank Test (mu = %g, Normality FAIL):\n", mu))
        res <- wilcox.test(y, mu = mu, exact = FALSE)
      }
    } else {
      cat("-> Executing Chi-Square Goodness-of-Fit Test:\n")
      res <- chisq.test(table(y))
    }
    print(res)
    return(invisible(res))
  }

  sub_df <- df_use[complete.cases(df_use[, c(response, predictors)]), ]
  y <- sub_df[[response]]
  x <- sub_df[[predictor]]

  if (prob_type == "indep_2group") {
    x_fac <- as.factor(x)
    g1 <- y[x_fac == levels(x_fac)[1]]
    g2 <- y[x_fac == levels(x_fac)[2]]
    s_pool <- sqrt(((length(g1) - 1)*var(g1) + (length(g2) - 1)*var(g2)) / (length(g1) + length(g2) - 2))
    cohens_d <- (mean(g1) - mean(g2)) / s_pool

    cat(sprintf("Group Summaries: %s (n=%d, Mean=%.3f, SD=%.3f, SEM=%.3f) vs. %s (n=%d, Mean=%.3f, SD=%.3f, SEM=%.3f)\n\n",
                levels(x_fac)[1], length(g1), mean(g1), sd(g1), sd(g1)/sqrt(length(g1)),
                levels(x_fac)[2], length(g2), mean(g2), sd(g2), sd(g2)/sqrt(length(g2))))

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
    cat(sprintf("Effect Size -> Cohen's d = %.3f\n", cohens_d))
    if (show_code) .biostat_last_audit$last_code <- biostat_code(df, response = response, predictor = predictor)
    return(invisible(res))

  } else if (prob_type == "paired_2group") {
    x_fac <- as.factor(x)
    g1 <- y[x_fac == levels(x_fac)[1]]
    g2 <- y[x_fac == levels(x_fac)[2]]
    diffs <- g1 - g2
    dz <- mean(diffs) / sd(diffs)

    cat(sprintf("Condition 1 (%s): n = %d, Mean = %.3f, SD = %.3f, Median = %.3f\n",
                levels(x_fac)[1], length(g1), mean(g1), sd(g1), median(g1)))
    cat(sprintf("Condition 2 (%s): n = %d, Mean = %.3f, SD = %.3f, Median = %.3f\n",
                levels(x_fac)[2], length(g2), mean(g2), sd(g2), median(g2)))
    cat(sprintf("Paired Diff (%s - %s): Mean = %.3f, SD = %.3f, SEM = %.3f\n\n",
                levels(x_fac)[1], levels(x_fac)[2], mean(diffs), sd(diffs), sd(diffs)/sqrt(length(diffs))))

    if (assump$normality_ok) {
      cat("-> Normality of Paired Differences PASS. Running Parametric Paired t-test:\n")
      res <- t.test(g1, g2, paired = TRUE)
      print(res)
      cat("(Non-parametric Wilcoxon Signed-Rank also shown for reference):\n")
      print(suppressWarnings(wilcox.test(g1, g2, paired = TRUE, exact = FALSE)))
    } else {
      cat("-> Normality of Paired Differences FAIL. Running Wilcoxon Signed-Rank Test:\n")
      res <- suppressWarnings(wilcox.test(g1, g2, paired = TRUE, exact = FALSE))
      print(res)
      cat("(Parametric Paired t-test also shown for reference):\n")
      print(t.test(g1, g2, paired = TRUE))
    }
    cat(sprintf("Paired Effect Size -> Cohen's dz = %.3f\n", dz))
    if (show_code) .biostat_last_audit$last_code <- biostat_code(df, response = response, predictor = predictor, paired = TRUE)
    return(invisible(res))

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
      res <- fit_aov
    } else if (assump$normality_ok && !assump$variance_ok) {
      cat("-> Normality PASS, Equal Variance FAIL. Running Welch's One-Way ANOVA:\n")
      res <- oneway.test(y ~ x_fac, var.equal = FALSE)
      print(res)
      cat("\n--- Post-Hoc: Pairwise Welch t-tests (Holm-adjusted) ---\n")
      print(pairwise.t.test(y, x_fac, pool.sd = FALSE, p.adjust.method = "holm"))
    } else {
      cat("-> Normality FAIL. Running Non-Parametric Kruskal-Wallis Test:\n")
      res <- kruskal.test(y ~ x_fac)
      print(res)
      if (res$p.value < alpha) {
        cat("\n--- Post-Hoc: Pairwise Wilcoxon Rank-Sum Tests (Holm-adjusted) ---\n")
        print(pairwise.wilcox.test(y, x_fac, p.adjust.method = "holm", exact = FALSE))
      }
    }
    if (show_code) .biostat_last_audit$last_code <- biostat_code(df, response = response, predictor = predictor)
    return(invisible(res))

  } else if (prob_type == "continuous_xy") {
    cat("-> Fitting Simple Linear Regression (OLS: y ~ x):\n")
    fit_lm <- lm(y ~ x)
    print(anova(fit_lm))
    print(summary(fit_lm))

    if (assump$bivariate_normal) {
      cat("\n-> Bivariate Normality PASS. Running Pearson Correlation (r):\n")
      print(cor.test(x, y, method = "pearson"))
    } else {
      cat("\n-> Bivariate Normality FAIL. Running Spearman Rank Correlation (rho):\n")
      print(suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)))
      cat("(Pearson Correlation also shown for reference):\n")
      print(cor.test(x, y, method = "pearson"))
    }
    if (show_code) .biostat_last_audit$last_code <- biostat_code(df, response = response, predictor = predictor)
    return(invisible(fit_lm))

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

message("Universal Gbiostats_toolkit.r (v4 - Superpower + Self-Checking Code & Math Engine) loaded!")
