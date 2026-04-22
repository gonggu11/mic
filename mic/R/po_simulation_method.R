library(icenReg)
library(survival)

prepare_po_data <- function(
    df,
    delta_col = "Delta",
    indicator_col = "Indicator",
    lower_bound = 0,
    upper_bound = NULL
) {
  
  # Check that the input is a data.frame
  if (!is.data.frame(df)) {
    stop("df must be a data.frame.")
  }
  
  # Check that the required columns exist
  required_cols <- c(delta_col, indicator_col)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Extract Delta and Indicator
  Delta <- df[[delta_col]]
  Indicator <- df[[indicator_col]]
  
  # Check that Delta is numeric
  if (!is.numeric(Delta)) {
    stop(delta_col, " must be numeric.")
  }
  
  # Check that Delta has no missing values
  if (any(is.na(Delta))) {
    stop(delta_col, " cannot contain NA.")
  }
  
  # Check that Indicator contains only 0 and 1
  if (!all(Indicator %in% c(0, 1))) {
    stop(indicator_col, " must contain only 0 and 1.")
  }
  
  # Check that lower_bound is a numeric scalar
  if (!is.numeric(lower_bound) || length(lower_bound) != 1 || is.na(lower_bound)) {
    stop("lower_bound must be a numeric scalar.")
  }
  
  # Validate upper_bound if it is provided
  if (!is.null(upper_bound)) {
    if (!is.numeric(upper_bound) || length(upper_bound) != 1 || is.na(upper_bound)) {
      stop("upper_bound must be NULL or a numeric scalar.")
    }
    
    if (upper_bound < max(Delta)) {
      stop("upper_bound must be at least as large as max(", delta_col, ").")
    }
  }
  
  # Build interval-censored representation
  # If Indicator = 1, then T <= Delta, so the interval is (lower_bound, Delta]
  # If Indicator = 0, then T > Delta, so the interval is (Delta, upper_bound]
  # or (Delta, Inf) if upper_bound is not provided
  interval_data <- data.frame(
    L = ifelse(Indicator == 1, lower_bound, Delta),
    R = ifelse(Indicator == 1, Delta, if (is.null(upper_bound)) Inf else upper_bound)
  )
  
  # Return the original data and the interval data separately
  list(
    original_data = df,
    interval_data = interval_data
  )
}

# Validate a user-supplied PO formula string
validate_po_formula <- function(
    df,
    user_formula
) {
  
  # Check that df is a data.frame
  if (!is.data.frame(df)) {
    stop("df must be a data.frame.")
  }
  
  # Check that user_formula is provided
  if (missing(user_formula) || is.null(user_formula)) {
    stop("user_formula must be provided.")
  }
  
  # Check that user_formula is a single character string
  if (!is.character(user_formula) || length(user_formula) != 1 || is.na(user_formula)) {
    stop("user_formula must be a single character string.")
  }
  
  # Remove leading and trailing whitespace
  user_formula <- trimws(user_formula)
  
  # Check whether the cleaned formula string is empty
  if (nchar(user_formula) == 0) {
    stop("user_formula cannot be empty.")
  }
  
  # Try to convert the user-supplied right-hand side into a valid R formula object
  rhs_formula <- tryCatch(
    as.formula(paste("~", user_formula)),
    error = function(e) {
      stop("Invalid user_formula. Unable to parse the formula string.")
    }
  )
  
  # Extract all variable names that appear in the parsed formula
  formula_vars <- all.vars(rhs_formula)
  
  # Identify variables used in the formula but missing from df
  missing_vars <- setdiff(formula_vars, names(df))
  
  # Stop if any variable in the formula is not found in df
  if (length(missing_vars) > 0) {
    stop(
      "The following variables in user_formula are not in df: ",
      paste(missing_vars, collapse = ", ")
    )
  }
  
  # Return the validated formula string
  user_formula
}




library(icenReg)
library(survival)

# Fit a PH interval-censored model using prepared data and a user-supplied formula
# Fit a PO interval-censored model using prepared data and a user-supplied formula
fit_po_model <- function(
    prepared_obj,
    user_formula,
    validate_formula_input = TRUE,
    drop_truth_columns = TRUE
) {
  
  # Check that prepared_obj is a list
  if (!is.list(prepared_obj)) {
    stop("prepared_obj must be a list.")
  }
  
  # Check that prepared_obj contains original_data and interval_data
  required_parts <- c("original_data", "interval_data")
  missing_parts <- setdiff(required_parts, names(prepared_obj))
  if (length(missing_parts) > 0) {
    stop(
      "prepared_obj is missing required components: ",
      paste(missing_parts, collapse = ", ")
    )
  }
  
  # Extract original data and interval data
  original_data <- prepared_obj$original_data
  interval_data <- prepared_obj$interval_data
  
  # Check that both are data.frames
  if (!is.data.frame(original_data)) {
    stop("prepared_obj$original_data must be a data.frame.")
  }
  if (!is.data.frame(interval_data)) {
    stop("prepared_obj$interval_data must be a data.frame.")
  }
  
  # Check that interval_data has exactly two columns
  if (ncol(interval_data) != 2) {
    stop("prepared_obj$interval_data must contain exactly two columns.")
  }
  
  # Check that original_data and interval_data have the same number of rows
  if (nrow(original_data) != nrow(interval_data)) {
    stop("original_data and interval_data must have the same number of rows.")
  }
  
  # Validate the user-supplied formula if requested
  if (validate_formula_input) {
    user_formula <- validate_po_formula(
      df = original_data,
      user_formula = user_formula
    )
  }
  
  # Start from original_data when building the fitting dataset
  fit_data <- original_data
  
  # Optionally remove latent truth columns such as T
  if (drop_truth_columns) {
    truth_cols <- intersect(c("T"), names(fit_data))
    if (length(truth_cols) > 0) {
      fit_data[truth_cols] <- NULL
    }
  }
  
  # Create a helper that generates a guaranteed non-conflicting internal column name
  make_safe_name <- function(base_name, existing_names) {
    candidate <- base_name
    counter <- 1
    
    while (candidate %in% existing_names) {
      candidate <- paste0(base_name, "_", counter)
      counter <- counter + 1
    }
    
    candidate
  }
  
  # Generate safe internal names for the two interval endpoints
  left_name <- make_safe_name(".po_left", names(fit_data))
  right_name <- make_safe_name(".po_right", c(names(fit_data), left_name))
  
  # Add the two interval endpoint columns to fit_data
  fit_data[[left_name]] <- interval_data[[1]]
  fit_data[[right_name]] <- interval_data[[2]]
  
  # Create an environment in which Surv is available
  formula_env <- new.env(parent = parent.frame())
  formula_env$Surv <- survival::Surv
  
  # Build the full PO model formula dynamically
  full_formula <- as.formula(
    paste0(
      "Surv(",
      left_name,
      ", ",
      right_name,
      ", type = 'interval2') ~ ",
      user_formula
    ),
    env = formula_env
  )
  
  # Fit the PO interval-censored model
  fitted_model <- icenReg::ic_sp(
    formula = full_formula,
    data = fit_data,
    model = "po"
  )
  
  # Return a structured fitted object
  out <- list(
    original_data = original_data,
    interval_data = interval_data,
    fit_data = fit_data,
    user_formula = user_formula,
    full_formula = full_formula,
    fitted_model = fitted_model,
    internal_left_name = left_name,
    internal_right_name = right_name,
    method = "po"
  )
  
  class(out) <- "po_model_fit"
  out
}


# Extract conditional and population MIC from a fitted PO model
extract_po_mic <- function(
    po_fit_obj,
    data_type = c("discrete", "continuous"),
    support = NULL,
    integration_upper = NULL
) {
  
  data_type <- match.arg(data_type)
  
  if (!is.list(po_fit_obj)) {
    stop("po_fit_obj must be a list.")
  }
  
  if (!("fitted_model" %in% names(po_fit_obj))) {
    stop("po_fit_obj must contain a fitted_model component.")
  }
  
  if (!("original_data" %in% names(po_fit_obj))) {
    stop("po_fit_obj must contain an original_data component.")
  }
  
  fitted_model <- po_fit_obj$fitted_model
  original_data <- po_fit_obj$original_data
  
  # Discrete case: estimate pmf on a finite support and then take expectation
  if (data_type == "discrete") {
    
    if (is.null(support)) {
      stop("For discrete data, support must be provided.")
    }
    
    if (!is.numeric(support) || length(support) < 2) {
      stop("support must be a numeric vector with length >= 2.")
    }
    
    # Placeholder:
    # Here you would later build a PO-based pmf extractor on the finite support.
    # For now, we keep the interface explicit and stop.
    stop("Discrete PO pmf extraction is not yet implemented. This should return p(t | X) on the supplied support.")
    
  } else {
    
    # Continuous case: use survival curves and integrate S(t | X) over t
    if (is.null(integration_upper)) {
      stop("For continuous data, integration_upper must be provided.")
    }
    
    if (!is.numeric(integration_upper) || length(integration_upper) != 1 || is.na(integration_upper)) {
      stop("integration_upper must be a numeric scalar.")
    }
    
    sc_obj <- icenReg::getSCurves(
      fit = fitted_model,
      newdata = original_data
    )
    
    Tbull_ints <- sc_obj$Tbull_ints
    S_curves <- sc_obj$S_curves
    
    if (is.null(Tbull_ints) || ncol(Tbull_ints) != 2) {
      stop("Could not extract valid Turnbull intervals from the fitted model.")
    }
    
    mean_from_surv_curve <- function(S_curve, intervals, upper) {
      left <- intervals[, 1]
      right <- intervals[, 2]
      
      right[is.infinite(right)] <- upper
      
      widths <- right - left
      sum(S_curve * widths)
    }
    
    conditional_mic <- vapply(
      S_curves,
      FUN = mean_from_surv_curve,
      FUN.VALUE = numeric(1),
      intervals = Tbull_ints,
      upper = integration_upper
    )
    
    population_mic <- mean(conditional_mic)
    
    out <- list(
      conditional_mic = conditional_mic,
      population_mic = population_mic,
      Tbull_ints = Tbull_ints,
      S_curves = S_curves,
      data_type = "continuous",
      integration_upper = integration_upper,
      method = "po"
    )
    
    class(out) <- "po_mic_estimate"
    return(out)
  }
}