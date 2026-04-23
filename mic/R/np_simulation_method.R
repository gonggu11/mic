library(icenReg)
library(survival)

prepare_np_data <- function(
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

# No formula validator is needed for the NP model,because the NPMLE fit does not use a user-supplied covariate formula.

fit_np_model <- function(
    prepared_obj,
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
  left_name <- make_safe_name(".np_left", names(fit_data))
  right_name <- make_safe_name(".np_right", c(names(fit_data), left_name))
  
  # Add the two interval endpoint columns to fit_data
  fit_data[[left_name]] <- interval_data[[1]]
  fit_data[[right_name]] <- interval_data[[2]]
  
  # Fit the NP interval-censored model
  fitted_model <- icenReg::ic_np(
    cbind(fit_data[[left_name]], fit_data[[right_name]])
  )
  
  # Return a structured fitted object
  out <- list(
    original_data = original_data,
    interval_data = interval_data,
    fit_data = fit_data,
    fitted_model = fitted_model,
    internal_left_name = left_name,
    internal_right_name = right_name,
    method = "np"
  )
  
  class(out) <- "np_model_fit"
  out
}

# Extract population MIC from a fitted NP / NPMLE model
extract_np_mic <- function(
    np_fit_obj,
    data_type = c("discrete", "continuous"),
    support = NULL,
    integration_upper = NULL
) {
  
  data_type <- match.arg(data_type)
  
  if (!is.list(np_fit_obj)) {
    stop("np_fit_obj must be a list.")
  }
  
  if (!("fitted_model" %in% names(np_fit_obj))) {
    stop("np_fit_obj must contain a fitted_model component.")
  }
  
  fitted_model <- np_fit_obj$fitted_model
  
  # -----------------------------
  # Discrete case
  # -----------------------------
  if (data_type == "discrete") {
    
    if (is.null(support)) {
      stop("For discrete data, support must be provided.")
    }
    
    if (!is.numeric(support) || length(support) < 2 || any(is.na(support))) {
      stop("support must be a numeric vector with length >= 2 and no NA.")
    }
    
    support <- sort(unique(support))
    
    if (is.null(fitted_model$p_hat)) {
      stop("fitted_model does not contain p_hat; cannot extract discrete pmf.")
    }
    
    if (is.null(fitted_model$T_bull_Intervals)) {
      stop("fitted_model does not contain T_bull_Intervals; cannot extract discrete pmf.")
    }
    
    p_hat <- fitted_model$p_hat
    tbull <- fitted_model$T_bull_Intervals
    
    if (!is.matrix(tbull) || nrow(tbull) != 2) {
      stop("fitted_model$T_bull_Intervals must be a 2 x k matrix.")
    }
    
    upper_endpoints <- tbull[2, ]
    
    if (length(p_hat) != length(upper_endpoints)) {
      stop("Lengths of p_hat and Turnbull intervals do not match.")
    }
    
    # For the discrete MIC setting, we map each NPMLE mass to the
    # right endpoint of its Turnbull interval.
    if (!all(upper_endpoints %in% support)) {
      stop("Some Turnbull interval upper endpoints are not contained in support.")
    }
    
    pmf <- numeric(length(support))
    names(pmf) <- as.character(support)
    
    for (j in seq_along(p_hat)) {
      key <- as.character(upper_endpoints[j])
      pmf[key] <- pmf[key] + p_hat[j]
    }
    
    # Numerical cleanup
    pmf <- pmf / sum(pmf)
    
    population_mic <- sum(support * pmf)
    
    out <- list(
      population_mic = population_mic,
      pmf = pmf,
      support = support,
      data_type = "discrete",
      method = "np"
    )
    
    class(out) <- "np_mic_estimate"
    return(out)
  }
  
  # -----------------------------
  # Continuous case
  # -----------------------------
  if (is.null(integration_upper)) {
    stop("For continuous data, integration_upper must be provided.")
  }
  
  if (!is.numeric(integration_upper) || length(integration_upper) != 1 || is.na(integration_upper)) {
    stop("integration_upper must be a numeric scalar.")
  }
  
  sc_obj <- icenReg::getSCurves(fit = fitted_model)
  
  Tbull_ints <- sc_obj$Tbull_ints
  S_curves <- sc_obj$S_curves
  
  if (is.null(Tbull_ints) || !is.matrix(Tbull_ints) || ncol(Tbull_ints) != 2) {
    stop("Could not extract valid Turnbull intervals from the fitted model.")
  }
  
  if (is.null(S_curves) || length(S_curves) < 1) {
    stop("Could not extract survival curve(s) from the fitted model.")
  }
  
  # ic_np without covariates should yield a single survival curve
  S_curve <- S_curves[[1]]
  
  if (length(S_curve) != nrow(Tbull_ints)) {
    stop("Length of survival curve does not match the number of Turnbull intervals.")
  }
  
  left <- Tbull_ints[, 1]
  right <- Tbull_ints[, 2]
  
  right[is.infinite(right)] <- integration_upper
  
  widths <- right - left
  if (any(widths < 0)) {
    stop("Invalid Turnbull intervals: some widths are negative.")
  }
  
  population_mic <- sum(S_curve * widths)
  
  out <- list(
    population_mic = population_mic,
    Tbull_ints = Tbull_ints,
    S_curve = S_curve,
    data_type = "continuous",
    integration_upper = integration_upper,
    method = "np"
  )
  
  class(out) <- "np_mic_estimate"
  return(out)
}
