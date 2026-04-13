# Convert current-status MIC data into interval-censored data
prepare_ph_data <- function(
    df,
    delta_col = "Delta",
    indicator_col = "Indicator",
    lower_bound = 0,
    upper_bound = NULL
) {
  
  # Check that df is a data.frame
  if (!is.data.frame(df)) {
    stop("df must be a data.frame.")
  }
  
  # Define the required columns
  required_cols <- c(delta_col, indicator_col)
  
  # Identify any missing required columns
  missing_cols <- setdiff(required_cols, names(df))
  
  # Stop if required columns are missing
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Extract the PRO change column
  Delta <- df[[delta_col]]
  
  # Extract the binary indicator column
  Indicator <- df[[indicator_col]]
  
  # Check that Delta is numeric
  if (!is.numeric(Delta)) {
    stop(delta_col, " must be numeric.")
  }
  
  # Check that Delta does not contain missing values
  if (any(is.na(Delta))) {
    stop(delta_col, " cannot contain NA.")
  }
  
  # Check that Indicator contains only 0 and 1
  if (!all(Indicator %in% c(0, 1))) {
    stop(indicator_col, " must contain only 0 and 1.")
  }
  
  # Check that lower_bound is a single numeric value
  if (!is.numeric(lower_bound) || length(lower_bound) != 1 || is.na(lower_bound)) {
    stop("lower_bound must be a numeric scalar.")
  }
  
  # If upper_bound is provided, check that it is valid
  if (!is.null(upper_bound)) {
    
    # upper_bound must be a single numeric value and not missing
    if (!is.numeric(upper_bound) || length(upper_bound) != 1 || is.na(upper_bound)) {
      stop("upper_bound must be NULL or a numeric scalar.")
    }
    
    # upper_bound must be at least as large as the largest observed Delta
    if (upper_bound < max(Delta)) {
      stop("upper_bound must be at least as large as max(", delta_col, ").")
    }
  }
  
  # Make a copy of the original data
  out <- df
  
  # Construct the left endpoint L
  # If Indicator = 1, then tau <= Delta, so L = lower_bound
  # If Indicator = 0, then tau > Delta, so L = Delta
  out$L <- ifelse(Indicator == 1, lower_bound, Delta)
  
  # Construct the right endpoint R
  # If Indicator = 1, then tau <= Delta, so R = Delta
  # If Indicator = 0, then tau > Delta:
  #   - use upper_bound if known
  #   - otherwise use Inf
  out$R <- ifelse(
    Indicator == 1,
    Delta,
    if (is.null(upper_bound)) Inf else upper_bound
  )
  
  # Return the transformed data frame
  out
}



# Build a PH model formula from user-specified covariates
build_ph_formula <- function(
    df,
    categorical_covariates = NULL,
    continuous_covariates = NULL,
    exclude_cols = c("Delta", "T", "Indicator", "L", "R")
) {
  
  # Check that df is a data.frame
  if (!is.data.frame(df)) {
    stop("df must be a data.frame.")
  }
  
  # Require at least one of the two covariate inputs
  if (is.null(categorical_covariates) && is.null(continuous_covariates)) {
    stop("You must specify at least one of categorical_covariates or continuous_covariates.")
  }
  
  # Replace NULL with empty character vectors
  if (is.null(categorical_covariates)) {
    categorical_covariates <- character(0)
  }
  if (is.null(continuous_covariates)) {
    continuous_covariates <- character(0)
  }
  
  # Check that both inputs are character vectors
  if (!is.character(categorical_covariates)) {
    stop("categorical_covariates must be NULL or a character vector.")
  }
  if (!is.character(continuous_covariates)) {
    stop("continuous_covariates must be NULL or a character vector.")
  }
  
  # Combine all requested covariates
  all_covariates <- c(categorical_covariates, continuous_covariates)
  
  # Remove duplicated names
  all_covariates <- unique(all_covariates)
  
  # Check that all requested covariates exist in df
  missing_covariates <- setdiff(all_covariates, names(df))
  if (length(missing_covariates) > 0) {
    stop(
      "The following covariates are not in df: ",
      paste(missing_covariates, collapse = ", ")
    )
  }
  
  # Check that excluded columns are not incorrectly used as covariates
  bad_covariates <- intersect(all_covariates, exclude_cols)
  if (length(bad_covariates) > 0) {
    stop(
      "The following columns cannot be used as covariates here: ",
      paste(bad_covariates, collapse = ", ")
    )
  }
  
  # Check that categorical and continuous covariates do not overlap
  overlap_covariates <- intersect(categorical_covariates, continuous_covariates)
  if (length(overlap_covariates) > 0) {
    stop(
      "These covariates were declared as both categorical and continuous: ",
      paste(overlap_covariates, collapse = ", ")
    )
  }
  
  # Build formula terms for categorical covariates
  categorical_terms <- if (length(categorical_covariates) > 0) {
    paste0("factor(", categorical_covariates, ")")
  } else {
    character(0)
  }
  
  # Continuous covariates are used as-is
  continuous_terms <- continuous_covariates
  
  # Combine all terms
  formula_terms <- c(continuous_terms, categorical_terms)
  
  # If no usable terms remain, return intercept-only formula
  if (length(formula_terms) == 0) {
    return("1")
  }
  
  # Collapse all terms into a single formula string
  paste(formula_terms, collapse = " + ")
}


