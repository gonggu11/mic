#' Generate MIC simulation data
#'
#' @param n Positive integer sample size.
#' @param type "discrete" or "continuous".
#' @param setup Optional numeric setup. Used only to provide defaults.
#' @param seed Optional random seed.
#' @param covariates Optional data.frame of covariates. Must have n rows.
#' @param covariate_generator Optional function taking n and returning a data.frame.
#' @param delta_support Optional support for discrete Delta.
#' @param t_support Optional support for discrete T.
#' @param coefs_D Optional named numeric vector for Delta model.
#' @param coefs_T Optional named numeric vector for T model.
#'        Special name "Delta" is allowed in coefs_T.
#' @param sd_D Optional positive numeric scalar.
#' @param sd_T Optional positive numeric scalar.
#' @param return_T Logical; if FALSE, drop T before returning.
#'
#' @return data.frame
#' @export
generate_mic_data <- function(
    n,
    type = c("discrete", "continuous"),
    setup = NULL,
    seed = NULL,
    covariates = NULL,
    covariate_generator = NULL,
    delta_support = NULL,
    t_support = NULL,
    coefs_D = NULL,
    coefs_T = NULL,
    sd_D = NULL,
    sd_T = NULL,
    return_T = TRUE
) {
  type <- match.arg(type)
  
  if (!is.null(seed)) set.seed(seed)
  
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n <= 0) {
    stop("n must be a positive integer.")
  }
  n <- as.integer(n)
  
  if (!is.null(setup) && (!is.numeric(setup) || length(setup) != 1 || is.na(setup))) {
    stop("setup must be NULL or a numeric scalar.")
  }
  
  if (!is.logical(return_T) || length(return_T) != 1 || is.na(return_T)) {
    stop("return_T must be TRUE or FALSE.")
  }
  
  # -------------------------
  # helper: default setup
  # -------------------------
  get_default_setup <- function(type, setup) {
    if (is.null(setup)) setup <- 0
    
    if (type == "discrete") {
      out <- list(
        delta_support = -5:10,
        t_support = 1:10,
        coefs_D = c(intercept = 1, sex = 0.5, race = 0, cont = 0.01),
        coefs_T = c(intercept = 1, Delta = 0, sex = 0.4, race = -0.2, cont = -0.01),
        sd_D = 5,
        sd_T = 2
      )
      
      if (setup == 0.1) {
        out$coefs_T["Delta"] <- 1
        out$sd_D <- 2
      } else if (setup == 0.2) {
        out$coefs_T["Delta"] <- 2
        out$sd_D <- 2
      } else if (setup == 0.3) {
        out$coefs_T["Delta"] <- 3
        out$sd_D <- 2
      } else if (setup == 0.4) {
        out$coefs_T["Delta"] <- 4
        out$sd_D <- 2
      } else if (setup == 0.5) {
        out$coefs_T["Delta"] <- 1.5
        out$sd_D <- 2
      } else if (setup == 1) {
        out$coefs_D["intercept"] <- 0
        out$coefs_T["Delta"] <- 0.8
      } else if (setup == 1.1) {
        out$coefs_T["sex"] <- 3
      } else if (setup != 0) {
        stop("Unsupported setup.")
      }
      
    } else {
      out <- list(
        delta_support = NULL,
        t_support = NULL,
        coefs_D = c(intercept = 1, sex = 0.5, race = 0, cont = 0.01),
        coefs_T = c(intercept = 1, Delta = 0, sex = 0.4, race = -0.2, cont = -0.01),
        sd_D = 5,
        sd_T = 2
      )
      
      if (setup == 0.1) {
        out$coefs_T["Delta"] <- 1
        out$sd_D <- 2
      } else if (setup == 0.2) {
        out$coefs_T["Delta"] <- 2
        out$sd_D <- 2
      } else if (setup == 0.3) {
        out$coefs_T["Delta"] <- 3
        out$sd_D <- 2
      } else if (setup == 0.4) {
        out$coefs_T["Delta"] <- 4
        out$sd_D <- 2
      } else if (setup == 0.5) {
        out$coefs_T["Delta"] <- 1.5
        out$sd_D <- 2
      } else if (setup == 1) {
        out$coefs_D["intercept"] <- 0
        out$coefs_T["Delta"] <- 0.8
      } else if (setup == 1.1) {
        out$coefs_T["sex"] <- 3
      } else if (setup != 0) {
        stop("Unsupported setup.")
      }
    }
    
    out
  }
  
  defaults <- get_default_setup(type = type, setup = setup)
  
  # -------------------------
  # helper: user overrides
  # -------------------------
  if (is.null(delta_support)) delta_support <- defaults$delta_support
  if (is.null(t_support)) t_support <- defaults$t_support
  if (is.null(coefs_D)) coefs_D <- defaults$coefs_D
  if (is.null(coefs_T)) coefs_T <- defaults$coefs_T
  if (is.null(sd_D)) sd_D <- defaults$sd_D
  if (is.null(sd_T)) sd_T <- defaults$sd_T
  
  # -------------------------
  # covariates
  # -------------------------
  if (!is.null(covariates) && !is.null(covariate_generator)) {
    stop("Provide only one of covariates and covariate_generator.")
  }
  
  if (!is.null(covariates)) {
    if (!is.data.frame(covariates)) stop("covariates must be a data.frame.")
    if (nrow(covariates) != n) stop("covariates must have n rows.")
    X <- covariates
  } else if (!is.null(covariate_generator)) {
    X <- covariate_generator(n)
    if (!is.data.frame(X)) stop("covariate_generator(n) must return a data.frame.")
    if (nrow(X) != n) stop("covariate_generator(n) must return n rows.")
  } else {
    X <- data.frame(
      sex = sample(0:1, n, replace = TRUE),
      race = sample(0:2, n, replace = TRUE),
      cont = rnorm(n, mean = 65, sd = 10)
    )
  }
  
  # -------------------------
  # helper: linear predictor
  # -------------------------
  eval_lp <- function(coefs, Xrow, Delta_value = NULL) {
    val <- 0
    
    if ("intercept" %in% names(coefs)) {
      val <- val + coefs["intercept"]
    }
    
    if ("Delta" %in% names(coefs)) {
      if (is.null(Delta_value)) stop("Delta_value is required for coefs containing 'Delta'.")
      val <- val + coefs["Delta"] * Delta_value
    }
    
    cov_names <- setdiff(names(coefs), c("intercept", "Delta"))
    
    for (nm in cov_names) {
      if (nm %in% names(Xrow)) {
        val <- val + coefs[nm] * Xrow[[nm]]
      }
    }
    
    as.numeric(val)
  }
  
  # -------------------------
  # checks
  # -------------------------
  if (!is.numeric(sd_D) || length(sd_D) != 1 || is.na(sd_D) || sd_D <= 0) {
    stop("sd_D must be a positive numeric scalar.")
  }
  if (!is.numeric(sd_T) || length(sd_T) != 1 || is.na(sd_T) || sd_T <= 0) {
    stop("sd_T must be a positive numeric scalar.")
  }
  
  if (type == "discrete") {
    if (!is.numeric(delta_support) || length(delta_support) < 2) {
      stop("delta_support must be a numeric vector with length >= 2.")
    }
    if (!is.numeric(t_support) || length(t_support) < 2) {
      stop("t_support must be a numeric vector with length >= 2.")
    }
  }
  
  # -------------------------
  # generate Delta
  # -------------------------
  Delta <- numeric(n)
  
  if (type == "discrete") {
    for (i in seq_len(n)) {
      mu_D <- eval_lp(coefs_D, X[i, , drop = FALSE])
      
      p <- dnorm(delta_support, mean = mu_D, sd = sd_D)
      p <- p / sum(p)
      
      Delta[i] <- sample(delta_support, size = 1, prob = p)
    }
  } else {
    for (i in seq_len(n)) {
      mu_D <- eval_lp(coefs_D, X[i, , drop = FALSE])
      Delta[i] <- rnorm(1, mean = mu_D, sd = sd_D)
    }
  }
  
  # -------------------------
  # generate T
  # -------------------------
  T <- numeric(n)
  
  if (type == "discrete") {
    for (i in seq_len(n)) {
      mu_T <- eval_lp(coefs_T, X[i, , drop = FALSE], Delta_value = Delta[i])
      
      p <- dnorm(t_support, mean = mu_T, sd = sd_T)
      p <- p / sum(p)
      
      T[i] <- sample(t_support, size = 1, prob = p)
    }
  } else {
    for (i in seq_len(n)) {
      mu_T <- eval_lp(coefs_T, X[i, , drop = FALSE], Delta_value = Delta[i])
      T[i] <- rnorm(1, mean = mu_T, sd = sd_T)
    }
  }
  
  Indicator <- as.integer(Delta >= T)
  
  df <- cbind(
    data.frame(Delta = Delta, T = T, Indicator = Indicator),
    X
  )
  
  if (!return_T) {
    df$T <- NULL
  }
  
  rownames(df) <- NULL
  return(df)
}