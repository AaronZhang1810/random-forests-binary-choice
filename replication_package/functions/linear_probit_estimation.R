##############################################################################

# Random Forests for Nonsmooth Binary Choice Models
# Simulation Study
# Author: Jihuan Zhang, Department of Economics, Indiana University Bloomington
# Email: zhangjih@iu.edu

##############################################################################


##############################################################################
######### Functions for linear probit #########
##############################################################################
make_linear_design <- function(X) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"

  if (nrow(X) < 1L || ncol(X) < 1L) {
    stop("X must contain at least one observation and one covariate.")
  }

  if (any(!is.finite(X))) {
    stop("X contains nonfinite values.")
  }

  n <- nrow(X)
  d <- ncol(X)

  variable_names <- colnames(X)
  if (is.null(variable_names)) {
    variable_names <- paste0("X", seq_len(d))
  }
  variable_names <- make.names(variable_names, unique = TRUE)

  intercept <- matrix(1, nrow = n, ncol = 1L)
  colnames(intercept) <- "(Intercept)"

  linear_terms <- X
  colnames(linear_terms) <- variable_names

  design_matrix <- cbind(
    intercept,
    linear_terms
  )

  storage.mode(design_matrix) <- "double"
  design_matrix
}


# Estimate the conventional linear probit model by Bernoulli maximum
# likelihood. glm.fit() uses iteratively reweighted least squares as the
# numerical algorithm, but the optimized criterion is the Bernoulli
# log-likelihood with a probit link.
fit_linear_probit_mle <- function(X,
                                  Y,
                                  max_iterations = 100L,
                                  tolerance = 1e-8) {
  Y <- as.numeric(Y)

  if (length(Y) != nrow(X)) {
    stop("The lengths of X and Y are incompatible.")
  }

  if (any(!Y %in% c(0, 1))) {
    stop("Y must be binary and coded as 0/1.")
  }

  if (length(unique(Y)) < 2L) {
    stop("The estimation sample contains only one outcome class.")
  }

  design_matrix <- make_linear_design(X)
  warning_messages <- character(0)

  fit <- withCallingHandlers(
    stats::glm.fit(
      x = design_matrix,
      y = Y,
      family = stats::binomial(link = "probit"),
      control = stats::glm.control(
        epsilon = tolerance,
        maxit = max_iterations,
        trace = FALSE
      ),
      intercept = TRUE,
      singular.ok = FALSE
    ),
    warning = function(w) {
      warning_messages <<- c(
        warning_messages,
        conditionMessage(w)
      )
      invokeRestart("muffleWarning")
    }
  )

  if (!isTRUE(fit$converged)) {
    stop("The conventional linear probit MLE did not converge.")
  }

  if (length(fit$coefficients) != ncol(design_matrix)) {
    stop("The fitted coefficient vector has the wrong dimension.")
  }

  if (any(!is.finite(fit$coefficients))) {
    stop("The conventional linear probit MLE produced nonfinite coefficients.")
  }

  list(
    coefficients = fit$coefficients,
    converged = isTRUE(fit$converged),
    boundary = isTRUE(fit$boundary),
    iterations = as.integer(fit$iter),
    rank = as.integer(fit$rank),
    n_parameters = as.integer(length(fit$coefficients)),
    deviance = as.numeric(fit$deviance),
    aic = as.numeric(fit$aic),
    warning_message = if (length(warning_messages) == 0L) {
      NA_character_
    } else {
      paste(unique(warning_messages), collapse = " | ")
    }
  )
}


# Estimate linear probit by nonlinear least squares on the probabilities:
#   sum_i (Y_i - Phi(z_i' beta))^2.
# A constant-variance quasi family gives this squared-error criterion in
# glm.fit(). The nonlinear criterion may have multiple local minima.
fit_linear_probit_ls <- function(X,
                                 Y,
                                 max_iterations = 1000L,
                                 tolerance = 1e-8) {
  design_matrix <- make_linear_design(X)
  Y <- as.numeric(Y)

  if (length(Y) != nrow(design_matrix)) {
    stop("The lengths of X and Y are incompatible.")
  }

  if (any(!Y %in% c(0, 1))) {
    stop("Y must be binary and coded as 0/1.")
  }

  if (length(unique(Y)) < 2L) {
    stop("The estimation sample contains only one outcome class.")
  }

  if (length(max_iterations) != 1L || !is.numeric(max_iterations) ||
      !is.finite(max_iterations) || max_iterations < 1 ||
      max_iterations != floor(max_iterations) ||
      max_iterations > .Machine$integer.max) {
    stop("max_iterations must be a positive integer <= .Machine$integer.max.")
  }

  if (length(tolerance) != 1L || !is.numeric(tolerance) ||
      !is.finite(tolerance) || tolerance <= 0 || tolerance >= 1) {
    stop("tolerance must be finite and strictly between zero and one.")
  }

  # Binary Y would give infinite probit indices under the default quasi
  # initialization. Start from the finite best constant probability instead.
  initial_coefficients <- setNames(
    numeric(ncol(design_matrix)),
    colnames(design_matrix)
  )
  initial_coefficients[1L] <- stats::qnorm(mean(Y))
  warning_messages <- character(0)

  fit <- withCallingHandlers(
    stats::glm.fit(
      x = design_matrix,
      y = Y,
      family = stats::quasi(link = "probit", variance = "constant"),
      start = initial_coefficients,
      control = stats::glm.control(
        epsilon = tolerance,
        maxit = as.integer(max_iterations),
        trace = FALSE
      ),
      intercept = TRUE,
      singular.ok = FALSE
    ),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (!isTRUE(fit$converged)) {
    stop("The linear probit least-squares fit did not converge.")
  }

  if (length(fit$coefficients) != ncol(design_matrix) ||
      any(!is.finite(fit$coefficients))) {
    stop("The linear probit least-squares fit produced invalid coefficients.")
  }

  # Evaluate the loss using the same pnorm predictions as the prediction helper.
  # glm.fit() clips extreme probit tails internally at machine precision.
  g_hat <- drop(design_matrix %*% fit$coefficients)
  p_hat <- stats::pnorm(g_hat)
  residual <- p_hat - Y
  loss <- sum(residual^2)
  gradient_norm <- max(abs(
    2 * crossprod(design_matrix, residual * stats::dnorm(g_hat))
  ))
  if (any(!is.finite(g_hat)) || !is.finite(loss) || !is.finite(gradient_norm)) {
    stop("The linear probit least-squares fit produced a nonfinite result.")
  }

  probability_tolerance <- 10 * .Machine$double.eps
  if (any(p_hat < probability_tolerance | p_hat > 1 - probability_tolerance)) {
    warning_messages <- c(
      warning_messages,
      paste0(
        "Fitted probabilities are numerically 0 or 1; ",
        "the least-squares criterion can be flat in these regions."
      )
    )
  }

  list(
    coefficients = fit$coefficients,
    converged = isTRUE(fit$converged),
    boundary = isTRUE(fit$boundary),
    iterations = as.integer(fit$iter),
    rank = as.integer(fit$rank),
    n_parameters = as.integer(length(fit$coefficients)),
    # Bernoulli likelihood diagnostics do not apply to this fit; report SSE
    # separately as loss, consistent with fit_spp_ls().
    deviance = NA_real_,
    aic = NA_real_,
    loss = as.numeric(loss),
    gradient_norm = as.numeric(gradient_norm),
    fit_method = "least_squares",
    warning_message = if (length(warning_messages) == 0L) {
      NA_character_
    } else {
      paste(unique(warning_messages), collapse = " | ")
    }
  )
}


predict_linear_probit <- function(object, newdata) {
  design_matrix <- make_linear_design(newdata)

  if (ncol(design_matrix) != length(object$coefficients)) {
    stop("The prediction design matrix is incompatible with the fitted model.")
  }

  g_hat <- drop(design_matrix %*% object$coefficients)
  p_hat <- stats::pnorm(g_hat)

  if (any(!is.finite(g_hat)) || any(!is.finite(p_hat))) {
    stop("The conventional linear probit estimator produced nonfinite predictions.")
  }

  list(
    g_hat = g_hat,
    p_hat = p_hat
  )
}



