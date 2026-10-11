##############################################################################

# Random Forests for Nonsmooth Binary Choice Models
# Simulation
# Author: Jihuan Zhang, Department of Economics, Indiana University Bloomington
# Email: zhangjih@iu.edu

##############################################################################



##############################################################################
######### Functions for second-order polynomial probit #########
##############################################################################
# Construct the complete polynomial basis of total degree at most two:
# intercept, all linear terms, all squared terms, and all pairwise interactions.
make_poly2_design <- function(X) {
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
  
  squared_terms <- X^2
  colnames(squared_terms) <- paste0(variable_names, "_sq")
  
  if (d >= 2L) {
    interaction_index <- utils::combn(d, 2L)
    
    interaction_terms <-
      X[, interaction_index[1L, ], drop = FALSE] *
      X[, interaction_index[2L, ], drop = FALSE]
    
    colnames(interaction_terms) <- paste0(
      variable_names[interaction_index[1L, ]],
      "_x_",
      variable_names[interaction_index[2L, ]]
    )
    
    design_matrix <- cbind(
      intercept,
      linear_terms,
      squared_terms,
      interaction_terms
    )
  } else {
    design_matrix <- cbind(
      intercept,
      linear_terms,
      squared_terms
    )
  }
  
  storage.mode(design_matrix) <- "double"
  design_matrix
}


# Fit SPP by MLE
fit_spp_mle <- function(X,
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
  
  design_matrix <- make_poly2_design(X)
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
    stop("The SPP probit MLE did not converge.")
  }
  
  if (length(fit$coefficients) != ncol(design_matrix)) {
    stop("The fitted coefficient vector has the wrong dimension.")
  }
  
  if (any(!is.finite(fit$coefficients))) {
    stop("The SPP probit MLE produced nonfinite coefficients.")
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


# Fit SPP by nonlinear least squares: sum_i (Y_i - Phi(z_i' beta))^2
# The analytic gradient is 2 * Z' ((Phi(Z beta) - Y) * phi(Z beta)).
# This criterion is nonconvex; nlminb() searches for a local minimum.
fit_spp_ls <- function(X,
                       Y,
                       max_iterations = 1000L,
                       tolerance = 1e-8) {
  design_matrix <- make_poly2_design(X)
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
      max_iterations > .Machine$integer.max / 2) {
    stop("max_iterations must be a positive integer <= .Machine$integer.max / 2.")
  }
  
  if (length(tolerance) != 1L || !is.numeric(tolerance) ||
      !is.finite(tolerance) || tolerance <= 0) {
    stop("tolerance must be positive and finite.")
  }
  
  if (any(!is.finite(design_matrix))) {
    stop("The polynomial design matrix contains nonfinite values.")
  }
  
  design_rank <- qr(design_matrix)$rank
  n_parameters <- ncol(design_matrix)
  if (design_rank < n_parameters) {
    stop("The SPP design matrix is rank deficient.")
  }
  
  # Start from the best constant probability, without requiring an MLE fit.
  initial_coefficients <- setNames(numeric(n_parameters), colnames(design_matrix))
  initial_coefficients[1L] <- stats::qnorm(mean(Y))
  
  squared_loss <- function(coefficients) {
    g_hat <- drop(design_matrix %*% coefficients)
    sum((Y - stats::pnorm(g_hat))^2)
  }
  
  loss_gradient <- function(coefficients) {
    g_hat <- drop(design_matrix %*% coefficients)
    residual <- stats::pnorm(g_hat) - Y
    drop(2 * crossprod(design_matrix, residual * stats::dnorm(g_hat)))
  }
  
  warning_messages <- character(0)
  fit <- withCallingHandlers(
    stats::nlminb(
      start = initial_coefficients,
      objective = squared_loss,
      gradient = loss_gradient,
      control = list(
        iter.max = as.integer(max_iterations),
        eval.max = as.integer(2 * max_iterations),
        rel.tol = tolerance
      )
    ),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  
  if (fit$convergence != 0L) {
    stop("The SPP least-squares fit did not converge: ", fit$message)
  }
  
  if (length(fit$par) != n_parameters || any(!is.finite(fit$par)) ||
      !is.finite(fit$objective)) {
    stop("The SPP least-squares fit produced a nonfinite or incompatible result.")
  }
  
  g_hat <- drop(design_matrix %*% fit$par)
  p_hat <- stats::pnorm(g_hat)
  gradient_norm <- max(abs(loss_gradient(fit$par)))
  if (any(!is.finite(g_hat)) || !is.finite(gradient_norm)) {
    stop("The SPP least-squares fit produced a nonfinite index or gradient.")
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
    coefficients = setNames(fit$par, colnames(design_matrix)),
    converged = TRUE,
    # These GLM diagnostics are not defined for the least-squares fit.
    boundary = NA,
    iterations = as.integer(fit$iterations),
    rank = as.integer(design_rank),
    n_parameters = as.integer(n_parameters),
    deviance = NA_real_,
    aic = NA_real_,
    loss = as.numeric(fit$objective),
    gradient_norm = gradient_norm,
    fit_method = "least_squares",
    optimizer_message = fit$message,
    warning_message = if (length(warning_messages) == 0L) {
      NA_character_
    } else {
      paste(unique(warning_messages), collapse = " | ")
    }
  )
}

# Make predictions on the test set
predict_spp <- function(object, newdata) {
  design_matrix <- make_poly2_design(newdata)
  
  if (ncol(design_matrix) != length(object$coefficients)) {
    stop("The prediction design matrix is incompatible with the fitted model.")
  }
  
  g_hat <- drop(design_matrix %*% object$coefficients)
  p_hat <- stats::pnorm(g_hat)
  
  if (any(!is.finite(g_hat)) || any(!is.finite(p_hat))) {
    stop("The SPP estimator produced nonfinite predictions.")
  }
  
  list(
    g_hat = g_hat,
    p_hat = p_hat
  )
}