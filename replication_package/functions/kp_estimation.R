##############################################################################

# Estimation for Binary Choice Models with Random Forests
# Author: Jihuan Zhang, Department of Economics, Indiana University Bloomington
# Email: zhangjih@iu.edu

##############################################################################



##############################################################################
######### Main functions for kernelized probit #########
##############################################################################
# Computing gram matrix using Gaussian kernel
gaussian_kernel <- function(X1, X2 = NULL, sigma = 1) {
  
  X1 <- as.matrix(X1)
  
  if (is.null(X2)) {
    X2 <- X1
  } else {
    X2 <- as.matrix(X2)
  }
  
  squared_norm_1 <- rowSums(X1^2)
  squared_norm_2 <- rowSums(X2^2)
  
  squared_distance <- outer(
    squared_norm_1,
    squared_norm_2,
    "+"
  ) - 2 * tcrossprod(X1, X2)
  
  # Remove small negative values caused by numerical rounding
  squared_distance <- pmax(squared_distance, 0)
  
  exp(-squared_distance / (2 * sigma^2))
}

# Computing leading eigenvalues for the gram matrix
leading_eigen <- function(K,
                          m_max,
                          eigenvalue_tolerance = 1e-10) {
  
  n <- nrow(K)
  m_max <- min(m_max, n)
  
  use_rspectra <- requireNamespace("RSpectra", quietly = TRUE) &&
    m_max < n - 1
  
  if (use_rspectra) {
    
    eig <- RSpectra::eigs_sym(
      K,
      k = m_max,
      which = "LA"
    )
    
    ordering <- order(eig$values, decreasing = TRUE)
    
    values <- eig$values[ordering]
    vectors <- eig$vectors[, ordering, drop = FALSE]
    
  } else {
    
    eig <- eigen(K, symmetric = TRUE)
    
    values <- eig$values[seq_len(m_max)]
    vectors <- eig$vectors[, seq_len(m_max), drop = FALSE]
  }
  
  # Remove eigenvalues that are numerically zero
  keep <- values >
    eigenvalue_tolerance * max(values)
  
  list(
    values = values[keep],
    vectors = vectors[, keep, drop = FALSE]
  )
}

# fit the kp model (intermediate; for cross-validation)
fit_kp_reduced <- function(Y, # binary outcome
                           eigenvectors, # n*m matrix, for spectral cutoff
                           eigenvalues, # m*1 vector, for spectral cutoff
                           B2, # the constraint
                           initial_value = NULL, # for gradient descent, default is zero
                           n_starts = 1,
                           seed = 1810,
                           max_iterations = 2000,
                           tolerance = 1e-8) {
  
  if (!requireNamespace("nloptr", quietly = TRUE)) {
    stop("The nloptr package is required.")
  }
  
  # create basic variables
  Y <- as.numeric(Y)
  U <- eigenvectors
  lambda <- eigenvalues
  n <- length(Y)
  m <- length(lambda)
  
  # set initial values
  if (is.null(initial_value)) {
    initial_value <- rep(0, m)
  }
  
  # Rescale the starting value if it does not satisfy the constraint
  # The function is used on the user-provided initial value
  # so that nloptr always receives a feasible starting point.
  project_to_constraint <- function(zeta) {
    
    quadratic_form <- sum(zeta^2 / lambda)
    
    if (!is.finite(quadratic_form)) {
      return(rep(0, m))
    }
    
    if (quadratic_form > B2) {
      zeta <- zeta * sqrt(B2 / quadratic_form)
    }
    
    zeta
  }
  
  initial_value <- project_to_constraint(initial_value)
  
  objective <- function(zeta) {
    
    latent_index <- drop(U %*% zeta)
    probability <- pnorm(latent_index)
    
    mean((Y - probability)^2)
  }
  
  objective_gradient <- function(zeta) {
    
    latent_index <- drop(U %*% zeta)
    probability <- pnorm(latent_index)
    
    drop(
      (2 / n) *
        crossprod(
          U,
          (probability - Y) * dnorm(latent_index)
        )
    )
  }
  
  inequality_constraint <- function(zeta) {
    
    sum(zeta^2 / lambda) - B2
  }
  
  inequality_jacobian <- function(zeta) {
    
    matrix(
      2 * zeta / lambda,
      nrow = 1
    )
  }
  
  starting_values <- vector("list", n_starts)
  starting_values[[1]] <- initial_value
  
  if (n_starts > 1) {
    
    set.seed(seed)
    
    for (s in 2:n_starts) {
      
      direction <- rnorm(m)
      direction <- direction / sqrt(sum(direction^2))
      
      # Generate a point strictly inside the feasible ellipsoid
      radius <- runif(1, min = 0.05, max = 0.75)
      
      starting_values[[s]] <-
        sqrt(B2) *
        sqrt(lambda) *
        radius *
        direction
    }
  }
  
  best_result <- NULL
  
  for (s in seq_len(n_starts)) {
    
    optimization <- tryCatch(
      nloptr::nloptr(
        x0 = starting_values[[s]],
        eval_f = objective,
        eval_grad_f = objective_gradient,
        eval_g_ineq = inequality_constraint,
        eval_jac_g_ineq = inequality_jacobian,
        opts = list(
          algorithm = "NLOPT_LD_SLSQP",
          xtol_rel = tolerance,
          maxeval = max_iterations,
          print_level = 0
        )
      ),
      error = function(e) NULL
    )
    
    if (is.null(optimization)) {
      next
    }
    
    zeta_hat <- project_to_constraint(optimization$solution)
    objective_value <- objective(zeta_hat)
    
    if (
      is.finite(objective_value) &&
      (
        is.null(best_result) ||
        objective_value < best_result$objective
      )
    ) {
      
      best_result <- list(
        zeta = zeta_hat,
        objective = objective_value,
        status = optimization$status,
        message = optimization$message,
        iterations = optimization$iterations
      )
    }
  }
  
  if (is.null(best_result)) {
    stop("The KP optimization failed for all starting values.")
  }
  
  best_result
}

# fit the kp model (final; after cross-validation)
fit_kp <- function(X,
                   Y,
                   B2,
                   m,
                   sigma = 1,
                   n_starts = 1,
                   seed = 1810,
                   max_iterations = 2000,
                   tolerance = 1e-8,
                   eigenvalue_tolerance = 1e-10) {
  
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  
  K <- gaussian_kernel(
    X1 = X,
    sigma = sigma
  )
  
  eig <- leading_eigen(
    K = K,
    m_max = m,
    eigenvalue_tolerance = eigenvalue_tolerance
  )
  
  if (length(eig$values) < m) {
    stop(
      "There are fewer than m numerically positive eigenvalues."
    )
  }
  
  eigenvalues <- eig$values[seq_len(m)]
  eigenvectors <- eig$vectors[, seq_len(m), drop = FALSE]
  
  reduced_fit <- fit_kp_reduced(
    Y = Y,
    eigenvectors = eigenvectors,
    eigenvalues = eigenvalues,
    B2 = B2,
    n_starts = n_starts,
    seed = seed,
    max_iterations = max_iterations,
    tolerance = tolerance
  )
  
  zeta_hat <- reduced_fit$zeta
  
  # delta_hat = U_m Lambda_m^{-1} zeta_hat
  delta_hat <- drop(
    eigenvectors %*% (zeta_hat / eigenvalues)
  )
  
  fitted_g <- drop(K %*% delta_hat)
  fitted_probability <- pnorm(fitted_g)
  
  model <- list(
    X_train = X,
    Y_train = Y,
    delta = delta_hat,
    zeta = zeta_hat,
    eigenvalues = eigenvalues,
    B2 = B2,
    B = sqrt(B2),
    m = m,
    sigma = sigma,
    fitted_g = fitted_g,
    fitted_probability = fitted_probability,
    objective = reduced_fit$objective,
    optimizer_status = reduced_fit$status,
    optimizer_message = reduced_fit$message
  )
  
  class(model) <- "kp_fit"
  
  model
}
# Prediction
predict.kp_fit <- function(object,
                           newdata,
                           type = c("probability", "g"),
                           ...) {
  
  type <- match.arg(type)
  
  X_new <- as.matrix(newdata)
  
  K_new <- gaussian_kernel(
    X1 = X_new,
    X2 = object$X_train,
    sigma = object$sigma
  )
  
  g_hat <- drop(K_new %*% object$delta)
  
  if (type == "g") {
    return(g_hat)
  }
  
  pnorm(g_hat)
}

# Make subsamples for cross-validation
make_stratified_folds <- function(Y,
                                  k_folds = 5,
                                  seed = 1810) {
  
  Y <- as.numeric(Y)
  fold_id <- integer(length(Y))
  
  set.seed(seed)
  
  for (class_value in sort(unique(Y))) {
    
    indices <- which(Y == class_value)
    indices <- sample(indices)
    
    fold_id[indices] <-
      rep(seq_len(k_folds), length.out = length(indices))
  }
  
  fold_id
}

# Cross-validation for selecting B^2 and m
cv_kp <- function(X,
                  Y,
                  B2_grid = c(
                    5, 10, 20, 30, 40,
                    50, 75, 100, 200, 500
                  ),
                  m_grid = c(10, 25, 50, 75, 100),
                  k_folds = 5,
                  sigma = 1,
                  seed = 1810,
                  n_starts = 1,
                  max_iterations = 2000,
                  tolerance = 1e-8,
                  eigenvalue_tolerance = 1e-10,
                  verbose = TRUE) {
  
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  
  if (!all(Y %in% c(0, 1))) {
    stop("Y must contain only zeros and ones.")
  }
  
  B2_grid <- sort(unique(B2_grid))
  m_grid <- sort(unique(as.integer(m_grid)))
  
  fold_id <- make_stratified_folds(
    Y = Y,
    k_folds = k_folds,
    seed = seed
  )
  
  # Since sigma is fixed, the full Gram matrix can be computed once.
  # Submatrices are then used within each fold.
  K_full <- gaussian_kernel(
    X1 = X,
    sigma = sigma
  )
  
  fold_results <- list()
  result_counter <- 1
  
  for (fold in seq_len(k_folds)) {
    
    if (verbose) {
      message(
        "Cross-validation fold ",
        fold,
        " of ",
        k_folds
      )
    }
    
    validation_index <- which(fold_id == fold)
    training_index <- which(fold_id != fold)
    
    Y_training <- Y[training_index]
    Y_validation <- Y[validation_index]
    
    K_training <- K_full[
      training_index,
      training_index,
      drop = FALSE
    ]
    
    K_validation <- K_full[
      validation_index,
      training_index,
      drop = FALSE
    ]
    
    eig <- leading_eigen(
      K = K_training,
      m_max = max(m_grid),
      eigenvalue_tolerance = eigenvalue_tolerance
    )
    
    available_m <- length(eig$values)
    valid_m_grid <- m_grid[m_grid <= available_m]
    
    if (length(valid_m_grid) == 0) {
      stop(
        "No candidate value of m is smaller than the numerical rank."
      )
    }
    
    for (m in valid_m_grid) {
      
      eigenvalues <- eig$values[seq_len(m)]
      eigenvectors <- eig$vectors[
        ,
        seq_len(m),
        drop = FALSE
      ]
      
      # Warm start when moving from a smaller B^2 to a larger B^2
      initial_value <- rep(0, m)
      
      for (B2 in B2_grid) {
        
        reduced_fit <- fit_kp_reduced(
          Y = Y_training,
          eigenvectors = eigenvectors,
          eigenvalues = eigenvalues,
          B2 = B2,
          initial_value = initial_value,
          n_starts = n_starts,
          seed = seed + 10000 * fold + 100 * m,
          max_iterations = max_iterations,
          tolerance = tolerance
        )
        
        initial_value <- reduced_fit$zeta
        
        delta_hat <- drop(
          eigenvectors %*%
            (reduced_fit$zeta / eigenvalues)
        )
        
        g_validation <- drop(
          K_validation %*% delta_hat
        )
        
        p_validation <- pnorm(g_validation)
        
        validation_loss <- mean(
          (Y_validation - p_validation)^2
        )
        
        fold_results[[result_counter]] <- data.frame(
          fold = fold,
          B2 = B2,
          B = sqrt(B2),
          m = m,
          validation_loss = validation_loss,
          optimizer_status = reduced_fit$status
        )
        
        result_counter <- result_counter + 1
      }
    }
  }
  
  fold_results <- do.call(
    rbind,
    fold_results
  )
  
  grouped_results <- split(
    fold_results,
    interaction(
      fold_results$B2,
      fold_results$m,
      drop = TRUE
    )
  )
  
  cv_results <- do.call(
    rbind,
    lapply(
      grouped_results,
      function(result) {
        
        data.frame(
          B2 = result$B2[1],
          B = result$B[1],
          m = result$m[1],
          mean_validation_loss =
            mean(result$validation_loss),
          sd_validation_loss =
            sd(result$validation_loss),
          number_of_folds =
            length(unique(result$fold))
        )
      }
    )
  )
  
  rownames(cv_results) <- NULL
  
  # Retain only combinations evaluated in every fold
  cv_results <- cv_results[
    cv_results$number_of_folds == k_folds,
    ,
    drop = FALSE
  ]
  
  cv_results <- cv_results[
    order(
      cv_results$mean_validation_loss,
      cv_results$m,
      cv_results$B2
    ),
    ,
    drop = FALSE
  ]
  
  best_B2 <- cv_results$B2[1]
  best_m <- cv_results$m[1]
  
  if (verbose) {
    message(
      "Selected B^2 = ",
      best_B2,
      ", B = ",
      round(sqrt(best_B2), 4),
      ", and m = ",
      best_m
    )
  }
  
  final_model <- fit_kp(
    X = X,
    Y = Y,
    B2 = best_B2,
    m = best_m,
    sigma = sigma,
    n_starts = n_starts,
    seed = seed,
    max_iterations = max_iterations,
    tolerance = tolerance,
    eigenvalue_tolerance = eigenvalue_tolerance
  )
  
  output <- list(
    model = final_model,
    best_B2 = best_B2,
    best_B = sqrt(best_B2),
    best_m = best_m,
    sigma = sigma,
    cv_results = cv_results,
    fold_results = fold_results,
    fold_id = fold_id
  )
  
  class(output) <- "kp_cv"
  
  output
}

# Prediction with cross-validation
predict.kp_cv <- function(object,
                          newdata,
                          type = c("probability", "g"),
                          ...) {
  
  type <- match.arg(type)
  
  predict(
    object$model,
    newdata = newdata,
    type = type
  )
}

# Batched KP prediction (to save memory during prediction stage)
predict_kp_g_batch <- function(object,
                               newdata,
                               batch_size = 500L) {
  
  X_new <- as.matrix(newdata)
  X_train <- as.matrix(object$X_train)
  
  if (ncol(X_new) != ncol(X_train)) {
    stop(
      "The training sample has ",
      ncol(X_train),
      " covariates, whereas the test sample has ",
      ncol(X_new),
      " covariates."
    )
  }
  
  n_new <- nrow(X_new)
  g_hat <- numeric(n_new)
  
  batch_start <- seq.int(
    from = 1L,
    to = n_new,
    by = batch_size
  )
  
  for (start_index in batch_start) {
    
    end_index <- min(
      start_index + batch_size - 1L,
      n_new
    )
    
    batch_index <- start_index:end_index
    
    K_batch <- gaussian_kernel(
      X1 = X_new[batch_index, , drop = FALSE],
      X2 = X_train,
      sigma = object$sigma
    )
    
    g_hat[batch_index] <- drop(
      K_batch %*% object$delta
    )
    
    rm(K_batch)
  }
  
  g_hat
}