##############################################################################

# Random Forests for Nonsmooth Binary Choice Models
# Simulation Study
# Author: Jihuan Zhang, Department of Economics, Indiana University Bloomington
# Email: zhangjih@iu.edu

##############################################################################



##############################################################################
######### Functions for DGPs #########
##############################################################################
generate_sin <- function(n,
                         beta = c(-0.47, -0.26, 0.15, 0.82, -0.60,
                                  0.80, 0.89, 0.32, 0.26, -0.88),
                         seed = 1810) {
  set.seed(seed)
  
  m <- 10
  X <- matrix(runif(n * m, -1, 1), nrow = n, ncol = m)
  colnames(X) <- paste0("X", 1:m)
  
  g <- rep(0, n)
  
  for (j in 1:5) {
    g <- g +
      beta[j] * sin(pi * X[, j] * X[, j + 5]) +
      beta[j + 5] * X[, j] * X[, j + 5]
  }
  
  U <- rnorm(n)
  Y <- as.integer(g - U > 0)
  
  data <- data.frame(Y = Y, X)
  
  return(list(
    data = data,
    g = g,
    U = U
  ))
}

generate_sigmoid <- function(n,
                             kappa = 1,
                             beta = c(-0.47, -0.26, 0.15, 0.82, -0.60,
                                      0.80, 0.89, 0.32, 0.26, -0.88),
                             seed = 1810) {
  set.seed(seed)
  
  m <- 10
  X <- matrix(runif(n * m, -1, 1), nrow = n, ncol = m)
  colnames(X) <- paste0("X", 1:m)
  
  alpha_kappa <- ( 0.25 * (1 - 2 * tanh(kappa / 2) / kappa) * sum(beta^2) )^(-0.5)
  
  g <- rep(0, n)
  
  for (j in 1:m) {
    g <- g +
      beta[j] * (plogis(kappa * X[, j]) - 0.5)
  }
  
  g <- alpha_kappa * g
  
  U <- rnorm(n)
  Y <- as.integer(g - U > 0)
  
  data <- data.frame(Y = Y, X)
  
  return(list(
    data = data,
    g = g,
    U = U
  ))
}

generate_threshold <- function(n,
                               beta = c(-0.47, -0.26, 0.15, 0.82, -0.60,
                                        0.80, 0.89, 0.32, 0.26, -0.88),
                               tau = rep(0, 10),
                               seed = 1810) {
  set.seed(seed)
  
  m <- 10
  X <- matrix(runif(n * m, -1, 1), nrow = n, ncol = m)
  colnames(X) <- paste0("X", 1:m)
  
  g <- rep(0, n)
  
  for (j in 1:10) {
    g <- g +
      beta[j] * (as.integer(X[, j] > tau[j]) - 0.5)
  }
  
  U <- rnorm(n)
  Y <- as.integer(g - U > 0)
  
  data <- data.frame(Y = Y, X)
  
  return(list(
    data = data,
    g = g,
    U = U
  ))
}

generate_linear <- function(n,
                            n_irrelevant = 10,
                            beta = c(-0.47, -0.26, 0.15, 0.82, -0.60,
                                     0.80, 0.89, 0.32, 0.26, -0.88),
                            seed = 1810) {
  set.seed(seed)
  
  m_relevant <- 10
  m <- m_relevant + n_irrelevant
  
  X <- matrix(runif(n * m, -1, 1), nrow = n, ncol = m)
  colnames(X) <- paste0("X", 1:m)
  
  g <- rep(0, n)
  
  for (j in 1:m_relevant) {
    g <- g + beta[j] * X[, j]
  }
  
  U <- rnorm(n)
  Y <- as.integer(g - U > 0)
  
  data <- data.frame(Y = Y, X)
  
  return(list(
    data = data,
    g = g,
    U = U
  ))
}

generate_binary_linear <- function(n,
                                   beta = c(-0.47, -0.26, 0.15, 0.82, -0.60,
                                            0.80, 0.89, 0.32, 0.26, -0.88),
                                   seed = 1810) {
  set.seed(seed)
  
  m <- 10
  
  X <- matrix(
    sample(
      c(-1 / sqrt(3), 1 / sqrt(3)),
      size = n * m,
      replace = TRUE
    ),
    nrow = n,
    ncol = m
  )
  
  colnames(X) <- paste0("X", 1:m)
  
  g <- as.numeric(X %*% beta)
  
  U <- rnorm(n)
  Y <- as.integer(g - U > 0)
  
  data <- data.frame(Y = Y, X)
  
  return(list(
    data = data,
    g = g,
    U = U
  ))
}

# Generate data for target DGP
generate_sample <- function(DGP, n, seed) {
  
  if (DGP == "sin") {
    
    return(
      generate_sin(
        n = n,
        seed = seed
      )
    )
  }
  
  if (DGP == "sigmoid_1") {
    
    return(
      generate_sigmoid(
        n = n,
        kappa = 1,
        seed = seed
      )
    )
  }
  
  if (DGP == "sigmoid_10") {
    
    return(
      generate_sigmoid(
        n = n,
        kappa = 10,
        seed = seed
      )
    )
  }
  
  if (DGP == "sigmoid_100") {
    
    return(
      generate_sigmoid(
        n = n,
        kappa = 100,
        seed = seed
      )
    )
  }
  
  if (DGP == "threshold") {
    
    return(
      generate_threshold(
        n = n,
        seed = seed
      )
    )
  }
  
  if (DGP == "linear") {
    
    return(
      generate_linear(
        n = n,
        n_irrelevant = 10,
        seed = seed
      )
    )
  }
  
  if (DGP == "binary") {
    
    return(
      generate_binary_linear(
        n = n,
        seed = seed
      )
    )
  }
  
  stop("Unknown DGP: ", DGP)
}
