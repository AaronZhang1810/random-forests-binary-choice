##############################################################################

# Random Forests for Nonsmooth Binary Choice Models
# Simulation Study: Linear Probit
# Author: Jihuan Zhang, Department of Economics, Indiana University Bloomington
# Email: zhangjih@iu.edu

##############################################################################



##############################################################################
######### Housekeeping #########
##############################################################################
# Set this script's folder as the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls()); cat("\f")

# Load packages and functions
library("openxlsx")

source("./functions/generate_data.R")
source("./functions/resume_simulation_results.R")



##############################################################################
######### Simulation settings #########
##############################################################################
# Specify output folder and names
xlsx_name  <- "./tables/simulation_results_LinearProbit.xlsx"
Rdata_name <- "./tables/simulation_results_LinearProbit.Rdata"

seed_grid <- seq(1810, 1809 + 1000, 1)
n_grid <- c(1000L, 2000L, 5000L)
DGP_grid <- c("sin","sigmoid_1","sigmoid_100","threshold","linear")
n_test <- 10000L
full_sample_size <- 10000L + n_test
excel_save_frequency <- 100L



##############################################################################
######### Initialize result table #########
##############################################################################
simulation_results <- expand.grid(
  DGP = DGP_grid,
  n = n_grid,
  seed = seed_grid,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# Test-sample approximation measures
simulation_results$g_mae <- NA_real_
simulation_results$g_mse <- NA_real_
simulation_results$p_mae <- NA_real_
simulation_results$p_mse <- NA_real_

# Classification performance
simulation_results$class_error <- NA_real_
simulation_results$class_error_true <- NA_real_

# Progress information
simulation_results$completed <- FALSE

simulation_results <- resume_simulation_results(
  simulation_results = simulation_results,
  Rdata_name = Rdata_name,
  required_result_columns = setdiff(
    names(simulation_results),
    c("DGP", "n", "seed", "completed")
  )
)



##############################################################################
######### Main Loop #########
##############################################################################
pending_rows <- which(!simulation_results$completed)
number_pending <- length(pending_rows)

for (replication_index in seq_along(pending_rows)) {
  row_index <- pending_rows[replication_index]
  
  DGP_current <- simulation_results$DGP[row_index]
  n_current <- as.integer(simulation_results$n[row_index])
  seed_current <- as.integer(simulation_results$seed[row_index])

  message(
    "Replication ",
    replication_index,
    " of ",
    number_pending,
    ": DGP = ",
    DGP_current,
    ", n = ",
    n_current,
    ", seed = ",
    seed_current
  )
  
  # Generate full sample
  full_sample <- generate_sample(DGP = DGP_current, n = full_sample_size, seed = seed_current)
  X_full <- as.matrix(full_sample$data[, -1, drop = FALSE])
  Y_full <- as.numeric(full_sample$data$Y)
  g_full <- as.numeric(full_sample$g)
  
  # Get training and test sets
  training_index <- seq_len(n_current)
  test_index <- 10001L:(10000L + n_test)
  
  X_train <- X_full[training_index, , drop = FALSE]
  Y_train <- Y_full[training_index]
  
  X_test <- X_full[test_index, , drop = FALSE]
  Y_test <- Y_full[test_index]
  
  # Get true g and p
  g_true <- g_full[test_index]
  p_true <- pnorm(g_true)
  
  # Fit linear probit
  probit_model <- glm(
    Y ~ .,
    data = data.frame(Y = Y_train, X_train),
    family = binomial(link = "probit"),
    singular.ok = FALSE
  )
  
  # Prediction
  g_hat <- predict(
    probit_model,
    newdata = as.data.frame(X_test),
    type = "link"
  )
  p_hat <- pnorm(g_hat)
  
  replication_result <- list(
    g_mae = mean(abs(g_hat - g_true)),
    g_mse = mean((g_hat - g_true)^2),
    p_mae = mean(abs(p_hat - p_true)),
    p_mse = mean((p_hat - p_true)^2),
    class_error = mean((p_hat > 0.5) != Y_test),
    class_error_true = mean((p_true > 0.5) != Y_test),
    converged = isTRUE(probit_model$converged)
  )
  
  # Write results
  simulation_results$g_mae[row_index] <- replication_result$g_mae
  simulation_results$g_mse[row_index] <- replication_result$g_mse
  simulation_results$p_mae[row_index] <- replication_result$p_mae
  simulation_results$p_mse[row_index] <- replication_result$p_mse
  simulation_results$class_error[row_index] <- replication_result$class_error
  simulation_results$class_error_true[row_index] <- replication_result$class_error_true
  simulation_results$completed[row_index] <- replication_result$converged
  
  save(simulation_results, file = Rdata_name)

  if (
    replication_index %% excel_save_frequency == 0L || replication_index == number_pending
  ) {
    openxlsx::write.xlsx(simulation_results, file = xlsx_name, overwrite = TRUE)
  }
}



##############################################################################
######### Final save #########
##############################################################################
save(simulation_results, file = Rdata_name)
openxlsx::write.xlsx(simulation_results, file = xlsx_name, overwrite = TRUE)

message("All linear probit replications are complete.")

