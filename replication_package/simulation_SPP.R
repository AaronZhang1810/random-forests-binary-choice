##############################################################################

# Random Forests for Nonsmooth Binary Choice Models
# Simulation Study: Second-Order Polynomial Probit (SPP)
# Author: Jihuan Zhang, Department of Economics, Indiana University Bloomington
# Email: zhangjih@iu.edu

##############################################################################



##############################################################################
######### Housekeeping #########
##############################################################################
# Set this script's folder as the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls()); cat("\f")

# Load packages
library("openxlsx")

# Load functions
source("./functions/generate_data.R")
source("./functions/spp_estimation.R")
source("./functions/resume_simulation_results.R")



##############################################################################
######### Simulation settings #########
##############################################################################
# Specify output folder and names
xlsx_name  <- "./tables/simulation_results_SPP.xlsx"
Rdata_name <- "./tables/simulation_results_SPP.Rdata"

seed_grid <- seq(1810, 1809 + 1000, 1)
n_grid <- c(1000L, 2000L, 5000L)
DGP_grid <- c("sin","sigmoid_1","sigmoid_100","threshold","linear")
n_test <- 10000L
full_sample_size <- 10000L + n_test
glm_max_iterations <- 100L
glm_tolerance <- 1e-8
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

# Fitted-model information
simulation_results$n_parameters <- NA_integer_

# Test-sample approximation measures
simulation_results$g_mae <- NA_real_
simulation_results$g_mse <- NA_real_
simulation_results$p_mae <- NA_real_
simulation_results$p_mse <- NA_real_

# Classification performance
simulation_results$class_error <- NA_real_
simulation_results$class_error_true <- NA_real_

# Progress and convergence information
simulation_results$completed <- FALSE
simulation_results$converged <- NA
simulation_results$warning_message <- NA_character_

simulation_results <- resume_simulation_results(
  simulation_results = simulation_results,
  Rdata_name = Rdata_name,
  required_result_columns = setdiff(
    names(simulation_results),
    c("DGP", "n", "seed", "completed", "warning_message")
  )
)
simulation_results$completed <- simulation_results$completed &
  (simulation_results$converged %in% TRUE)



##############################################################################
######### Main Monte Carlo loop #########
##############################################################################

pending_rows <- which(!simulation_results$completed)
number_pending <- length(pending_rows)

if (number_pending == 0L) {
  message("All SPP Monte Carlo replications have already been completed.")
}

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

  replication_result <- tryCatch(
    local({
      full_sample <- generate_sample(
        DGP = DGP_current,
        n = full_sample_size,
        seed = seed_current
      )

      X_full <- as.matrix(full_sample$data[, -1, drop = FALSE])
      Y_full <- as.numeric(full_sample$data$Y)
      g_full <- as.numeric(full_sample$g)

      training_index <- seq_len(n_current)
      test_index <- 10001L:(10000L + n_test)

      X_train <- X_full[training_index, , drop = FALSE]
      Y_train <- Y_full[training_index]

      X_test <- X_full[test_index, , drop = FALSE]
      Y_test <- Y_full[test_index]
      g_true <- g_full[test_index]
      p_true <- stats::pnorm(g_true)
      
      # Fit SPP
      spp_model <- fit_spp_mle(
        X = X_train,
        Y = Y_train,
        max_iterations = glm_max_iterations,
        tolerance = glm_tolerance
      )
      
      # Prediction on the test set
      prediction <- predict_spp(
        object = spp_model,
        newdata = X_test
      )

      g_hat <- prediction$g_hat
      p_hat <- prediction$p_hat

      list(
        successful = TRUE,
        n_parameters = spp_model$n_parameters,
        converged = spp_model$converged,
        warning_message = spp_model$warning_message,
        g_mae = mean(abs(g_hat - g_true)),
        g_mse = mean((g_hat - g_true)^2),
        p_mae = mean(abs(p_hat - p_true)),
        p_mse = mean((p_hat - p_true)^2),
        class_error = mean((p_hat > 0.5) != Y_test),
        class_error_true = mean((p_true > 0.5) != Y_test)
      )
    }),
    error = function(e) {
      list(
        successful = FALSE,
        warning_message = paste0("Error: ", conditionMessage(e))
      )
    }
  )

  simulation_results$warning_message[row_index] <-
    replication_result$warning_message

  if (replication_result$successful) {
    simulation_results$n_parameters[row_index] <-
      replication_result$n_parameters
    simulation_results$converged[row_index] <-
      replication_result$converged

    simulation_results$g_mae[row_index] <- replication_result$g_mae
    simulation_results$g_mse[row_index] <- replication_result$g_mse
    simulation_results$p_mae[row_index] <- replication_result$p_mae
    simulation_results$p_mse[row_index] <- replication_result$p_mse
    simulation_results$class_error[row_index] <-
      replication_result$class_error
    simulation_results$class_error_true[row_index] <-
      replication_result$class_error_true

    simulation_results$completed[row_index] <- TRUE
  } else {
    simulation_results$completed[row_index] <- FALSE
    simulation_results$converged[row_index] <- NA

    warning(
      "Replication failed for DGP = ",
      DGP_current,
      ", n = ",
      n_current,
      ", seed = ",
      seed_current,
      ": ",
      replication_result$warning_message,
      call. = FALSE
    )
  }

  save(simulation_results, file = Rdata_name)

  if (
    replication_index %% excel_save_frequency == 0L ||
    replication_index == number_pending
  ) {
    openxlsx::write.xlsx(
      simulation_results,
      file = xlsx_name,
      overwrite = TRUE
    )
  }

  rm(replication_result)
  invisible(gc())
}



##############################################################################
######### Final save #########
##############################################################################

save(simulation_results, file = Rdata_name)

openxlsx::write.xlsx(
  simulation_results,
  file = xlsx_name,
  overwrite = TRUE
)

number_incomplete <- sum(!simulation_results$completed)

if (number_incomplete == 0L) {
  message("The SPP Monte Carlo simulation is complete.")
} else {
  warning(
    number_incomplete,
    " SPP replications remain incomplete. Restart the script to retry them.",
    call. = FALSE
  )
}
