##############################################################################

# Random Forests for Nonsmooth Binary Choice Models
# Simulation Study: Kernelized Probit -- Parallel Version
# Author: Jihuan Zhang, Department of Economics, Indiana University Bloomington
# Email: zhangjih@iu.edu

##############################################################################



##############################################################################
######### Housekeeping #########
##############################################################################
# Set this script's folder as the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls()); cat("\f")

# load packages
library("openxlsx")
library("RSpectra") # for large scale spectral decomposition
library("nloptr")   # for nonlinear optimization
library("parallel")

# load functions
source("./functions/generate_data.R")
source("./functions/kp_estimation.R")
source("./functions/resume_simulation_results.R")



##############################################################################
######### Simulation settings #########
##############################################################################
# Specify output folder and names
xlsx_name  <- "./tables/simulation_results_KP.xlsx"
Rdata_name <- "./tables/simulation_results_KP.Rdata"

seed_grid <- seq(1810, 1809 + 1000, 1) # for MC replications
n_grid <- c(1000, 2000, 5000) # sample size
DGP_grid <- c("sin", "sigmoid_1", "sigmoid_100", "threshold", "linear")

# Yan (2026) p.16 footnote 22
B2_grid <- c(5, 10, 20, 30, 40, 50, 75, 100, 200, 500) # for constraint
m_grid <- c(10, 25, 50, 75, 100, 150, 200) # for spectral cutoff

k_folds <- 5L # Number of CV folds
sigma <- 1 # Gaussian-kernel bandwidth
n_starts <- 1L # Number of optimizer starting values
prediction_batch_size <- 500L # Number of test observations in each prediction batch
excel_save_frequency <- 32L # attempted replications between Excel writes

# Parallel settings
requested_workers <- 8L # Number of cores
detected_workers <- parallel::detectCores(logical = TRUE)
if (is.na(detected_workers)) {
  detected_workers <- requested_workers + 1L
}
available_workers <- max(1L, detected_workers - 1L)
number_workers <- min(requested_workers, available_workers)

# Tasks are submitted in batches so that the master process can checkpoint
# completed results regularly. Only the master process writes output files.
parallel_batch_size <- 16L



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

simulation_results$m_cv <- NA_integer_
simulation_results$B2_cv <- NA_real_
simulation_results$cv_error <- NA_real_

# Test-sample approximation measures
simulation_results$g_mae <- NA_real_
simulation_results$g_mse <- NA_real_
simulation_results$p_mae <- NA_real_
simulation_results$p_mse <- NA_real_

# Numerical and computational information
simulation_results$optimizer_status <- NA_integer_
simulation_results$runtime_seconds <- NA_real_

# Progress information
simulation_results$completed <- FALSE
simulation_results$error_message <- NA_character_

# Resume previous simulation results
simulation_results <- resume_simulation_results(
  simulation_results = simulation_results,
  Rdata_name = Rdata_name
)



##############################################################################
######### One Monte Carlo replication #########
##############################################################################
run_kp_replication <- function(task) {
  row_index <- task$row_index
  DGP_current <- task$DGP
  n_current <- task$n
  seed_current <- task$seed

  start_time <- proc.time()[["elapsed"]]

  replication_result <- tryCatch(
    {
      # This task-specific seed makes results reproducible independently of
      # worker assignment and load-balancing order.
      set.seed(seed_current)

      # Generate training and test samples
      full_sample <- generate_sample(
        DGP = DGP_current,
        n = 20000,
        seed = seed_current
      )

      X_full <- as.matrix(full_sample$data[, -1, drop = FALSE])
      Y_full <- full_sample$data$Y
      g_full <- full_sample$g

      training_index <- seq_len(n_current)
      test_index <- 10001:20000

      X_train <- X_full[training_index, , drop = FALSE]
      Y_train <- Y_full[training_index]

      X_test <- X_full[test_index, , drop = FALSE]
      g_true <- g_full[test_index]
      p_true <- pnorm(g_true)

      # Cross-validation and full-sample estimation
      kp_model <- cv_kp(
        X = X_train,
        Y = Y_train,
        B2_grid = B2_grid,
        m_grid = m_grid,
        k_folds = k_folds,
        sigma = sigma,
        seed = seed_current + 10000L,
        n_starts = n_starts,
        max_iterations = 2000,
        tolerance = 1e-8,
        eigenvalue_tolerance = 1e-10,
        verbose = FALSE
      )

      # Reject unsuccessful optimization statuses
      if (
        any(kp_model$fold_results$optimizer_status < 0) ||
        kp_model$model$optimizer_status < 0
      ) {
        stop("At least one KP optimization returned a negative status.")
      }

      # KP predictions
      g_hat <- predict_kp_g_batch(
        object = kp_model$model,
        newdata = X_test,
        batch_size = prediction_batch_size
      )

      p_hat <- pnorm(g_hat)

      if (
        any(!is.finite(g_hat)) ||
        any(!is.finite(p_hat))
      ) {
        stop("The KP estimator produced nonfinite test predictions.")
      }

      list(
        successful = TRUE,
        row_index = row_index,
        DGP = DGP_current,
        n = n_current,
        seed = seed_current,
        m_cv = as.integer(kp_model$best_m),
        B2_cv = kp_model$best_B2,
        cv_error = kp_model$cv_results$mean_validation_loss[1],
        g_mae = mean(abs(g_hat - g_true)),
        g_mse = mean((g_hat - g_true)^2),
        p_mae = mean(abs(p_hat - p_true)),
        p_mse = mean((p_hat - p_true)^2),
        optimizer_status = as.integer(kp_model$model$optimizer_status),
        error_message = NA_character_
      )
    },

    error = function(e) {
      list(
        successful = FALSE,
        row_index = row_index,
        DGP = DGP_current,
        n = n_current,
        seed = seed_current,
        error_message = conditionMessage(e)
      )
    }
  )

  replication_result$runtime_seconds <-
    proc.time()[["elapsed"]] - start_time

  replication_result
}



##############################################################################
######### Master-process update function #########
##############################################################################
update_simulation_results <- function(simulation_results, replication_result) {
  row_index <- replication_result$row_index

  simulation_results$runtime_seconds[row_index] <-
    replication_result$runtime_seconds

  if (replication_result$successful) {
    simulation_results$m_cv[row_index] <-
      replication_result$m_cv

    simulation_results$B2_cv[row_index] <-
      replication_result$B2_cv

    simulation_results$cv_error[row_index] <-
      replication_result$cv_error

    simulation_results$g_mae[row_index] <-
      replication_result$g_mae

    simulation_results$g_mse[row_index] <-
      replication_result$g_mse

    simulation_results$p_mae[row_index] <-
      replication_result$p_mae

    simulation_results$p_mse[row_index] <-
      replication_result$p_mse

    simulation_results$optimizer_status[row_index] <-
      replication_result$optimizer_status

    simulation_results$completed[row_index] <- TRUE
    simulation_results$error_message[row_index] <- NA_character_
  } else {
    simulation_results$completed[row_index] <- FALSE
    simulation_results$error_message[row_index] <-
      replication_result$error_message
  }

  simulation_results
}



##############################################################################
######### Parallel Monte Carlo loop #########
##############################################################################
pending_rows <- which(!simulation_results$completed)
number_pending <- length(pending_rows)
project_directory <- getwd()

if (number_pending == 0L) {
  message("All Monte Carlo replications have already been completed.")
} else {
  message(
    "Starting ",
    number_pending,
    " pending replications with ",
    number_workers,
    " PSOCK workers."
  )

  # A PSOCK cluster is used rather than a fork cluster for portability and to
  # isolate the numerical optimizer in each worker process.
  cluster <- parallel::makePSOCKcluster(number_workers)

  tryCatch(
    {
      # Export the project path before initializing the workers.
      parallel::clusterExport(
        cluster,
        varlist = "project_directory",
        envir = globalenv()
      )

      # Prevent each worker from using additional BLAS/OpenMP threads. This
      # avoids severe oversubscription when several spectral decompositions and
      # optimizations run simultaneously.
      parallel::clusterEvalQ(
        cluster,
        {
          Sys.setenv(
            OMP_NUM_THREADS = "1",
            OPENBLAS_NUM_THREADS = "1",
            MKL_NUM_THREADS = "1",
            VECLIB_MAXIMUM_THREADS = "1",
            NUMEXPR_NUM_THREADS = "1"
          )

          library("RSpectra")
          library("nloptr")

          setwd(project_directory)
          source("./functions/generate_data.R")
          source("./functions/kp_estimation.R")

          NULL
        }
      )

      parallel::clusterExport(
        cluster,
        varlist = c(
          "B2_grid",
          "m_grid",
          "k_folds",
          "sigma",
          "n_starts",
          "prediction_batch_size",
          "run_kp_replication"
        ),
        envir = globalenv()
      )

      task_table <- data.frame(
        row_index = pending_rows,
        DGP = simulation_results$DGP[pending_rows],
        n = simulation_results$n[pending_rows],
        seed = simulation_results$seed[pending_rows],
        stringsAsFactors = FALSE
      )

      task_batches <- split(
        seq_len(nrow(task_table)),
        ceiling(seq_len(nrow(task_table)) / parallel_batch_size)
      )

      number_batches <- length(task_batches)
      number_attempted <- 0L
      last_excel_write <- 0L

      for (batch_index in seq_along(task_batches)) {
        batch_positions <- task_batches[[batch_index]]

        batch_tasks <- lapply(
          batch_positions,
          function(position) {
            list(
              row_index = task_table$row_index[position],
              DGP = task_table$DGP[position],
              n = task_table$n[position],
              seed = task_table$seed[position]
            )
          }
        )

        message(
          "Submitting batch ",
          batch_index,
          " of ",
          number_batches,
          " (",
          length(batch_tasks),
          " replications)."
        )

        # Load balancing is useful because optimization time can vary greatly
        # across DGPs, sample sizes, and tuning parameters.
        batch_results <- parallel::parLapplyLB(
          cluster,
          batch_tasks,
          run_kp_replication
        )

        for (replication_result in batch_results) {
          simulation_results <- update_simulation_results(
            simulation_results = simulation_results,
            replication_result = replication_result
          )

          if (replication_result$successful) {
            message(
              "Completed: DGP = ",
              replication_result$DGP,
              ", n = ",
              replication_result$n,
              ", seed = ",
              replication_result$seed,
              ", runtime = ",
              round(replication_result$runtime_seconds, 2),
              " seconds."
            )
          } else {
            warning(
              "Replication failed: DGP = ",
              replication_result$DGP,
              ", n = ",
              replication_result$n,
              ", seed = ",
              replication_result$seed,
              ". ",
              replication_result$error_message,
              call. = FALSE
            )
          }
        }

        number_attempted <- number_attempted + length(batch_results)

        # Checkpoint after every batch. Failed replications remain incomplete and
        # will be retried automatically when the script is restarted.
        save(simulation_results, file = Rdata_name)

        if (
          number_attempted - last_excel_write >= excel_save_frequency ||
          batch_index == number_batches
        ) {
          openxlsx::write.xlsx(
            simulation_results,
            file = xlsx_name,
            overwrite = TRUE
          )

          last_excel_write <- number_attempted
        }

        rm(batch_results, batch_tasks)
        invisible(gc())
      }
    },

    finally = {
      parallel::stopCluster(cluster)
    }
  )
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
  message("The parallel KP Monte Carlo simulation is complete.")
} else {
  warning(
    number_incomplete,
    " replications remain incomplete. Restart the script to retry them.",
    call. = FALSE
  )
}
