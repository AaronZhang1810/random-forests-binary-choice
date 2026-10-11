##############################################################################

# Monte Carlo Simulation for Binary Choice Models
# Simulation Study: Standard Regression Forest Regression -- Parallel Version
# Author: Jihuan Zhang, Department of Economics, Indiana University Bloomington
# Email: zhangjih@iu.edu

##############################################################################



##############################################################################
######### Housekeeping #########
##############################################################################
# Set this script's folder as the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls()); cat("\f")

# load packages and functions
library("ranger")
library("openxlsx")
library("parallel")

# load functions
source("./functions/generate_data.R")
source("./functions/resume_simulation_results.R")



##############################################################################
######### Parameter Setting #########
##############################################################################
# Specify output folder and names
xlsx_name <- "./tables/simulation_results_RF.xlsx"
Rdata_name <- "./tables/simulation_results_RF.Rdata"

DGP_grid <- c("sin", "sigmoid_1", "sigmoid_100", "threshold", "linear")
n_test <- 10000L # test sample size
B <- 500L # number of trees to grow
seed_grid <- 1810:(1809 + 1000L) # for MC replications
k_folds <- 5L
n_grid <- c(1000L, 2000L, 5000L) # sample size
min_node_size_grid <- c(5L, 25L, 50L, 75L, 100L)
excel_save_frequency <- 10L # Write Excel after this many replications
num_threads <- 1L # Keep ranger single-threaded inside each parallel worker

# Parallelize independent Monte Carlo replications across R processes.
physical_cores <- parallel::detectCores(logical = FALSE)
if (is.na(physical_cores)) {
  physical_cores <- parallel::detectCores(logical = TRUE)
}
if (is.na(physical_cores)) {
  physical_cores <- 2L
}

n_workers <- max(1L, min(8L, as.integer(physical_cores) - 2L))

# Process two tasks per worker before checkpointing to improve load balancing.
parallel_batch_size <- max(n_workers, 2L * n_workers)



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

# Cross-validated forest hyperparameters.
simulation_results$mtry_cv <- NA_integer_
simulation_results$depth_cv <- NA_integer_
simulation_results$min_node_size_cv <- NA_integer_
simulation_results$cv_error <- NA_real_

# Boundary diagnostics for Phi^{-1}(p_hat_RF).
simulation_results$clip_epsilon <- NA_real_
simulation_results$n_clipped_test <- NA_integer_
simulation_results$share_clipped_test <- NA_real_
simulation_results$p_hat_min <- NA_real_
simulation_results$p_hat_max <- NA_real_

# Test-sample approximation measures.
simulation_results$g_mae <- NA_real_
simulation_results$g_mse <- NA_real_
simulation_results$p_mae <- NA_real_
simulation_results$p_mse <- NA_real_

# Predictive measures against the realized binary test outcome.
simulation_results$brier_score <- NA_real_
simulation_results$log_loss_clipped <- NA_real_
simulation_results$class_error <- NA_real_
simulation_results$class_error_true <- NA_real_

# Numerical and progress information.
simulation_results$runtime_seconds <- NA_real_
simulation_results$completed <- FALSE
simulation_results$error_message <- NA_character_

# Resume previous simulation results #########
simulation_results <- resume_simulation_results(
  simulation_results = simulation_results,
  Rdata_name = Rdata_name,
  required_result_columns = setdiff(
    names(simulation_results),
    c("DGP", "n", "seed", "completed", "error_message")
  )
)




##############################################################################
######### Helper functions #########
##############################################################################

validate_probability_vector <- function(probability,
                                        expected_length,
                                        label) {
  probability <- as.numeric(probability)

  if (length(probability) != expected_length) {
    stop(label, " has an incorrect length.")
  }
  if (any(!is.finite(probability))) {
    stop(label, " contains nonfinite values.")
  }

  numerical_tolerance <- sqrt(.Machine$double.eps)
  if (any(
    probability < -numerical_tolerance |
    probability > 1 + numerical_tolerance
  )) {
    stop(label, " contains values outside [0,1].")
  }

  # Correct only negligible floating-point excursions.
  pmin(1, pmax(0, probability))
}



##############################################################################
######### One Monte Carlo replication #########
##############################################################################

# Each worker receives only the identifying information for one replication.
# All large objects are created and destroyed inside the worker process.
run_rfplugin_replication <- function(task) {
  row_index <- as.integer(task$row_index)
  DGP_current <- as.character(task$DGP)
  n_current <- as.integer(task$n)
  seed_current <- as.integer(task$seed)

  replication_start_time <- proc.time()[["elapsed"]]

  replication_result <- tryCatch(
    local({
      ########################################################################
      # Generate the same common training and test samples as the RFP script
      ########################################################################

      full_sample <- generate_sample(
        DGP = DGP_current,
        n = 20000L,
        seed = seed_current
      )

      X_full <- full_sample$data[, -1, drop = FALSE]
      Y_full <- as.numeric(full_sample$data$Y)
      g_full <- as.numeric(full_sample$g)

      training_index <- seq_len(n_current)
      test_index <- 10001:20000

      X_train <- X_full[training_index, , drop = FALSE]
      Y_train <- Y_full[training_index]
      X_test <- X_full[test_index, , drop = FALSE]
      Y_test <- Y_full[test_index]

      g_true <- g_full[test_index]
      p_true <- stats::pnorm(g_true)

      d <- ncol(X_train)

      ########################################################################
      # Candidate hyperparameters: identical tree grid to the RFP simulation
      ########################################################################

      mtry_grid <- sort(unique(pmin(
        pmax(
          as.integer(round(c(sqrt(d), d / 3, d))),
          1L
        ),
        d
      )))

      # In ranger, max.depth = 0 denotes unlimited depth.
      depth_grid <- sort(unique(as.integer(c(
        0L,
        floor(0.5 * log2(n_current)),
        floor(log2(n_current)),
        floor(1.5 * log2(n_current)),
        floor(2 * log2(n_current))
      ))))

      forest_tuning_grid <- expand.grid(
        mtry = mtry_grid,
        depth = depth_grid,
        min.node.size = min_node_size_grid,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )

      number_forest_candidates <- nrow(forest_tuning_grid)

      ########################################################################
      # Construct the same stratified fold assignment rule as the RFP script
      ########################################################################

      fold_id <- integer(n_current)
      set.seed(seed_current + 10000L)

      for (class_value in c(0, 1)) {
        class_indices <- which(Y_train == class_value)

        if (length(class_indices) < k_folds) {
          stop("One outcome class contains fewer observations than k_folds.")
        }

        class_indices <- sample(class_indices)
        fold_id[class_indices] <- rep(
          seq_len(k_folds),
          length.out = length(class_indices)
        )
      }

      ########################################################################
      # Tune the ordinary forest using direct-probability Brier loss
      ########################################################################

      cv_sse <- matrix(
        NA_real_,
        nrow = number_forest_candidates,
        ncol = k_folds
      )

      for (fold in seq_len(k_folds)) {
        validation_index <- which(fold_id == fold)
        estimation_index <- which(fold_id != fold)

        X_estimation <- X_train[estimation_index, , drop = FALSE]
        Y_estimation <- Y_train[estimation_index]
        X_validation <- X_train[validation_index, , drop = FALSE]
        Y_validation <- Y_train[validation_index]

        for (candidate_index in seq_len(number_forest_candidates)) {
          mtry_candidate <- forest_tuning_grid$mtry[candidate_index]
          depth_candidate <- forest_tuning_grid$depth[candidate_index]
          min_node_size_candidate <-
            forest_tuning_grid$min.node.size[candidate_index]

          # Fit the standard random forest regression
          forest_cv <- ranger(
            y = Y_estimation,
            x = X_estimation,
            classification = FALSE,
            probability = FALSE,
            num.trees = B,
            mtry = mtry_candidate,
            splitrule = "variance",
            max.depth = depth_candidate,
            min.node.size = min_node_size_candidate,
            replace = TRUE,
            sample.fraction = 1,
            keep.inbag = FALSE,
            write.forest = TRUE,
            oob.error = FALSE,
            importance = "none",
            seed = seed_current + 100000L * fold,
            num.threads = num_threads,
            verbose = FALSE
          )

          p_validation <- predict(
            forest_cv,
            data = X_validation,
            type = "response",
            num.threads = num_threads,
            verbose = FALSE
          )$predictions

          p_validation <- validate_probability_vector(
            probability = p_validation,
            expected_length = length(validation_index),
            label = "A CV direct-forest prediction"
          )

          cv_sse[candidate_index, fold] <- sum(
            (Y_validation - p_validation)^2
          )

          rm(forest_cv, p_validation)
        }

        rm(
          X_estimation,
          Y_estimation,
          X_validation,
          Y_validation,
          validation_index,
          estimation_index
        )
        invisible(gc())
      }

      cv_error <- rowSums(cv_sse) / n_current

      if (any(!is.finite(cv_error))) {
        stop("At least one forest candidate has nonfinite CV error.")
      }

      # Break exact ties in favor of a shallower tree, a larger leaf-size
      # constraint, and a smaller mtry. Treat depth = 0 as unlimited depth.
      depth_complexity <- ifelse(
        forest_tuning_grid$depth == 0L,
        Inf,
        forest_tuning_grid$depth
      )

      best_index <- order(
        cv_error,
        depth_complexity,
        -forest_tuning_grid$min.node.size,
        forest_tuning_grid$mtry
      )[1L]

      best_mtry <- as.integer(forest_tuning_grid$mtry[best_index])
      best_depth <- as.integer(forest_tuning_grid$depth[best_index])
      best_min_node_size <- as.integer(
        forest_tuning_grid$min.node.size[best_index]
      )
      best_cv_error <- cv_error[best_index]

      ########################################################################
      # Refit the selected ordinary forest on the complete training sample
      ########################################################################

      forest_full <- ranger(
        y = Y_train,
        x = X_train,
        classification = FALSE,
        probability = FALSE,
        num.trees = B,
        mtry = best_mtry,
        splitrule = "variance",
        max.depth = best_depth,
        min.node.size = best_min_node_size,
        replace = TRUE,
        sample.fraction = 1,
        keep.inbag = FALSE,
        write.forest = TRUE,
        oob.error = FALSE,
        importance = "none",
        seed = seed_current + 900000L,
        num.threads = num_threads,
        verbose = FALSE
      )

      p_hat <- predict(
        forest_full,
        data = X_test,
        type = "response",
        num.threads = num_threads,
        verbose = FALSE
      )$predictions

      p_hat <- validate_probability_vector(
        probability = p_hat,
        expected_length = n_test,
        label = "The selected direct-forest prediction"
      )

      ########################################################################
      # Naive inverse-link estimate of the latent index
      ########################################################################

      # The raw p_hat is retained for every probability metric. Clipping is
      # used only because qnorm(0) and qnorm(1) are infinite. This prespecified
      # epsilon vanishes as the training sample grows.
      clip_epsilon <- 0.5 / (n_current + 1)
      p_for_inverse <- pmin(
        1 - clip_epsilon,
        pmax(clip_epsilon, p_hat)
      )
      clipped_test <- p_for_inverse != p_hat
      g_hat <- stats::qnorm(p_for_inverse)

      if (any(!is.finite(g_hat))) {
        stop("The inverse-probit plug-in estimate contains nonfinite values.")
      }

      ########################################################################
      # Test-sample performance
      ########################################################################

      list(
        successful = TRUE,
        best_mtry = best_mtry,
        best_depth = best_depth,
        best_min_node_size = best_min_node_size,
        best_cv_error = best_cv_error,
        clip_epsilon = clip_epsilon,
        n_clipped_test = as.integer(sum(clipped_test)),
        share_clipped_test = mean(clipped_test),
        p_hat_min = min(p_hat),
        p_hat_max = max(p_hat),
        g_mae = mean(abs(g_hat - g_true)),
        g_mse = mean((g_hat - g_true)^2),
        p_mae = mean(abs(p_hat - p_true)),
        p_mse = mean((p_hat - p_true)^2),
        brier_score = mean((Y_test - p_hat)^2),
        log_loss_clipped = -mean(
          Y_test * log(p_for_inverse) +
            (1 - Y_test) * log1p(-p_for_inverse)
        ),
        class_error = mean((p_hat > 0.5) != Y_test),
        class_error_true = mean((p_true > 0.5) != Y_test),
        error_message = NA_character_
      )
    }),
    error = function(e) {
      list(
        successful = FALSE,
        error_message = conditionMessage(e)
      )
    }
  )

  replication_result$row_index <- row_index
  replication_result$runtime_seconds <-
    proc.time()[["elapsed"]] - replication_start_time

  replication_result
}



##############################################################################
######### Parallel Monte Carlo loop #########
##############################################################################

pending_rows <- which(!simulation_results$completed)
number_pending <- length(pending_rows)

if (number_pending == 0L) {
  message(
    "All ordinary-RF plug-in Monte Carlo replications have already ",
    "been completed."
  )
} else {
  message(
    "Starting ",
    number_pending,
    " pending ordinary-RF plug-in replications with ",
    n_workers,
    " parallel workers. Each ranger fit uses ",
    num_threads,
    " thread."
  )

  project_dir <- getwd()

  cluster <- parallel::makePSOCKcluster(n_workers)
  pending_batches <- split(
    pending_rows,
    ceiling(seq_along(pending_rows) / parallel_batch_size)
  )

  processed_pending <- 0L
  results_since_excel <- 0L
  parallel_loop_error <- NULL
  start_time <- Sys.time()

  tryCatch(
    {
      parallel::clusterExport(
        cluster,
        varlist = c(
          "project_dir",
          "validate_probability_vector",
          "run_rfplugin_replication",
          "n_test",
          "B",
          "k_folds",
          "min_node_size_grid",
          "num_threads"
        ),
        envir = environment()
      )

      parallel::clusterEvalQ(
        cluster,
        {
          setwd(project_dir)
          library("ranger")
          source("./functions/generate_data.R")
          NULL
        }
      )

      successful_result_fields <- c(
        mtry_cv = "best_mtry",
        depth_cv = "best_depth",
        min_node_size_cv = "best_min_node_size",
        cv_error = "best_cv_error",
        clip_epsilon = "clip_epsilon",
        n_clipped_test = "n_clipped_test",
        share_clipped_test = "share_clipped_test",
        p_hat_min = "p_hat_min",
        p_hat_max = "p_hat_max",
        g_mae = "g_mae",
        g_mse = "g_mse",
        p_mae = "p_mae",
        p_mse = "p_mse",
        brier_score = "brier_score",
        log_loss_clipped = "log_loss_clipped",
        class_error = "class_error",
        class_error_true = "class_error_true"
      )

      for (batch_index in seq_along(pending_batches)) {
        batch_rows <- pending_batches[[batch_index]]
        batch_tasks <- lapply(
          batch_rows,
          function(row_index) {
            list(
              row_index = row_index,
              DGP = simulation_results$DGP[row_index],
              n = simulation_results$n[row_index],
              seed = simulation_results$seed[row_index]
            )
          }
        )

        batch_results <- parallel::parLapplyLB(
          cluster,
          batch_tasks,
          run_rfplugin_replication
        )

        # Only the master R process updates and writes the result table.
        for (replication_result in batch_results) {
          row_index <- replication_result$row_index
          simulation_results$runtime_seconds[row_index] <-
            replication_result$runtime_seconds

          if (replication_result$successful) {
            for (result_name in names(successful_result_fields)) {
              source_name <- successful_result_fields[[result_name]]
              simulation_results[[result_name]][row_index] <-
                replication_result[[source_name]]
            }

            simulation_results$completed[row_index] <- TRUE
            simulation_results$error_message[row_index] <- NA_character_
          } else {
            simulation_results$completed[row_index] <- FALSE
            simulation_results$error_message[row_index] <-
              replication_result$error_message

            warning(
              "Replication failed for DGP = ",
              simulation_results$DGP[row_index],
              ", n = ",
              simulation_results$n[row_index],
              ", seed = ",
              simulation_results$seed[row_index],
              ": ",
              replication_result$error_message,
              call. = FALSE
            )
          }
        }

        processed_pending <- processed_pending + length(batch_results)
        results_since_excel <- results_since_excel + length(batch_results)

        save(
          simulation_results,
          file = Rdata_name
        )

        if (
          results_since_excel >= excel_save_frequency ||
          batch_index == length(pending_batches)
        ) {
          openxlsx::write.xlsx(
            simulation_results,
            file = xlsx_name,
            overwrite = TRUE
          )
          results_since_excel <- 0L
        }

        rm(batch_tasks, batch_results)
        invisible(gc())

        message(sprintf(
          paste0(
            "Completed %d of %d pending replications; ",
            "total hours elapsed: %.4f"
          ),
          processed_pending,
          number_pending,
          as.numeric(difftime(
            Sys.time(),
            start_time,
            units = "hours"
          ))
        ))
      }
    },
    error = function(e) {
      parallel_loop_error <<- e
    },
    finally = {
      parallel::stopCluster(cluster)
    }
  )

  if (!is.null(parallel_loop_error)) {
    stop(parallel_loop_error)
  }
}



##############################################################################
######### Final save #########
##############################################################################
save(simulation_results, file = Rdata_name)
openxlsx::write.xlsx(simulation_results, file = xlsx_name, overwrite = TRUE)

number_incomplete <- sum(!simulation_results$completed)

if (number_incomplete == 0L) {
  message("The Standard random forest simulation is complete.")
} else {
  warning(
    number_incomplete,
    " standard random forest replications remain incomplete. ",
    "Restart the script to retry them.",
    call. = FALSE
  )
}
