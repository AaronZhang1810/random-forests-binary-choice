##############################################################################

# Random Forests for Nonsmooth Binary Choice Models
# Simulation Study: Random Forest Probit
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
library("ranger")
library("openxlsx")

source("./functions/generate_data.R")
source("./functions/resume_simulation_results.R")



##############################################################################
######### Simulation setting #########
##############################################################################
# Specify output folder and names
xlsx_name  = "./tables/simulation_results_RFP.xlsx"
Rdata_name = "./tables/simulation_results_RFP.Rdata"

DGP_grid <- c("sin","sigmoid_1","sigmoid_100","threshold","linear")
n_test <- 10000 # test sample size
B <- 500 # number of trees to grow
seed_grid <- 1810:(1809 + 1000L) # for MC replications
k_folds <- 5L
n_grid <- c(1000,2000,5000) # sample size
min_node_size_grid <- c(5,25,50,75,100)
varphi_grid <- c(0.001,0.1,1) # parameter for pseudo count
excel_save_frequency <- 10L # Write Excel after this many forest-parameter groups
num_threads <- 1L # Use one ranger thread for reproducibility



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

# Cross-validated hyperparameters
simulation_results$mtry_cv <- NA_integer_
simulation_results$depth_cv <- NA_integer_
simulation_results$min_node_size_cv <- NA_integer_
simulation_results$varphi_cv <- NA_real_
simulation_results$cv_error <- NA_real_

# Test-sample performance of the RFP estimator
simulation_results$g_mae <- NA_real_
simulation_results$g_mse <- NA_real_

# Plug-in CCP implied by the RFP estimator: F(g_hat)
simulation_results$p_mae <- NA_real_
simulation_results$p_mse <- NA_real_

# Direct regression-forest CCP estimator, retained as a diagnostic benchmark
simulation_results$p_rf_mae <- NA_real_
simulation_results$p_rf_mse <- NA_real_

# Classification performance
simulation_results$class_error <- NA_real_
simulation_results$class_error_rf <- NA_real_
simulation_results$class_error_true <- NA_real_

# Numerical and progress information
simulation_results$runtime_seconds <- NA_real_
simulation_results$completed <- FALSE
simulation_results$error_message <- NA_character_

# Resume previous simulation results
simulation_results <- resume_simulation_results(
  simulation_results = simulation_results,
  Rdata_name = Rdata_name,
  required_result_columns = setdiff(
    names(simulation_results),
    c("DGP", "n", "seed", "completed", "error_message")
  )
)



##############################################################################
######### Main Monte Carlo loop #########
##############################################################################
# Check what's left to be done
pending_rows <- which(!simulation_results$completed)
number_pending <- length(pending_rows)

if (number_pending == 0L) {
  message("All RFP Monte Carlo replications have already been completed.")
}

start_time <- Sys.time()

for (replication_index in seq_along(pending_rows)) {
  row_index <- pending_rows[replication_index]
  
  DGP_current <- simulation_results$DGP[row_index]
  n_current <- simulation_results$n[row_index]
  seed_current <- simulation_results$seed[row_index]
  
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
  
  replication_start_time <- proc.time()[["elapsed"]]
  
  replication_result <- tryCatch(
    local({ # use "local" to save memory
      
      
      
      ########################################################################
      # Generate training and test samples
      ########################################################################
      full_sample <- generate_sample(DGP = DGP_current, n = 20000L, seed = seed_current)
      
      X_full <- full_sample$data[,-1,drop = FALSE]
      Y_full <- as.numeric(full_sample$data$Y)
      
      g_full <- full_sample$g
      
      training_index <- seq_len(n_current)
      
      test_index <- 10001:20000
      
      X_train <- X_full[training_index,,drop = FALSE]
      Y_train <- Y_full[training_index]
      
      X_test <- X_full[test_index,,drop = FALSE]
      Y_test <- Y_full[test_index]
      
      g_true <- g_full[test_index]
      p_true <- pnorm(g_true)
      
      d <- ncol(X_train)
      
      
      
      ########################################################################
      # Candidate hyperparameters
      ########################################################################
      # mtry_grid <- unique(floor(c(sqrt(d),d/3,d)))
      mtry_grid <- sort(unique(pmin(pmax(as.integer(round(c(sqrt(d),d/3,d))),1L),d)))
      
      # In ranger, max.depth = 0 denotes unlimited depth.
      depth_grid <- sort(
        unique(
          as.integer(
            c(
              0L,
              floor(0.5 * log2(n_current)),
              floor(log2(n_current)),
              floor(1.5 * log2(n_current)),
              floor(2 * log2(n_current))
            )
          )
        )
      )
      
      forest_tuning_grid <- expand.grid(
        mtry = mtry_grid,
        depth = depth_grid,
        min.node.size = min_node_size_grid,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      
      number_forest_candidates <- nrow(forest_tuning_grid)
      
      # a different varphi does not need to grow a new forest
      number_varphi <- length(varphi_grid)
      
      
      
      ########################################################################
      # Construct stratified K-fold assignments
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
      # K-fold cross-validation
      ########################################################################
      
      # cv_sse[candidate, varphi, fold] stores squared-error sums using
      # Phi(g_hat), where g_hat averages the treewise inverse-probit estimates.
      cv_sse <- array(
        NA_real_,
        dim = c(
          number_forest_candidates,
          number_varphi,
          k_folds
        )
      )
      
      # Calculate CV error for all folds
      for (fold in seq_len(k_folds)) {
        
        # Get out-of-fold training sample and in-fold validation sample
        validation_index <- which(fold_id == fold)
        
        estimation_index <- which(fold_id != fold)
        
        X_estimation <- X_train[estimation_index,,drop = FALSE]
        Y_estimation <- Y_train[estimation_index]
        
        X_validation <- X_train[validation_index,,drop = FALSE]
        Y_validation <- Y_train[validation_index]
        
        # Calculate CV error for all hyperparameter choices for the current fold
        for (candidate_index in seq_len(number_forest_candidates)) {
          
          mtry_candidate <- forest_tuning_grid$mtry[candidate_index]
          depth_candidate <- forest_tuning_grid$depth[candidate_index]
          min_node_size_candidate <- forest_tuning_grid$min.node.size[candidate_index]
          
          # Numeric Y with splitrule = "variance" gives regression trees
          # grown using the squared-loss CART criterion.
          md_cv <- ranger(
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
            keep.inbag = TRUE,
            write.forest = TRUE,
            oob.error = FALSE,
            importance = "none",
            seed = seed_current + 100000L * fold,
            num.threads = num_threads,
            verbose = FALSE
          )
          
          # Identifies the leaf reached by each observation in every tree
          # The result is an n_estimation*B matrix
          terminal_estimation <- predict(
            md_cv,
            data = X_estimation,
            type = "terminalNodes",
            num.threads = num_threads,
            verbose = FALSE
          )$predictions
          
          terminal_validation <- predict(
            md_cv,
            data = X_validation,
            type = "terminalNodes",
            num.threads = num_threads,
            verbose = FALSE
          )$predictions
          
          # Get number of times an observation is selected in a bootstrap sample
          if (is.list(md_cv$inbag.counts)) {
            
            inbag_counts_cv <- do.call(
              cbind,
              md_cv$inbag.counts
            )
            
          } else {
            
            inbag_counts_cv <- as.matrix(
              md_cv$inbag.counts
            )
          }
          
          # Correct the orientation for inbag_counts_cv if necessary
          if (
            nrow(inbag_counts_cv) != nrow(X_estimation) &&
            ncol(inbag_counts_cv) == nrow(X_estimation)
          ) {
            inbag_counts_cv <- t(
              inbag_counts_cv
            )
          }
          
          if (
            nrow(inbag_counts_cv) != nrow(X_estimation) ||
            ncol(inbag_counts_cv) != B
          ) {
            stop(
              "The dimensions of a CV in-bag count matrix are incorrect."
            )
          }
          
          g_sum_validation <- matrix(
            0,
            nrow = nrow(X_validation),
            ncol = number_varphi
          )
          
          # Calculate g_hat for all validation observations 
          for (tree in seq_len(B)) {
            
            # Calculate within-leaf sample mean
            # An OOB observation has inbag_counts_cv[i, tree] == 0,
            # so its contribution is zero to both quantities
            leaf_count <- rowsum(
              inbag_counts_cv[, tree],
              group = terminal_estimation[, tree],
              reorder = FALSE
            )
            
            leaf_success <- rowsum(
              inbag_counts_cv[, tree] * Y_estimation,
              group = terminal_estimation[, tree],
              reorder = FALSE
            )
            
            leaf_count_vector <- drop(leaf_count)
            leaf_success_vector <- drop(leaf_success)
            
            if (any(leaf_count_vector <= 0)) {
              stop(
                "A CV terminal leaf has a nonpositive in-bag sample size."
              )
            }
            
            leaf_ids <- as.integer(rownames(leaf_count))
            
            validation_leaf_index <- match(
              terminal_validation[, tree],
              leaf_ids
            )
            
            if (anyNA(validation_leaf_index)) {
              stop(
                paste0(
                  "A validation terminal node could not be matched ",
                  "to a fitted leaf."
                )
              )
            }
            
            leaf_probability_raw <- leaf_success_vector / leaf_count_vector
            
            pure_leaf <-
              leaf_success_vector == 0 |
              leaf_success_vector == leaf_count_vector
            
            for (varphi_index in seq_len(number_varphi)) {
              
              varphi_candidate <- varphi_grid[varphi_index]
              
              # Impure leaves keep their ordinary sample proportions.
              leaf_probability_for_inverse <- leaf_probability_raw
              
              # The pseudo-count is applied only to pure leaves.
              leaf_probability_for_inverse[pure_leaf] <-
                (
                  leaf_success_vector[pure_leaf] +
                    varphi_candidate
                ) /
                (
                  leaf_count_vector[pure_leaf] +
                    2 * varphi_candidate
                )
              
              # within-leaf sample mean for every validation data point
              p_tree_for_inverse <-
                leaf_probability_for_inverse[validation_leaf_index]
              
              if (
                any(p_tree_for_inverse <= 0) ||
                any(p_tree_for_inverse >= 1)
              ) {
                stop(
                  paste0(
                    "A probability used in the CV inverse probit ",
                    "is not in (0,1)."
                  )
                )
              }
              
              g_sum_validation[,
                varphi_index
              ] <- g_sum_validation[
                ,
                varphi_index
              ] + qnorm(
                p_tree_for_inverse
              )
            }
          }
          
          g_validation <- g_sum_validation / B
          
          # Calculate structured p_hat for each validation data point
          p_validation_rfp <- pnorm(g_validation)
          
          if (
            any(!is.finite(g_validation)) ||
            any(!is.finite(p_validation_rfp))
          ) {
            stop(
              "A CV forest produced nonfinite RFP predictions."
            )
          }
          
          # Calculate cross-validation error for the current hyperparameters & fold
          for (varphi_index in seq_len(number_varphi)) {
            
            cv_sse[
              candidate_index,
              varphi_index,
              fold
            ] <- sum(
              (
                Y_validation -
                  p_validation_rfp[, varphi_index]
              )^2
            )
          }
          
          rm(
            md_cv,
            terminal_estimation,
            terminal_validation,
            inbag_counts_cv,
            g_sum_validation,
            g_validation,
            p_validation_rfp
          )
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
      
      # Calculate average CV error across folds
      cv_error_matrix <- apply(cv_sse, MARGIN = c(1, 2), FUN = sum) / n_current
      
      # Create a joint set of hyperparameters with CV error
      joint_tuning_grid <- expand.grid(
        forest_candidate = seq_len(number_forest_candidates),
        varphi_index = seq_len(number_varphi),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      
      joint_tuning_grid$mtry <- forest_tuning_grid$mtry[
        joint_tuning_grid$forest_candidate
      ]
      
      joint_tuning_grid$depth <- forest_tuning_grid$depth[
        joint_tuning_grid$forest_candidate
      ]
      
      joint_tuning_grid$min.node.size <-
        forest_tuning_grid$min.node.size[
          joint_tuning_grid$forest_candidate
        ]
      
      joint_tuning_grid$varphi <- varphi_grid[
        joint_tuning_grid$varphi_index
      ]
      
      joint_tuning_grid$cv_error <- cv_error_matrix[
        cbind(
          joint_tuning_grid$forest_candidate,
          joint_tuning_grid$varphi_index
        )
      ]
      
      if (any(!is.finite(joint_tuning_grid$cv_error))) {
        stop(
          "At least one joint hyperparameter combination has nonfinite CV error."
        )
      }
      
      # Obtain the best tuning parameters
      # Break exact ties in favor of a shallower tree, a larger
      # min.node.size, a smaller mtry, and finally a smaller varphi.
      depth_complexity <- ifelse(
        joint_tuning_grid$depth == 0L,
        Inf,
        joint_tuning_grid$depth
      )
      
      best_order <- order(
        joint_tuning_grid$cv_error,
        depth_complexity,
        -joint_tuning_grid$min.node.size,
        joint_tuning_grid$mtry,
        joint_tuning_grid$varphi
      )
      
      best_index <- best_order[1]
      
      # Extract optimal hyperparameters
      best_mtry <- as.integer(joint_tuning_grid$mtry[best_index])
      best_depth <- as.integer(joint_tuning_grid$depth[best_index])
      best_min_node_size <- as.integer(joint_tuning_grid$min.node.size[best_index])
      best_varphi <- joint_tuning_grid$varphi[best_index]
      
      # Extract min CV error
      best_cv_error <- joint_tuning_grid$cv_error[best_index]
      
      
      
      ########################################################################
      # Fit the selected forest on the full training sample
      ########################################################################
      md_full <- ranger(
        y = Y_train,
        x = X_train,
        classification = FALSE,
        probability = FALSE,
        num.trees = B,
        mtry = best_mtry, # optimal
        splitrule = "variance",
        max.depth = best_depth, # optimal
        min.node.size = best_min_node_size, # optimal
        replace = TRUE,
        sample.fraction = 1,
        keep.inbag = TRUE,
        write.forest = TRUE,
        oob.error = FALSE,
        importance = "none",
        seed = seed_current + 900000L,
        num.threads = num_threads,
        verbose = FALSE
      )
      
      terminal_train <- predict(
        md_full,
        data = X_train,
        type = "terminalNodes",
        num.threads = num_threads,
        verbose = FALSE
      )$predictions
      
      terminal_test <- predict(
        md_full,
        data = X_test,
        type = "terminalNodes",
        num.threads = num_threads,
        verbose = FALSE
      )$predictions
      
      if (is.list(md_full$inbag.counts)) {
        
        inbag_counts <- do.call(
          cbind,
          md_full$inbag.counts
        )
        
      } else {
        
        inbag_counts <- as.matrix(
          md_full$inbag.counts
        )
      }
      
      if (
        nrow(inbag_counts) != nrow(X_train) &&
        ncol(inbag_counts) == nrow(X_train)
      ) {
        inbag_counts <- t(
          inbag_counts
        )
      }
      
      if (
        nrow(inbag_counts) != nrow(X_train) ||
        ncol(inbag_counts) != B
      ) {
        stop(
          "The dimensions of the full-sample in-bag matrix are incorrect."
        )
      }
      
      
      
      ########################################################################
      # Construct test-sample estimates
      ########################################################################
      
      # Direct regression-forest CCP estimator
      p_rf_sum_test <- numeric(n_test)
      
      # RFP systematic-component estimator
      g_sum_test <- numeric(n_test)
      
      for (tree in seq_len(B)) {
        
        leaf_count <- rowsum(
          inbag_counts[, tree],
          group = terminal_train[, tree],
          reorder = FALSE
        )
        
        leaf_success <- rowsum(
          inbag_counts[, tree] * Y_train,
          group = terminal_train[, tree],
          reorder = FALSE
        )
        
        leaf_count_vector <- drop(
          leaf_count
        )
        
        leaf_success_vector <- drop(
          leaf_success
        )
        
        if (any(leaf_count_vector <= 0)) {
          stop(
            "A fitted terminal leaf has a nonpositive in-bag sample size."
          )
        }
        
        leaf_ids <- as.integer(
          rownames(leaf_count)
        )
        
        test_leaf_index <- match(
          terminal_test[, tree],
          leaf_ids
        )
        
        if (anyNA(test_leaf_index)) {
          stop(
            "A test terminal node could not be matched to a fitted leaf."
          )
        }
        
        leaf_probability_raw <-
          leaf_success_vector /
          leaf_count_vector
        
        p_rf_sum_test <- p_rf_sum_test +
          leaf_probability_raw[
            test_leaf_index
          ]
        
        pure_leaf <-
          leaf_success_vector == 0 |
          leaf_success_vector == leaf_count_vector
        
        leaf_probability_for_inverse <-
          leaf_probability_raw
        
        leaf_probability_for_inverse[pure_leaf] <-
          (
            leaf_success_vector[pure_leaf] +
              best_varphi
          ) /
          (
            leaf_count_vector[pure_leaf] +
              2 * best_varphi
          )
        
        p_tree_for_inverse <-
          leaf_probability_for_inverse[
            test_leaf_index
          ]
        
        if (
          any(p_tree_for_inverse <= 0) ||
          any(p_tree_for_inverse >= 1)
        ) {
          stop(
            "A test probability used in the inverse probit is not in (0,1)."
          )
        }
        
        g_sum_test <- g_sum_test + qnorm(
          p_tree_for_inverse
        )
      }
      
      p_rf_hat <- p_rf_sum_test / B
      g_hat <- g_sum_test / B
      p_hat <- pnorm(
        g_hat
      )
      
      if (
        any(!is.finite(p_rf_hat)) ||
        any(!is.finite(g_hat)) ||
        any(!is.finite(p_hat))
      ) {
        stop(
          "The selected RFP estimator produced nonfinite test predictions."
        )
      }
      
      
      
      ########################################################################
      # Test-sample performance
      ########################################################################
      
      g_mae <- mean(abs(g_hat - g_true))
      g_mse <- mean((g_hat -g_true)^2)
      
      # Main probability measure: plug-in probability implied by RFP.
      p_mae <- mean(abs(p_hat -p_true))
      p_mse <- mean((p_hat -p_true)^2)
      
      # Diagnostic probability measure: direct regression-forest CCP.
      p_rf_mae <- mean(abs(p_rf_hat -p_true))
      p_rf_mse <- mean((p_rf_hat -p_true)^2)
      
      class_error <- mean((p_hat > 0.5) != Y_test)
      class_error_rf <- mean((p_rf_hat > 0.5) != Y_test)
      
      class_error_true <- mean((p_true > 0.5) != Y_test)
      
      output <- list(
        successful = TRUE,
        best_mtry = best_mtry,
        best_depth = best_depth,
        best_min_node_size = best_min_node_size,
        best_varphi = best_varphi,
        best_cv_error = best_cv_error,
        g_mae = g_mae,
        g_mse = g_mse,
        p_mae = p_mae,
        p_mse = p_mse,
        p_rf_mae = p_rf_mae,
        p_rf_mse = p_rf_mse,
        class_error = class_error,
        class_error_rf = class_error_rf,
        class_error_true = class_error_true,
        error_message = NA_character_
      )
      
      rm(
        full_sample,
        X_full,
        Y_full,
        g_full,
        X_train,
        Y_train,
        X_test,
        Y_test,
        g_true,
        p_true,
        forest_tuning_grid,
        joint_tuning_grid,
        fold_id,
        cv_sse,
        cv_error_matrix,
        md_full,
        terminal_train,
        terminal_test,
        inbag_counts,
        p_rf_sum_test,
        g_sum_test,
        p_rf_hat,
        g_hat,
        p_hat
      )
      
      output
    }),
    
    error = function(e) {
      
      list(
        successful = FALSE,
        error_message = conditionMessage(e)
      )
    }
  )
  
  replication_runtime <- proc.time()[["elapsed"]] - replication_start_time
  
  
  
  ##########################################################################
  # Update result table
  ##########################################################################
  
  simulation_results$runtime_seconds[row_index] <- replication_runtime
  
  if (replication_result$successful) {
    
    simulation_results$mtry_cv[row_index] <- replication_result$best_mtry
    simulation_results$depth_cv[row_index] <- replication_result$best_depth
    simulation_results$min_node_size_cv[row_index] <-replication_result$best_min_node_size
    simulation_results$varphi_cv[row_index] <-replication_result$best_varphi
    simulation_results$cv_error[row_index] <-replication_result$best_cv_error
    simulation_results$g_mae[row_index] <-replication_result$g_mae
    simulation_results$g_mse[row_index] <-replication_result$g_mse
    simulation_results$p_mae[row_index] <-replication_result$p_mae
    simulation_results$p_mse[row_index] <-replication_result$p_mse
    simulation_results$p_rf_mae[row_index] <-replication_result$p_rf_mae
    simulation_results$p_rf_mse[row_index] <-replication_result$p_rf_mse
    simulation_results$class_error[row_index] <-replication_result$class_error
    simulation_results$class_error_rf[row_index] <-replication_result$class_error_rf
    simulation_results$class_error_true[row_index] <-replication_result$class_error_true
    simulation_results$completed[row_index] <- TRUE
    simulation_results$error_message[row_index] <- NA_character_
    
  } else {
    
    simulation_results$completed[row_index] <- FALSE
    
    simulation_results$error_message[row_index] <- replication_result$error_message
    
    warning(
      "Replication failed: ",
      replication_result$error_message,
      call. = FALSE
    )
  }
  
  
  
  ##########################################################################
  # Save after every replication
  ##########################################################################
  
  save(simulation_results,file = Rdata_name)
  
  # Update Excel file after every several replications
  if (
    replication_index %% excel_save_frequency == 0L || replication_index == number_pending
  ) {
    openxlsx::write.xlsx(simulation_results, file = xlsx_name, overwrite = TRUE)
  }
  
  rm(replication_result)
  invisible(gc())
  
  message(
    sprintf(
      "Current replication: %d; total hours elapsed: %.4f",
      replication_index,
      as.numeric(difftime(Sys.time(), start_time, units = "hours"))
    )
  )
}



##############################################################################
######### Final save #########
##############################################################################
save(simulation_results, file = Rdata_name)
openxlsx::write.xlsx(simulation_results, file = xlsx_name, overwrite = TRUE)

number_incomplete <- sum(!simulation_results$completed)

if (number_incomplete == 0L) {
  
  message("The RFP Monte Carlo simulation is complete.")
  
} else {
  
  warning(
    number_incomplete,
    " RFP replications remain incomplete. ",
    "Restart the script to retry them.",
    call. = FALSE
  )
}
