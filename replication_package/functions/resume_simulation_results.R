##############################################################################

# Random Forests for Nonsmooth Binary Choice Models
# Simulation Study
# Author: Jihuan Zhang, Department of Economics, Indiana University Bloomington
# Email: zhangjih@iu.edu

##############################################################################



##############################################################################
######### Checkpoint helpers #########
##############################################################################
# Resume rows from an existing checkpoint
resume_simulation_results <- function(
    simulation_results,
    Rdata_name,
    required_result_columns = character(0),
    simulation_configuration = NULL,
    allow_legacy_configuration = TRUE) {

  required_key_columns <- c("DGP", "n", "seed")

  validate_result_table <- function(results, label) {
    if (!is.data.frame(results)) {
      stop(label, " simulation_results must be a data frame.")
    }

    missing_columns <- setdiff(
      c(required_key_columns, "completed"),
      colnames(results)
    )

    if (length(missing_columns) > 0L) {
      stop(
        label,
        " simulation_results is missing: ",
        paste(missing_columns, collapse = ", "),
        "."
      )
    }

    if (!is.logical(results$completed)) {
      stop(label, " simulation_results$completed must be logical.")
    }

    if (anyNA(results[, required_key_columns, drop = FALSE])) {
      stop(label, " simulation keys must not contain missing values.")
    }

    keys <- do.call(
      paste,
      c(results[required_key_columns], list(sep = "|"))
    )

    if (anyDuplicated(keys)) {
      stop(label, " simulation_results contains duplicate DGP-n-seed rows.")
    }

    keys
  }

  values_are_valid <- function(values) {
    valid <- !is.na(values)

    if (is.numeric(values)) {
      valid <- valid & is.finite(values)
    } else if (is.character(values)) {
      valid <- valid & nzchar(values)
    }

    valid
  }

  current_key <- validate_result_table(simulation_results, "Current")

  missing_required <- setdiff(
    required_result_columns,
    colnames(simulation_results)
  )
  if (length(missing_required) > 0L) {
    stop(
      "The current simulation table is missing required result columns: ",
      paste(missing_required, collapse = ", "),
      "."
    )
  }

  simulation_results$completed[is.na(simulation_results$completed)] <- FALSE

  if (!file.exists(Rdata_name)) {
    return(simulation_results)
  }

  previous_environment <- new.env(parent = emptyenv())
  loaded_objects <- load(Rdata_name, envir = previous_environment)

  if (!("simulation_results" %in% loaded_objects)) {
    stop("The existing RData file does not contain simulation_results.")
  }

  previous_results <- get(
    "simulation_results",
    envir = previous_environment,
    inherits = FALSE
  )
  previous_key <- validate_result_table(previous_results, "Existing")

  if (!is.null(simulation_configuration)) {
    if ("simulation_configuration" %in% loaded_objects) {
      previous_configuration <- get(
        "simulation_configuration",
        envir = previous_environment,
        inherits = FALSE
      )

      if (!identical(previous_configuration, simulation_configuration)) {
        stop(
          paste0(
            "The checkpoint was produced with a different simulation ",
            "configuration. Use a new checkpoint name or restore the old ",
            "settings instead of mixing results."
          )
        )
      }
    } else if (!isTRUE(allow_legacy_configuration)) {
      stop(
        "The checkpoint has no simulation_configuration metadata. ",
        "Use a new checkpoint or explicitly allow this legacy checkpoint."
      )
    } else {
      warning(
        paste0(
          "The existing checkpoint predates configuration metadata. ",
          "Compatible completed rows will be retained; the next checkpoint ",
          "write will record the current configuration."
        ),
        call. = FALSE
      )
    }
  }

  previous_index <- match(current_key, previous_key)
  matched_rows <- which(!is.na(previous_index))

  columns_to_copy <- intersect(
    setdiff(
      colnames(simulation_results),
      c(required_key_columns, "completed")
    ),
    colnames(previous_results)
  )

  for (column_name in columns_to_copy) {
    current_column <- simulation_results[[column_name]]
    previous_column <- previous_results[[column_name]]

    compatible_type <-
      (is.numeric(current_column) && is.numeric(previous_column)) ||
      identical(typeof(current_column), typeof(previous_column))

    if (!compatible_type) {
      warning(
        "Skipping incompatible checkpoint column: ",
        column_name,
        call. = FALSE
      )
      next
    }

    simulation_results[[column_name]][matched_rows] <-
      previous_column[previous_index[matched_rows]]
  }

  simulation_results$completed[] <- FALSE

  if (length(matched_rows) > 0L) {
    candidate_complete <- previous_results$completed[
      previous_index[matched_rows]
    ]
    candidate_complete[is.na(candidate_complete)] <- FALSE

    for (column_name in required_result_columns) {
      if (!(column_name %in% colnames(previous_results))) {
        candidate_complete[] <- FALSE
        break
      }

      candidate_complete <- candidate_complete & values_are_valid(
        simulation_results[[column_name]][matched_rows]
      )
    }

    simulation_results$completed[matched_rows] <- candidate_complete
  }

  message(length(matched_rows), " previous simulation rows were matched.")
  message(sum(simulation_results$completed), " replications are already complete.")

  simulation_results
}