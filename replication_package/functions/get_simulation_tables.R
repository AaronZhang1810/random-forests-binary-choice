##############################################################################

# Random Forests for Nonsmooth Binary Choice Models
# Simulation: Generate LaTeX RMSE and MAD Tables from Simulation Results
# Methods:
#   RFP = Random Forest Probit
#   KP  = Kernelized Probit
#   SPP = Second-order Polynomial Probit
#   LP  = Conventional Linear Probit
#   RF  = Random Forest Regression: direct p estimate, inverse-link g estimate
# Writes three versions of each table: all methods, no RF, and RFP versus RF.
#
# Author: Jihuan Zhang, Department of Economics, Indiana University Bloomington
# Email: zhangjih@iu.edu

##############################################################################



##############################################################################
######### Read and check one simulation file #########
##############################################################################
library("openxlsx")
read_simulation_results <- function(file_name,
                                    method_name,
                                    optional = FALSE) {
  
  # RF may be absent or have no completed runs while its simulation is running.
  # A typed empty result still allows every requested table to be generated.
  empty_results <- data.frame(
    DGP = character(), n = integer(), seed = integer(),
    g_mae = numeric(), g_mse = numeric(),
    p_mae = numeric(), p_mse = numeric(),
    completed = logical(), method = character(),
    stringsAsFactors = FALSE
  )
  
  if (!file.exists(file_name)) {
    if (optional) {
      warning("No result file for ", method_name, ": ", file_name,
              ". Its table entries will be --.", call. = FALSE)
      return(empty_results)
    }
    stop(
      "The following simulation result file does not exist: ",
      file_name
    )
  }
  
  results <- openxlsx::read.xlsx(
    file_name
  )
  
  required_columns <- c(
    "DGP",
    "n",
    "seed",
    "g_mae",
    "g_mse",
    "p_mae",
    "p_mse",
    "completed"
  )
  
  if (!all(
    required_columns %in% colnames(results)
  )) {
    
    missing_columns <- setdiff(
      required_columns,
      colnames(results)
    )
    
    stop(
      "The file ",
      file_name,
      " is missing the following required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
  
  # Keep only completed Monte Carlo replications
  completed_indicator <-
    results$completed %in% c(
      TRUE,
      1,
      "TRUE",
      "True",
      "true",
      "1"
    )
  
  results <- results[
    completed_indicator,
    required_columns,
    drop = FALSE
  ]
  
  if (nrow(results) == 0L) {
    if (optional) {
      warning("No completed replications in ", file_name,
              ". Its table entries will be --.", call. = FALSE)
      return(empty_results)
    }
    stop(
      "No completed replications were found in ",
      file_name,
      "."
    )
  }
  
  metric_columns <- c(
    "g_mae",
    "g_mse",
    "p_mae",
    "p_mse"
  )
  
  if (any(
    !vapply(
      results[metric_columns],
      function(x) all(is.finite(x)),
      logical(1)
    )
  )) {
    stop(
      "The completed rows in ",
      file_name,
      " contain nonfinite MAE or MSE values."
    )
  }
  
  # Each DGP-n-seed combination should identify one Monte Carlo replication.
  simulation_key <- paste(
    results$DGP,
    results$n,
    results$seed,
    sep = "|"
  )
  
  if (anyDuplicated(simulation_key)) {
    stop(
      "The file ",
      file_name,
      " contains duplicate DGP-n-seed rows."
    )
  }
  
  results$method <- method_name
  
  results
}



##############################################################################
######### Read all five estimators #########
##############################################################################

RFP_results <- read_simulation_results(
  file_name = RFP_file,
  method_name = "RFP"
)

KP_results <- read_simulation_results(
  file_name = KP_file,
  method_name = "KP"
)

SPP_results <- read_simulation_results(
  file_name = SPP_file,
  method_name = "SPP"
)

LP_results <- read_simulation_results(
  file_name = LP_file,
  method_name = "LP"
)

RF_results <- read_simulation_results(
  file_name = RF_file,
  method_name = "RF",
  optional = TRUE
)


##############################################################################
######### Combine simulation results #########
##############################################################################

simulation_results <- rbind(
  RFP_results,
  KP_results,
  SPP_results,
  LP_results,
  RF_results
)

simulation_results$method <- factor(
  simulation_results$method,
  levels = method_order
)

simulation_results$DGP <- factor(
  simulation_results$DGP,
  levels = DGP_order
)



##############################################################################
######### Calculate Monte Carlo summaries #########
##############################################################################
summary_observed <- aggregate(
  cbind(
    g_mae,
    g_mse,
    p_mae,
    p_mse
  ) ~ DGP + n + method,
  data = simulation_results,
  FUN = mean
)

summary_observed$g_rmse <- sqrt(
  summary_observed$g_mse
)

summary_observed$p_rmse <- sqrt(
  summary_observed$p_mse
)

summary_observed$g_mad <- summary_observed$g_mae
summary_observed$p_mad <- summary_observed$p_mae


# Number of completed Monte Carlo replications used in each cell.
replication_count <- aggregate(
  seed ~ DGP + n + method,
  data = simulation_results,
  FUN = length
)

colnames(replication_count)[
  colnames(replication_count) == "seed"
] <- "n_replications"

summary_observed <- merge(
  summary_observed,
  replication_count,
  by = c(
    "DGP",
    "n",
    "method"
  ),
  all.x = TRUE,
  sort = FALSE
)


##############################################################################
######### Construct the complete table grid #########
##############################################################################

# Construct all table cells so the LaTeX output can still be generated if a
# particular DGP-n-method cell currently has zero completed replications.
# Such a cell is printed as "--".

complete_grid <- expand.grid(
  DGP = DGP_order,
  n = n_order,
  method = method_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

summary_table <- merge(
  complete_grid,
  summary_observed,
  by = c(
    "DGP",
    "n",
    "method"
  ),
  all.x = TRUE,
  sort = FALSE
)

summary_table$n_replications[is.na(summary_table$n_replications)] <- 0L

summary_table$method <- factor(
  summary_table$method,
  levels = method_order
)

summary_table$DGP <- factor(
  summary_table$DGP,
  levels = DGP_order
)


##############################################################################
######### Display interim replication counts #########
##############################################################################

message(
  "Completed Monte Carlo replications currently used in each table cell:"
)

print(
  summary_table[
    order(
      summary_table$DGP,
      summary_table$n,
      summary_table$method
    ),
    c(
      "DGP",
      "n",
      "method",
      "n_replications"
    )
  ],
  row.names = FALSE
)

rf_comparison_counts <- merge(
  summary_table[as.character(summary_table$method) == "RFP",
                c("DGP", "n", "n_replications")],
  summary_table[as.character(summary_table$method) == "RF",
                c("DGP", "n", "n_replications")],
  by = c("DGP", "n"), suffixes = c("_RFP", "_RF")
)
if (any(rf_comparison_counts$n_replications_RFP !=
        rf_comparison_counts$n_replications_RF)) {
  message(
    "NOTE: RFP and RF have unequal completed replication counts in some cells. ",
    "Every table uses all completed replications for each method, not matched seeds."
  )
}



##############################################################################
######### LaTeX formatting functions #########
##############################################################################

format_number <- function(x,
                          digits = number_digits) {
  
  formatted <- rep(
    "--",
    length(x)
  )
  
  finite_values <- is.finite(x)
  
  formatted[finite_values] <- sprintf(
    paste0(
      "%.",
      digits,
      "f"
    ),
    x[finite_values]
  )
  
  formatted
}


format_sample_size <- function(n) {
  
  n_formatted <- format(
    as.integer(n),
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
  
  # In math mode, use {,} so LaTeX does not add punctuation spacing.
  n_formatted <- gsub(
    ",",
    "{,}",
    n_formatted,
    fixed = TRUE
  )
  
  paste0(
    "$",
    n_formatted,
    "$"
  )
}



##############################################################################
######### General function to generate one LaTeX table #########
##############################################################################

generate_latex_table <- function(metric = c("RMSE", "MAD"),
                                 tex_name,
                                 methods = method_order) {
  
  metric <- match.arg(metric)
  
  if (missing(tex_name) || length(tex_name) != 1L || !nzchar(tex_name)) {
    stop("tex_name must be a nonempty file name.")
  }
  
  if (!is.character(methods) || length(methods) == 0L ||
      anyNA(methods) || anyDuplicated(methods) ||
      !all(methods %in% method_order)) {
    stop("methods must be a nonempty, unique subset of method_order.")
  }
  n_methods <- length(methods)
  n_columns <- 1L + 2L * n_methods
  method_header <- paste0("& ", paste(methods, collapse = " & "))
  column_block <- paste(rep("c", n_methods), collapse = "")
  
  if (metric == "RMSE") {
    g_column <- "g_rmse"
    p_column <- "p_rmse"
  } else {
    g_column <- "g_mad"
    p_column <- "p_mad"
  }
  
  
  extract_table_row <- function(DGP_current,
                                n_current) {
    
    current_results <- summary_table[
      as.character(summary_table$DGP) == DGP_current &
        summary_table$n == n_current,
      ,
      drop = FALSE
    ]
    
    method_index <- match(
      methods,
      as.character(current_results$method)
    )
    
    if (anyNA(method_index)) {
      stop(
        "At least one method is missing for DGP = ",
        DGP_current,
        " and n = ",
        n_current,
        "."
      )
    }
    
    current_results <- current_results[
      method_index,
      ,
      drop = FALSE
    ]
    
    g_values <- format_number(
      current_results[[g_column]]
    )
    
    p_values <- format_number(
      current_results[[p_column]]
    )
    
    paste0(
      "\t",
      format_sample_size(n_current),
      "\n",
      "\t& ",
      paste(
        g_values,
        collapse = " & "
      ),
      "\n",
      "\t& ",
      paste(
        p_values,
        collapse = " & "
      ),
      "\n",
      "\t\\\\"
    )
  }
  
  
  
  ##########################################################################
  ######### Generate LaTeX tabular content only #########
  ##########################################################################
  
  # Caption, label, font size, table environment, and notes are deliberately
  # specified manually in the LaTeX manuscript.
  
  table_lines <- c(
    paste0("\\begin{tabular}{l@{\\hspace{2em}}", column_block,
           "@{\\hspace{2em}}", column_block, "}"),
    "\\toprule",
    "&",
    paste0(
      "\\multicolumn{", n_methods, "}{c}{",
      metric,
      " of $\\widehat{g}$}"
    ),
    "&",
    paste0(
      "\\multicolumn{", n_methods, "}{c}{",
      metric,
      " of $\\widehat{p}$}"
    ),
    "\\\\",
    paste0("\\cmidrule(lr){2-", n_methods + 1L, "}"),
    paste0("\\cmidrule(lr){", n_methods + 2L, "-", n_columns, "}"),
    "$n$",
    method_header,
    method_header,
    "\\\\",
    "\\midrule",
    ""
  )
  
  for (DGP_current in DGP_order) {
    
    table_lines <- c(
      table_lines,
      paste0(
        "\t\\multicolumn{", n_columns, "}{l}{\\textit{",
        panel_titles[[DGP_current]],
        "}}\\\\"
      )
    )
    
    for (n_current in n_order) {
      
      table_lines <- c(
        table_lines,
        extract_table_row(
          DGP_current = DGP_current,
          n_current = n_current
        )
      )
    }
    
    if (DGP_current != tail(DGP_order, 1L)) {
      table_lines <- c(
        table_lines,
        "\t\\addlinespace",
        ""
      )
    }
  }
  
  table_lines <- c(
    table_lines,
    "\t\\bottomrule",
    "\\end{tabular}"
  )
  
  
  
  ##########################################################################
  ######### Save and display LaTeX table #########
  ##########################################################################
  
  writeLines(
    table_lines,
    con = tex_name
  )
  
  cat(
    paste(
      table_lines,
      collapse = "\n"
    ),
    "\n"
  )
  
  message(
    "The LaTeX ",
    metric,
    " tabular content was saved to: ",
    tex_name
  )
  
  invisible(table_lines)
}



##############################################################################
######### Generate all six tables #########
##############################################################################

RMSE_table <- generate_latex_table(
  metric = "RMSE",
  tex_name = RMSE_tex_name
)

MAD_table <- generate_latex_table(
  metric = "MAD",
  tex_name = MAD_tex_name
)

RMSE_no_RF_table <- generate_latex_table(
  metric = "RMSE",
  tex_name = RMSE_no_RF_tex_name,
  methods = method_order_no_RF
)

MAD_no_RF_table <- generate_latex_table(
  metric = "MAD",
  tex_name = MAD_no_RF_tex_name,
  methods = method_order_no_RF
)

RMSE_RFP_vs_RF_table <- generate_latex_table(
  metric = "RMSE",
  tex_name = RMSE_RFP_vs_RF_tex_name,
  methods = method_order_RFP_vs_RF
)

MAD_RFP_vs_RF_table <- generate_latex_table(
  metric = "MAD",
  tex_name = MAD_RFP_vs_RF_tex_name,
  methods = method_order_RFP_vs_RF
)
