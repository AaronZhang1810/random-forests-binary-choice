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
######### Housekeeping #########
##############################################################################
# Set this script's folder as the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls()); cat("\f")

# Load packages
library("openxlsx")



##############################################################################
######### File names #########
##############################################################################

# Simulation output files
RFP_file <- "./tables/simulation_results_RFP.xlsx"
KP_file  <- "./tables/simulation_results_KP.xlsx"
SPP_file <- "./tables/simulation_results_SPP.xlsx"
LP_file  <- "./tables/simulation_results_LinearProbit.xlsx"
RF_file  <- "./tables/simulation_results_RF.xlsx"

# Full tables, including RF immediately after LP.
RMSE_tex_name <- "./tables/simulation_table_RMSE.tex"
MAD_tex_name  <- "./tables/simulation_table_MAD.tex"

# Original four-method tables and two-method comparison tables.
RMSE_no_RF_tex_name <- "./tables/simulation_table_RMSE_no_RF.tex"
MAD_no_RF_tex_name  <- "./tables/simulation_table_MAD_no_RF.tex"
RMSE_RFP_vs_RF_tex_name <- "./tables/simulation_table_RMSE_RFP_vs_RF.tex"
MAD_RFP_vs_RF_tex_name  <- "./tables/simulation_table_MAD_RFP_vs_RF.tex"



##############################################################################
######### Table Settings #########
##############################################################################
# Specify orders
method_order <- c("RFP","KP","SPP","LP","RF")
method_order_no_RF <- setdiff(method_order, "RF")
method_order_RFP_vs_RF <- c("RFP", "RF")
DGP_order <- c("sin","sigmoid_1","sigmoid_100","threshold","linear")
n_order <- c(1000L,2000L,5000L)

# Specify panel titles
panel_titles <- c(
  sin = "Panel A: Sinusoidal",
  sigmoid_1 = "Panel B: Sigmoid, $\\kappa=1$",
  sigmoid_100 = "Panel C: Sigmoid, $\\kappa=100$",
  threshold = "Panel D: Threshold",
  linear = "Panel E: Sparse linear"
)

# Specify number of decimals
number_digits <- 4L



##############################################################################
######### Generate tables #########
##############################################################################
source("./functions/get_simulation_tables.R")


