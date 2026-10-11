##############################################################################

# Random Forests for Nonsmooth Binary Choice Models
# Simulation Study
# The program replicates Tables 1-4
# Author: Jihuan Zhang, Department of Economics, Indiana University Bloomington
# Email: zhangjih@iu.edu

##############################################################################



##############################################################################
######### Run simulation #########
##############################################################################
# Linear probit (LP)
source("./simulation_LinearProbit.R")

# Second-order polynomial probit (SPP)
source("./simulation_SPP.R")

# Kernalization-based probit (KP)
source("./simulation_KP_parallel.R")

# Random forest probit (RFP)
source("./simulation_RFP_parallel.R")

# Standard random forest regression (RF)
source("./simulation_RF_parallel.R")



##############################################################################
######### Generate tables #########
##############################################################################
source("./generate_simulation_tables.R")