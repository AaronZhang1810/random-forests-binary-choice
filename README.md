# Random Forests for Nonsmooth Binary Choice Models

**Jihuan Zhang**  
Department of Economics, Indiana University Bloomington

[Read the job market paper (PDF)](jmp.pdf) — October 11, 2026.

This repository contains the paper, its simulation replication code, and saved simulation results. The supplied package supports Tables 1–4 in the paper. It does not include the data or code for the LendingClub empirical application.

## Repository contents

- `jmp.pdf`: the job market paper.
- `replication_package/main.R`: runs the simulation drivers and generates tables.
- `replication_package/simulation_*.R`: drivers for linear probit (LP), second-order polynomial probit (SPP), kernelized probit (KP), random forest probit (RFP), and standard random forest regression (RF).
- `replication_package/functions/`: data generators, estimators, checkpoint recovery, and table generation.
- `replication_package/tables/`: saved `.Rdata` checkpoints, Excel results, and generated LaTeX tables.
- `replication_package/README.txt`: the original package documentation, dated October 9, 2026.

`main.R` uses the parallel KP, RFP, and RF drivers. Serial alternatives are included for KP and RFP.

## Requirements

The original package documentation reports testing with **R 4.5.2 on Windows 11**. The simulation and table drivers use **RStudio** to determine their working directory.

Install the required packages in R:

```r
install.packages(c("ranger", "openxlsx", "RSpectra", "nloptr", "rstudioapi"))
```

The scripts also use `parallel`, `stats`, and `utils`, which are included with R. Package versions are not pinned.

## Generate tables from the supplied results

1. Download or clone the repository.
2. Open `replication_package/generate_simulation_tables.R` in RStudio and keep it as the active editor tab.
3. Source that script.

The script reads the supplied Excel results and writes six LaTeX table fragments to `replication_package/tables/`: RMSE and MAD tables for all five methods, for the four methods excluding RF, and for RFP versus RF.

## Run or resume simulations

1. Start a separate RStudio session; the simulation drivers clear the R workspace.
2. Open `replication_package/main.R` and keep it as the active editor tab.
3. Set the working directory to `replication_package/`, then source `main.R`.

The current drivers require RStudio. Running `Rscript main.R` without adapting their working-directory setup is not supported.

The default design uses five data-generating processes, sample sizes of 1,000, 2,000, and 5,000, and 1,000 Monte Carlo replications per combination. The parallel drivers use up to eight workers. Worker counts and other settings can be adjusted in the individual drivers.

Existing `.Rdata` checkpoints are loaded automatically, and completed replications are retained. To recompute from scratch, use a separate working copy and move its supplied checkpoints out of `tables/` before running. Use separate output filenames when changing simulation settings, and do not run serial and parallel versions of the same method simultaneously against the same output files.

The scripts overwrite their result workbooks and generated tables. Table generation uses completed replications and prints the number used in each cell; check these counts before treating an interrupted run as complete.

## Contact

Questions or error reports: [zhangjih@iu.edu](mailto:zhangjih@iu.edu).

