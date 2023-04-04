# Bayesian viral phylogenetic source attribution from consensus sequences and time since infection estimates

This repository includes code and partial data for the analyses in Blenkinsop, Sofocleous, di Lauro et al. aRxiv (2023)

## Data
The data folder contains the following input files:
* Anonymised transmission pairs including covariates:
    * age of transmitter and recipient at estimated time of transmission of the incident case
    * time elapsed
    * patristic distance
* Data from Belgian study for fitting the evolutionary clock meta-analysis model [cite]

## Code
The scripts folder contains the following scripts to run the simulation study:
1) make-stan-data.R - simulate data and run stan model using cmdstan. The input arguments at the top of the file can be modified to run the vanilla, covariate or HSGP BMM accordingly. The average number of phylogenetically possible sources per incident case is modified by changing the argument p_pairs, with the relationship p=1/c.
2) post-processing.R - for a continuous source covariate (age). Produce convergence diagnostics, summarise generated quantities including transmission flows
2) post-processing-binary.R - for a binary source covariate. Produce convergence diagnostics, summarise generated quantities including transmission flows


The scripts folder contains the following scripts to run the application to Amsterdam MSM are as follows:
1) tree-distances.R - extract and save patristic distances by HIV subtype from Amsterdam MSM phylogenies
1) formulate-Amsterdam-pairs.R - formulate phylogenetically possible pairs using epidemiological and clinical data. This script cannot be run due to data confidentiality.
2) make-stan-data-Amsterdam.R - prepare stan data from prepared anonymised phylogenetically possible transmission pairs in /data/ directory
2) post-processing-Amsterdam.R - produce convergence diagnostics, summarise generated quantities including transmission flows
