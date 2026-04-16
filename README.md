# DSE4231 Group 3 Project

This project replicates and extends a causal mediation analysis of the Job Corps program using IPW, Random Forest, and Double Machine Learning approaches, following Huber (2014).

## Code Structure
- `data/` – dataset used in the analysis  
- `src/` – R scripts for replication and extensions  
  - `replication.R` – baseline IPW replication using probit models   
  - `rf_trim.R` – IPW replication using Random Forest-based propensity score models
  - `double_ml.R` – DML estimation using the `causalweight` package
  - `rf_probit_comparison.R` – covariate balance comparison between the RF and probit-based propensity score models

All results and plots are stored in `src/results` for ease of access and reproducibility, given the computational intensity of the bootstrap procedure.
