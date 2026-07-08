# PCLake_LD_CL

Repository for Lake District pathways modelling using PCLake+ as part of the PLURALAKES project.

Repository is organised as follows:

`data/` includes the raw data files - not on GH

`R/` contains helper functions for extracting data (met, sagis etc.), settling up files, and running the models

`scripts/` contains the workflow scripts:

`scaling_SAGIS.R` extracts, scales, and archives the raw SAGIS data (from data request) to a time series to be used as inputs in PCLake

`run_flake.R` runs the FLake model for all the lakes in the uklakes dataset to extract water temperature time series. FLake is run using two different light attenuation values to model a "turbid" and a "clear" lake. Then obtain an average of these.

`validate_flake.R` script to validate the output of the FLake model runs using observational data.

`baseline_runs.R` runs the PCLake model for all the uklakes dataset

`model_validation.R` script to validate the output of the PCLake model runs using observational data.

Other helper files:
- `restart_states.txt` list of the states needed to restart PCLake
