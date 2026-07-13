# PCLake Lake District Modelling

Repository for Lake District pathways modelling using PCLake+ as part of the PLURALAKES project.

Repository is organised as follows:

`data/` includes the raw data files - not on GH

`R/` contains helper functions for extracting data (met, sagis etc.), settling up files, and running the models

`scripts/` contains the workflow scripts:

-   `00_download_era5.py` a _python_ script to download and unpack the era-land dataset for the region. Get data from each grid point. 

-   `01_run_flake.R` runs the FLake model for all the lakes in the uklakes dataset to extract water temperature time series. FLake is run using two different light attenuation values to model a "turbid" and a "clear" lake. Then obtain an average of these.

-   `02_scaling_SAGIS.R` extracts, scales, and archives the raw SAGIS data (from data request) to a time series to be used as inputs in PCLake

-   `03_run_PCLake.R` runs the PCLake model for all the uklakes dataset

and the validation steps (a, b,..):

-   `a_validate_flake.R` script to validate the output of the FLake model runs using observational data.

-   `b_model_validation.R` script to validate the output of the PCLake model runs using observational data.

Other helper files: `restart_states.txt` list of the states needed to restart PCLake

## PCLake set-up:

This implementation of PCLake uses the DATM file saved in the main PCLake directory (see the set up in the main [PCLake repository](https://github.com/pcmodel/PCModel/) to show where this file sits). The file name should be ./PCModel/Licence_agreement/I_accept/PCModel1350/PCModel/3.00/Models/PCLake+/6.13.16/.

To run PCLake in R uses the scripts from ./PCModel/Licence_agreement/I_accept/PCModel1350/PCModel/3.00/Models/PCLake+/6.13.16/PCShell/scripts

## Setting up the project:

-   This directory only includes the optimisation project files. To actually run the workflow you need a lot more files and set up!

-   Start with cloning the PCLake repository or open a docker container that contains a fixed version of PCLake

-   Then clone this project repository into the *workcases* subdirectory for PCShell (R implementation) (./PCModel/Licence_agreement/I_accept/PCModel1350/PCModel/3.00/Models/PCLake+/6.13.16/PCShell/work_cases)

-   Then open the R project (PCLake_LD_CL.Rproj).

-   Once you've run some of the scripts it will generate the `model_code`, `source_cpp`, and `source_cpp_adjusted` directories which include the C++ files and the compiled executables.
