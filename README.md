# GGSM-4
This repository contains MATLAB model code, parameterization details and results of the Generalized Global Sustainability Model (GGSM). <br>
New features: Water stress based price model and feedback model

## Main model

The main model is written in MATLAB.
This code has been tested to function smoothly when run with MATLAB Online R2022a (Tested on 2022 Jun 19).

 - Run a_run_this_file_scenario_popex_consinc.m for running the model. Select whether Scenario 4 to be activated from this file.
 - b_main_model_file.m contains the main model. Select the water price model related parameters from this file.


## Parameterization
 -  Parameters used for the water price model are showed in `Water price parameter matrix.xlsx` along with the detailed computation process.
 -  Input data used for the parameterization is also included in the same file. 
 -  Sources and processing used to obtain this input data are included `Parameterization-input-data` directory. 

## Results

Results include data files for Scenario 1, 4 and sensitivity of threshold price for agricultural elastic regime transition.
Scenarios 1 and 4 defined by Nisal et al (2022) [https://doi.org/10.1371/journal.pone.0267403].

   - `results-GGSM-4-Scenario-1.zip` contains results of Scenario 1: business as ususal
   - `results-GGSM-4-Scenario-4.zip` contains results of Scenario 4: population explosion with consumption increase
   - `results-GGSM-4-Sensitivity.zip` contains results of Scenario 1 with threshold values 700% and 2100% (scenario 1 default 1500%)

Generalized directory structure is as follows:
   - Scenario directories
     - zip of results
          Various result files with read me file which includes the scenario information: The .csv result files contain a header row. The .mat files can be read using MATLAB
