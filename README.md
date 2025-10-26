# Repository containing paper results
## paper "Bayesian Information Fusion for State Estimation in AC/DC Networks with Unbalanced HVDC Grids" subbmitted to the PSCC 2026 (under review)
[![Build Status](https://github.com/MartaVanin/ACDCStateEstimation.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/MartaVanin/ACDCStateEstimation.jl/actions/workflows/CI.yml?query=branch%3Amaster)



# Instructions:

## Network data and measurement sets 

- The network data are avaible in the  sub-directory test/data/matacdc_scripts
- The measurement sets employed for both test systems used in the paper are avaible in the folder test/data/meas_set. The folder also contains a file explaining the measurement sets files

## Paper results 

### Section A
The scripts used to generate of Section A can be reproduced by executing the following scripts  

- "examples/stationary_MC5_wls.jl"
- "examples/stationary_MC5_wlav.jl"

Both scripts will execute 100 Monte Carlo simulations, in the "case5" test case, the first one emplying the WLS objective, and the second one the WLAV, for all formulations employed in the paper (Prior, Bay and Hyb). 
The excution of both scripts will result in the files :

- "stationary_MC5_wls_errors.csv"
- "stationary_MC5_wlav_errors.csv" 

These files contains the reference values for all state variables and the estimated values for these quantities in all simulations, by using these values it is possible to reproduce TAB. I, and Fig.4. 

Note: To generate Tab I was generated without considering DC SMU, these masurements need to be removed in  "test/data/meas_set/meas_set_case5.csv". 

### Section B

The script used to generate the results in section B is "examples/time_series_simul_case118.jl", this script will execute a time series for the methods tested in the paper (Prior, Bay and Hyb). Note that the time series is repeated 100 times in the Monte Carlo simulation. 

Note: This execution can take time!!! 

The measurement residuals objective is defined in the line 67 of the script 

-obj="rwls" executes the Weighted Least Squares
-obj="rwlav" executes the Weighted Least Absolute Values

The excution of both scripts will result in the files:

df_errors_total_time_series_simul_case118_{obj}.csv

With the results in this file, the Fig. 5 can be generated.

### Section C

The scripts used to generate the results in this section are:

- "examples/stationary_MC5_wlav_bellow_bound.jl"
- "examples/stationary_MC5_wlav_bellow_bound.jl"

The script is similar to the ones used in Section A, to obtain the results in Fig. 6. The value of the bound used to the Pmax of the Conv 3 is altered in line 64 of the script. The scripts have to be runned for each different bound value.

Running the scripts will generate the following files:  

- "stationary_MC5_wls_errors.csv"
- "stationary_MC5_wlav_errors.csv" 

via the result in these files it is possible to generate Fig. 6.

### Section D

The scripts used to generate the results in this section are:

- "examples/stationary_MC5_wlav_time_bench.jl"
- "examples/stationary_MC5_wls_time_bench.jl"
- "examples/stationary_MC118_wlav_time_bench.jl"
- "examples/stationary_MC118_wls_time_bench.jl"

Each script will run a monte carlo simulation and the resulting execution times will be avaible in the file:

- "stationary_MC5_wls_conv.csv"
- "stationary_MC5_wlav_conv.csv"
- "stationary_MC118_wls_conv.csv"
- "stationary_MC118_wlav_conv.csv"
With the execution times, and convergence information. The content of these files can be used to generate Table IV. Note that this table depends on the computer used. 

