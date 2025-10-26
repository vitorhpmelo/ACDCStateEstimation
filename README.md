# Repository containing the data used in the paper "Bayesian Information Fusion for State Estimation in AC/DC Networks with Unbalanced HVDC Grids" subbmitted to the PSCC 2026
[![Build Status](https://github.com/MartaVanin/ACDCStateEstimation.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/MartaVanin/ACDCStateEstimation.jl/actions/workflows/CI.yml?query=branch%3Amaster)



## Instructions:

### Network data and measurement sets 

- The network data are avaible in the  sub-directory test/data/matacdc_scripts
- The measurement sets employed for both test systems used in the paper are avaible in the folder test/data/meas_set

### Paper results 

The scripts used to generate of Section A can be reproduced by executing the following scripts  
    [1] examples/stationary_MC5_wls.jl
    [2] examples/stationary_MC5_wlav.jl
Both scripts will execute Monte Carlo simulations, in the "case5" test case, the first one emplying the WLS objective, and the second one the WLAV, for all formulations employed in the paper (Prior, Bay and Hyb). 
The excution of both scripts will result in two files :"stationary_MC5_wls_errors.csv" and "stationary_MC5_wlav_errors.csv". These files contains the reference values for all state variables and the estimated values for these quantities in all simulations, by using these values it is possible to reproduce TAB. I, and Fig.6. 
Note: To generate Tab I was generated without considering DC SMU, these masurements need to be removed in  test/data/meas_set/meas_set_case5.csv. 


