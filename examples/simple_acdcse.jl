#####################################################
#####################################################
#### A SIMPLE EXAMPLE OF HOW TO USE THE PACKAGE ####
#####################################################
#####################################################

import ACDCStateEstimation as _ACDCSE

import Ipopt
import PowerModels as _PM
import PowerModelsMCDC as _PMMCDC

nlp_optimizer = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 0, "sb" => "yes"
)

####################################################################
# PART 1: NOISELESS STATE ESTIMATION
#
# IF THE NETWORK MODEL USED BY THE STATE ESTIMATOR MATCHES
# THE GROUND-TRUTH (REAL NETWORK, OR POWER FLOW MODEL IF SYNTHETIC),
# THEN THERE IS NO NOISE FROM MODEL ERRORS. THE ONLY NOISE IS FROM 
# THE MEASUREMENTS. BUT HERE, WE SET THAT TO ZERO, TOO.
#
# SYNTHETIC MEASUREMENTS ARE BUILT BY RUNNING A POWER FLOW (OPF, 
# ACTUALLY) THAT WILL GIVE US VOLTAGES AND CURRENTS, TO WHICH 
# WE CAN ADD NOISE... OR NOT
# NOISE IS GAUSSIAN AND HAS A CERTAIN σ THAT WE NEED TO INDICATE.
# NOTE THAT WE CAN CREATE (NOISY) MEASUREMENTS FROM THE PF INPUTS,
# TOO (AS A GENERAL RULE, THE MORE THE MEASUREMENTS, 
# THE MERRIER THE SE) 
####################################################################

# get network data (note it is the same for PF and SE)
data_se = _ACDCSE.quickget_case5()
data_pf = _ACDCSE.quickget_case5()

# run power flow
result = _PMMCDC.solve_mcdcopf(data_pf, _PM.ACPPowerModel, nlp_optimizer)
_ACDCSE.get_dc_power!(result, data_pf) # -> the `result` dictionary only has dc currents and voltages, here we posprocess to get the powers

# below, the supposed variances for each measurement type, needed to build synthetic measurements
σ_dict = Dict{String, Float64}(
    "vm"       => 0.01, # (AC) voltage magnitude
    "va"       => 0.01, # (AC) voltage angle
    "pg"       => 0.01, # Active power injection from (AC) generator
    "qg"       => 0.01, # Reactive power injection from (AC) generator
    "pd"       => 0.01, # Active power injection from (AC) load
    "qd"       => 0.01, # Rective power injection from (AC) load
    "p_ac"     => 0.01, # Active power flow on an AC branch   (either injection or flow, really)
    "q_ac"     => 0.01, # Reactive power flow on an AC branch (either injection or flow, really)
    "vdcm"     => 0.01, # DC voltage magnitude
    "p_dc"     => 0.01, # DC power (either injection or flow)
    "i_dcgrid" => 0.01,  # DC current (either injection or flow)
    "cm"       => 0.01
)

_ACDCSE.add_current_flows_to_se_result!(result)

# let's build the synthetic measurements for the SE using the variances and the power flow results
# ⚠⚠⚠ IMPORTANT! ⚠⚠⚠ `sample_error` = false means NO noise is added to the PF inputs and outputs
_ACDCSE.powerflow2measurements!(data_se, result, σ_dict, sample_error = false, measurements = ["vm", "va", "p_to", "q_to", "vdcm", "pg", 
                                "pd", "qg", "qd", "p_dc_to", "p_dc_fr"])

# this function adds some settings to the state estimation, e.g., says if we minimize WLS or WLAV. 
# moreover, it remove PV buses from the data (not compatible with conventional SE).
# finally, it removes slack buses (except one), if there are more than one. Here, bus 1 is kept as slack
_ACDCSE.prepare_data_for_se_default!(data_se, exceptions = [1]) 

# and now, let's run state estimation.
# ⚠⚠⚠ because there is no noise, the output should be identical to the powerflow inputs/outputs,
# ⚠⚠⚠ and the objective value should be zero!
se_res = _ACDCSE.solve_acdcse(data_se, _PM.ACPPowerModel, nlp_optimizer)

####################################################################
# PART 2: STATE ESTIMATION WITH MEASUREMENT NOISE
# 
# VIRTUALLY AS THE ABOVE, BUT WE ADD NOISE BY SETTING
# `sample_error = true`
####################################################################

_ACDCSE.powerflow2measurements!(data_se, result, σ_dict, sample_error = true, measurements = ["vm", "va", "p_to", "q_to", "vdcm", "pg", 
                                "pd", "qg", "qd", "p_dc_to", "p_dc_fr"])

# you can appreciate how the objective is not zero anymore, and the estimated variable values are a bit different from the OPF results
# in general, higher σ -> higher noise -> higher difference between SE and PF results
se_res = _ACDCSE.solve_acdcse(data_se, _PM.ACPPowerModel, nlp_optimizer) 

####################################################################
# PART 3: STATE ESTIMATION WITH NETWORK MODEL ERRORS
#
# IF THE NETWORK MODEL USED BY THE STATE ESTIMATOR DOES *NOT* MATCH 
# THE GROUND-TRUTH, WE MIGHT BE IN TROUBLE. BUT THE ESTIMATOR CAN 
# HELP FIND OUT, E.G., WHICH LINE HAS BEEN DE-ENERGIZED 
# IN REAL-LIFE/SYNTHETIC DATA, WITHOUT INFORMING THE CONTROL ROOM, 
# I.E., THE ESTIMATION PEOPLE.
####################################################################

# let's just deenergize one of the branches...
_ACDCSE.deenergize_ac_branches!(data_pf, 11)

# ...and re-run the (optimal) power flow on the new data
result = _PMMCDC.solve_mcdcopf(data_pf, _PM.ACPPowerModel, nlp_optimizer)
_ACDCSE.get_dc_power!(result, data_pf) # -> the `result` dictionary only has dc currents and voltages, here we posprocess to get the powers

# now let's load synthetic measurements from the new power flow results
# you can play with the `sample_error` boolean to see how the errors/residuals behave 😃
_ACDCSE.powerflow2measurements!(data_se, result, σ_dict, sample_error = false, measurements = ["vm", "va", "p_to", "q_to", "vdcm", "pg", 
                                "pd", "qg", "qd", "p_dc_to", "p_dc_fr"])

# but let's run the state estimation on the old data (no de-energized branch)
# you can see that the residuals now are quite different, as the SE expect a different model than what we have
# so the estimates can be quite far from the measured values! 
se_res = _ACDCSE.solve_acdcse(data_se, _PM.ACPPowerModel, nlp_optimizer)

# btw, the state estimation result dictionary also does not report DC power,
# but you can get it from currents and voltages with the below
_ACDCSE.get_dc_power!(se_res, data_se)