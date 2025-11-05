#%%
import ACDCStateEstimation as _ACDCSE #ACDC State estimation package

import Ipopt # Interior point solver
using JuMP # JuMP package for mathematical optimization
import PowerModels as _PM # ACpower models package from ELECTA
import PowerModelsMCDC as _PMMCDC # MCDC power models from ELECTA
import Plots as plt# plotting package
using CSV
using DataFrames
using Printf
using LinearAlgebra
import Statistics as st 
using Enzyme
include("build_data.jl")

include("meas_equations.jl")

include("utils.jl")

include("meas_set_tests.jl")

include("fase.jl")
include("data_analysis.jl")
include("simulation_functions.jl")


nlp_optimizer_pf = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-8, "print_level" => 5, "sb" => "yes"
)



data_pf= _ACDCSE.quickget_cigre_B4()
data_se= _ACDCSE.quickget_cigre_B4()

reference=[1]
#%%

set_fixed_bus_voltages!(data_pf)
set_fixed_gen_pg!(data_pf)

result, σ_dict, data_wls_sr = generate_data_basic_acdcse(data_pf, data_se, nlp_optimizer_pf,"no_branch",reference, sample_error = false);
