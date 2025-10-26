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


using LaTeXStrings # latex strings package

include("build_data.jl")

include("meas_equations.jl")

include("utils.jl")

include("meas_set_tests.jl")

include("fase.jl")
include("data_analysis.jl")
include("simulation_functions.jl")




JuMP._CONSTRAINT_LIMIT_FOR_PRINTING[] = 1000000


nlp_optimizer_pf = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-8, "print_level" => 5, "sb" => "yes"
)



ipot_file_fase="ipopt_fase.out"
ipot_file_wls="ipopt_wls.out"
ipot_file_scada="ipopt_scada.out"

nlp_optimizer_pmu = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-5, "print_level" => 5,  "sb" => "yes", "max_iter" => 10000, "output_file" => ipot_file_fase
)


nlp_optimizer_hyb = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-5, "print_level" => 5, "sb" => "yes", "max_iter" => 10000, "output_file" => ipot_file_wls
)

nlp_optimizer_scada = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-5, "print_level" => 5, "sb" => "yes", "max_iter" => 10000, "output_file" => ipot_file_scada
)


min_sigma = 1e-5

name_experiment="stationary_MC5_wlav"

data_wls_sr = _ACDCSE.quickget_case5_paper() #loads network data
data_pf = _ACDCSE.quickget_case5_paper() #loads network data




reference=[1,6,11] #slack buses

set_fixed_bus_voltages!(data_pf)
result, σ_dict, data_wls_sr = generate_data_basic_acdcse(data_pf, data_wls_sr, nlp_optimizer_pf,"all", sample_error = false); # solves powerflow and generates7 SE data

N=100
n0=1

meas_set=CSV.read(joinpath(_ACDCSE.ACDCSE_dir(), "test/data/meas_set/meas_set_case5.csv"),DataFrame; stringtype=String);

d_meas_set=create_dmeas_set(meas_set,data_pf)

d_keys,d_prec=filtermeas_set!(data_wls_sr,d_meas_set,sample_error=false,min_sig=min_sigma)


ipopt_out=Dict("pmu"=>ipot_file_fase,"hyb"=>ipot_file_wls,"scada"=>ipot_file_scada)
se_objective=Dict("scada"=>"rwlav","pmu"=>"rwlav" ,"hyb"=>"rwlav")



dfs_errors = DataFrame[]
dfs_conv = DataFrame[]
dfs_priors = DataFrame[]
set_converter_constraints!(data_wls_sr)

for a in [1e-4]
    df_errors, df_conv,df_prior = run_fase_experiment_stationary_timebench(N, n0, a, min_sigma, result, data_wls_sr, d_prec, d_keys, d_meas_set, nlp_optimizer_hyb, nlp_optimizer_scada, nlp_optimizer_pmu, ipopt_out;cutoff=1e0,se_objective=se_objective)
    df_errors[!, :a] .= a
    df_conv[!, :a] .= a
    push!(dfs_errors, df_errors)
    push!(dfs_conv, df_conv)
    push!(dfs_priors, df_prior) 
end


df_errors_all =  vcat(dfs_errors...)
df_conv_all =  vcat(dfs_conv...)
df_priors_all = vcat(dfs_priors...)

CSV.write("$(name_experiment)_priors.csv", df_priors_all)
CSV.write("$(name_experiment)_errors.csv", df_errors_all)
CSV.write("$(name_experiment)_conv.csv", df_conv_all)