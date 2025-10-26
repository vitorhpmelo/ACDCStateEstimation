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
    Ipopt.Optimizer, "tol" => 1e-5, "print_level" => 5, "sb" => "yes", "max_iter" => 20000, "output_file" => ipot_file_wls
)

nlp_optimizer_scada = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-5, "print_level" => 5, "sb" => "yes", "max_iter" => 10000, "output_file" => ipot_file_scada
)


min_sigma = 1e-5
a=1e-3
n=1 #seed

#read input data
data_wls_sr = _ACDCSE.quickget_case118_paper() #loads network data
data_pf = _ACDCSE.quickget_case118_paper() #loads network data


ipopt_out=Dict("pmu"=>ipot_file_fase,"hyb"=>ipot_file_wls,"scada"=>ipot_file_scada)

reference=[69,119,120,121,122,123] #slack buses

#%%
#run pf and generate measurements
result, σ_dict, data_wls_sr = generate_data_basic_acdcse(data_pf, data_wls_sr, nlp_optimizer_pf,"all",reference, sample_error = false);

N=1
d_conv_scada = Dict(
    "termination_status" => Vector{Any}(undef, N),
    "n_iters" => Vector{Any}(undef, N),
)

d_conv_pmu = Dict(
    "termination_status" => Vector{Any}(undef, N),
    "n_iters" => Vector{Any}(undef, N),
)

d_conv_hyb = Dict(
    "termination_status" => Vector{Any}(undef, N),
    "n_iters" => Vector{Any}(undef, N),
)

df_errors_total = DataFrame()




meas_set=CSV.read("meas_set_case118.csv",DataFrame; stringtype=String);
d_meas_set=create_dmeas_set(meas_set,data_pf)


#filter measurement set
d_keys,d_prec=filtermeas_set!(data_wls_sr,d_meas_set,sample_error=false,min_sig=min_sigma)
ipopt_out=Dict("pmu"=>ipot_file_fase,"hyb"=>ipot_file_wls,"scada"=>ipot_file_scada)
se_objective=Dict("scada"=>"rwls","pmu"=>"rwls" ,"hyb"=>"rwls")

scada_measurements = ["scada", "dc", "conv"]
d_meas_set_scada = Dict()
for key in scada_measurements
    d_meas_set_scada[key] = deepcopy(d_meas_set[key])
end

#generate data for fase, wls and scada
data_hyb = deepcopy(data_wls_sr) #copy data to avoid modifying original data
introduce_noise!(data_hyb, d_prec, d_keys, seed=n, min_sig=min_sigma, bound=true)
data_pmu = deepcopy(data_hyb)
data_scada = deepcopy(data_hyb)

split_meas_set!(d_prec, d_keys, data_scada, ["scada", "dc", "conv"])
split_meas_set!(d_prec, d_keys, data_pmu, ["pmu", "dtu"]) #splits the measurements into SCADA and PMU sets;

data_scada["se_settings"]["criterion"] = se_objective["scada"]
se_res_scada = _ACDCSE.solve_acdcse(data_scada, _PM.ACPPowerModel, nlp_optimizer_scada) #Solve SE
n_iters = add_iterations_from_log!(ipopt_out["scada"])
df_error_scada = df_error_var(result, se_res_scada)
df_error_scada[!, :se] .= "scada"
df_error_scada[!, :n] .= n


d_conv_scada["termination_status"][n] = string(se_res_scada["termination_status"])
d_conv_scada["n_iters"][n] = n_iters


data_hyb["se_settings"]["criterion"] = se_objective["hyb"]
se_res_hyb = _ACDCSE.solve_acdcse(data_hyb, _PM.ACPPowerModel, nlp_optimizer_hyb) #Solve WLS
n_iters = add_iterations_from_log!(ipopt_out["hyb"])
d_conv_hyb["n_iters"][n] = n_iters
df_error_hyb = df_error_var(result, se_res_hyb)
df_error_hyb[!, :se] .= "hyb"
df_error_hyb[!, :n] .= n
d_conv_hyb["termination_status"][n] = string(se_res_hyb["termination_status"])


CSV.write("meas_hyb.csv",meas_dict_to_dataframe(data_hyb["meas"]))
CSV.write("meas_pmu.csv",meas_dict_to_dataframe(data_pmu["meas"]))
CSV.write("meas_scada.csv",meas_dict_to_dataframe(data_scada["meas"]))


bus_PQ_lst=[]
for (key, meas) in data_hyb["meas"]
    if meas["cmp"] == :load
        push!(bus_PQ_lst, data_hyb["load"]["$(meas["cmp_id"])"]["load_bus"])
    elseif meas["cmp"] == :gen
        push!(bus_PQ_lst, data_hyb["gen"]["$(meas["cmp_id"])"]["gen_bus"])
    end
end
CSV.write("bus_PQ_lst.csv", DataFrame(bus_PQ_lst = bus_PQ_lst))
#%%
net = generate_data_prior!(se_res_scada, data_scada, data_pmu, d_meas_set_scada, d=a, min_sig=min_sigma, print_info=true) #Generate prior data for FASE
        


data_pmu["se_settings"]["criterion"] = se_objective["pmu"]
se_res_pmu = _ACDCSE.solve_acdcse(data_pmu, _PM.ACPPowerModel, nlp_optimizer_pmu) #Solve FASE

n_iters = add_iterations_from_log!(ipopt_out["pmu"])

df_error_pmu = df_error_var(result, se_res_pmu)
df_error_pmu[!, :se] .= "pmu"
df_error_pmu[!, :n] .= n

df_errors_total = vcat(df_errors_total, df_error_pmu, df_error_hyb, df_error_scada)
        
d_conv_pmu["termination_status"][n-n0] = string(se_res_pmu["termination_status"])
d_conv_pmu["n_iters"][n-n0] = n_iters



