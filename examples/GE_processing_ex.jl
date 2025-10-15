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
    Ipopt.Optimizer, "tol" => 1e-5, "print_level" => 5, "sb" => "yes", "max_iter" => 20000, "output_file" => ipot_file_wls
)

nlp_optimizer_scada = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-5, "print_level" => 5, "sb" => "yes", "max_iter" => 10000, "output_file" => ipot_file_scada
)


min_sigma = 1e-5
a=1e-3
n=2 #seed

#read input data
data_wls_sr = _ACDCSE.quickget_case5alt() #loads network data
data_pf = _ACDCSE.quickget_case5alt() #loads network data


ipopt_out=Dict("pmu"=>ipot_file_fase,"hyb"=>ipot_file_wls,"scada"=>ipot_file_scada)

reference=[1,6,11] #slack buses


#run pf and generate measurement
set_fixed_bus_voltages!(data_pf)
result, σ_dict, data_wls_sr = generate_data_basic_acdcse(data_pf, data_wls_sr, nlp_optimizer_pf,"all",reference, sample_error = false);



d_conv_scada = Dict(
    "termination_status" => "any string"
)

d_conv_pmu = Dict(
    "termination_status" => "any string"
)

d_conv_hyb = Dict(
    "termination_status" => "any string"
)

df_errors_total = DataFrame()



meas_set=CSV.read("meas_set_case5.csv",DataFrame; stringtype=String);
d_meas_set=create_dmeas_set(meas_set,data_pf)


#filter measurement set
d_keys,d_prec=filtermeas_set!(data_wls_sr,d_meas_set,sample_error=false,min_sig=min_sigma)

ipopt_out=Dict("pmu"=>ipot_file_fase,"hyb"=>ipot_file_wls,"scada"=>ipot_file_scada)
se_objective=Dict("scada"=>"rwls","pmu"=>"rwls" ,"hyb"=>"rwls")

scada_measurements = ["scada", "dc", "conv"]
d_meas_set_scada = Dict()
for key in keys(d_meas_set)
    if key in scada_measurements
        d_meas_set_scada[key] = deepcopy(d_meas_set[key])
    end
end

d_meas_set_pmu = Dict()
pmu_measurements = ["pmu", "dtu"]
for key in pmu_measurements
    d_meas_set_pmu[key] = deepcopy(d_meas_set[key])
end


#generate data for fase, wls and scada
data_hyb = deepcopy(data_wls_sr) #copy data to avoid modifying original data


introduce_noise!(data_hyb, d_prec, d_keys, seed=n, min_sig=min_sigma, bound=true, bound_factor=2) #introduce noise to measurements
data_pmu = deepcopy(data_hyb)
data_scada = deepcopy(data_hyb)



split_meas_set!(d_keys, data_scada, ["scada", "dc", "conv"])
split_meas_set!(d_keys, data_pmu, ["pmu", "dtu"]) #splits the measurements into SCADA and PMU sets;

data_scada["se_settings"]["criterion"] = se_objective["scada"]

set_converter_constraints!(data_scada, Pmax_factor=2.0, Pmin_factor=2.0, Vmax_factor=1.3, Vmin_factor=0.7)
μ_original = data_scada["meas"]["192"]["dst"][2].μ
σ_original = data_scada["meas"]["192"]["dst"][2].σ
μ_eg=μ_original+10*σ_original
σ_eg=abs(μ_eg)*0.02/3
data_scada["meas"]["192"]["dst"][2]=_ACDCSE._DST.Normal(μ_eg, σ_eg)


#set a very low error for the DC current measurement to avoid numerical issues
se_res_scada = _ACDCSE.solve_acdcse(data_scada, _PM.ACPPowerModel, nlp_optimizer_scada) #Solve SE
df_error_scada = df_error_var(result, se_res_scada)
df_error_scada[!, :se] .= "scada"
df_error_scada[!, :n] .= n

#%%



#%%


# net_posterior = generate_data_posterior_cov!(se_res_pmu, data_pmu,d_keys)

Ω_1,r_1,W_1,G1,network_info_1 = generate_res_cov_new!(se_res_scada, data_scada, d_keys,min_sig=1e-5,print_info=true)


Ω_2,r_2,W_2,G2,network_info_2 = generate_res_cov_new!(se_res_scada, data_scada, d_keys,min_sig=1e-8,print_info=true)
# Ωqr,rqr,Wqr= generate_res_cov_QR!(se_res_scada, data_scada, d_keys;min_sig=1e-4,print_info=true)



rn_2= similar(r_2)
for (i,r_i) in enumerate(r_2)
    if Ω_2[i,i] < 1e-10
        rn_2[i]=0.0
    else
        rn_2[i]=r_i/sqrt(Ω_2[i,i])
    end
end


rn_1= similar(r_1)
k=0
for (i,r_i) in enumerate(r_1)
    if Ω_1[i,i] < 1e-10
        rn_1[i]=0.0
    else
        rn_1[i]=r_i/sqrt(Ω_1[i,i])

       global k
       k=k+1
    end
end

max_rn_1, idx_rn_1 = findmax(rn_1)
println("Maximum value of rn_1: ", max_rn_1, " at index ", idx_rn_1)

max_rn_2, idx_rn_2 = findmax(rn_2)
println("Maximum value of rn_2: ", max_rn_2, " at index ", idx_rn_2)



d_1 = abs.(diag(Ω_1))
d_1 = filter(!=(0.0), d_1)
log10d_1 = log10.(d_1)

d_2 = abs.(diag(Ω_2))
d_2 = filter(!=(0.0), d_2)
log10d_2 = log10.(d_2)

plt.plot(log10d_1, label="log10(|diag(Ω_1)|)")
plt.plot!(log10d_2, label="log10(|diag(Ω_2)|)")
plt.xlabel!("Index")
plt.ylabel!("Value")
plt.title!("Absolute Values of diag(Ω)")




max_val, max_idx = findmax(rn)


sqrt(W[50,50])