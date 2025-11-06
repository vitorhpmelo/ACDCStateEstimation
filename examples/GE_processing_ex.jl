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





nlp_optimizer_pf = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-8, "print_level" => 5, "sb" => "yes"
)




ipot_file_wls="ipopt_wls.out"



nlp_optimizer_scada = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-5, "print_level" => 5, "sb" => "yes", "max_iter" => 10000
)


min_sigma = 1e-5
a=1e-3
n=2 #seed

#read input data
data_wls_sr = _ACDCSE.quickget_case5alt() #loads network data
data_pf = _ACDCSE.quickget_case5alt() #loads network data



reference=[1,6,11] #slack buses


#run pf and generate measurement
set_fixed_bus_voltages!(data_pf)
result, σ_dict, data_wls_sr = generate_data_basic_acdcse(data_pf, data_wls_sr, nlp_optimizer_pf,"all",reference, sample_error = false);




meas_set=CSV.read(joinpath(_ACDCSE.ACDCSE_dir(), "test/data/meas_set/meas_set_case5.csv"),DataFrame; stringtype=String);
d_meas_set=create_dmeas_set(meas_set,data_pf)


#filter measurement set
d_keys,d_prec=filtermeas_set!(data_wls_sr,d_meas_set,sample_error=false,min_sig=min_sigma)




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
data_scada = deepcopy(data_hyb)



split_meas_set!(d_keys, data_scada, ["scada", "dc", "conv"])


#set a very low error for the DC current measurement to avoid numerical issues
se_res_scada = _ACDCSE.solve_acdcse(data_scada, _PM.ACPPowerModel, nlp_optimizer_scada) #Solve SE


# net_posterior = generate_data_posterior_cov!(se_res_pmu, data_pmu,d_keys)

Ω_1,r_1,W_1,G1,network_info_1 = generate_res_cov_QR!(se_res_scada, data_scada, d_keys,min_sig=1e-5,print_info=true)


# Ωqr,rqr,Wqr= generate_res_cov_QR!(se_res_scada, data_scada, d_keys;min_sig=1e-4,print_info=true)





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
# find second highest distinct value and its index
if length(rn_1) < 2
    println("Array too short to find a second highest value.")
else
    sorted_idx = sortperm(rn_1, rev=true)  # indices sorted by value desc
    first_idx = sorted_idx[1]
    first_val = rn_1[first_idx]

    second_idx = nothing
    second_val = nothing
    for idx in Iterators.drop(sorted_idx, 1)
        if rn_1[idx] != first_val
            second_idx = idx
            second_val = rn_1[idx]
            break
        end
    end

    if second_idx === nothing
        println("No second distinct maximum found; all values equal to ", first_val)
    else
        println("Second highest value of rn_1: ", second_val, " at index ", second_idx)
    end
end

# plot rn_1 and mark the maximum
inds = collect(1:length(rn_1))
p = plt.plot(inds, rn_1, label="rn_1", lw=2)
plt.scatter!(p, [idx_rn_1], [max_rn_1], color=:red, label="max", ms=6)
plt.xlabel!(p, "Measurement index")
plt.ylabel!(p, "Normalized residual (rn_1)")
plt.title!(p, "Normalized residuals rn_1")
display(p)


d_1 = abs.(diag(Ω_1))
d_1 = filter(!=(0.0), d_1)
log10d_1 = log10.(d_1)





plt.plot(log10d_1, label="log10(|diag(Ω_1)|)")
plt.xlabel!("Index")
plt.ylabel!("Value")
plt.title!("Absolute Values of diag(Ω)")




xxmax_val, max_idx = findmax(rn)


sqrt(W[50,50])