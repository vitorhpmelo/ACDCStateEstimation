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

include("simulation_functions.jl")

include("meas_equations.jl")

include("utils.jl")

include("meas_set_tests.jl")

include("fase.jl")
include("data_analysis.jl")







JuMP._CONSTRAINT_LIMIT_FOR_PRINTING[] = 1000000



nlp_optimizer = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-8, "print_level" => 5, "sb" => "yes"
)

ipot_file_fase="ipopt_fase.out"
ipot_file_wls="ipopt_wls.out"
ipot_file_scada="ipopt_scada.out"

nlp_optimizer_fase = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-5, "print_level" => 5, "sb" => "yes", "max_iter" => 10000, "output_file" => ipot_file_fase
)


nlp_optimizer_wls = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-5, "print_level" => 5, "sb" => "yes", "max_iter" => 10000, "output_file" => ipot_file_wls
)


nlp_optimizer_scada = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-5, "print_level" => 5, "sb" => "yes", "max_iter" => 10000, "output_file" => ipot_file_scada
)

min_sigma = 1e-5

a=1e-4
obj="rwls" #defines the objective function used in the simulation (rwls = Weighted Least Squares , rwlav Weighted Least Absolute Values)
name_simulation = "time_series_simul_case118_$obj"

data_wls_sr = _ACDCSE.quickget_case118_paper() #loads network data
data_pf = _ACDCSE.quickget_case118_paper() #loads network data
reference=[69,119,120,121,122,123] #slack buses


print("Generating data...\n")




T=50
results=Vector{Dict}(undef, T) # to store results with size T
data_pf_all=Vector{Dict}(undef, T) # to store all data with size T
data_se_all=Vector{Dict}(undef, T) # to store all data with size T
data_scada_all=Vector{Dict}(undef, T) # to store all data with size T
data_pmu_all=Vector{Dict}(undef, T) # to store all data with size T
data_wls_all=Vector{Dict}(undef, T) # to store all data with size T


# create data
for t in 1:T
    data_pf_all[t]=deepcopy(data_pf)
end


modify_all_loads_randomly!(data_pf_all,118,2,0.5) # modify all loads randomly over time

# Draw a set of random integers in a predefined range

# Example usage:
buses_2_mod,instant_load_mod,percent_mod = modify_loads_sudden_change!(data_pf_all, 5, 118, 110, T)



# Print buses_2_mod as a LaTeX table
lineone = join(["$bus" for bus in buses_2_mod], " & ")
line_two = join(["$(instant_load/10)" for instant_load in instant_load_mod], " & ")
line_three = join([@sprintf("%.2f", percent) for percent in percent_mod], " & ")

# Sort lineone (buses_2_mod) in ascending order and reorder line_two and line_three accordingly
sorted_indices = sortperm(buses_2_mod)
lineone_sorted = join(["$(buses_2_mod[i])" for i in sorted_indices], " & ")
line_two_sorted = join(["$(instant_load_mod[i]/10)" for i in sorted_indices], " & ")
line_three_sorted = join([@sprintf("%.2f", percent_mod[i]) for i in sorted_indices], " & ")

print(join([lineone_sorted, line_two_sorted, line_three_sorted], " \\\\ \n"))

bus2mod=[6]

modify_bus_v_over_time_ramp!(data_pf_all,25,3, bus2mod,2.0)


bus2mod=[8]

modify_bus_v_over_time_ramp!(data_pf_all,45,3, bus2mod,2.0)





wind_gen=[56,57,58] #

create_se_data!(data_pf_all, data_se_all, results, nlp_optimizer, reference; wind_gen=wind_gen,wind_gen_var=1.0); # create se data for all time steps

meas_set=CSV.read(joinpath(_ACDCSE.ACDCSE_dir(), "test/data/meas_set/meas_set_case118.csv"),DataFrame; stringtype=String);




N=100
n0=1

sampling_scada=2000 #in ms
sampling_pmu=100 #in ms


dfs_meas_all_out=Vector(undef, N)
dfs_meas_hyb_out=Vector(undef, N)
dfs_meas_pmu_out=Vector(undef, N)
dfs_meas_scada_out=Vector(undef, N)



df_errors_total=DataFrame()
df_conv_total=DataFrame()
df_conv_prior_total=DataFrame()

for data in data_se_all
    set_converter_constraints!(data)
end


for n in 1:N
    print("Running simulation for n=$n...\n")
    df_errors, df_conv, df_conv_prior = run_se_simulation_map!(
        T, n, data_se_all, data_scada_all, data_pmu_all, data_wls_all, results,
        sampling_scada, sampling_pmu,
        meas_set, min_sigma, a, nlp_optimizer_scada, nlp_optimizer_wls, nlp_optimizer_fase,
        ipot_file_fase, ipot_file_wls; se_objective=Dict("scada"=>"$obj","pmu"=>"$obj" ,"hyb"=>"$obj"), with_noise=true, noise_bound=true, noise_bound_factor=2.5
    )
    global df_errors_total, df_conv_total, df_conv_prior_total
    df_errors_total = vcat(df_errors_total, df_errors)
    df_conv_total = vcat(df_conv_total, df_conv)
    df_conv_prior_total = vcat(df_conv_prior_total, df_conv_prior)

end




CSV.write(
    "d_conv$name_simulation.csv",   
    df_conv_total
)



CSV.write(
    "df_errors_total_$name_simulation.csv",
    df_errors_total
)