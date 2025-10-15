#%%
import ACDCStateEstimation as _ACDCSE #ACDC State estimation package

import Ipopt # Interior point solver
import PowerModels as _PM # ACpower models package from ELECTA
import PowerModelsMCDC as _PMMCDC # MCDC power models from ELECTA
import Plots # plotting package
using CSV
using DataFrames
using LinearAlgebra
using Enzyme
using BenchmarkTools
import Statistics as st 

using LaTeXStrings # latex strings package

include("utils.jl")
include("build_data.jl")

include("meas_equations.jl")
include("fase.jl")
include("meas_set_tests.jl")
include("data_analysis.jl")
include("SE_autodiff.jl")


nlp_optimizer = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 5, "mu_strategy" => "adaptive", "sb" => "yes"
)



data_se = _ACDCSE.quickget_case118_paper()
data_pf = _ACDCSE.quickget_case118_paper() #loads network data



result, σ_dict, data_se = generate_data_basic_acdcse(data_pf, data_se, nlp_optimizer, sample_error = false);

#%%    



meas_set=CSV.read(joinpath(_ACDCSE.ACDCSE_dir(), "test/data/meas_set/meas_set_case118.csv"),DataFrame; stringtype=String);

d_meas_set=create_dmeas_set(meas_set,data_pf);



filtermeas_set!(data_se,d_meas_set,sample_error=false,min_sig=1e-4) 

data=data_se
reference_buses=[69,119,120,121,122,123]
for bus in reference_buses
    data["bus"]["$bus"]["bus_type"]=3
end

data_se["bus"]
se_res = _ACDCSE.solve_acdcse(data, _PM.ACPPowerModel, nlp_optimizer) #Solve SE





network_info = build_network_info(data, result)

network_info["reference_buses"] = [69,119,120,121,122,123]


 [213, 146, 148, 149, 150, 151]

x,dx_norms=SE_WLS_timebench(network_info, d_meas_set; min_sig=1e-7, itmax=10)

for (i, x_i) in enumerate(network_info["x"])
    
    println("x[$i] = $x_i", "dx =", x_i-x[i])

end

x_info = Vector{String}()

for x in network_info["x_info"]
    push!(x_info, x.label)
end


 df_x = DataFrame(
    x_ref = network_info["x"],
    x_est = x,
    x_label = x_info,
    dx = abs.(x - network_info["x"])
)
println(df_x)
CSV.write("df_x_case118.csv", df_x)

println(norm(x-network_info["x"]), " norm of difference between x and x_ref")


Plots.plot(dx_norms, xlabel="Iteration", ylabel="‖dx‖",legend=false, yscale=:log10)

z_val = [m.z for m in network_info["z"]]
g_val = [0.0 for _ in network_info["g"]]
z_label = Vector{String}()

h_val=vcat(z_val, g_val)
hx=network_info["h_vec"](x)
maximum(abs.(h_val-hx))
#calcular com x verdadeiro e x convergido


df=df_error_var(result, se_res)
df_sorted = sort(df, :error, rev=true)
CSV.write("sorted_errors_case118.csv", df_sorted)
# Filter df to select certain variables
