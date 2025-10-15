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

import Statistics as st 

using LaTeXStrings # latex strings package

include("utils.jl")
include("build_data.jl")

include("meas_equations.jl")
include("fase.jl")
include("meas_set_tests.jl")


function flat_start(x, network_info)
    for (key, value) in network_info["vm_ruler"]
        if !occursin("-c", key)
            x[value] = 1.0
        else
            x[value] = 1.0 + 0.01 # perturb voltage magnitudes
        end
    end
    for (key, value) in network_info["vdc_ruler"]
        if occursin("dc-1", key)
            x[value] = 1.0
        elseif occursin("dc-2", key)
            x[value] = -1.0
        else
            x[value] = 0.0
        end
    end
end


function determine_reference_buses(reference_buses, network_info)
    cols_to_remove = Vector{Int}()
    for bus in reference_buses
        push!(cols_to_remove, network_info["va_ruler"][string(bus)])
    end
    return cols_to_remove
end

function update_x!(x, dx, cols_to_remove)
    j = 1
    for i in eachindex(x)
        if i in cols_to_remove
            continue
        else
            x[i] = x[i] + dx[j]
            j += 1
        end
    end
end



nlp_optimizer = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 5, "sb" => "yes"
)



data_se = _ACDCSE.quickget_case5()
data_pf = _ACDCSE.quickget_case5() #loads network data




result, σ_dict, data_se = generate_data_basic_acdcse(data_pf, data_se, nlp_optimizer, sample_error = false);

data=data_se





meas_set=CSV.read("meas_set.csv",DataFrame; stringtype=String);

d_meas_set=create_dmeas_set(meas_set,data_pf);



filtermeas_set!(data_se,d_meas_set,sample_error=false,min_sig=1e-9) 



se_res = _ACDCSE.solve_acdcse(data_se, _PM.ACPPowerModel, nlp_optimizer) #Solve SE



network_info = build_network_info(data, result)

weight_z_autodiffse!(network_info, d_meas_set; min_sig=1e-7 )

h=create_h(data,network_info)

g=create_g(network_info)

W = create_W_autodiff(network_info, h, g)





h_vec(x) = [f(x) for f in vcat(h, g)]



#%%
x=zeros(length(network_info["x"]))
mode_rt = Enzyme.set_runtime_activity(Enzyme.Forward)
H = Enzyme.jacobian(mode_rt, h_vec, x)[1]
cols_to_remove = determine_reference_buses(reference_buses, network_info)
flat_start(x, network_info)
z1_val = [m.z for m in network_info["z"]]
z2_val = [0.0 for _ in network_info["g"]]
itmax = 10
z_val = vcat(z1_val, z2_val)

it=0
while it<itmax
    flat_start(x, network_info)
    it += 1
    H = Enzyme.jacobian(mode_rt, h_vec, x)[1]
    Hx = H[:, Not(cols_to_remove)]
    Gain = H' * W * H
    dz=z_val - h_vec(x)
    b=H'*W*dz
    dx= Gain\b
    update_x!(x, dx, cols_to_remove)    
    println("Norm of dx: ", norm(dx))
    if norm(dx) < 1e-6
        break
    end
end

 












