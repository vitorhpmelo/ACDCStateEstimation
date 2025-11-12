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
include("forecast_functions.jl")


nlp_optimizer_pf = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-8, "print_level" => 5, "sb" => "yes"
)




function populate_meas_from_forecast!(data::Dict{String,Any}, res_t::Dict{String,Any}, load_and_gen_data_fore_t; prec_meas::Float64=0.05, prec_fore::Float64=0.05, prec_virtual::Float64=1e-5)
    data["meas"] = Dict{String,Any}() 
    create_vm!(data, res_t; prec=prec_meas)
    create_va!(data, res_t; prec=prec_meas)
    create_vmf!(data, res_t; prec=prec_meas)
    create_vaf!(data, res_t; prec=prec_meas)
    create_vac!(data, res_t; prec=prec_meas)
    create_vmc!(data, res_t; prec=prec_meas)
    create_vmdc!(data, res_t; prec=prec_meas)

    for row in eachrow(load_and_gen_data_fore_t)
        if row[:type] == 0
            n = length(data["meas"]) + 1
            data["meas"][string(n)] = Dict{String,Any}()
            load = determine_loads_to_modify([row[:bus]], data)[1]
            m = data["meas"][string(n)]
            m["crit"]   = "rwls"
            m["cmp_id"] = parse(Int64, load)
            m["cmp"]    = :load
            m["var"]    = :pd
            μ = row[:var] * data["load"][string(load)]["pd"]
            σ = prec_fore * μ / 3
            if σ < prec_virtual
                σ = prec_virtual
            end
            m["dst"] = [_ACDCSE._DST.Normal(μ, σ)]

            gen = determine_gens_to_modify([row[:bus]], data)
            if length(gen) == 0
                continue
            end
            gen = gen[1]
            n = length(data["meas"]) + 1
            data["meas"][string(n)] = Dict{String,Any}()
            m = data["meas"][string(n)]
            m["crit"]   = "rwls"
            m["cmp_id"] = parse(Int64, gen)
            m["cmp"]    = :gen
            m["var"]    = :pg
            μ = res_t["solution"]["gen"][string(gen)]["pg"]
            σ = prec_fore * μ / 3
            if σ < prec_virtual
                σ = prec_virtual
            end
            m["dst"] = [_ACDCSE._DST.Normal(μ, σ)]

        elseif row[:type] == 1
            gen = determine_gens_to_modify([row[:bus]], data)[1]
            n = length(data["meas"]) + 1
            data["meas"][string(n)] = Dict{String,Any}()
            m = data["meas"][string(n)]
            m["crit"]   = "rwls"
            m["cmp_id"] = parse(Int64, gen)
            m["cmp"]    = :gen
            m["var"]    = :pg
            μ = row[:var] * data["gen"][string(gen)]["pg"]
            σ = prec_fore * μ / 3
            if σ < prec_virtual
                σ = prec_virtual
            end
            m["dst"] = [_ACDCSE._DST.Normal(μ, σ)]

            load = determine_loads_to_modify([row[:bus]], data)
            if length(load) == 0
                continue
            end
            load = load[1]
            n = length(data["meas"]) + 1
            data["meas"][string(n)] = Dict{String,Any}()
            m = data["meas"][string(n)]
            m["crit"]   = "rwls"
            m["cmp_id"] = parse(Int64, load)
            m["cmp"]    = :load
            m["var"]    = :pd
            μ = res_t["solution"]["load"][string(load)]["pd"]
            σ = prec_fore * μ / 3
            if σ < prec_virtual
                σ = prec_virtual
            end
            m["dst"] = [_ACDCSE._DST.Normal(μ, σ)]
        end
    end

    return data
end



data_pf= _ACDCSE.quickget_cigre_B4()
data_se= _ACDCSE.quickget_cigre_B4()

reference=[1,2,3,4]

load_and_gen_data = CSV.read(joinpath(_ACDCSE.ACDCSE_dir(),"test/data/load_and_gen/CIGRE_B4_measured.csv"),DataFrame; stringtype=String);
load_and_gen_data_fore_mean = CSV.read(joinpath(_ACDCSE.ACDCSE_dir(),"test/data/load_and_gen/CIGRE_B4_forecasted.csv"),DataFrame; stringtype=String);
load_and_gen_data_fore_P10 = CSV.read(joinpath(_ACDCSE.ACDCSE_dir(),"test/data/load_and_gen/CIGRE_B4_forecasted_P10.csv"),DataFrame; stringtype=String);
load_and_gen_data_fore_P90 = CSV.read(joinpath(_ACDCSE.ACDCSE_dir(),"test/data/load_and_gen/CIGRE_B4_forecasted_P90.csv"),DataFrame; stringtype=String);


time_steps=sort(unique(load_and_gen_data[!,:t]))
buses=sort(unique(load_and_gen_data[!,:bus]))
load_data=load_and_gen_data[load_and_gen_data[!,:type].==0,:]
gen_data=load_and_gen_data[load_and_gen_data[!,:type].==1,:]

data_pfs=Vector{Dict{String,Any}}(undef,length(time_steps))
data_ses=Vector{Dict{String,Any}}(undef,length(time_steps))


modfify_loads_fp!(time_steps, load_data, data_pf, data_pfs, data_se, data_ses)
modify_gen_fp(time_steps, gen_data, data_pf, data_pfs, data_se, data_ses)


for data in data_pfs
    set_fixed_bus_voltages!(data)
    set_fixed_busdc_voltages!(data)
    set_fixed_gen_pg_wind!(data,[3,4])
end

results_fp=Vector{Any}(undef,length(time_steps))
for t in time_steps
    println(reference)
    results_fp[t], σ_dict, data_ses[t] = generate_data_basic_acdcse(data_pfs[t], data_ses[t], nlp_optimizer_pf,"no_branch",reference, sample_error = true);
end

#%%

se_res=Vector{Any}(undef,length(time_steps))
for t in time_steps
    se_res[t]=_ACDCSE.solve_acdcse(data_ses[t], _PM.ACPPowerModel, nlp_optimizer_pf)
end

data_forecast_mean=Vector{Dict{String,Any}}(undef,length(time_steps))
data_forecast_P10=Vector{Dict{String,Any}}(undef,length(time_steps))
data_forecast_P90=Vector{Dict{String,Any}}(undef,length(time_steps))

for t in time_steps
    data_forecast_mean[t]=deepcopy(data_ses[t])
    data_forecast_P10[t]=deepcopy(data_ses[t])
    data_forecast_P90[t]=deepcopy(data_ses[t])
end





for t in time_steps
    load_and_gen_data_fore_mean_t=load_and_gen_data_fore_mean[load_and_gen_data_fore_mean[!,:t].==t,:]
    data_forecast_mean[t]=populate_meas_from_forecast!(data_forecast_mean[t], se_res[t], load_and_gen_data_fore_mean_t; prec_meas=0.05, prec_fore=0.05, prec_virtual=1e-5)

    load_and_gen_data_fore_P10_t=load_and_gen_data_fore_P10[load_and_gen_data_fore_P10[!,:t].==t,:]
    data_forecast_P10[t]=populate_meas_from_forecast!(data_forecast_P10[t], se_res[t], load_and_gen_data_fore_P10_t; prec_meas=0.05, prec_fore=0.05, prec_virtual=1e-5)

    load_and_gen_data_fore_P90_t=load_and_gen_data_fore_P90[load_and_gen_data_fore_P90[!,:t].==t,:]
    data_forecast_P90[t]=populate_meas_from_forecast!(data_forecast_P90[t], se_res[t], load_and_gen_data_fore_P90_t; prec_meas=0.05, prec_fore=0.05, prec_virtual=1e-5)

end

res_fore_mean=Vector{Any}(undef,length(time_steps))
res_fore_P10=Vector{Any}(undef,length(time_steps))
res_fore_P90=Vector{Any}(undef,length(time_steps))


for t in time_steps
    res_fore_mean[t]=_ACDCSE.solve_acdcse(data_forecast_mean[t], _PM.ACPPowerModel, nlp_optimizer_pf)
    res_fore_P10[t]=_ACDCSE.solve_acdcse(data_forecast_P10[t], _PM.ACPPowerModel, nlp_optimizer_pf)
    res_fore_P90[t]=_ACDCSE.solve_acdcse(data_forecast_P90[t], _PM.ACPPowerModel, nlp_optimizer_pf)
end

v_real=[]
v_P10=[]
v_P90=[]
v_mean=[]
for t in time_steps
    push!(v_mean, res_fore_mean[t]["solution"]["busdc"]["4"]["vm"][1])
    push!(v_P10, res_fore_P10[t]["solution"]["busdc"]["4"]["vm"][1])
    push!(v_P90, res_fore_P90[t]["solution"]["busdc"]["4"]["vm"][1])
    push!(v_real, results_fp[t]["solution"]["busdc"]["4"]["vm"][1])
end
p = plt.plot(time_steps, v_real, label="real", lw=2, marker=:circle)

plt.plot!(time_steps, v_mean, label="mean", lw=2, marker=:diamond)
plt.plot!(time_steps, v_P10, label="P10", lw=1, ls=:dash, marker=:utriangle)
plt.plot!(time_steps, v_P90, label="P90", lw=1, ls=:dash, marker=:dtriangle)
plt.xlabel!("Time step")
plt.ylabel!("Vm (busdc 4)")
plt.title!("Forecast vs Real voltages")

plt.savefig("v_forecast_vs_real.png")

