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

nlp_optimizer_fore = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-5, "print_level" => 5,  "sb" => "yes", "max_iter" => 10000
)




function include_set_point_meas!(data;prec::Float64 = 0.05, prec_virtual =1e-5)
    create_vm_meas_set!(data; prec=prec, prec_virtual=prec_virtual)
    create_vmdc_meas_set!(data; prec=prec, prec_virtual=prec_virtual)
end

function include_forecast_measv2!(data::Dict{String,Any},res_t::Dict{String,Any},load_and_gen_data_fore_t;prec_fore::Float64=0.05, prec_virtual::Float64=1e-5)
    
    data["meas"] = Dict{String,Any}() 

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


function generate_data_forecastsv2!(time_steps, data_forecast_mean, data_forecast_P10, data_forecast_P90,
                                  se_res, load_and_gen_data_fore_mean, load_and_gen_data_fore_P10, load_and_gen_data_fore_P90;
                                  prec_fore::Float64=0.01, prec_meas::Float64=0.05)

    t_se = 1
    for t in time_steps
        (t-1) % 4 == 0 ? t_se = t : nothing

        load_and_gen_data_fore_mean_t = load_and_gen_data_fore_mean[load_and_gen_data_fore_mean[!,:t].==t, :]
        
        
        include_forecast_measv2!(data_forecast_mean[t], se_res[t], load_and_gen_data_fore_mean_t; prec_fore=prec_fore, prec_virtual=1e-5)
        include_set_point_meas!(data_forecast_mean[t]; prec=prec_meas, prec_virtual=1e-5)
        # include_set_point_meas!(data_forecast_mean[t]; prec=prec_meas, prec_virtual=1e-5)
        
        load_and_gen_data_fore_P10_t = load_and_gen_data_fore_P10[load_and_gen_data_fore_P10[!,:t].==t, :]
        # include_se_estimates_meas!(data_forecast_P10[t], se_res[t_se]; prec_meas=prec_meas, prec_virtual=1e-5)
        include_forecast_measv2!(data_forecast_P10[t], se_res[t], load_and_gen_data_fore_P10_t; prec_fore=prec_fore, prec_virtual=1e-5)
        include_set_point_meas!(data_forecast_P10[t]; prec=prec_meas, prec_virtual=1e-5)
        # include_set_point_meas!(data_forecast_P10[t]; prec=prec_meas, prec_virtual=1e-5)

        load_and_gen_data_fore_P90_t = load_and_gen_data_fore_P90[load_and_gen_data_fore_P90[!,:t].==t, :]
        # include_se_estimates_meas!(data_forecast_P90[t], se_res[t_se]; prec_meas=prec_meas, prec_virtual=1e-5)
        include_forecast_measv2!(data_forecast_P90[t], se_res[t], load_and_gen_data_fore_P90_t; prec_fore=prec_fore, prec_virtual=1e-5)
        include_set_point_meas!(data_forecast_P90[t]; prec=prec_meas, prec_virtual=1e-5)
        # include_set_point_meas!(data_forecast_P90[t]; prec=prec_meas, prec_virtual=1e-5)
    end

end


function make_forecast_comparison_df(lst, time_steps, res_fore_mean, res_fore_P10, res_fore_P90, results_fp)
    df_out = DataFrame(time=Int[], real=Float64[], P10=Float64[], P90=Float64[], mean=Float64[], variable=String[], component=String[],cmp_id= String[], c = String[])
    for var_plot in lst
        v_real = Float64[]
        v_P10  = Float64[]
        v_P90  = Float64[]
        v_mean = Float64[]
        cmp_id = var_plot[3]
        var    = var_plot[2]
        cmp    = var_plot[1]
        c      = var_plot[4]
        for t in time_steps
            mean_val = res_fore_mean[t]["solution"][cmp][string(cmp_id)][var][c]
            max_val  = res_fore_P90[t]["solution"][cmp][string(cmp_id)][var][c]
            min_val  = res_fore_P10[t]["solution"][cmp][string(cmp_id)][var][c]

            push!(v_mean, st.median([mean_val, max_val, min_val]))
            push!(v_P10, minimum([mean_val, max_val, min_val]))
            push!(v_P90, maximum([mean_val, max_val, min_val]))
            push!(v_real, results_fp[t]["solution"][cmp][string(cmp_id)][var][c])
        end
        var_col = fill(string(var), length(v_real))
        cmp_col = fill(string(cmp),  length(v_real))
        cmp_id_col = fill(string(cmp_id), length(v_real))
        c_col = fill(string(c), length(v_real))
        df = DataFrame(time = time_steps, real = v_real, P10 = v_P10, P90 = v_P90, mean = v_mean, variable = var_col, component = cmp_col,cmp_id= cmp_id_col, c = c_col)
        append!(df_out, df)
    end
    return df_out
end


data_pf= _ACDCSE.quickget_cigre_B4() # Load the reference files se
data_se= _ACDCSE.quickget_cigre_B4() # Load the reference files se

reference=[1,2,3,4]

load_and_gen_data = CSV.read(joinpath(_ACDCSE.ACDCSE_dir(),"test/data/load_and_gen/day_ahead_3/CIGRE_B4_measured.csv"),DataFrame; stringtype=String); #load the measured data
load_and_gen_data_fore_mean = CSV.read(joinpath(_ACDCSE.ACDCSE_dir(),"test/data/load_and_gen/day_ahead_3/CIGRE_B4_forecasted.csv"),DataFrame; stringtype=String); #load the forecast data
load_and_gen_data_fore_P10 = CSV.read(joinpath(_ACDCSE.ACDCSE_dir(),"test/data/load_and_gen/day_ahead_3/CIGRE_B4_forecasted_P10.csv"),DataFrame; stringtype=String); #load the forecast data
load_and_gen_data_fore_P90 = CSV.read(joinpath(_ACDCSE.ACDCSE_dir(),"test/data/load_and_gen/day_ahead_3/CIGRE_B4_forecasted_P90.csv"),DataFrame; stringtype=String); #load the forecast data


time_steps=sort(unique(load_and_gen_data[!,:t]))

buses=sort(unique(load_and_gen_data[!,:bus]))
load_data=load_and_gen_data[load_and_gen_data[!,:type].==0,:]
gen_data=load_and_gen_data[load_and_gen_data[!,:type].==1,:]



data_pfs=Vector{Dict{String,Any}}(undef,length(time_steps))
data_ses=Vector{Dict{String,Any}}(undef,length(time_steps))


modfify_loads_fp!(time_steps, load_data, data_pf, data_pfs, data_se, data_ses) # set for every time step the loads according to the load_data
modify_gen_fp(time_steps, gen_data, data_pf, data_pfs, data_se, data_ses) # set for every time step the generators according to the gen_data (wind gen) info


#set constraints
for data in data_pfs
    set_fixed_bus_voltages!(data) # set fixed voltages
    set_fixed_busdc_voltages!(data) # set fixed voltages
    set_fixed_gen_pg_wind_conv!(data,[3,4]) # set fixed power 
end


# run the power flow for every time step
results_fp=Vector{Any}(undef,length(time_steps)) 
solution_status_fp=Vector{String}(undef,length(time_steps))
for t in time_steps
    results_fp[t], σ_dict, data_ses[t] = generate_data_basic_acdcse(data_pfs[t], data_ses[t], nlp_optimizer_pf,"no_branch",reference, sample_error = false);
    solution_status_fp[t] = string(results_fp[t]["termination_status"])
end

#% run a se to estimate the gen 
se_res=Vector{Any}(undef,length(time_steps))
solution_status_se=Vector{String}(undef,length(time_steps))
for t in time_steps
    se_res[t]=_ACDCSE.solve_acdcse(data_ses[t], _PM.ACPPowerModel, nlp_optimizer_pf)
    solution_status_se[t] = string(se_res[t]["termination_status"])
end


# create the data for the forecasts
data_forecast_mean=Vector{Dict{String,Any}}(undef,length(time_steps)) 
data_forecast_P10=Vector{Dict{String,Any}}(undef,length(time_steps))
data_forecast_P90=Vector{Dict{String,Any}}(undef,length(time_steps))

for t in time_steps
    data_forecast_mean[t]=deepcopy(data_ses[t])
    data_forecast_P10[t]=deepcopy(data_ses[t])
    data_forecast_P90[t]=deepcopy(data_ses[t])
end

#generate the data
generate_data_forecastsv2!(time_steps, data_forecast_mean, data_forecast_P10, data_forecast_P90,
                              se_res, load_and_gen_data_fore_mean, load_and_gen_data_fore_P10, load_and_gen_data_fore_P90;
                              prec_fore=0.01, prec_meas=0.05)


                              

res_fore_mean=Vector{Any}(undef,length(time_steps))
res_fore_P10=Vector{Any}(undef,length(time_steps))
res_fore_P90=Vector{Any}(undef,length(time_steps))

solution_fore_mean=Vector{String}(undef,length(time_steps))
solution_fore_P10=Vector{String}(undef,length(time_steps))
solution_fore_P90=Vector{String}(undef,length(time_steps))  


for t in time_steps
    # set_fixed_busdc_voltages!(data_forecast_mean[t])
    # set_fixed_busdc_voltages!(data_forecast_P10[t])
    # set_fixed_busdc_voltages!(data_forecast_P90[t])
    data_forecast_mean[t]["se_settings"]["criterion"]="rwls"
    data_forecast_P10[t]["se_settings"]["criterion"]="rwls"
    data_forecast_P90[t]["se_settings"]["criterion"]="rwls"

    res_fore_mean[t]=_ACDCSE.solve_acdcse(data_forecast_mean[t], _PM.ACPPowerModel, nlp_optimizer_fore)
    res_fore_P10[t]=_ACDCSE.solve_acdcse(data_forecast_P10[t], _PM.ACPPowerModel, nlp_optimizer_fore)
    res_fore_P90[t]=_ACDCSE.solve_acdcse(data_forecast_P90[t], _PM.ACPPowerModel, nlp_optimizer_fore)
    solution_fore_mean[t] = string(res_fore_mean[t]["termination_status"])
    solution_fore_P10[t] = string(res_fore_P10[t]["termination_status"])
    solution_fore_P90[t] = string(res_fore_P90[t]["termination_status"])
end 


mean=[]
P10=[]
P90=[]
for i in time_steps
    push!(mean, res_fore_mean[i]["solve_time"])
    push!(P10, res_fore_P10[i]["solve_time"])
    push!(P90, res_fore_P90[i]["solve_time"])
end

mean_solve_time = st.median(mean)
P10_solve_time  = st.median(P10)
P90_solve_time  = st.median(P90)

max_mean_solve_time = st.maximum(mean)
max_P10_solve_time  = st.maximum(P10)
max_P90_solve_time  = st.maximum(P90)

min_mean_solve_time = st.minimum(mean)
min_P10_solve_time  = st.minimum(P10)
min_P90_solve_time  = st.minimum(P90)    



lst=[]
for (key,bus) in data_pfs[1]["busdc"]
    push!(lst, ["busdc","vm", parse(Int64,key), 1])
    push!(lst, ["busdc","vm", parse(Int64, key), 2])
end

for (key,branchdc) in data_pfs[1]["branchdc"]
    push!(lst, ["branchdc","i_from", parse(Int64,key), 1])
    push!(lst, ["branchdc","i_from", parse(Int64, key), 2])
end


for (key,convdc) in data_pfs[1]["convdc"]
    push!(lst, ["convdc","vmfilt", parse(Int64,key), 1])
    push!(lst, ["convdc","vmfilt", parse(Int64, key), 2])
end



df_forecast_comparison = make_forecast_comparison_df(lst, time_steps, res_fore_mean, res_fore_P10, res_fore_P90, results_fp)

CSV.write("test/data/acdcpaper_2026_long_term_forecast_comparison.csv", df_forecast_comparison)