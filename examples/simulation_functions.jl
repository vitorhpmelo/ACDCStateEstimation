using BenchmarkTools


function silence(f, args...; kwargs...)
    redirect_stdout(devnull) do
        redirect_stderr(devnull) do
            return f(args...; kwargs...)
        end
    end
end

macro silence(expr)
    return quote
        redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                $(esc(expr))
            end
        end
    end
end


function add_iterations_from_log!(txt)
    logtxt = read(txt, String)
    return match(r"Number of Iterations\.*:\s*(\d+)", logtxt)[1]
end

function set_fixed_bus_voltages!(data_pf)
    for (k, bus) in data_pf["bus"]
        if bus["bus_type"] != 1
            bus["vmax"] = bus["vm"]*1.001
            bus["vmin"] = bus["vm"]*0.999
        end
    end
end

function set_fixed_busdc_voltages!(data_pf)
    relax=0.01
    p=1.0+relax
    n=1.0-relax
    for (k, conv) in data_pf["convdc"]
        if conv["type_dc"] == 2
            key_busdc=conv["busdc_i"]
            vset=data_pf["busdc"]["$key_busdc"]["Vdc"]
            data_pf["busdc"]["$key_busdc"]["Vdcmax"] = [vset[1],vset[2],0.1]
            data_pf["busdc"]["$key_busdc"]["Vdcmin"] = [vset[1],vset[2],-0.1]
        end
    end
end


function set_fixed_bus_conv_voltages!(data_pf)
    for (k, conv) in data_pf["convdc"]
        if conv["type_ac"] != 1
            busac=conv["busac_i"]
            data_pf["bus"]["$busac"]["vmax"] = data_pf["bus"]["$busac"]["vm"]*1.01
            data_pf["bus"]["$busac"]["vmin"] = data_pf["bus"]["$busac"]["vm"]*0.99
        end
    end
end

function set_fixed_gen_pg!(data_pf)
    for (k, gen) in data_pf["gen"]
        bus=gen["gen_bus"]
        if data_pf["bus"]["$bus"]["bus_type"] == 2
            gen["pmax"] = gen["pg"]
            gen["pmin"] = gen["pg"]
        end
    end
end

function set_fixed_gen_pg_wind!(data_pf,buslist=[])
    for bus in buslist
        for (key, gen) in data_pf["gen"]
            if gen["gen_bus"]==bus
                pg= gen["pg"]
                gen["pmax"] = pg
                gen["pmin"] = pg
            end
        end
    end
end

function set_fixed_gen_pg_result(data_pf,result)
    for (k, gen) in data_pf["gen"]
        bus=gen["gen_bus"]

        if data_pf["bus"]["$bus"]["bus_type"] == 2
            pg= abs(result["solution"]["gen"]["$k"]["pg"])
            gen["pmax"] = pg*1.05
            gen["pmin"] = pg*0.95
        end
    end
end

function set_fixed_gen_pg_result_wg(data_pf,result,wg,wind_var=0.1)

    wind_var=wind_var/100    
    pg= abs(result["solution"]["gen"]["$wg"]["pg"])
    pg_sample = rand(_ACDCSE._DST.Truncated(_ACDCSE._DST.Normal(pg, pg*wind_var), pg-pg*wind_var, pg+pg*wind_var))
    data_pf["gen"]["$wg"]["pmax"] = pg_sample*1.01
    data_pf["gen"]["$wg"]["pmin"] = pg_sample*0.99

end

function modify_loads_over_time_rand!(data_pf_all, load2mod; ratio=1,seed=1,limit=2)
    T = length(data_pf_all)
    _ACDCSE._RAN.seed!(seed)
    N_loads=length(load2mod)
    for t in 2:T
        χ_p = _ACDCSE._RAN.rand(_ACDCSE._DST.Truncated(_ACDCSE._DST.Normal(0, 1), -limit, limit), N_loads)
        χ_q = _ACDCSE._RAN.rand(_ACDCSE._DST.Truncated(_ACDCSE._DST.Normal(0, 1), -limit, limit), N_loads)

        for (i, load) in enumerate(load2mod)

            μp = data_pf_all[t-1]["load"][load]["pd"]
            μq = data_pf_all[t-1]["load"][load]["qd"]
            σp = (abs(μp) * (ratio/100))/3
            σq = (abs(μq) * (ratio/100))/3
            print("Modifying load $load at time $t: P from $μp to $(μp + σp * χ_p[i]), Q from $μq to $(μq + σq * χ_q[i])\n")

            data_pf_all[t]["load"][load]["pd"] = μp + σp * χ_p[i]
            data_pf_all[t]["load"][load]["qd"] = μq + σq * χ_q[i]
        end
    end
end


function modify_loads_over_time_ramp!(data_pf_all, tstart,Tlen, load2mod,inc=5)
    
    inc=inc/100
    kp=Dict()
    kq=Dict()
    for load in load2mod
        kp[load] = data_pf_all[tstart]["load"][load]["pd"]
        kq[load] = data_pf_all[tstart]["load"][load]["qd"]
    end
    Δt=Tlen
    for t in tstart:tstart+Tlen
        for load in load2mod
            data_pf_all[t]["load"][load]["pd"] = data_pf_all[t]["load"][load]["pd"]+kp[load]*inc*(t-tstart)/Δt
            data_pf_all[t]["load"][load]["qd"] = data_pf_all[t]["load"][load]["qd"]+kq[load]*inc*(t-tstart)/Δt
        end
    end

    for t in tstart+Tlen+1:length(data_pf_all)
        for load in load2mod
            data_pf_all[t]["load"][load]["pd"] = data_pf_all[t]["load"][load]["pd"]+kp[load]*inc
            data_pf_all[t]["load"][load]["qd"] = data_pf_all[t]["load"][load]["qd"]+kq[load]*inc
        end
    end
end


function modify_bus_v_over_time_ramp!(data_pf_all, tstart,Tlen, bus2mod,inc=5.0)
    
    inc=inc/100
    kv=Dict()
    
    for bus in bus2mod
        kv["$bus"] = data_pf_all[tstart]["bus"]["$bus"]["vm"]
    end
    Δt=Tlen
    for t in tstart:tstart+Tlen
        for bus in bus2mod
            println(bus)
            data_pf_all[t]["bus"]["$bus"]["vm"] = data_pf_all[t]["bus"]["$bus"]["vm"]+kv["$bus"]*inc*(t-tstart)/Δt
        end
    end

    for t in tstart+Tlen+1:length(data_pf_all)
        for bus in bus2mod
            data_pf_all[t]["bus"]["$bus"]["vm"] = data_pf_all[t]["bus"]["$bus"]["vm"]+kv["$bus"]*inc
        end
    end
end




function determine_loads_to_modify(bus2mod,data_pf)
    loads=keys(data_pf["load"])
    load2mod=[]
    for load in loads
        if data_pf["load"][load]["load_bus"] in bus2mod
            push!(load2mod, load)
        end
    end
    return load2mod
end



function determine_gens_to_modify(bus2mod,data_pf)
    gens=keys(data_pf["gen"])
    gens2mod=[]
    for gen in gens
        if data_pf["gen"][gen]["gen_bus"] in bus2mod
            push!(gens2mod, gen)
        end
    end
    return gens2mod
end

function create_se_data!(data_pf_all, data_se_all, results, nlp_optimizer,reference;wind_gen=[],wind_gen_var=0.0)
    T = length(data_pf_all)
    for t in 1:T
        data_se_all[t] = deepcopy(data_pf_all[t])
    end

    for t in 1:T
        set_fixed_bus_voltages!(data_pf_all[t])
        set_fixed_bus_conv_voltages!(data_pf_all[t])

        if 1 < t

            set_fixed_gen_pg_result(data_pf_all[t],results[t-1])
            for wg in wind_gen
                set_fixed_gen_pg_result_wg(data_pf_all[t],results[t-1],wg,wind_gen_var)
            end
        end

        result, σ_dict, data_se_t = generate_data_basic_acdcse(
            data_pf_all[t], data_se_all[t], nlp_optimizer, "all", reference, sample_error=false
        )
        results[t] = result
        data_se_all[t] = data_se_t
    end

end




function run_se_simulation_map!(
    T::Int, n::Int, data_se_all::Vector{Dict}, data_scada_all::Vector{Dict}, data_pmu_all::Vector{Dict}, data_hyb_all::Vector{Dict}, results_pf::Vector{Dict},
    sampling_scada::Int, sampling_pmu::Int,
    meas_set::DataFrame, min_sigma::Float64, a::Float64, nlp_optimizer, nlp_optimizer_hyb, nlp_optimizer_pmu,
    ipot_file_pmu, ipot_file_hyb;se_objective=Dict("scada"=>"rwls","pmu"=>"rwls" ,"hyb"=>"rwls"), with_noise=true, noise_bound=true, noise_bound_factor=2.5
)
    scada_measurements = ["scada", "dc", "conv"]
    pmu_measurements = ["pmu", "dtu"]

    ratio_scada_pmu = Int(sampling_scada / sampling_pmu)
    Tscada = Vector(1:ratio_scada_pmu:T)
    d_meas_set = create_dmeas_set(meas_set, data_se_all[1])

    d_meas_set_scada = Dict()
    for key in scada_measurements
        d_meas_set_scada[key] = deepcopy(d_meas_set[key])
    end

    d_meas_set_pmu = Dict()
    for key in pmu_measurements
        d_meas_set_pmu[key] = deepcopy(d_meas_set[key])
    end

    df_errors = DataFrame(var=[], comp=[], cmp_id=[], c=[], val_pf=[], val_se=[], error=[], se=[], t=[], n=[])
    df_conv = DataFrame(n=[], t=[], termination_status_hyb=[], termination_status_pmu=[], n_iters_hyb=[],n_iters_pmu=[])
    df_conv_prior = DataFrame(n=[], t=[], termination_status_scada=[])

    se_res_prior = Dict()

    for t in 1:T
        data_se_t = deepcopy(data_se_all[t])
        d_keys, d_prec = filtermeas_set!(data_se_t, d_meas_set, sample_error=false, min_sig=min_sigma) #filters the measurement set

        if with_noise==true
            introduce_noise!(data_se_t, d_prec, d_keys, seed=(t+(n-1)*T), min_sig=min_sigma, bound=true)
        end

        data_scada_all[t] = deepcopy(data_se_t)
        data_pmu_all[t] = deepcopy(data_se_t)

        split_meas_set!(d_keys, data_scada_all[t], scada_measurements)
        split_meas_set!(d_keys, data_pmu_all[t], pmu_measurements) #splits the measurements into SCADA and PMU sets;

        data_hyb_all[t] = deepcopy(data_pmu_all[t])

        t_scada = maximum(filter(x -> x <= t, Tscada))

        for (key, meas) in data_scada_all[t_scada]["meas"]
            if haskey(data_hyb_all[t]["meas"], key)
                max_key = maximum(parse.(Int, keys(data_hyb_all[t]["meas"])))
                new_key = max_key + 1
                data_hyb_all[t]["meas"]["$new_key"] = deepcopy(meas)
            else
                data_hyb_all[t]["meas"][key] = deepcopy(meas)
            end
        end

        if t in Tscada
            data_scada_all[t]["se_settings"]["criterion"] = se_objective["scada"]
            se_res_scada = _ACDCSE.solve_acdcse(data_scada_all[t], _PM.ACPPowerModel, nlp_optimizer) #Solve SE
            df_error_scada = df_error_var(results_pf[t], se_res_scada)
            df_error_scada[!, :se] .= "scada"
            df_error_scada[!, :t] .= t
            df_error_scada[!, :n] .= n

            df_conv_prior = vcat(df_conv_prior, DataFrame(n=n, t=t, termination_status_scada=se_res_scada["termination_status"]))

            se_res_prior = deepcopy(se_res_scada)
            net = generate_data_prior!(se_res_prior, data_scada_all[t], data_pmu_all[t], d_keys, d=a, min_sig=min_sigma, print_info=false)
            df_errors = vcat(df_errors, df_error_scada)
        else
            print("Using prior from t=$t_scada for t=$t\n")
            print(se_res_prior)
            net = generate_data_prior!(se_res_prior, data_scada_all[t_scada], data_pmu_all[t], d_keys, d=a, min_sig=min_sigma, print_info=false)
        end

        data_hyb_all[t]["se_settings"]["criterion"] = se_objective["hyb"]
        se_res_hyb = _ACDCSE.solve_acdcse(data_hyb_all[t], _PM.ACPPowerModel, nlp_optimizer_hyb) #Solve WLS
        df_error_hyb = df_error_var(results_pf[t], se_res_hyb)
        df_error_hyb[!, :se] .= "hyb"
        df_error_hyb[!, :t] .= t
        df_error_hyb[!, :n] .= n
        df_errors = vcat(df_errors, df_error_hyb)

        data_pmu_all[t]["se_settings"]["criterion"] = se_objective["pmu"]
        se_res_pmu = _ACDCSE.solve_acdcse(data_pmu_all[t], _PM.ACPPowerModel, nlp_optimizer_pmu) #Solve pmu

        df_conv = vcat(df_conv, DataFrame(
            n=n, t=t,
            termination_status_hyb=se_res_hyb["termination_status"],
            termination_status_pmu=se_res_pmu["termination_status"],
            n_iters_pmu=add_iterations_from_log!(ipot_file_pmu),
            n_iters_hyb=add_iterations_from_log!(ipot_file_hyb)
        ))
        df_error_pmu = df_error_var(results_pf[t], se_res_pmu)
        df_error_pmu[!, :se] .= "pmu"
        df_error_pmu[!, :t] .= t
        df_error_pmu[!, :n] .= n
        df_errors = vcat(df_errors, df_error_pmu)
    end

    return df_errors, df_conv, df_conv_prior
end


function run_se_simulation_map_test!(
    T::Int, n::Int, data_se_all::Vector{Dict}, data_scada_all::Vector{Dict}, data_pmu_all::Vector{Dict}, data_hyb_all::Vector{Dict}, results_pf::Vector{Dict},
    sampling_scada::Int, sampling_pmu::Int,
    meas_set::DataFrame, min_sigma::Float64, a::Float64, nlp_optimizer, nlp_optimizer_hyb, nlp_optimizer_pmu,
    ipot_file_pmu, ipot_file_hyb;se_objective=Dict("scada"=>"rwls","pmu"=>"rwls" ,"hyb"=>"rwls"), with_noise=true, noise_bound=true, noise_bound_factor=2.5
)
    scada_measurements = ["scada", "dc", "conv"]
    pmu_measurements = ["pmu", "dtu"]

    ratio_scada_pmu = Int(sampling_scada / sampling_pmu)
    Tscada = Vector(1:ratio_scada_pmu:T)
    d_meas_set = create_dmeas_set(meas_set, data_se_all[1])

    d_meas_set_scada = Dict()
    for key in scada_measurements
        d_meas_set_scada[key] = deepcopy(d_meas_set[key])
    end

    d_meas_set_pmu = Dict()
    for key in pmu_measurements
        d_meas_set_pmu[key] = deepcopy(d_meas_set[key])
    end

    
    dfs_scada=DataFrame[]
    dfs_pmu=DataFrame[]
    dfs_hyb=DataFrame[]
    
    d_keys_all=[]


    for t in 1:T
        data_se_t = deepcopy(data_se_all[t])
        d_keys, d_prec = filtermeas_set!(data_se_t, d_meas_set, sample_error=false, min_sig=min_sigma) #filters the measurement set


        if with_noise==true
            introduce_noise!(data_se_t, d_prec, d_keys, seed=(t+(n-1)*T), min_sig=min_sigma, bound=true)
        end

        data_scada_all[t] = deepcopy(data_se_t)
        data_pmu_all[t] = deepcopy(data_se_t)

        split_meas_set!(d_keys, data_scada_all[t], scada_measurements)
        split_meas_set!(d_keys, data_pmu_all[t], pmu_measurements) #splits the measurements into SCADA and PMU sets;

        data_hyb_all[t] = deepcopy(data_pmu_all[t])

        t_scada = maximum(filter(x -> x <= t, Tscada))

        for (key, meas) in data_scada_all[t_scada]["meas"]
            if haskey(data_hyb_all[t]["meas"], key)
                max_key = maximum(parse.(Int, keys(data_hyb_all[t]["meas"])))
                new_key = max_key + 1
                data_hyb_all[t]["meas"]["$new_key"] = deepcopy(meas)
            else
                data_hyb_all[t]["meas"][key] = deepcopy(meas)
            end
        end
        
        if t in Tscada
            println("t=$t, t_scada=$t_scada")
            
            df_meas_scada= meas_dict_to_dataframe(data_scada_all[t]["meas"])

            df_meas_scada[!, :se] .= "scada"
            df_meas_scada[!, :t] .= t
            df_meas_scada[!, :n] .= n
            dfs_scada = vcat(dfs_scada, df_meas_scada)
            print(t)

            # se_res_prior = deepcopy(se_res_scada)
            # net = generate_data_prior!(se_res_prior, data_scada_all[t], data_pmu_all[t], d_keys, d=a, min_sig=min_sigma, print_info=false)
            # df_errors = vcat(df_errors, df_error_scada)
        else
            # print("Using prior from t=$t_scada for t=$t\n")
            # print(se_res_prior)
            # net = generate_data_prior!(se_res_prior, data_scada_all[t_scada], data_pmu_all[t], d_keys, d=a, min_sig=min_sigma, print_info=false)
        end

        df_meas_hyb = meas_dict_to_dataframe(data_hyb_all[t]["meas"])
        df_meas_hyb[!, :se] .= "hyb"
        df_meas_hyb[!, :t] .= t
        df_meas_hyb[!, :n] .= n
        dfs_hyb = vcat(dfs_hyb, df_meas_hyb)
        df_meas_pmu = meas_dict_to_dataframe(data_pmu_all[t]["meas"])
        df_meas_pmu[!, :se] .= "pmu"
        df_meas_pmu[!, :t] .= t
        df_meas_pmu[!, :n] .= n
        dfs_pmu = vcat(dfs_pmu, df_meas_pmu)

        
    end
    df_meas_scada=vcat(dfs_scada...)
    df_meas_hyb=vcat(dfs_hyb...)
    df_meas_pmu=vcat(dfs_pmu...)
    return df_meas_scada, df_meas_hyb, df_meas_pmu,d_keys_all

end


function run_se_simulation_map_only_hyb!(
    T::Int, n::Int, data_se_all::Vector{Dict}, data_scada_all::Vector{Dict}, data_pmu_all::Vector{Dict}, data_hyb_all::Vector{Dict}, results_pf::Vector{Dict},
    sampling_scada::Int, sampling_pmu::Int,
    meas_set::DataFrame, min_sigma::Float64, a::Float64, nlp_optimizer, nlp_optimizer_hyb, nlp_optimizer_pmu,
    ipot_file_pmu, ipot_file_hyb;se_objective=Dict("scada"=>"rwls","pmu"=>"rwls" ,"hyb"=>"rwls"), with_noise=true, noise_bound=true, noise_bound_factor=2.5
)
    scada_measurements = ["scada", "dc", "conv"]
    pmu_measurements = ["pmu", "dtu"]

    ratio_scada_pmu = Int(sampling_scada / sampling_pmu)
    Tscada = Vector(1:ratio_scada_pmu:T)
    d_meas_set = create_dmeas_set(meas_set, data_se_all[1])

    d_meas_set_scada = Dict()
    for key in scada_measurements
        d_meas_set_scada[key] = deepcopy(d_meas_set[key])
    end

    d_meas_set_pmu = Dict()
    for key in pmu_measurements
        d_meas_set_pmu[key] = deepcopy(d_meas_set[key])
    end

    df_errors = DataFrame(var=[], comp=[], cmp_id=[], c=[], val_pf=[], val_se=[], error=[], se=[], t=[], n=[])
    df_conv = DataFrame(n=[], t=[], termination_status_hyb=[], n_iters_hyb=[],)
    df_conv_prior = DataFrame(n=[], t=[], termination_status_scada=[])

    se_res_prior = Dict()

    for t in 1:T
        data_se_t = deepcopy(data_se_all[t])
        d_keys, d_prec = filtermeas_set!(data_se_t, d_meas_set, sample_error=false, min_sig=min_sigma) #filters the measurement set

        if with_noise==true
            introduce_noise!(data_se_t, d_prec, d_keys, seed=(t+(n-1)*T), min_sig=min_sigma, bound=true,bound_factor=noise_bound_factor)
        end

        data_scada_all[t] = deepcopy(data_se_t)
        data_pmu_all[t] = deepcopy(data_se_t)

        split_meas_set!( d_keys, data_scada_all[t], scada_measurements)
        split_meas_set!( d_keys, data_pmu_all[t], pmu_measurements) #splits the measurements into SCADA and PMU sets;

        data_hyb_all[t] = deepcopy(data_pmu_all[t])

        t_scada = maximum(filter(x -> x <= t, Tscada))

        for (key, meas) in data_scada_all[t_scada]["meas"]
            if haskey(data_hyb_all[t]["meas"], key)
                max_key = maximum(parse.(Int, keys(data_hyb_all[t]["meas"])))
                new_key = max_key + 1
                data_hyb_all[t]["meas"]["$new_key"] = deepcopy(meas)
            else
                data_hyb_all[t]["meas"][key] = deepcopy(meas)
            end
        end

        if t in Tscada
            data_scada_all[t]["se_settings"]["criterion"] = se_objective["scada"]
            se_res_scada = _ACDCSE.solve_acdcse(data_scada_all[t], _PM.ACPPowerModel, nlp_optimizer) #Solve SE
            df_error_scada = df_error_var(results_pf[t], se_res_scada)
            df_error_scada[!, :se] .= "scada"
            df_error_scada[!, :t] .= t
            df_error_scada[!, :n] .= n

            df_conv_prior = vcat(df_conv_prior, DataFrame(n=n, t=t, termination_status_scada=se_res_scada["termination_status"]))

            #se_res_prior = deepcopy(se_res_scada)
            #net = generate_data_prior!(se_res_prior, data_scada_all[t], data_pmu_all[t], d_meas_set_scada, d=a, min_sig=min_sigma, print_info=false)
            df_errors = vcat(df_errors, df_error_scada)
        else
            #print("Using prior from t=$t_scada for t=$t\n")
            #print(se_res_prior)
            #net = generate_data_prior!(se_res_prior, data_scada_all[t_scada], data_pmu_all[t], d_meas_set_scada, d=a, min_sig=min_sigma, print_info=false)
        end
        
        data_hyb_all[t]["se_settings"]["criterion"] = se_objective["hyb"]
        se_res_hyb = _ACDCSE.solve_acdcse(data_hyb_all[t], _PM.ACPPowerModel, nlp_optimizer_hyb) #Solve WLS
        df_error_hyb = df_error_var(results_pf[t], se_res_hyb)
        df_error_hyb[!, :se] .= "hyb"
        df_error_hyb[!, :t] .= t
        df_error_hyb[!, :n] .= n
        df_errors = vcat(df_errors, df_error_hyb)

        #se_res_pmu = _ACDCSE.solve_acdcse(data_pmu_all[t], _PM.ACPPowerModel, nlp_optimizer_pmu) #Solve pmu

        df_conv = vcat(df_conv, DataFrame(
            n=n, t=t,
            termination_status_hyb=se_res_hyb["termination_status"],
            n_iters_hyb=add_iterations_from_log!(ipot_file_hyb)
        ))

    end

    return df_errors, df_conv, df_conv_prior
end




function run_fase_experiment_stationary(N, n0, a, min_sigma, result, data_hyb_sr, d_prec, d_keys, d_meas_set, nlp_optimizer_hyb, nlp_optimizer_scada, nlp_optimizer_pmu, ipopt_files;
    se_objective=Dict("scada"=>"rwls","pmu"=>"rwls" ,"hyb"=>"rwls"),cutoff::Float64=1e-6)

    scada_measurements = ["scada", "dc", "conv"]
    d_meas_set_scada = Dict()
    for key in scada_measurements
        d_meas_set_scada[key] = deepcopy(d_meas_set[key])
    end

    pmu_measurements = ["pmu", "dtu"]
    d_meas_set_pmu = Dict()
    for key in pmu_measurements
        d_meas_set_pmu[key] = deepcopy(d_meas_set[key])
    end



    d_conv_scada = Dict(
        "termination_status" => Vector{Any}(undef, N),
        "n_iters" => Vector{Any}(undef, N),
        "solve_time" => Vector{Any}(undef, N),
    )

    d_conv_pmu = Dict(
        "termination_status" => Vector{Any}(undef, N),
        "n_iters" => Vector{Any}(undef, N),
        "solve_time" => Vector{Any}(undef, N),
    )

    d_conv_hyb = Dict(
        "termination_status" => Vector{Any}(undef, N),
        "n_iters" => Vector{Any}(undef, N),
        "solve_time" => Vector{Any}(undef, N),
    )

    df_errors_total = DataFrame()

    dfs_priors_all = DataFrame[]

    for n in n0+1:n0+N
        data_hyb = deepcopy(data_hyb_sr) #copy data to avoid modifying original data

        introduce_noise!(data_hyb, d_prec, d_keys, seed=n, min_sig=min_sigma, bound=true)

        data_pmu = deepcopy(data_hyb)
        data_scada = deepcopy(data_hyb)


        split_meas_set!(d_keys, data_scada, ["scada", "dc", "conv"])
        split_meas_set!(d_keys, data_pmu, ["pmu", "dtu"]) #splits the measurements into SCADA and PMU sets;]


        data_scada["se_settings"]["criterion"] = se_objective["scada"]
        se_res_scada = _ACDCSE.solve_acdcse(data_scada, _PM.ACPPowerModel, nlp_optimizer_scada,false) #Solve SE
        n_iters = add_iterations_from_log!(ipopt_files["scada"])
        df_error_scada = df_error_var(result, se_res_scada)
        df_error_scada[!, :se] .= "scada"
        df_error_scada[!, :n] .= n


        d_conv_scada["termination_status"][n-n0] = string(se_res_scada["termination_status"])
        d_conv_scada["n_iters"][n-n0] = n_iters
        d_conv_scada["solve_time"][n-n0] = se_res_scada["solve_time"]


        data_hyb["se_settings"]["criterion"] = se_objective["hyb"]
        se_res_hyb = _ACDCSE.solve_acdcse(data_hyb, _PM.ACPPowerModel, nlp_optimizer_hyb,false) #Solve WLS
        n_iters = add_iterations_from_log!(ipopt_files["hyb"])
        d_conv_hyb["n_iters"][n-n0] = n_iters
        df_error_hyb = df_error_var(result, se_res_hyb)
        df_error_hyb[!, :se] .= "hyb"
        df_error_hyb[!, :n] .= n
        d_conv_hyb["termination_status"][n-n0] = string(se_res_hyb["termination_status"])
        d_conv_hyb["solve_time"][n-n0] = se_res_hyb["solve_time"]



        net = generate_data_prior!(se_res_scada, data_scada, data_pmu, d_keys, d=a, min_sig=min_sigma, print_info=true, cutoff=cutoff) #Generate prior data for FASE



        data_pmu["se_settings"]["criterion"] = se_objective["pmu"]
        se_res_pmu = _ACDCSE.solve_acdcse(data_pmu, _PM.ACPPowerModel, nlp_optimizer_pmu,false) #Solve FASE
        df_prior_internal = extract_prior_dataframe(data_pmu)
        df_prior_internal[!, :n] .= n

        dfs_priors_all = push!(dfs_priors_all, df_prior_internal)

        n_iters = add_iterations_from_log!(ipopt_files["pmu"])
        
        df_error_pmu = df_error_var(result, se_res_pmu)
        df_error_pmu[!, :se] .= "pmu"
        df_error_pmu[!, :n] .= n
        
        df_errors_total = vcat(df_errors_total, df_error_pmu, df_error_hyb, df_error_scada)
        
        d_conv_pmu["termination_status"][n-n0] = string(se_res_pmu["termination_status"])
        d_conv_pmu["n_iters"][n-n0] = n_iters
        d_conv_pmu["solve_time"][n-n0] = se_res_pmu["solve_time"]

    end

    df_conv_pmu = DataFrame(d_conv_pmu)
    df_conv_hyb = DataFrame(d_conv_hyb)
    df_conv_scada = DataFrame(d_conv_scada)
    df_conv_pmu[!, :se].= "pmu"
    df_conv_hyb[!, :se].= "hyb"
    df_conv_scada[!, :se].= "scada"
    df_conv = vcat(df_conv_pmu, df_conv_hyb, df_conv_scada)
    df_priors_all=vcat(dfs_priors_all...)
    return df_errors_total, df_conv, df_priors_all
end

function run_fase_experiment_stationary_timebench(N, n0, a, min_sigma, result, data_hyb_sr, d_prec, d_keys, d_meas_set, nlp_optimizer_hyb, nlp_optimizer_scada, nlp_optimizer_pmu, ipopt_files;
    se_objective=Dict("scada"=>"rwls","pmu"=>"rwls" ,"hyb"=>"rwls"),cutoff::Float64=1e-6)

    scada_measurements = ["scada", "dc", "conv"]
    d_meas_set_scada = Dict()
    for key in scada_measurements
        d_meas_set_scada[key] = deepcopy(d_meas_set[key])
    end

    pmu_measurements = ["pmu", "dtu"]
    d_meas_set_pmu = Dict()
    for key in pmu_measurements
        d_meas_set_pmu[key] = deepcopy(d_meas_set[key])
    end



    d_conv_scada = Dict(
        "termination_status" => Vector{Any}(undef, N),
        "n_iters" => Vector{Any}(undef, N),
        "solve_time" => Vector{Any}(undef, N),
        "P_time" => Vector{Any}(undef, N),
    )

    d_conv_pmu = Dict(
        "termination_status" => Vector{Any}(undef, N),
        "n_iters" => Vector{Any}(undef, N),
        "solve_time" => Vector{Any}(undef, N),
        "P_time" => Vector{Any}(undef, N),
    )

    d_conv_hyb = Dict(
        "termination_status" => Vector{Any}(undef, N),
        "n_iters" => Vector{Any}(undef, N),
        "solve_time" => Vector{Any}(undef, N),
        "P_time" => Vector{Any}(undef, N),
    )

    df_errors_total = DataFrame()

    dfs_priors_all = DataFrame[]

    for n in n0+1:n0+N
        data_hyb = deepcopy(data_hyb_sr) #copy data to avoid modifying original data

        introduce_noise!(data_hyb, d_prec, d_keys, seed=n, min_sig=min_sigma, bound=true)

        data_pmu = deepcopy(data_hyb)
        data_scada = deepcopy(data_hyb)


        split_meas_set!(d_keys, data_scada, ["scada", "dc", "conv"])
        split_meas_set!(d_keys, data_pmu, ["pmu", "dtu"]) #splits the measurements into SCADA and PMU sets;]


        data_scada["se_settings"]["criterion"] = se_objective["scada"]
        se_res_scada = _ACDCSE.solve_acdcse(data_scada, _PM.ACPPowerModel, nlp_optimizer_scada) #Solve SE
        n_iters = add_iterations_from_log!(ipopt_files["scada"])
        df_error_scada = df_error_var(result, se_res_scada)
        df_error_scada[!, :se] .= "scada"
        df_error_scada[!, :n] .= n


        d_conv_scada["termination_status"][n-n0] = string(se_res_scada["termination_status"])
        d_conv_scada["n_iters"][n-n0] = n_iters
        d_conv_scada["solve_time"][n-n0] = se_res_scada["solve_time"]


        data_hyb["se_settings"]["criterion"] = se_objective["hyb"]
        se_res_hyb = _ACDCSE.solve_acdcse(data_hyb, _PM.ACPPowerModel, nlp_optimizer_hyb) #Solve WLS
        n_iters = add_iterations_from_log!(ipopt_files["hyb"])
        d_conv_hyb["n_iters"][n-n0] = n_iters
        df_error_hyb = df_error_var(result, se_res_hyb)
        df_error_hyb[!, :se] .= "hyb"
        df_error_hyb[!, :n] .= n
        d_conv_hyb["termination_status"][n-n0] = string(se_res_hyb["termination_status"])
        d_conv_hyb["solve_time"][n-n0] = se_res_hyb["solve_time"]



        t = @belapsed generate_data_prior!($se_res_scada, $data_scada, $data_pmu, $d_keys, d=$a, min_sig=$min_sigma, print_info=true, cutoff=$cutoff) #Generate prior data for FASE
        d_conv_pmu["P_time"][n-n0] = t
        d_conv_hyb["P_time"][n-n0] = t
        d_conv_scada["P_time"][n-n0] = t


        data_pmu["se_settings"]["criterion"] = se_objective["pmu"]
        se_res_pmu = _ACDCSE.solve_acdcse(data_pmu, _PM.ACPPowerModel, nlp_optimizer_pmu) #Solve FASE
        df_prior_internal = extract_prior_dataframe(data_pmu)
        df_prior_internal[!, :n] .= n

        dfs_priors_all = push!(dfs_priors_all, df_prior_internal)

        n_iters = add_iterations_from_log!(ipopt_files["pmu"])
        
        df_error_pmu = df_error_var(result, se_res_pmu)
        df_error_pmu[!, :se] .= "pmu"
        df_error_pmu[!, :n] .= n
        
        df_errors_total = vcat(df_errors_total, df_error_pmu, df_error_hyb, df_error_scada)
        
        d_conv_pmu["termination_status"][n-n0] = string(se_res_pmu["termination_status"])
        d_conv_pmu["n_iters"][n-n0] = n_iters
        d_conv_pmu["solve_time"][n-n0] = se_res_pmu["solve_time"]

    end

    df_conv_pmu = DataFrame(d_conv_pmu)
    df_conv_hyb = DataFrame(d_conv_hyb)
    df_conv_scada = DataFrame(d_conv_scada)
    df_conv_pmu[!, :se].= "pmu"
    df_conv_hyb[!, :se].= "hyb"
    df_conv_scada[!, :se].= "scada"
    df_conv = vcat(df_conv_pmu, df_conv_hyb, df_conv_scada)
    df_priors_all=vcat(dfs_priors_all...)
    return df_errors_total, df_conv, df_priors_all
end



function run_fase_experiment_stationary_print_meas(N, n0, a, min_sigma, result, data_hyb_sr, d_prec, d_keys, d_meas_set, nlp_optimizer_hyb, nlp_optimizer_scada, nlp_optimizer_pmu, ipopt_files;
    se_objective=Dict("scada"=>"rwls","pmu"=>"rwls" ,"hyb"=>"rwls"),cutoff::Float64=1e-6)

    scada_measurements = ["scada", "dc", "conv"]
    d_meas_set_scada = Dict()
    for key in scada_measurements
        d_meas_set_scada[key] = deepcopy(d_meas_set[key])
    end

    pmu_measurements = ["pmu", "dtu"]
    d_meas_set_pmu = Dict()
    for key in pmu_measurements
        d_meas_set_pmu[key] = deepcopy(d_meas_set[key])
    end





    dfs_meas_set_scada = DataFrame[]
    dfs_meas_set_hyb = DataFrame[]
    dfs_meas_set_pmu = DataFrame[]

    for n in n0+1:n0+N
        data_hyb = deepcopy(data_hyb_sr) #copy data to avoid modifying original data

        introduce_noise!(data_hyb, d_prec, d_keys, seed=n, min_sig=min_sigma, bound=true)

        data_pmu = deepcopy(data_hyb)
        data_scada = deepcopy(data_hyb)


        split_meas_set!(d_keys, data_scada, ["scada", "dc", "conv"])
        split_meas_set!(d_keys, data_pmu, ["pmu", "dtu"]) #splits the measurements into SCADA and PMU sets;]


        data_scada["se_settings"]["criterion"] = se_objective["scada"]
        df_meas_set_scada = meas_dict_to_dataframe(data_scada["meas"])
        df_meas_set_scada[!, :n] .= n
        push!(dfs_meas_set_scada, df_meas_set_scada)


        df_meas_set_hyb = meas_dict_to_dataframe(data_hyb["meas"])
        df_meas_set_hyb[!, :n] .= n
        push!(dfs_meas_set_hyb, df_meas_set_hyb)

       

       
        data_pmu["se_settings"]["criterion"] = se_objective["pmu"]
        df_meas_set_pmu = meas_dict_to_dataframe(data_pmu["meas"])
        df_meas_set_pmu[!, :n] .= n
        push!(dfs_meas_set_pmu, df_meas_set_pmu)

    end

    df_scada = vcat(dfs_meas_set_scada...)
    df_hyb = vcat(dfs_meas_set_hyb...)
    df_pmu = vcat(dfs_meas_set_pmu...)  
    return df_scada, df_hyb, df_pmu
end





function add_iterations_from_log!()
    logtxt = read("ipopt_fase.out", String)
    return match(r"Number of Iterations\.*:\s*(\d+)", logtxt)[1]
end

function add_iterations_from_log!(txt)
    logtxt = read(txt, String)
    return match(r"Number of Iterations\.*:\s*(\d+)", logtxt)[1]
end


function modify_all_loads_randomly!(data_pf_all,num_buses,seed,ratio)
    bus2mod = range(1, stop=num_buses, step=1)
    load2mod = determine_loads_to_modify(bus2mod, data_pf_all[1])
    modify_loads_over_time_rand!(data_pf_all, load2mod; seed=seed, ratio=ratio)
end


function modify_loads_sudden_change!(data_pf_all, percent_of_buses_to_modify, num_buses, seed, T)
    num_integers = Int(round(num_buses * percent_of_buses_to_modify / 100))
    range_start = 1
    range_end = num_buses
    _ACDCSE._RAN.seed!(seed)  # Set seed for reproducibility
    buses_2_mod = rand(range_start:range_end, num_integers)
    instant_load_mod = rand(0:T, num_integers)
    percent_mod = rand(-15:15, num_integers)

    for (i, bus) in enumerate(buses_2_mod)
        bus2mod = [bus]
        t_i = instant_load_mod[i]
        P = percent_mod[i]
        load2mod = determine_loads_to_modify(bus2mod, data_pf_all[1])
        modify_loads_over_time_ramp!(data_pf_all, t_i, 3, load2mod, P)
    end
    return buses_2_mod,instant_load_mod,percent_mod
end


function set_converter_constraints!(data;Pmax_factor::Float64=1.2,Pmin_factor::Float64=1.2, Vmax_factor::Float64=1.05,Vmin_factor::Float64=0.95)
    for (key, conv) in data["convdc"]
        for (i, Pacrated) in enumerate(data["convdc"][key]["Pacrated"])
            data["convdc"][key]["Pacrated"][i] = Pmax_factor * Pacrated
        end
        for (i, Qacrated) in enumerate(data["convdc"][key]["Qacrated"])
            data["convdc"][key]["Qacrated"][i] = Pmax_factor * Qacrated
        end
        for (i, Pacmax) in enumerate(data["convdc"][key]["Pacmax"])
            data["convdc"][key]["Pacmax"][i] = Pmax_factor * Pacmax
        end
        for (i, Qacmax) in enumerate(data["convdc"][key]["Qacmax"])
            data["convdc"][key]["Qacmax"][i] = Pmax_factor * Qacmax
        end
        for (i, Pacmin) in enumerate(data["convdc"][key]["Pacmin"])
            data["convdc"][key]["Pacmin"][i] = Pmin_factor * Pacmin
        end
        for (i, Qacmin) in enumerate(data["convdc"][key]["Qacmin"])
            data["convdc"][key]["Qacmin"][i] = Pmin_factor * Qacmin
        end
        for (i, Vmmax) in enumerate(data["convdc"][key]["Vmmax"])
            data["convdc"][key]["Vmmax"][i] = Vmax_factor * Vmmax
        end
        for (i, Vmmin) in enumerate(data["convdc"][key]["Vmmin"])
            data["convdc"][key]["Vmmin"][i] = Vmin_factor * Vmmin
        end
    end
end
