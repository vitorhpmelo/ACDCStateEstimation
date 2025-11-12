
function modfify_loads_fp!(time_steps, load_data, data_pf, data_pfs, data_se, data_ses)
    buses=sort(unique(load_data[!,:bus]))
    for t in time_steps
        data_pfs[t] = deepcopy(data_pf)
        data_ses[t] = deepcopy(data_se)
        for bus in buses
            load_data_t_bus = load_data[(load_data[!,:t].==t) .& (load_data[!,:bus].==bus),:]
            if nrow(load_data_t_bus) > 0
                load = determine_loads_to_modify([bus], data_pf)[1]
                data_pfs[t]["load"][string(load)]["pd"] = load_data_t_bus[1,:var] * data_pfs[t]["load"][string(load)]["pd"]
                data_ses[t]["load"][string(load)]["qd"] = load_data_t_bus[1,:var] * data_ses[t]["load"][string(load)]["qd"]
            end
        end
    end
end

function modify_gen_fp(time_steps,gen_data,data_pf,data_pfs,data_se,data_ses)
    buses=sort(unique(gen_data[!,:bus]))   
    for t in time_steps
        data_pfs[t] = deepcopy(data_pf)
        data_ses[t] = deepcopy(data_se)
        for bus in buses
            gen_data_t_bus=gen_data[(gen_data[!,:t].==t) .& (gen_data[!,:bus].==bus),:]
            if nrow(gen_data_t_bus)>0
                gen=determine_gens_to_modify([bus], data_pf)[1]
                data_pfs[t]["gen"][string(gen)]["pg"] = gen_data_t_bus[1,:var]*data_pfs[t]["gen"][string(gen)]["pg"]
            end
        end
    end
end


function create_vm!(data::Dict{String,Any}, res_t::Dict{String,Any}; prec::Float64=0.05, perc_virtual::Float64=1e-5)
    n=length(data["meas"])+1
    for (key, _) in data["bus"]
        data["meas"][string(n)] = Dict{String,Any}()
        m = data["meas"][string(n)]
        m["crit"]   = "rwls"
        m["cmp_id"] = parse(Int64, key)
        m["cmp"]    = :bus
        m["var"]    = :vm
        μ = res_t["solution"]["bus"][key]["vm"]
        σ = prec * μ / 3
        if σ < perc_virtual
            σ = perc_virtual
        end
        m["dst"]    = [_ACDCSE._DST.Normal(μ, σ)]
        n += 1
    end

end

function create_va!(data::Dict{String,Any}, res_t::Dict{String,Any}; prec::Float64=0.05, perc_virtual::Float64=1e-5)
    n = length(data["meas"])+1
    for (key, _) in data["bus"]
        data["meas"][string(n)] = Dict{String,Any}()
        m = data["meas"][string(n)]
        m["crit"]   = "rwls"
        m["cmp_id"] = parse(Int64, key)
        m["cmp"]    = :bus
        m["var"]    = :va
        μ = res_t["solution"]["bus"][key]["va"]
        σ = prec * μ / 3
        if σ < perc_virtual
            σ = perc_virtual
        end
        m["dst"]    = [_ACDCSE._DST.Normal(μ, σ)]
        n += 1
    end
end







function create_va!(data::Dict{String,Any}, res_t::Dict{String,Any}; prec::Float64=0.05, perc_virtual::Float64=1e-5)
    n = length(data["meas"])+1
    for (key, _) in data["bus"]
        data["meas"][string(n)] = Dict{String,Any}()
        m = data["meas"][string(n)]
        m["crit"]   = "rwls"
        m["cmp_id"] = parse(Int64, key)
        m["cmp"]    = :bus
        m["var"]    = :va
        μ = res_t["solution"]["bus"][key]["va"]
        σ = prec * μ / 3
        if σ < perc_virtual
            σ = perc_virtual
        end
        m["dst"]    = [_ACDCSE._DST.Normal(μ, σ)]
        n += 1
    end
end


function create_vmf!(data::Dict{String,Any}, res_t::Dict{String,Any}; prec::Float64=0.05, perc_virtual::Float64=1e-5)
    n=length(data["meas"])+1
    for (key, _) in data["convdc"]
        data["meas"][string(n)] = Dict{String,Any}()
        m = data["meas"][string(n)]
        m["crit"]   = "rwls"
        m["cmp_id"] = parse(Int64, key)
        m["cmp"]    = :convdc
        m["var"]    = :vmf
        μ_s = res_t["solution"]["convdc"][key]["vmfilt"]
        σ_s = Float64[]
        for μ in μ_s
            σ = prec * μ / 3
            if σ < perc_virtual
                σ = perc_virtual
            end
            push!(σ_s, σ)
        end
        m["dst"]    = [_ACDCSE._DST.Normal(μ, σ) for (μ, σ) in zip(μ_s, σ_s)]
        n += 1
    end
end

function create_vaf!(data::Dict{String,Any}, res_t::Dict{String,Any}; prec::Float64=0.05, perc_virtual::Float64=1e-5)
    n=length(data["meas"])+1
    for (key, _) in data["convdc"]
        data["meas"][string(n)] = Dict{String,Any}()
        m = data["meas"][string(n)]
        m["crit"]   = "rwls"
        m["cmp_id"] = parse(Int64, key)
        m["cmp"]    = :convdc
        m["var"]    = :vaf
        μ_s = res_t["solution"]["convdc"][key]["vafilt"]
        σ_s = Float64[]
        for μ in μ_s
            σ = prec * μ / 3
            if σ < perc_virtual
                σ = perc_virtual
            end
            push!(σ_s, σ)
        end
        m["dst"]    = [_ACDCSE._DST.Normal(μ, σ) for (μ, σ) in zip(μ_s, σ_s)]
        n += 1
    end
end

function create_vmc!(data::Dict{String,Any}, res_t::Dict{String,Any}; prec::Float64=0.05, perc_virtual::Float64=1e-5)
    n=length(data["meas"])+1
    for (key, _) in data["convdc"]
        data["meas"][string(n)] = Dict{String,Any}()
        m = data["meas"][string(n)]
        m["crit"]   = "rwls"
        m["cmp_id"] = parse(Int64, key)
        m["cmp"]    = :convdc
        m["var"]    = :vmc
        μ_s = res_t["solution"]["convdc"][key]["vmconv"]
        σ_s = Float64[]
        for μ in μ_s
            σ = prec * μ / 3
            if σ < perc_virtual
                σ = perc_virtual
            end
            push!(σ_s, σ)
        end
        m["dst"]    = [_ACDCSE._DST.Normal(μ, σ) for (μ, σ) in zip(μ_s, σ_s)]
        n += 1
    end
end

function create_vac!(data::Dict{String,Any}, res_t::Dict{String,Any}; prec::Float64=0.05, perc_virtual::Float64=1e-5)
    n=length(data["meas"])+1
    for (key, _) in data["convdc"]
        data["meas"][string(n)] = Dict{String,Any}()
        m = data["meas"][string(n)]
        m["crit"]   = "rwls"
        m["cmp_id"] = parse(Int64, key)
        m["cmp"]    = :convdc
        m["var"]    = :vac
        μ_s = res_t["solution"]["convdc"][key]["vaconv"]
        σ_s = Float64[]
        for μ in μ_s
            σ = prec * μ / 3
            if σ < perc_virtual
                σ = perc_virtual
            end
            push!(σ_s, σ)
        end
        m["dst"]    = [_ACDCSE._DST.Normal(μ, σ) for (μ, σ) in zip(μ_s, σ_s)]
        n += 1
    end
end


function create_vmdc!(data::Dict{String,Any}, res_t::Dict{String,Any}; prec::Float64=0.05, perc_virtual::Float64=1e-5)
    n=length(data["meas"])+1
    for (key, _) in data["busdc"]
        data["meas"][string(n)] = Dict{String,Any}()
        m = data["meas"][string(n)]
        m["crit"]   = "rwls"
        m["cmp_id"] = parse(Int64, key)
        m["cmp"]    = :busdc
        m["var"]    = :vdcm
        μ_s = res_t["solution"]["busdc"][key]["vm"]
        σ_s = Float64[]
        for μ in μ_s
            σ = prec * μ / 3
            if σ < perc_virtual
                σ = perc_virtual
            end
            push!(σ_s, σ)
        end
        m["dst"]    = [_ACDCSE._DST.Normal(μ, σ) for (μ, σ) in zip(μ_s, σ_s)]
        n += 1
    end
end