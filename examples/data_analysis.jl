function df_error_var(result_pf,result_se; filter_cmp_var=true,vars=["va","vm","vafilt","vmfilt","vaconv","vmconv"])


    dfout=DataFrame(var=[],comp=[],cmp_id=[],c=[],val_pf=[],val_se=[],error=[])



    "convdc" ∈ keys(result_se["solution"]) && for (key_conv, conv) in result_se["solution"]["convdc"]
        for (key_var,var) in conv
            if filter_cmp_var && key_var ∈ vars
                for i in range(1,length(var))
                    d_dummy=Dict()
                    d_dummy["var"]=key_var
                    d_dummy["comp"]="convdc"
                    d_dummy["cmp_id"]=key_conv
                    d_dummy["c"]=i
                    d_dummy["val_se"]=var[i]
                    d_dummy["val_pf"]=result_pf["solution"]["convdc"][key_conv][key_var][i]
                    d_dummy["error"]=abs(var[i]-result_pf["solution"]["convdc"][key_conv][key_var][i])
                    push!(dfout,d_dummy)
                end
            end
        end
    end



    for (key_bus, bus) in result_pf["solution"]["bus"]
        for (key_var,var) in bus
            if filter_cmp_var && key_var ∈ vars
                for i in range(1,length(var))
                    d_dummy=Dict()
                    d_dummy["var"]=key_var
                    d_dummy["comp"]="bus"
                    d_dummy["cmp_id"]=key_bus
                    d_dummy["c"]=i
                    d_dummy["val_pf"]=var[i]
                    d_dummy["val_se"]=result_se["solution"]["bus"][key_bus][key_var][i]
                    d_dummy["error"]=abs(var[i]-result_se["solution"]["bus"][key_bus][key_var][i])
                    push!(dfout,d_dummy)
                end
            end
        end
    end



    "busdc" ∈ keys(result_se["solution"]) && for (key_busdc, busdc) in result_pf["solution"]["busdc"]
        for (key_var,var) in busdc
            if filter_cmp_var && key_var ∈ vars
                for i in range(1,length(var))
                    d_dummy=Dict()
                    d_dummy["var"]=key_var
                    d_dummy["comp"]="busdc"
                    d_dummy["cmp_id"]=key_busdc
                    d_dummy["c"]=i
                    d_dummy["val_pf"]=var[i]
                    d_dummy["val_se"]=result_se["solution"]["busdc"][key_busdc][key_var][i]
                    d_dummy["error"]=abs(var[i]-result_se["solution"]["busdc"][key_busdc][key_var][i])
                    push!(dfout,d_dummy)
                end
            end
        end
    end


    return dfout
end


function filter_cpmandvar(cmp,var,cmps=["bus"],vars=["va","vm"]) ::Bool
    interesting_comp = cmp in cmps
    interesting_vars = var in vars
    interesting_comp && interesting_vars
end


function get_true_prior(cmp,var,cmp_id,conn,result)
    if cmp == :bus
        val = var==:vm ? result["solution"]["bus"]["$cmp_id"]["vm"] : result["solution"]["bus"]["$cmp_id"]["va"]
    elseif cmp == :convdc
        sub_cmp = occursin("f","$var") ? "filt" : "conv"
        sub_var = occursin("a","$var") ? "a" : "m"
        key_var="v"*sub_var*sub_cmp
        val=result["solution"]["convdc"]["$cmp_id"][key_var][conn]
    elseif cmp == :busdc
        val=result["solution"]["busdc"]["$cmp_id"]["vm"][conn]
    end
    return val
end




function print_network_info(network_info)
    open("x_info.csv", "w") do io
    println(io, "var,cmp,cmp_id,x_label")
        for (i, x_i) in enumerate(network_info["x_info"])
            println(io, "$(x_i.var),$(x_i.cmp),$(x_i.cmp_id),$(x_i.label)")
        end
    end

    open("z_info.csv", "w") do io
    println(io, "index,label,z,sigma")
        for (i, z_i) in enumerate(network_info["z"])
            println(io, "$i,$(z_i.z_label),$(z_i.z),$(z_i.σ)")
        end
    end


    # open("zvirtual_info.csv", "w") do io
    # println(io, "index,label,z,sigma")
    #     for (i, z_i) in enumerate(network_info["zvirtual"])
    #         println(io, "$i,$(z_i.z_label),$(z_i.z),$(z_i.σ)")
    #     end
    # end
    
    open("g_info.csv", "w") do io
    println(io, "index,label,sigma")
        for (i, g_i) in enumerate(network_info["g"])
            println(io, "$i,$(g_i.g_label),$(g_i.σ)")
        end
    end
    open("w_info.csv", "w") do io
    println(io, "index,sigma")
        for (i, g_i) in enumerate(network_info["σ"])
            println(io, "$i,$(g_i)")
        end
    end
    open("Jacobian.csv", "w") do io
        for i in 1:size(network_info["H"], 1)
            println(io, join(network_info["H"][i, :], ", "))
        end
    end

end


function get_true_prior(cmp,var,cmp_id,conn,result)
    if cmp == :bus
        val = var==:vm ? result["solution"]["bus"]["$cmp_id"]["vm"] : result["solution"]["bus"]["$cmp_id"]["va"]
    elseif cmp == :convdc
        sub_cmp = occursin("f","$var") ? "filt" : "conv"
        sub_var = occursin("a","$var") ? "a" : "m"
        key_var="v"*sub_var*sub_cmp
        val=result["solution"]["convdc"]["$cmp_id"][key_var][conn]
    elseif cmp == :busdc
        val=result["solution"]["busdc"]["$cmp_id"]["vm"][conn]
    end
    return val
end

function extract_prior_dataframe(data)
    y_vals = Dict()
    n = length(data["prior"])

    y_vals["aij"]=Vector{Any}(undef, n)
    y_vals["var_i"]=Vector{Any}(undef, n)
    y_vals["var_j"]=Vector{Any}(undef, n)
    y_vals["π_i"]=Vector{Float64}(undef, n)
    y_vals["π_j"]=Vector{Float64}(undef, n)
    y_vals["cmp_id_i"]=Vector{Any}(undef, n)
    y_vals["cmp_id_j"]=Vector{Any}(undef, n)
    y_vals["conn_i"]=Vector{Any}(undef, n)
    y_vals["conn_j"]=Vector{Any}(undef, n)
    y_vals["cmp_i"]=Vector{Any}(undef, n)
    y_vals["cmp_j"]=Vector{Any}(undef, n)
    y_vals["i_i"]=Vector{Any}(undef, n)
    y_vals["i_j"]=Vector{Any}(undef, n)
    y_vals["key"]=Vector{Any}(undef, n)

    i=1 
    for (key,prior) in data["prior"]
        y_vals["key"][i]=key
        y_vals["aij"][i] = data["prior"][key]["a_ij"]
        y_vals["var_i"][i] = data["prior"][key]["var_i"]
        y_vals["var_j"][i] = data["prior"][key]["var_j"]
        y_vals["π_i"][i] = data["prior"][key]["π_i"]
        y_vals["π_j"][i] = data["prior"][key]["π_j"]
        y_vals["cmp_id_i"][i] = data["prior"][key]["cmp_id_i"]
        y_vals["cmp_id_j"][i] = data["prior"][key]["cmp_id_j"]
        y_vals["i_i"][i] = data["prior"][key]["i_i"]
        y_vals["i_j"][i] = data["prior"][key]["i_j"]
        y_vals["conn_i"][i] = data["prior"][key]["conn_i"]
        y_vals["conn_j"][i] = data["prior"][key]["conn_j"]
        y_vals["cmp_i"][i] = data["prior"][key]["cmp_i"]
        y_vals["cmp_j"][i] = data["prior"][key]["cmp_j"]

        i += 1
    end
    return DataFrame(y_vals)
end

function add_iterations_from_log!()
    logtxt = read("ipopt_fase.out", String)
    return match(r"Number of Iterations\.*:\s*(\d+)", logtxt)
end


function meas_dict_to_dataframe(d_meas_set_se)
    df_meas_set_se = DataFrame(key=[], cmp=[], cmp_id=[], var=[], c=[], value=[], sigma=[])
    for (key, item) in pairs(d_meas_set_se)
        c = length(item["dst"])
        for i in 1:c
            df = DataFrame(
                key = key,
                cmp = item["cmp"],
                cmp_id = item["cmp_id"],
                var = item["var"],
                c = i,
                value = item["dst"][i].μ,
                sigma = item["dst"][i].σ
            )
            append!(df_meas_set_se, df)
        end
    end
    return df_meas_set_se
end
