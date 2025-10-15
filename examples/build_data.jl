mutable struct Meas_se
    var::Symbol
    cmp::Symbol
    cmp_id::Int
    direction::Symbol
    conn::Int
    z::Float64
    z_label::String
    σ::Float64
    prec::Float64
    key::String
end



mutable struct Convcons_se
    cmp_id::String
    conn::Int
    g_label::String
    σ::Float64
end

struct State_variable
    var::Symbol
    cmp::Symbol
    cmp_id::Int64
    conn::Int
    label::String
end

"""
this function returns the bus for load (d) or generator (g)
"""
function get_bus_d_o_g(data, var, cmp_id)

    if var == :pd || var == :qd
        return data["load"]["$cmp_id"]["load_bus"]
    else
        return data["gen"]["$cmp_id"]["gen_bus"]
    end
end

"""
This function gets the dc bus conection of a converter or a branch  
"""
function get_ivdcs(config,connected_at)
    if config == 2
        i_vdcs = [1,2,3]
    elseif connected_at == 0
        i_vdcs = [1,2]
    elseif connected_at == 1
        i_vdcs = [1,3]
    elseif connected_at == 2
        i_vdcs = [2,3]
    end
    return i_vdcs
end

"""
calculates the liquid active power injection for a bus"
"""
function calc_liq_P(data, respc_bus)
    P = 0
    for (key, meas) in data["meas"]
        if meas["var"] == :pd || meas["var"] == :pg 
            cmp_id = meas["cmp_id"]
            bus = get_bus_d_o_g(data, meas["var"], cmp_id)
            if bus == respc_bus
                if meas["var"] == :pd
                    P -= meas["dst"][1].μ
                else # pg
                    P += meas["dst"][1].μ   
                end
            end
        end
    end
    return P
end

"""
"calculates the liquid reactive power injection for a bus"
"""
function calc_liq_Q(data, respc_bus)

    Q = 0
    for (key, meas) in data["meas"]
        if meas["var"] == :qd || meas["var"] == :qg
            cmp_id = meas["cmp_id"]
            bus = get_bus_d_o_g(data, meas["var"], cmp_id)
            if bus == respc_bus
                if meas["var"] == :qd
                    Q -= meas["dst"][1].μ
                else # qg
                    Q += meas["dst"][1].μ
                end
            end
        end
    end
    return Q
end


"""
This function builds a dictionary with bus information for calculation of the literal Jacobian
"""
function build_bus_calc(data)

    bus_calc = Dict{String, Any}()

    for (key, bus) in data["bus"]
        bus_calc[key] = Dict{String, Any}(
        "bus_i" => bus["bus_i"],
        "branch_fr" => Vector{String}(),
        "branch_to" => Vector{String}(),
        "shunt" => Vector{String}(),
        "gs" => 0,
        "bs" => 0,
        "load" => Vector{String}(),
        "gen" => Vector{String}(),
        "convdc" => Vector{String}(),
        )
    end
    for (key, bus) in data["bus"]
        bus_calc[key]["bus_i"] = bus["bus_i"]
    end
    for (key, load) in data["load"]
        push!(bus_calc["$(load["load_bus"])"]["load"], key) 
    end
    for (key, gen) in data["gen"]
        push!(bus_calc["$(gen["gen_bus"])"]["gen"], key)
    end

    for (key, branch) in data["branch"]
        push!(bus_calc["$(branch["f_bus"])"]["branch_fr"], key)
        push!(bus_calc["$(branch["t_bus"])"]["branch_to"], key)
    end
    if "convdc" ∈ keys(data)
        for (key, convdc) in data["convdc"]
            push!(bus_calc["$(convdc["busac_i"])"]["convdc"], key)
        end
    end

    for (key, shunt) in data["shunt"]
        push!(bus_calc["$(shunt["shunt_bus"])"]["shunt"], key)
        bus_calc["$(shunt["shunt_bus"])"]["gs"] += shunt["gs"]
        bus_calc["$(shunt["shunt_bus"])"]["bs"] += shunt["bs"]
    end

    return bus_calc
end


"""
This function builds a dictionary with branch information for calculation of the literal Jacobian
"""
function build_branch_calc(data)

    branch_calc = Dict{String, Any}()
    for (key, branch) in data["branch"]
        branch_calc[key] = Dict{String, Any}()
        br_r = branch["br_r"]
        br_x = branch["br_x"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        br_g = real(1 / complex(br_r, br_x))
        br_b = imag(1 / complex(br_r, br_x))
        tap = 1 / branch["tap"]
        branch_calc[key]["f_bus"] = "$(branch["f_bus"])"
        branch_calc[key]["t_bus"] = "$(branch["t_bus"])"
        branch_calc[key]["gii"] = g_to + (tap^2) * br_g
        branch_calc[key]["gij"] = -(tap) * br_g
        branch_calc[key]["gji"] = -(tap) * br_g
        branch_calc[key]["gjj"] = g_fr + br_g
        branch_calc[key]["bii"] = b_to + (tap^2) * br_b
        branch_calc[key]["bij"] = -(tap) * br_b
        branch_calc[key]["bji"] = -(tap) * br_b
        branch_calc[key]["bjj"] = b_fr + br_b
    end
    return branch_calc
end

"""
This function builds a dictionary with converter information for calculation of the literal Jacobian
"""
function build_conv_calc(data)

    conv_calc = Dict{String, Any}() # Dictionary to hold converter calculations
    if !("convdc" ∈ keys(data))
        return conv_calc
    end
    for (key, conv) in data["convdc"]
        conv_calc[key] = Dict{String, Any}()
        conv_calc[key]["conv_i"]=key
        conv_calc[key]["busac_i"] = conv["busac_i"]
        conv_calc[key]["busdc_i"] = conv["busdc_i"]

        c = conv["conv_confi"] == 2 ? 2 : 1
        conv_calc[key]["conn"] = c

        conv_calc[key]["vdcp"] = Vector{Int}(undef, c)
        conv_calc[key]["vdcn"] = Vector{Int}(undef, c)
        conv_calc[key]["LossA"] = conv["LossA"]
        conv_calc[key]["LossB"] = conv["LossB"]
        conv_calc[key]["LossC"] = conv["LossCinv"] # TODO check with Marta, the code uses the inverter losses when they are different 
        
        i_vdcs=get_ivdcs(conv["conv_confi"], conv["connect_at"])
        
        if conv["conv_confi"]==2
            conv_calc[key]["vdcp"][1] = i_vdcs[1]
            conv_calc[key]["vdcn"][1] = i_vdcs[3]
            conv_calc[key]["vdcp"][2] = i_vdcs[3]
            conv_calc[key]["vdcn"][2] = i_vdcs[2]
        elseif 1 in i_vdcs
            conv_calc[key]["vdcp"][1] = i_vdcs[1]
            conv_calc[key]["vdcn"][1] = i_vdcs[2]
        else
            conv_calc[key]["vdcp"][1] = i_vdcs[2]
            conv_calc[key]["vdcn"][1] = i_vdcs[1]
        end 

        if conv["transformer"] == 1
            conv_calc[key]["tf"]= Vector{Dict{String, Any}}(undef, c)
            for i in 1:c
                tap=1/conv["tm"][i]
                br_r = conv["rtf"][i]
                br_x = conv["xtf"][i]
                br_g = real(1 / complex(br_r, br_x))
                br_b = imag(1 / complex(br_r, br_x))
                conv_calc[key]["tf"][i] = Dict{String, Any}()
                conv_calc[key]["tf"][i]["f_bus"] = "$(conv["busac_i"])"
                if conv["reactor"]==1
                    conv_calc[key]["tf"][i]["t_bus"] = "$key-f-$i"
                else
                    conv_calc[key]["tf"][i]["t_bus"] = "$key-c-$i"
                end
                conv_calc[key]["tf"][i]["gii"] = (tap^2) * br_g
                conv_calc[key]["tf"][i]["gij"] = -(tap) * br_g
                conv_calc[key]["tf"][i]["gji"] = -(tap) * br_g
                conv_calc[key]["tf"][i]["gjj"] = br_g
                conv_calc[key]["tf"][i]["bii"] = (tap^2) * br_b
                conv_calc[key]["tf"][i]["bij"] = -(tap) * br_b
                conv_calc[key]["tf"][i]["bji"] = -(tap) * br_b
                conv_calc[key]["tf"][i]["bjj"] = br_b
            end

        end
        if conv["reactor"]==1
            conv_calc[key]["pr"]= Vector{Dict{String, Any}}(undef, c)
            for i in 1:c
                br_r = conv["rc"][i]
                br_x = conv["xc"][i]
                br_g = real(1 / complex(br_r, br_x))
                br_b = imag(1 / complex(br_r, br_x))
                conv_calc[key]["pr"][i] = Dict{String, Any}()
                if conv["transformer"]==1
                    conv_calc[key]["pr"][i]["f_bus"] = "$key-f-$i"
                else
                    conv_calc[key]["pr"][i]["f_bus"] = "$(conv["busac_i"])"
                end
                conv_calc[key]["pr"][i]["t_bus"] = "$key-c-$i"
                conv_calc[key]["pr"][i]["gii"] =  br_g
                conv_calc[key]["pr"][i]["gij"] = - br_g
                conv_calc[key]["pr"][i]["gji"] = - br_g
                conv_calc[key]["pr"][i]["gjj"] = br_g
                conv_calc[key]["pr"][i]["bii"] =  br_b
                conv_calc[key]["pr"][i]["bij"] = - br_b
                conv_calc[key]["pr"][i]["bji"] = - br_b
                conv_calc[key]["pr"][i]["bjj"] = br_b
            end        
        end
        if conv["filter"]==1
            conv_calc[key]["f"] = Vector{Dict{String, Number}}(undef, c)
            for i in 1:c
                conv_calc[key]["f"][i] = Dict{String, Number}()
                conv_calc[key]["f"][i]["bf"] = conv["bf"][i]
            end
        end
    end
    return conv_calc
end




 """
 This function adds the bus voltages in the state vector `x` and creates rulers for the bus voltages.
 It also adds the bus voltage labels to `x_labels`, for mapping the x vector
 it modifies the `x` and `x_labels` vectors in place.
 It returns two dictionaries: `vm_ruler` and `va_ruler`, which map bus IDs to their respective indices in `x`.
 """
function add_x_and_rulers_bus!(x, x_info, data, res)

    vm_ruler = Dict{String, Int}()
    va_ruler = Dict{String, Int}()
    

    for key in sort(collect(keys(data["bus"])))
        push!(x, res["solution"]["bus"][key]["vm"])
        i = length(x)
        b=data["bus"][key]["bus_i"]
        vm_ruler["$b"] = i
        push!(x_info, State_variable(:vm,:bus,data["bus"][key]["index"],1,"vm_$(key)")) 
    end

    for key in sort(collect(keys(data["bus"])))
        push!(x, res["solution"]["bus"][key]["va"])
        i = length(x)
        b=data["bus"][key]["bus_i"]
        va_ruler["$b"] = i
        push!(x_info, State_variable(:va,:bus,data["bus"][key]["index"],1,"va_$(key)"))
    end

    return vm_ruler, va_ruler
end



"""
This function adds the converter bus voltages in the state vector `x` and creates rulers for the bus voltages.
It also adds the bus voltage labels to `x_labels`, for mapping the x vector
it modifies the `x` and `x_labels` vectors in place.
It returns two dictionaries: `vm_ruler` and `va_ruler`, which map converter bus IDs to their respective indices in `x`.
Converter buses ids are a string composed by three fields: 'converter_id-(internal bus_type)-(converter station)
internal bus_type can be 'f' for filter, 'c' for converter
converter station can be 1 or 2, depending on the converter configuration
"""
function add_x_and_rulers_conv!(x, x_info, data, res)



    vm_ruler = Dict{String, Int}()
    va_ruler = Dict{String, Int}() 


    !("convdc" ∈ keys(data)) && return vm_ruler, va_ruler

    for key in sort(collect(keys(data["convdc"])))
        conn = data["convdc"][key]["conv_confi"] == 2 ? 2 : 1
        !haskey(res["solution"]["convdc"][key],"vmfilt") ? continue : nothing
        for c in 1:conn
            push!(x, res["solution"]["convdc"][key]["vmfilt"][c])
            i = length(x)
            vm_ruler["$key-f-$c"] = i
            push!(x_info, State_variable(:vmf, :convdc, data["convdc"][key]["index"], c, "vmf_$(key)-$(c)"))
        end
    end

    for key in sort(collect(keys(data["convdc"])))
        conn = data["convdc"][key]["conv_confi"] == 2 ? 2 : 1
        !haskey(res["solution"]["convdc"][key],"vafilt") ? continue : nothing
        for c in 1:conn
            push!(x, res["solution"]["convdc"][key]["vafilt"][c])
            i = length(x)
            va_ruler["$key-f-$c"] = i
            push!(x_info, State_variable(:vaf, :convdc, data["convdc"][key]["index"], c, "vaf_$(key)-$(c)"))
        end
    end

    
    for key in sort(collect(keys(data["convdc"])))
        conn = data["convdc"][key]["conv_confi"] == 2 ? 2 : 1
        !haskey(res["solution"]["convdc"][key],"vmconv") ? continue : nothing
        for c in 1:conn
            push!(x, res["solution"]["convdc"][key]["vmconv"][c])
            i = length(x)
            vm_ruler["$key-c-$c"] = i
            push!(x_info, State_variable(:vmc, :convdc, data["convdc"][key]["index"], c, "vmc_$(key)-$(c)"))
        end
    end

    for key in sort(collect(keys(data["convdc"])))
        conn = data["convdc"][key]["conv_confi"] == 2 ? 2 : 1
        !haskey(res["solution"]["convdc"][key],"vaconv") ? continue : nothing
        for c in 1:conn
            push!(x, res["solution"]["convdc"][key]["vaconv"][c])
            i = length(x)
            va_ruler["$key-c-$c"] = i
            push!(x_info, State_variable(:vac, :convdc, data["convdc"][key]["index"], c, "vac_$(key)-$(c)"))
        end
    end

    return vm_ruler, va_ruler
end



function build_branchdc_calc(data)
    branchdc_calc = Dict{String, Any}()

    if !("branchdc" ∈ keys(data))
        return branchdc_calc
    end

    for (key, branchdc) in data["branchdc"]
        branchdc_calc[key] = Dict{String, Any}()

        conn = branchdc["line_confi"] == 2 ? 3 : 2
        branchdc_calc[key]["conn"] = conn

        i_vdcs = get_ivdcs(branchdc["line_confi"], branchdc["connect_at"])
        branchdc_calc[key]["i_vdcs"] = i_vdcs
        branchdc_calc[key]["conn"] = conn
        branchdc_calc[key]["conductor"] = Vector{Dict{String, Any}}(undef, conn)

        for c in 1:conn
            branchdc_calc[key]["conductor"][c] = Dict{String, Any}()
            branchdc_calc[key]["conductor"][c]["f_busdc"] = "$(branchdc["fbusdc"])-dc-$(i_vdcs[c])"
            branchdc_calc[key]["conductor"][c]["t_busdc"] = "$(branchdc["tbusdc"])-dc-$(i_vdcs[c])"
            branchdc_calc[key]["conductor"][c]["r"] = branchdc["r"][c]
        end
        branchdc_calc[key]["return_z"] = branchdc["return_z"]
    end

    return branchdc_calc
end

function build_busdc_calc(data)

    busdc_calc = Dict{String, Any}()

    if !("busdc" ∈ keys(data))
        return busdc_calc
    end

    for (key, busdc) in data["busdc"]
        busdc_calc[key] = Dict{String, Any}(
        "busdc_i" => busdc["busdc_i"],
        "branchdc_fr" => Vector{String}(),
        "branchdc_to" => Vector{String}(),
        "convdc" => Vector{String}(),
        "conn" => 3,
        "i_vdcs" => [1,2,3] #TODO implement this for buses without the three voltages connected
        )
    end
    for (key, busdc) in data["busdc"]
        busdc_calc[key]["busdc_i"] = busdc["busdc_i"]
    end
    for (key, branchdc) in data["branchdc"]
        push!(busdc_calc["$(branchdc["fbusdc"])"]["branchdc_fr"], key)
        push!(busdc_calc["$(branchdc["tbusdc"])"]["branchdc_to"], key)
    end

    for (key, convdc) in data["convdc"]
        push!(busdc_calc["$(convdc["busdc_i"])"]["convdc"], key)
    end

    return busdc_calc
end

function add_x_and_ruler_dc!(x, x_info, data, res)

    vdc_ruler = Dict{String, Int}()
    !("busdc" ∈ keys(data)) && return vdc_ruler

    for key in sort(collect(keys(data["busdc"])))
        conn = 3
        for c in 1:conn
            push!(x, res["solution"]["busdc"][key]["vm"][c])
            i = length(x)
            vdc_ruler["$key-dc-$c"] = i
            push!(x_info, State_variable(:vdcm, :busdc, data["busdc"][key]["index"], c, "vdc_$(key)-$(c)"))
        end
    end

    return vdc_ruler

end

function build_z(data)

    z = Vector{Meas_se}()
    bus_P_injc_added = Vector{Int64}()
    bus_Q_injc_added = Vector{Int64}()

    for i in sort(parse.(Int,keys(data["meas"])))
        key="$(i)"
        meas=data["meas"][key]
        if meas["var"] in [:p,:q,:cr,:ci]
            push!(z, Meas_se(meas["var"],meas["cmp"],meas["cmp_id"][1],meas["direction"],1,meas["dst"][1].μ,"$(meas["var"])_$(meas["cmp_id"][2])-$(meas["cmp_id"][3])",meas["dst"][1].σ,1.0,key))
        elseif meas["var"] in [:vm,:va]
            push!(z, Meas_se(meas["var"],meas["cmp"],meas["cmp_id"][1],:na,1,meas["dst"][1].μ ,"$(meas["var"])_$(meas["cmp_id"][1])",meas["dst"][1].σ,1.0,key))
        elseif  meas["var"] in [:pd, :pg]  # add injection measurements
            cmp_id = meas["cmp_id"] # gets the load or gen id
            bus = get_bus_d_o_g(data, meas["var"], cmp_id) #gets the bus
            (bus in bus_P_injc_added) ? continue : nothing 
            #loads and gens meas need to be condensed into a single measurement
            push!(z, Meas_se(:pinj, :bus, bus, :na, 1, calc_liq_P(data, bus), "pinj_$(bus)",meas["dst"][1].σ,1.0,key))
            push!(bus_P_injc_added, bus)
        elseif meas["var"] in [:qd, :qg] # reactive power injections
            cmp_id = meas["cmp_id"] # gets the load or gen id
            bus = get_bus_d_o_g(data, meas["var"], cmp_id) #gets the bus
            (bus in bus_Q_injc_added) ? continue : nothing
            #loads and gens meas need to be condensed into a single measurement
            push!(z, Meas_se(:qinj, :bus, bus, :na, 1, calc_liq_Q(data, bus), "qinj_$(bus)",meas["dst"][1].σ,1.0,key))
            push!(bus_Q_injc_added, bus)
            # P is now a function that sums the results of all functions in h_int
        elseif meas["var"] in [:pconv_ac,:qconv_ac] #inverted sing
            for c in 1:length(meas["dst"])
                push!(z, Meas_se(meas["var"], :convdc, meas["cmp_id"], :na,c ,-meas["dst"][c].μ , "$(meas["var"])_$(meas["cmp_id"])-$(c)",meas["dst"][c].σ,1.0,key))
            end
        elseif meas["var"] in [:pconv_tf_to,:qconv_tf_to,:pconv_tf_fr,:qconv_tf_fr]
            for c in 1:length(meas["dst"])
                push!(z, Meas_se(meas["var"], :convdc, meas["cmp_id"], :na,c ,meas["dst"][c].μ , "$(meas["var"])_$(meas["cmp_id"])-$(c)",meas["dst"][c].σ,1.0,key))
            end
        elseif meas["var"] in [:pconv_pr_fr,:qconv_pr_fr]
            for c in 1:length(meas["dst"])
                push!(z, Meas_se(meas["var"], :convdc, meas["cmp_id"], :na,c ,meas["dst"][c].μ , "$(meas["var"])_$(meas["cmp_id"])-$(c)",meas["dst"][c].σ,1.0,key))
            end
        elseif meas["var"] in [:vmc,:vac,:vmf,:vaf]
            for c in 1:length(meas["dst"])
                push!(z, Meas_se(meas["var"], :convdc, meas["cmp_id"], :na,c ,meas["dst"][c].μ , "$(meas["var"])_$(meas["cmp_id"])-$(c)",meas["dst"][c].σ,1.0,key))
            end
        elseif meas["var"] in [:i_dcgrid, :p_dcgrid]
            for c in 1:length(meas["dst"])
                push!(z, Meas_se(meas["var"], :branchdc, meas["cmp_id"][1], meas["direction"],c ,meas["dst"][c].μ , "$(meas["var"])_$(meas["cmp_id"][2])-$(meas["cmp_id"][3])-$(c)",meas["dst"][c].σ,1.0,key))
            end
        elseif meas["var"] in [:vdcm]
            for c in 1:length(meas["dst"])
                push!(z, Meas_se(meas["var"], :busdc, meas["cmp_id"], :na,c ,meas["dst"][c].μ , "$(meas["var"])_$(meas["cmp_id"])-$(c)",meas["dst"][c].σ,1.0,key))
            end
        elseif meas["var"] == :mconv
            for c in 1:length(meas["dst"])
                push!(z, Meas_se(meas["var"], :convdc, meas["cmp_id"], :na,c ,meas["dst"][c].μ , "$(meas["var"])_$(meas["cmp_id"])-$(c)",meas["dst"][c].σ,1.0,key))
            end
        end

    end

    for (key,bus) in data["bus"]
        if bus["bus_type"] == 3
            push!(z, Meas_se(:va,:bus,bus["bus_i"],:na,1,0.0 ,"ref_$(key)",1e-4,1.0,"-1"))
        end
    end
    return z
end


function build_g(conv_calc)
    g = Vector{Convcons_se}()

    for conv in keys(conv_calc)
        for conn in 1:conv_calc[conv]["conn"]
            push!(g, Convcons_se(conv, conn, "conv_$(conv)_conn_$(conn)", 1.0))
        end
    end

    return g
end

function build_g_inj(conv_calc)
    g = Vector{Convcons_se}()

    for conv in keys(conv_calc)
        if haskey(conv_calc[conv],"pr") && haskey(conv_calc[conv],"tf")
            for conn in 1:conv_calc[conv]["conn"]
                push!(g, Convcons_se(conv, conn, "conv_$(conv)_conn_$(conn)", 1.0))
            end
        end
    end

    return g
end


function build_g_Pinj(conv_calc)
    g = Vector{Convcons_se}()

    for conv in keys(conv_calc)
        if haskey(conv_calc[conv],"pr") && haskey(conv_calc[conv],"tf")
            for conn in 1:conv_calc[conv]["conn"]
                push!(g, Convcons_se(conv, conn, "Pconv_$(conv)_conn_$(conn)", 1.0))
            end
        end
    end

    return g
end



function build_g_Qinj(conv_calc)
    g = Vector{Convcons_se}()

    for conv in keys(conv_calc)
        if haskey(conv_calc[conv],"pr") && haskey(conv_calc[conv],"tf")
            for conn in 1:conv_calc[conv]["conn"]
                push!(g, Convcons_se(conv, conn, "Qconv_$(conv)_conn_$(conn)", 1.0))
            end
        end
    end

    return g
end




function build_network_info(data, result)
    x = Vector{Float64}()
    x_info = Vector{State_variable}()

    vm_ruler_ac, va_ruler_ac = add_x_and_rulers_bus!(x,x_info , data, result)

    vm_ruler_conv, va_ruler_conv = add_x_and_rulers_conv!(x, x_info, data, result)
    vdc_ruler = add_x_and_ruler_dc!(x, x_info, data, result)

    # Merge bus and converter rulers into unified dictionaries
    vm_ruler = merge(vm_ruler_ac, vm_ruler_conv)
    va_ruler = merge(va_ruler_ac, va_ruler_conv)

    branch_calc = build_branch_calc(data)
    bus_calc = build_bus_calc(data)
    conv_calc = build_conv_calc(data)
    branchdc_calc = build_branchdc_calc(data)
    busdc_calc = build_busdc_calc(data)

    z = build_z(data)

    # z_virtual = build_z_virtual(data)    
    g = build_g(conv_calc)

    #c = create_cons(data)

    return Dict(
        "x" => x,
        "x_info" => x_info,
        "vm_ruler" => vm_ruler,
        "va_ruler" => va_ruler,
        "vdc_ruler" => vdc_ruler,
        "branch_calc" => branch_calc,
        "bus_calc" => bus_calc,
        "conv_calc" => conv_calc,
        "branchdc_calc" => branchdc_calc,
        "busdc_calc" => busdc_calc,
        "z" => z,
        # "z_virtual" => z_virtual,
        "g" => g
    )
end



function build_network_info_bad_data(data, result)
    x = Vector{Float64}()
    x_info = Vector{State_variable}()

    vm_ruler_ac, va_ruler_ac = add_x_and_rulers_bus!(x,x_info , data, result)

    vm_ruler_conv, va_ruler_conv = add_x_and_rulers_conv!(x, x_info, data, result)
    vdc_ruler = add_x_and_ruler_dc!(x, x_info, data, result)

    # Merge bus and converter rulers into unified dictionaries
    vm_ruler = merge(vm_ruler_ac, vm_ruler_conv)
    va_ruler = merge(va_ruler_ac, va_ruler_conv)

    branch_calc = build_branch_calc(data)
    bus_calc = build_bus_calc(data)
    conv_calc = build_conv_calc(data)
    branchdc_calc = build_branchdc_calc(data)
    busdc_calc = build_busdc_calc(data)

    z = build_z(data)

    # z_virtual = build_z_virtual(data)    
    g = build_g(conv_calc)

    g_Pinj = build_g_Pinj(conv_calc)
    g_Qinj = build_g_Qinj(conv_calc)

    return Dict(
        "x" => x,
        "x_info" => x_info,
        "vm_ruler" => vm_ruler,
        "va_ruler" => va_ruler,
        "vdc_ruler" => vdc_ruler,
        "branch_calc" => branch_calc,
        "bus_calc" => bus_calc,
        "conv_calc" => conv_calc,
        "branchdc_calc" => branchdc_calc,
        "busdc_calc" => busdc_calc,
        "z" => z,
        # "z_virtual" => z_virtual,
        "g" => g,
        "g_Pinj" => g_Pinj,
        "g_Qinj" => g_Qinj
    )
end


function build_network_info_posterior(data, result)
    x = Vector{Float64}()
    x_info = Vector{State_variable}()

    vm_ruler_ac, va_ruler_ac = add_x_and_rulers_bus!(x,x_info , data, result)

    vm_ruler_conv, va_ruler_conv = add_x_and_rulers_conv!(x, x_info, data, result)
    vdc_ruler = add_x_and_ruler_dc!(x, x_info, data, result)

    # Merge bus and converter rulers into unified dictionaries
    vm_ruler = merge(vm_ruler_ac, vm_ruler_conv)
    va_ruler = merge(va_ruler_ac, va_ruler_conv)

    branch_calc = build_branch_calc(data)
    bus_calc = build_bus_calc(data)
    conv_calc = build_conv_calc(data)
    branchdc_calc = build_branchdc_calc(data)
    busdc_calc = build_busdc_calc(data)

    z, z_keys = build_z_posterior_with_keys(data)


    return Dict(
        "x" => x,
        "x_info" => x_info,
        "vm_ruler" => vm_ruler,
        "va_ruler" => va_ruler,
        "vdc_ruler" => vdc_ruler,
        "branch_calc" => branch_calc,
        "bus_calc" => bus_calc,
        "conv_calc" => conv_calc,
        "branchdc_calc" => branchdc_calc,
        "busdc_calc" => busdc_calc,
        "z" => z,
        "z_keys" => z_keys,
        # "z_virtual" => z_virtual,

    )
end



caculate_fx(h,x) = [f(x) for f in h]


function create_W_autodiff(network_info)
    N_z=length(network_info["z"])
    haskey(network_info,"z_virtual") ? N_zv=length(network_info["z_virtual"]) : N_zv=0
    haskey(network_info,"g") ? N_g=length(network_info["g"]) : N_g=0
    N=N_z+N_zv+N_g
    
    σ = Vector{Float64}(undef, N)
    
    
    for i in 1:N_z
        σ[i] = network_info["z"][i].σ
    end
    if N_zv > 0
        for i in 1:N_zv
            σ[i + N_z] = network_info["z_virtual"][i].σ
        end
    end
    if N_g > 0
        for i in 1:N_g
            σ[i + N_z + N_zv] = network_info["g"][i].σ
        end
    end

    network_info["σ"] = σ
    return diagm(σ.^(-2))
end


function create_W_autodiff_eg(network_info)
    N_z=length(network_info["z"])
    haskey(network_info,"z_virtual") ? N_zv=length(network_info["z_virtual"]) : N_zv=0
    haskey(network_info,"g") ? N_g=length(network_info["g"]) : N_g=0
    haskey(network_info,"g_Pinj") ? N_gPinj=length(network_info["g_Pinj"]) : N_gPinj=0
    haskey(network_info,"g_Qinj") ? N_gQinj=length(network_info["g_Qinj"]) : N_gQinj=0
     
    N=N_z+N_zv+N_g +N_gPinj + N_gQinj
    
    σ = Vector{Float64}(undef, N)
    
    
    for i in 1:N_z
        σ[i] = network_info["z"][i].σ
    end
    if N_zv > 0
        for i in 1:N_zv
            σ[i + N_z] = network_info["z_virtual"][i].σ
        end
    end
    if N_g > 0
        for i in 1:N_g
            σ[i + N_z + N_zv] = network_info["g"][i].σ
        end
    end
    if N_gPinj > 0
        for i in 1:N_gPinj
            σ[i + N_z + N_zv + N_g] = network_info["g_Pinj"][i].σ
        end
    end
    if N_gQinj > 0
        for i in 1:N_gQinj
            σ[i + N_z + N_zv + N_g + N_gPinj] = network_info["g_Qinj"][i].σ
        end
    end 
        
    network_info["σ"] = σ
    return diagm(σ.^(-2))
end

function determine_virtual_buses(data)
    load_buses = Int64[]
    for (key, load) in data["load"]
        push!(load_buses, load["load_bus"])
    end

    gen_buses = Int64[]
    for (key, gen) in data["gen"]
        push!(gen_buses, gen["gen_bus"])
    end
    all_buses = collect(parse.(Int64, keys(data["bus"])))

    virtual_buses = setdiff(all_buses,union(load_buses,gen_buses))
    return virtual_buses
end

function build_z_virtual(data)
    virtual_buses = determine_virtual_buses(data)
    z_virtual_p = Vector{Meas_se}(undef, length(virtual_buses))
    z_virtual_q = Vector{Meas_se}(undef, length(virtual_buses))
    for (n, bus) in enumerate(virtual_buses)
        z_virtual_p[n]=Meas_se(:pinj, :bus, bus, :na, 1, 0.0, "pinj_$(bus)",1e-5,1.0,"-1")
        z_virtual_q[n]=Meas_se(:qinj, :bus, bus, :na, 1, 0.0, "qinj_$(bus)",1e-5,1.0,"-1")
    end
    return vcat(z_virtual_p, z_virtual_q)
end
  
