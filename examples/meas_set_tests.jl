
function calc_sigma!(z_i::Meas_se, prec::Float64; min_sig::Float64=1e-4)
    z_i.prec=prec
    σ = abs(z_i.z) * z_i.prec/3

    z_i.σ = σ>min_sig ? σ : min_sig
end

function calc_sigma!(z_i::Meas_se; min_sig::Float64=1e-4)
    z_i.σ = min_sig
end


function calc_sigma!(g_i::Convcons_se; prec::Float64=1e-8)
    g_i.σ = prec
end



function create_dmeas_set(meas_set,data_pf,branch_types_ac=["p-scada","c-pmu"],branch_types_dc=["i_dcgrid-dc","p_dcgrid-dc"])

    d_meas_set=Dict()
    types = unique([split(type,"-")[2] for type in names(meas_set)])

    
    for type in types
        d_meas_set[type] = Dict()
    end

    for col in names(meas_set)
        (meas,type) = split(col, "-")
        d_meas_set[type][meas]=collect(skipmissing(meas_set[:,col]))
    end



    for element in branch_types_ac
        (meas,type) = split(element,"-")
        for i in range(1,length(d_meas_set[type][meas]))
            bran=d_meas_set[type][meas][i]
            k,m = split(bran,",")
            k=parse(Int64,k)
            m=parse(Int64,m)
            for (key,value) in data_pf["branch"]
                if (value["f_bus"] == k) & (value["t_bus"] == m)
                    d_meas_set[type][meas][i]= key*",from"
                    break
                elseif (value["f_bus"] == m) & (value["t_bus"] == k)
                    d_meas_set[type][meas][i]= key*",to"
                    break
                end
            end
        end
    end

    "branchdc" ∈ keys(data_pf) ? nothing : return d_meas_set #if there is no DC branch, return the meas_set as it is
    for element in branch_types_dc
        (meas,type) = split(element,"-")
        for i in range(1,length(d_meas_set[type][meas]))
            bran=d_meas_set[type][meas][i]
            k,m = split(bran,",")
            k=parse(Int64,k)
            m=parse(Int64,m)
            for (key,value) in data_pf["branchdc"]
                if (value["fbusdc"] == k) & (value["tbusdc"] == m)
                    d_meas_set[type][meas][i]= key*",from"
                    break
                elseif (value["fbusdc"] == m) & (value["tbusdc"] == k)
                    d_meas_set[type][meas][i]= key*",to"
                    break
                end
            end
        end
    end




    return d_meas_set
end


function determine_meas_id(value,data, current_rect=true)


    if value["var"] in [:pd, :qd, :pg, :qg,:p,:q]
        type = ["scada"]
        if value["var"] in [:pd, :qd, :pg, :qg]
            meas=["pinj"]
            cmp = occursin("d", string(value["var"])) ? "load" : "gen"
            s = data[cmp][string(value["cmp_id"])][cmp * "_bus"] 
        elseif value["var"] in [:p, :q]
            meas=["p"]
            branch=value["cmp_id"]
            direction=value["direction"]
            s="$(branch[1]),$direction"
        end
    elseif value["var"] in [:vm]
        type = ["scada", "pmu"]
        meas=["vm", "v"]
        s = value["cmp_id"]
    elseif value["var"] in [:va,:cr,:ci,:cm,:ca]
        currents= current_rect ? [:cr,:ci] : [:cm,:ca]
        type = ["pmu"]
        if value["var"] in [:va]
            meas=["v"]
            s=value["cmp_id"]
        elseif value["var"] in currents
            meas=["c"]
            branch=value["cmp_id"]
            direction=value["direction"]
            s="$(branch[1]),$direction"
        end
    elseif value["var"] in [:vdcm,:i_dcgrid,:p_dcgrid]
        type = ["dc","dtu"]
        if value["var"] == :vdcm
            meas=["vdcm","vdcm"]
            s=value["cmp_id"]
        elseif value["var"] == :i_dcgrid
            meas=["i_dcgrid","i_dcgrid"]
            branch=value["cmp_id"]
            direction=value["direction"]
            s="$(branch[1]),$direction"
        elseif value["var"] == :p_dcgrid
            meas=["p_dcgrid","p_dcgrid"]
            branch=value["cmp_id"]
            direction=value["direction"]
            s="$(branch[1]),$direction"
        end
    elseif value["var"] in [:pconv_dc,:pconv_ac, :qconv_ac, :pconv_pr_fr, :qconv_pr_fr, :pconv_tf_fr, :qconv_tf_fr, :pconv_tf_to, :qconv_tf_to, :vmc, :vmf, :mconv]
        type = ["conv"]
        if value["var"] in [:pconv_ac, :qconv_ac]
            meas=["pconv_ac"]
            s=value["cmp_id"]
        elseif value["var"] in [:pconv_pr_fr, :qconv_pr_fr]
            meas=["pconv_pr_fr"]
            s=value["cmp_id"]
        elseif value["var"] in [:pconv_tf_fr, :qconv_tf_fr]
            meas=["pconv_tf_fr"]
            s=value["cmp_id"]
        elseif value["var"] in [:pconv_tf_to, :qconv_tf_to]
            meas=["pconv_tf_to"]
            s=value["cmp_id"]
        elseif value["var"] in [:vmc]
            meas=["vmc"]
            s=value["cmp_id"]
        elseif value["var"] in [:vmf]
            meas=["vmf"]
            s=value["cmp_id"]
        elseif value["var"] in [:mconv]
            meas=["mconv"]
            s=value["cmp_id"]
        elseif value["var"] in [:pconv_dc]
            meas=["pconvdc"]
            s=value["cmp_id"]
        end
    else
        @warn "Measurement type $(value["var"]) not recognized. No measurement ID will be assigned."
        type = []
        meas = []
        s = ""
    end
    return type, meas, s
end

function filtermeas_set!(data_se,d_meas_set;sample_error=false,min_sig=1e-6,seed=1)

    _ACDCSE._RAN.seed!(seed)
    prec_type = Dict(
        "scada"=> 0.01,
        "pmu"=> 0.001,
        "dc"=> 0.01,
        "conv"=> 0.01,
        "dtu"=> 0.001
    )
    prec_var = Dict(
        'v'=>1,
        'c'=>1,
        'i'=>1,
        'p'=>2,
        'q'=>2,
        'm'=>1
    )

    d_prec=Dict()
    d_keys=Dict()
    for (key, value) in d_meas_set
        d_prec[key] = Dict()
        d_keys[key] = Dict()
        for (meas, value2) in value
            d_prec[key][meas] = prec_type[key] * prec_var[meas[1]]
            d_keys[key][meas] = []
        end
    end

    selec_keys=Vector()
    meas_types=Dict()
    for (key, value) in d_meas_set  
        meas_types[key]=[]
    end




    for (key,value) in data_se["meas"]

        type, meas, s = determine_meas_id(value,data_se)

        for (i,type_i) in enumerate(type)
            
            meas_i = meas[i]

            type_i in keys(d_meas_set) ? nothing : continue
            if s in d_meas_set[type_i][meas_i]

                push!(meas_types[type_i],meas_i)
                push!(d_keys[type_i][meas_i],key)
                push!(selec_keys,key)
            end
        end
    end


    for (key,value) in data_se["meas"]
        if !(key in selec_keys) 
            delete!(data_se["meas"],key)
        end
    end

   if "pmu" in keys(d_meas_set) && "scada" in keys(d_meas_set)
        copy_pmus_vm!(d_keys,data_se)
   end 
   
    if "dc" in keys(d_meas_set) && "dtu" in keys(d_meas_set)
        copy_dtu!(d_keys,data_se)
   end 
    
    for (type,types_meas) in d_keys
        for (meas, meas_type) in types_meas
            for keys in meas_type
                N=length(data_se["meas"][keys]["dst"])
                σ_true =[maximum([abs(data_se["meas"][keys]["dst"][i].μ* d_prec[type][meas]/3), min_sig]) for i in 1:N]
                μ_sampled = sample_error ? [_ACDCSE._RAN.rand(_ACDCSE._DST.Normal(data_se["meas"][keys]["dst"][i].μ, σ_true[i]), ) for i in 1:N] : [data_se["meas"][keys]["dst"][i].μ for i ∈ 1:N]
                σ_sampled= [maximum([abs(μ_sampled[i]* d_prec[type][meas]/3), min_sig]) for i in 1:N]
                dst=[_ACDCSE._DST.Normal(μ_sampled[i], σ_sampled[i]) for i in 1:N]
                data_se["meas"][keys]["dst"]=dst
            end
        end
    end

    return d_keys, d_prec

end






function prepare_hyb_meas_set!(data_se,d_meas_set;sample_error=false,min_sig=1e-6,seed=1)

    prec_type = Dict(
        "scada"=> 0.01,
        "pmu"=> 0.001,
        "dc"=> 0.01,
        "conv"=> 0.01,
        "dtu"=> 0.001
    )
    prec_var = Dict(
        'v'=>1,
        'c'=>1,
        'i'=>1,
        'p'=>2,
        'q'=>2,
        'm'=>1
    )

    d_prec=Dict()
    d_keys=Dict()
    for (key, value) in d_meas_set
        d_prec[key] = Dict()
        d_keys[key] = Dict()
        for (meas, value2) in value
            d_prec[key][meas] = prec_type[key] * prec_var[meas[1]]
            d_keys[key][meas] = []
        end
    end

    selec_keys=Vector()
    meas_types=Dict()
    for (key, value) in d_meas_set  
        meas_types[key]=[]
    end




    for (key,value) in data_se["meas"]

        type, meas, s = determine_meas_id(value,data_se)

        for (i,type_i) in enumerate(type)
            
            println("type_i: $type_i, meas_i: $(meas), s: $s")
            meas_i = meas[i]

            type_i in keys(d_meas_set) ? nothing : continue
            if s in d_meas_set[type_i][meas_i]
                println("type_i: $type_i, meas_i: $meas_i, s: $s")
                push!(meas_types[type_i],meas_i)
                push!(d_keys[type_i][meas_i],key)
                push!(selec_keys,key)
            end
        end
    end


    for (key,value) in data_se["meas"]
        if !(key in selec_keys) 
            delete!(data_se["meas"],key)
        end
    end

   if "pmu" in keys(d_meas_set) && "scada" in keys(d_meas_set)
        copy_pmus_vm!(d_keys,data_se)
   end 
   
    if "dc" in keys(d_meas_set) && "dtu" in keys(d_meas_set)
        copy_dtu!(d_keys,data_se)
   end

    return d_keys, d_prec

end




function split_meas_set!(d_keys,data_se,meas_types=["scada","pmu","dc","conv","dtu"])


    select_keys=Vector()
    for type in meas_types
        for var in keys(d_keys[type])
            append!(select_keys, d_keys[type][var])
        end
    end


    for (key,value) in data_se["meas"]
        if !(key in select_keys) 
            delete!(data_se["meas"],key)
        end
    end


end



function introduce_noise!(data_se, d_prec,d_keys; seed=1,min_sig=1e-6,bound=true,bound_factor=3)

    _ACDCSE._RAN.seed!(seed)

    for (type,types_meas) in d_keys
        for (meas, meas_type) in types_meas
            for keys in meas_type
                N=length(data_se["meas"][keys]["dst"])
                
                μ_true = [data_se["meas"][keys]["dst"][i].μ for i in 1:N]
                σ_true =[abs(μ_true[i]* d_prec[type][meas]/3) for i in 1:N]

                μ_sampled = Vector{Float64}(undef, N)
                for i in 1:N
                    if σ_true[i] > min_sig/100
                        if bound
                            μ_sampled[i] = _ACDCSE._RAN.rand(_ACDCSE._DST.Normal(μ_true[i], σ_true[i]))
                            if abs(μ_sampled[i]-μ_true[i]) > bound_factor*σ_true[i]
                                μ_sampled[i] = μ_true[i] + sign(μ_sampled[i]-μ_true[i]) * bound_factor*σ_true[i]
                            end
                        else
                            μ_sampled[i] = _ACDCSE._RAN.rand(_ACDCSE._DST.Normal(μ_true[i], σ_true[i]))
                        end
                    else
                        μ_sampled[i] = μ_true[i]
                    end
                end



                σ_sampled= [maximum([abs(μ_sampled[i]* d_prec[type][meas]/3), min_sig]) for i in 1:N]
                dst=[_ACDCSE._DST.Normal(μ_sampled[i], σ_sampled[i]) for i in 1:N]
                data_se["meas"][keys]["dst"]=dst
            end
        end
    end


end

function filtermeas_set_basic!(data_se,d_meas_set)

    selec_keys=Vector()

    


    for (key,value) in data_se["meas"]
        if (value["var"]==:pd) | (value["var"]==:qd) | (value["var"]==:pg) | (value["var"]==:qg)
            bus=value["cmp_id"]
            if bus in d_meas_set["p_ac_scada"]
                push!(selec_keys,key)
            end
        elseif (value["var"]==:vm)
            bus=value["cmp_id"]
            if bus in d_meas_set["vm_scada"]
                push!(selec_keys,key)
            end
            if bus in d_meas_set["v_pmu"]
                push!(selec_keys,key)
            end
        elseif (value["var"]==:va)
            bus=value["cmp_id"]
            if bus in d_meas_set["v_pmu"]
                push!(selec_keys,key)
            end
        elseif (value["var"]==:p) || (value["var"]==:q)
            branch=value["cmp_id"]
            direction=value["direction"]
            s="$(branch[1]),$direction"
            if s in d_meas_set["pf_ac_scada"]
                push!(selec_keys,key)
            end
        elseif (value["var"]==:vdcm)
            bus=value["cmp_id"]
            if bus in d_meas_set["vdcm"]
                push!(selec_keys,key)
            end
        elseif (value["var"]==:p_dcgrid)
            branch=value["cmp_id"]
            direction=value["direction"]
            s="$(branch[1]),$direction"
            if s in d_meas_set["p_dc"]
                push!(selec_keys,key)
            end
        end
    end


    for (key,value) in data_se["meas"]
        if !(key in selec_keys) 
            delete!(data_se["meas"],key)
        end
    end

end




function weight_z_autodiffse_old!(network_info,d_meas_set;min_sig=1e-5,prec_dict=Dict(),filter=false)

    if prec_dict == Dict()
        prec_dict=Dict(
            "p_scada" => 0.02, #ac network power measurements
            "vm_scada" => 0.01, #ac network voltage measurements
            "vm_dc"=> 0.01, #dc network voltage measurements
            "p_dc" => 0.01,  # dc network measurements
            "i_dc" => 0.01, # dc network measurements
            "v_pmu" => 0.001, # pmu ac network measurements
            "Ipol_pmu"=>0.001, # pmu ac network measurements
            "Irec_pmu"=>0.001, # pmu ac network measurements
            "p_conv"=> 0.01,  # power measurements at converters
            "pdc_conv" => 0.01, # power measurements at dc converters
            "vm_conv"=> 0.01, # power measurements at converters
            "m_conv" => 0.01, # voltage ratio converter
        )
    end

    meas_types= keys(d_meas_set)


    selec_keys=Vector()
    meas_types=Dict()

    for (key, value) in d_meas_set  
        meas_types[key]=[]
    end

    i=0
    meas_exclude=[]
    for z_i in network_info["z"]
        if z_i.var in [:pd, :qd, :pg, :qg, :pinj, :qinj, :p,:q] #SCADA measurements
            type="scada"
            calc_sigma!(z_i,prec_dict["p_scada"]; min_sig=min_sig)
            if (z_i.var == :p) || (z_i.var == :q)
                branch=z_i.cmp_id
                direction=z_i.direction
                "$(z_i.cmp_id),$direction" in d_meas_set[type]["p"] ? nothing : push!(meas_exclude,i)    
            else
                bus=z_i.cmp_id
                bus in d_meas_set[type]["pinj"] ? nothing : push!(meas_exclude,i)
            end
        elseif z_i.var in [:vm]
            type="scada"
            if z_i.cmp_id in d_meas_set[type]["vm"]
                calc_sigma!(z_i,prec_dict["vm_scada"]; min_sig=min_sig)
            elseif z_i.cmp_id in d_meas_set["pmu"]["v"]
                calc_sigma!(z_i,prec_dict["v_pmu"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:va]
            bus=z_i.cmp_id
            type="pmu"
            if ("pmu" in keys(d_meas_set)) && (bus in d_meas_set[type]["v"])
                calc_sigma!(z_i,prec_dict["v_pmu"]; min_sig=min_sig)
            elseif occursin("ref", z_i.z_label)
                calc_sigma!(z_i,prec_dict["v_pmu"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:cr,:ci,:cm,:ca]
            branch=z_i.cmp_id
            direction=z_i.direction
            s="$branch,$direction"
            type="pmu"
            if (s in d_meas_set[type]["c"]) || (s in d_meas_set[type]["c"])
                calc_sigma!(z_i,prec_dict["Irec_pmu"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:vdcm]
            bus=z_i.cmp_id
            if bus in d_meas_set["dc"]["vdcm"]
                calc_sigma!(z_i,prec_dict["vm_dc"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:p_dcgrid]
            branch=z_i.cmp_id
            direction=z_i.direction
            s="$branch,$direction"
            if s in d_meas_set["dc"]["p_dcgrid"]
                calc_sigma!(z_i,prec_dict["p_dc"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:i_dcgrid]
            branch=z_i.cmp_id
            direction=z_i.direction
            s="$branch,$direction"
            if s in d_meas_set["dc"]["i_dcgrid"]
                calc_sigma!(z_i,prec_dict["i_dc"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:pconv_ac, :qconv_ac,:pconv_pr_fr,:qconv_pr_fr,:pconv_tf_fr,:qconv_tf_fr,:pconv_tf_to,:qconv_tf_to]
            conv=z_i.cmp_id
            type="conv"
            if (z_i.var in [:pconv_ac, :qconv_ac]) && (conv in d_meas_set[type]["pconv_ac"])
                calc_sigma!(z_i,prec_dict["p_conv"]; min_sig=min_sig)
            elseif (z_i.var in [:pconv_pr_fr, :qconv_pr_fr]) && (conv in d_meas_set[type]["pconv_pr_fr"])
                calc_sigma!(z_i,prec_dict["p_conv"]; min_sig=min_sig)
            elseif (z_i.var in [:pconv_tf_fr, :qconv_tf_fr]) && (conv in d_meas_set[type]["pconv_tf_fr"])
                calc_sigma!(z_i,prec_dict["p_conv"]; min_sig=min_sig)
            elseif (z_i.var in [:pconv_tf_to, :qconv_tf_to]) && (conv in d_meas_set[type]["pconv_tf_to"])
                calc_sigma!(z_i,prec_dict["p_conv"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:vmc,:vmf]
            conv=z_i.cmp_id
            if (z_i.var == :vmc  && (conv in d_meas_set["conv"]["vmc"])) || (z_i.var == :vmf && (conv in d_meas_set["conv"]["vmf"]))
                calc_sigma!(z_i,prec_dict["vm_conv"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var == :mconv
            conv=z_i.cmp_id
            if conv in d_meas_set["conv"]["mconv"]
                calc_sigma!(z_i,prec_dict["m_conv"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)   
            end
        elseif z_i.var in [:pconv_dc,:pconv_dcg_shunt,:pconv_dcg]
            conv=z_i.cmp_id
            if ((conv in d_meas_set["conv"]["pconv_dc"]) && (z_i.var == :pconv_dc)) || ((conv in d_meas_set["conv"]["pconv_dcg_shunt"]) && (z_i.var == :pconv_dcg_shunt)) || ((conv in d_meas_set["conv"]["pconv_dcg"]) && (z_i.var == :pconv_dcg))
                calc_sigma!(z_i,prec_dict["p_dc"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        end
        i= i + 1
    end
        
    filter==true ? deleteat!(network_info["z"], meas_exclude) : nothing

    # for zv_i in network_info["z_virtual"]
    #     calc_sigma!(zv_i,min_sig=min_sig)
    # end

    for g_i in network_info["g"]
        calc_sigma!(g_i,prec=min_sig)
    end

end


function weight_z_autodiffse!(network_info,d_keys;min_sig=1e-5,prec_dict=Dict(),filter=false)

    if prec_dict == Dict()
        prec_dict=Dict(
            "p_scada" => 0.02, #ac network power measurements
            "vm_scada" => 0.01, #ac network voltage measurements
            "vm_dc"=> 0.01, #dc network voltage measurements
            "vm_dtu"=> 0.001, #dc network voltage measurements
            "p_dc" => 0.02,  # dc network measurements
            "i_dc" => 0.01, # dc network measurements
            "i_dtu" => 0.001, # dc network measurements
            "v_pmu" => 0.001, # pmu ac network measurements
            "Ipol_pmu"=>0.001, # pmu ac network measurements
            "Irec_pmu"=>0.001, # pmu ac network measurements
            "p_conv"=> 0.02,  # power measurements at converters
            "pdc_conv" => 0.01, # power measurements at dc converters
            "vm_conv"=> 0.01, # power measurements at converters
            "m_conv" => 0.01, # voltage ratio converter
        )
    end





    i=0
    meas_exclude=[]

    for z_i in network_info["z"]
        if z_i.var in [:pd, :qd, :pg, :qg, :pinj, :qinj, :p,:q] #SCADA measurements
            type="scada"
            calc_sigma!(z_i,prec_dict["p_scada"]; min_sig=min_sig)
            if (z_i.var == :p) || (z_i.var == :q)
                z_i.key in d_keys["scada"]["p"] ? nothing : push!(meas_exclude,i)
            else
                z_i.key in d_keys[type]["pinj"] ? nothing : push!(meas_exclude,i)
            end
        elseif z_i.var in [:vm]
            type="scada"
            if haskey(d_keys, "scada") && (z_i.key in d_keys[type]["vm"])
                calc_sigma!(z_i,prec_dict["vm_scada"]; min_sig=min_sig)
            elseif haskey(d_keys, "pmu") && (z_i.key in d_keys["pmu"]["v"])
                calc_sigma!(z_i,prec_dict["v_pmu"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:va]
            bus=z_i.cmp_id
            type="pmu"
            if ("pmu" in keys(d_keys)) && (z_i.key in d_keys[type]["v"])
                calc_sigma!(z_i,prec_dict["v_pmu"]; min_sig=min_sig)
            elseif occursin("ref", z_i.z_label)
                calc_sigma!(z_i,prec_dict["v_pmu"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:cr,:ci,:cm,:ca]
            branch=z_i.cmp_id
            type="pmu"
            if z_i.key in d_keys[type]["c"]
                calc_sigma!(z_i,prec_dict["Irec_pmu"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:vdcm]
            if z_i.key in d_keys["dc"]["vdcm"]
                calc_sigma!(z_i,prec_dict["vm_dc"]; min_sig=min_sig)
            elseif z_i.key in d_keys["dtu"]["vdcm"]
                calc_sigma!(z_i,prec_dict["vm_dtu"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:p_dcgrid]
            if z_i.key in d_keys["dc"]["p_dcgrid"]
                calc_sigma!(z_i,prec_dict["p_dc"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:i_dcgrid]
            if z_i.key in d_keys["dc"]["i_dcgrid"]
                calc_sigma!(z_i,prec_dict["i_dc"]; min_sig=min_sig)
            elseif z_i.key in d_keys["dtu"]["i_dcgrid"]
                calc_sigma!(z_i,prec_dict["i_dtu"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:pconv_ac, :qconv_ac,:pconv_pr_fr,:qconv_pr_fr,:pconv_tf_fr,:qconv_tf_fr,:pconv_tf_to,:qconv_tf_to]
            conv=z_i.cmp_id
            type="conv"
            if (z_i.var in [:pconv_ac, :qconv_ac]) && (z_i.key in d_keys[type]["pconv_ac"])
                calc_sigma!(z_i,prec_dict["p_conv"]; min_sig=min_sig)
            elseif (z_i.var in [:pconv_pr_fr, :qconv_pr_fr]) && (z_i.key in d_keys[type]["pconv_pr_fr"])
                calc_sigma!(z_i,prec_dict["p_conv"]; min_sig=min_sig)
            elseif (z_i.var in [:pconv_tf_fr, :qconv_tf_fr]) && (z_i.key in d_keys[type]["pconv_tf_fr"])
                calc_sigma!(z_i,prec_dict["p_conv"]; min_sig=min_sig)
            elseif (z_i.var in [:pconv_tf_to, :qconv_tf_to]) && (z_i.key in d_keys[type]["pconv_tf_to"])
                calc_sigma!(z_i,prec_dict["p_conv"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:vmc,:vmf]
            conv=z_i.cmp_id
            if (z_i.var == :vmc  && (z_i.key in d_keys["conv"]["vmc"])) || (z_i.var == :vmf && (z_i.key in d_keys["conv"]["vmf"]))
                calc_sigma!(z_i,prec_dict["vm_conv"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var == :mconv
            if z_i.key in d_keys["conv"]["mconv"]
                calc_sigma!(z_i,prec_dict["m_conv"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)   
            end
        elseif z_i.var in [:pconv_dc,:pconv_dcg_shunt,:pconv_dcg]
            if ((z_i.key in d_keys["conv"]["pconv_dc"]) && (z_i.var == :pconv_dc)) || ((z_i.key in d_keys["conv"]["pconv_dcg_shunt"]) && (z_i.var == :pconv_dcg_shunt)) || ((z_i.key in d_keys["conv"]["pconv_dcg"]) && (z_i.var == :pconv_dcg))
                calc_sigma!(z_i,prec_dict["p_dc"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        end
        i= i + 1
    end
        
    filter==true ? deleteat!(network_info["z"], meas_exclude) : nothing

    # for zv_i in network_info["z_virtual"]
    #     calc_sigma!(zv_i,min_sig=min_sig)
    # end

    for g_i in network_info["g"]
        calc_sigma!(g_i,prec=min_sig)
    end
    
    if haskey(network_info, "g_Pinj")
        for g_i in network_info["g_Pinj"]
            calc_sigma!(g_i,prec=min_sig)
        end
        if haskey(network_info, "g_Qinj")
            for g_i in network_info["g_Qinj"]
                calc_sigma!(g_i,prec=min_sig)
            end
        end
    end



end



function copy_pmus_vm!(d_keys,data_se)

    keys2copy=intersect(d_keys["pmu"]["v"], d_keys["scada"]["vm"])
   
    for key in keys2copy
       
        max_key=maximum(parse.(Int,keys(data_se["meas"])))
        new_key=max_key+1
        data_se["meas"]["$new_key"]=deepcopy(data_se["meas"][key])
            
        if key in d_keys["pmu"]["v"]
            push!(d_keys["pmu"]["v"],"$new_key")
            d_keys["pmu"]["v"]=setdiff(d_keys["pmu"]["v"], ["$key"])
        end

    end
end



function copy_dtu!(d_keys,data_se)

    keys2copy=intersect(d_keys["dtu"]["i_dcgrid"], d_keys["dc"]["i_dcgrid"])

    for key in keys2copy
        max_key=maximum(parse.(Int,keys(data_se["meas"])))
        new_key=max_key+1
        data_se["meas"]["$new_key"]=deepcopy(data_se["meas"][key])

        if key in d_keys["dtu"]["i_dcgrid"]
            push!(d_keys["dtu"]["i_dcgrid"],"$new_key")
            d_keys["dtu"]["i_dcgrid"]=setdiff(d_keys["dtu"]["i_dcgrid"], ["$key"])
        end
    end

    keys2copy=intersect(d_keys["dtu"]["vdcm"], d_keys["dc"]["vdcm"])

    for key in keys2copy
        max_key=maximum(parse.(Int,keys(data_se["meas"])))
        new_key=max_key+1
        data_se["meas"]["$new_key"]=deepcopy(data_se["meas"][key])

        if key in d_keys["dtu"]["vdcm"]
            push!(d_keys["dtu"]["vdcm"],"$new_key")
            d_keys["dtu"]["vdcm"]=setdiff(d_keys["dtu"]["vdcm"], ["$key"])
        end
    end
end



function weight_z_autodiffse2!(network_info,d_keys;min_sig=1e-5,prec_dict=Dict())

    if prec_dict == Dict()
        prec_dict=Dict(
            "p_scada" => 0.02, #ac network power measurements
            "vm_scada" => 0.01, #ac network voltage measurements
            "vm_dc"=> 0.01, #dc network voltage measurements
            "p_dc" => 0.01,  # dc network measurements
            "i_dc" => 0.01, # dc network measurements
            "v_pmu" => 0.001, # pmu ac network measurements
            "Ipol_pmu"=>0.001, # pmu ac network measurements
            "Irec_pmu"=>0.001, # pmu ac network measurements
            "p_conv"=> 0.01,  # power measurements at converters
            "pdc_conv" => 0.01, # power measurements at dc converters
            "vm_conv"=> 0.01, # power measurements at converters
            "m_conv" => 0.01, # voltage ratio converter
        )
    end




    i=0
    meas_exclude=[]
    for (i, z_i) in enumerate(network_info["z"])
        z_key = network_info["z_keys"][i]

        if z_i.var in [:pd, :qd, :pg, :qg, :pinj, :qinj, :p,:q] #SCADA measurements
            type="scada"
            calc_sigma!(z_i,prec_dict["p_scada"]; min_sig=min_sig)
            if (z_i.var == :p) || (z_i.var == :q)
                branch=z_i.cmp_id
                direction=z_i.direction
                "$(z_i.cmp_id),$direction" in d_meas_set[type]["p"] ? nothing : push!(meas_exclude,i)    
            else
                bus=z_i.cmp_id
                bus in d_meas_set[type]["pinj"] ? nothing : push!(meas_exclude,i)
            end
        elseif z_i.var in [:vm]
            type="scada"
            if ("scada" in z_key) && (z_key in d_keys["scada"]["vm"])
                calc_sigma!(z_i,prec_dict["vm_scada"]; min_sig=min_sig)
            elseif z_key in d_meas_set["pmu"]["v"]
                calc_sigma!(z_i,prec_dict["v_pmu"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:va]
            bus=z_i.cmp_id
            type="pmu"
            if ("pmu" in keys(d_meas_set)) && (bus in d_meas_set[type]["v"])
                calc_sigma!(z_i,prec_dict["v_pmu"]; min_sig=min_sig)
            elseif occursin("ref", z_i.z_label)
                calc_sigma!(z_i,prec_dict["v_pmu"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:cr,:ci,:cm,:ca]
            branch=z_i.cmp_id
            direction=z_i.direction
            s="$branch,$direction"
            type="pmu"
            if (s in d_meas_set[type]["c"]) || (s in d_meas_set[type]["c"])
                calc_sigma!(z_i,prec_dict["Irec_pmu"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:vdcm]
            bus=z_i.cmp_id
            if bus in d_meas_set["dc"]["vdcm"]
                calc_sigma!(z_i,prec_dict["vm_dc"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:p_dcgrid]
            branch=z_i.cmp_id
            direction=z_i.direction
            s="$branch,$direction"
            if s in d_meas_set["dc"]["p_dcgrid"]
                calc_sigma!(z_i,prec_dict["p_dc"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:i_dcgrid]
            branch=z_i.cmp_id
            direction=z_i.direction
            s="$branch,$direction"
            if s in d_meas_set["dc"]["i_dcgrid"]
                calc_sigma!(z_i,prec_dict["i_dc"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:pconv_ac, :qconv_ac,:pconv_pr_fr,:qconv_pr_fr,:pconv_tf_fr,:qconv_tf_fr,:pconv_tf_to,:qconv_tf_to]
            conv=z_i.cmp_id
            type="conv"
            if (z_i.var in [:pconv_ac, :qconv_ac]) && (conv in d_meas_set[type]["pconv_ac"])
                calc_sigma!(z_i,prec_dict["p_conv"]; min_sig=min_sig)
            elseif (z_i.var in [:pconv_pr_fr, :qconv_pr_fr]) && (conv in d_meas_set[type]["pconv_pr_fr"])
                calc_sigma!(z_i,prec_dict["p_conv"]; min_sig=min_sig)
            elseif (z_i.var in [:pconv_tf_fr, :qconv_tf_fr]) && (conv in d_meas_set[type]["pconv_tf_fr"])
                calc_sigma!(z_i,prec_dict["p_conv"]; min_sig=min_sig)
            elseif (z_i.var in [:pconv_tf_to, :qconv_tf_to]) && (conv in d_meas_set[type]["pconv_tf_to"])
                calc_sigma!(z_i,prec_dict["p_conv"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var in [:vmc,:vmf]
            conv=z_i.cmp_id
            if (z_i.var == :vmc  && (conv in d_meas_set["conv"]["vmc"])) || (z_i.var == :vmf && (conv in d_meas_set["conv"]["vmf"]))
                calc_sigma!(z_i,prec_dict["vm_conv"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        elseif z_i.var == :mconv
            conv=z_i.cmp_id
            if conv in d_meas_set["conv"]["mconv"]
                calc_sigma!(z_i,prec_dict["m_conv"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)   
            end
        elseif z_i.var in [:pconv_dc,:pconv_dcg_shunt,:pconv_dcg]
            conv=z_i.cmp_id
            if ((conv in d_meas_set["conv"]["pconv_dc"]) && (z_i.var == :pconv_dc)) || ((conv in d_meas_set["conv"]["pconv_dcg_shunt"]) && (z_i.var == :pconv_dcg_shunt)) || ((conv in d_meas_set["conv"]["pconv_dcg"]) && (z_i.var == :pconv_dcg))
                calc_sigma!(z_i,prec_dict["p_dc"]; min_sig=min_sig)
            else
                push!(meas_exclude,i)
            end
        end
        i= i + 1
    end
        
    filter==true ? deleteat!(network_info["z"], meas_exclude) : nothing

    # for zv_i in network_info["z_virtual"]
    #     calc_sigma!(zv_i,min_sig=min_sig)
    # end

    for g_i in network_info["g"]
        calc_sigma!(g_i,prec=min_sig)
    end




end
