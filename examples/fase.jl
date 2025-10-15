using SparseArrays 

function h_vec(x, h::Vector{Function})
    return [f(x) for f in h]
end


function compute_gain(h::Vector{Function}, x::Vector{Float64}, W::Matrix{Float64})


    J = zeros(length(h), length(x))

    compute_Jacobian!(h, x, J, length(h))

 
    Jsparse = sparse(J)
    Wsparse = sparse(W)
    Gain_sparse = Jsparse' * Wsparse * Jsparse
    Gain= Matrix(Gain_sparse)

    return Gain,J
end


function update_gain!(network_info)

    h = network_info["h"]
    x = network_info["x"]
    W = network_info["W"]
    J = network_info["H"]

    compute_Jacobian!(h, x, H, length(h))
    Jsparse = sparse(H)
    Wsparse = sparse(W)
    Gain_sparse = Jsparse' * Wsparse * Jsparse

    Gain= Matrix(Gain_sparse)
    return Gain,H
end



function compute_Jacobian!(h::Vector{Function}, x::Vector{Float64},J::Matrix{Float64}, m::Int)

    for i in 1:m
        J[i, :] = Enzyme.gradient(Enzyme.Forward,
                                # pass the function itself as a *constant* pointer
                                Const(h[i]),
                                x)[1]
    end

end

function generate_data_prior_diag!(result,data,d::Float64=1e-2)


   data["prior"]=Dict{String, Any}()

    m_idx=isempty(data["prior"]) ? 1 : maximum(parse.(Int, keys(data["prior"])))
    for (b,bus) in data["bus"] 
        c=1
        data["prior"]["$m_idx"]=Dict(
            "a_ij"=>d,
            "cmp_i"=> :bus,
            "cmp_j"=> :bus,
            "var_i"=> :vm,
            "var_j"=> :vm,
            "conn_i" => c,
            "conn_j" => c,
            "cmp_id_i" => parse(Int,b),
            "cmp_id_j"=> parse(Int,b),
            "π_i"=>result["solution"]["bus"][b]["vm"][c],
            "π_j"=>result["solution"]["bus"][b]["vm"][c]
        )
        m_idx+=1
    end


    for (b,bus) in data["bus"] 
        c=1
        data["prior"]["$m_idx"]=Dict(
            "a_ij"=>d,
            "cmp_i"=> :bus,
            "cmp_j"=> :bus,
            "var_i"=> :va,
            "var_j"=> :va,
            "conn_i" => c,
            "conn_j" => c,
            "cmp_id_i" => parse(Int,b),
            "cmp_id_j"=> parse(Int,b),
            "π_i"=>result["solution"]["bus"][b]["va"][c],
            "π_j"=>result["solution"]["bus"][b]["va"][c]
        )
        m_idx+=1
    end


    for (bdc,busdc) in data["busdc"]
        for c in 1:length(result["solution"]["busdc"][bdc]["vm"]) 
            data["prior"]["$m_idx"]=Dict(
                "a_ij"=>d,
                "cmp_i"=> :busdc,
                "cmp_j"=> :busdc,
                "var_i"=> :vdcm,
                "var_j"=> :vdcm,
                "conn_i" => c,
                "conn_j" => c,
                "cmp_id_i" => parse(Int,bdc),
                "cmp_id_j"=> parse(Int,bdc),
                "π_i"=>result["solution"]["busdc"][bdc]["vm"][c],
                "π_j"=>result["solution"]["busdc"][bdc]["vm"][c]
            )
            m_idx+=1
        end
    end


    for (cv,conv) in data["convdc"]
        for c in 1:length(result["solution"]["convdc"][cv]["vmconv"]) 
            data["prior"]["$m_idx"]=Dict(
                "a_ij"=>d,
                "cmp_i"=> :convdc,
                "cmp_j"=> :convdc,
                "var_i"=> :vmc,
                "var_j"=> :vmc,
                "conn_i" => c,
                "conn_j" => c,
                "cmp_id_i" => parse(Int,cv),
                "cmp_id_j"=> parse(Int,cv),
                "π_i"=>result["solution"]["convdc"][cv]["vmconv"][c],
                "π_j"=>result["solution"]["convdc"][cv]["vmconv"][c]
            )
            m_idx+=1
        end
    end


    for (cv,conv) in data["convdc"]
        for c in 1:length(result["solution"]["convdc"][cv]["vaconv"]) 
            data["prior"]["$m_idx"]=Dict(
                "a_ij"=>d,
                "cmp_i"=> :convdc,
                "cmp_j"=> :convdc,
                "var_i"=> :vac,
                "var_j"=> :vac,
                "conn_i" => c,
                "conn_j" => c,
                "cmp_id_i" => parse(Int,cv),
                "cmp_id_j"=> parse(Int,cv),
                "π_i"=>result["solution"]["convdc"][cv]["vaconv"][c],
                "π_j"=>result["solution"]["convdc"][cv]["vaconv"][c]
            )
            m_idx+=1
        end
    end

    
    for (cv,conv) in data["convdc"]
        for c in 1:length(result["solution"]["convdc"][cv]["vmfilt"]) 
            data["prior"]["$m_idx"]=Dict(
                "a_ij"=>d,
                "cmp_i"=> :convdc,
                "cmp_j"=> :convdc,
                "var_i"=> :vmf,
                "var_j"=> :vmf,
                "conn_i" => c,
                "conn_j" => c,
                "cmp_id_i" => parse(Int,cv),
                "cmp_id_j"=> parse(Int,cv),
                "π_i"=>result["solution"]["convdc"][cv]["vmfilt"][c],
                "π_j"=>result["solution"]["convdc"][cv]["vmfilt"][c]
            )
            m_idx+=1
        end
    end


    for (cv,conv) in data["convdc"]
        for c in 1:length(result["solution"]["convdc"][cv]["vafilt"]) 
            data["prior"]["$m_idx"]=Dict(
                "a_ij"=>d,
                "cmp_i"=> :convdc,
                "cmp_j"=> :convdc,
                "var_i"=> :vaf,
                "var_j"=> :vaf,
                "conn_i" => c,
                "conn_j" => c,
                "cmp_id_i" => parse(Int,cv),
                "cmp_id_j"=> parse(Int,cv),
                "π_i"=>result["solution"]["convdc"][cv]["vafilt"][c],
                "π_j"=>result["solution"]["convdc"][cv]["vafilt"][c]
            )
            m_idx+=1
        end
    end
end 





function generate_data_prior!(result,data_prior,data_posterior,d_keys_prior;d::Float64=1e-2,min_sig::Float64=1e-4,print_info=false,cutoff::Float64=1e0)


    network_info = build_network_info(data_prior, result)
    weight_z_autodiffse!(network_info, d_keys_prior; min_sig=min_sig )

    x_info=network_info["x_info"] 
    x=network_info["x"]

    network_info["f_h"]=create_h(network_info)
    # network_info["f_hv"]=create_h_virtual(network_info)
    
    network_info["f_g"]=create_g(network_info)

    
    # network_info["f_ht"] = vcat(network_info["f_h"], network_info["f_hv"], network_info["f_g"])
    network_info["f_ht"] = vcat(network_info["f_h"], network_info["f_g"])
    network_info["W"]=create_W_autodiff(network_info)

    

    

    for (key, bus) in result["solution"]["bus"]
    data_posterior["bus"][key]["vm_start"] = bus["vm"]
    data_posterior["bus"][key]["va_start"] = bus["va"]
    end

    for (key, conv) in result["solution"]["convdc"]
        data_posterior["convdc"][key]["vmf_start"] = conv["vmfilt"]
        data_posterior["convdc"][key]["vaf_start"] = conv["vafilt"]
        data_posterior["convdc"][key]["vmc_start"] = conv["vmconv"]
        data_posterior["convdc"][key]["vac_start"] = conv["vaconv"]
    end

    for (key,budc) in result["solution"]["busdc"]
        data_posterior["busdc"][key]["vm_start"] = [budc["vm"][1],budc["vm"][2],0.0]
    end

    if "H" in keys(network_info)
        update_gain!(network_info)
    else
        Gain,H=compute_gain(network_info["f_ht"], x, network_info["W"])
        network_info["H"] = H
    end

    print_info && print_network_info(network_info)
    
    data_posterior["prior"]=Dict{String, Any}()
    m_idx=1
    length_x = length(x_info)
    
    for (i, x_i) in enumerate(network_info["x_info"])
        for (j, x_j) in enumerate(network_info["x_info"])
            if abs(Gain[i, j]) > cutoff
                data_posterior["prior"]["$m_idx"]=Dict(
                    "a_ij"=>1/(d*Gain[i, j]), 
                    "cmp_i"=>x_i.cmp,
                    "cmp_j"=>x_j.cmp,
                    "var_i"=>x_i.var,
                    "var_j"=>x_j.var,
                    "i_i"=> i,
                    "i_j"=> j, 
                    "conn_i" => x_i.conn,
                    "conn_j" => x_j.conn,
                    "cmp_id_i" => x_i.cmp_id,
                    "cmp_id_j"=> x_j.cmp_id,
                    "π_i"=>network_info["x"][i],
                    "π_j"=>network_info["x"][j]
                )  
                m_idx += 1          
            end
        end
    end

    return network_info
    
end 




function generate_data_posterior_cov!(result,data_posterior,d_keys_posterior;min_sig::Float64=1e-4,print_info=false,cutoff::Float64=1e0)


    network_info = build_network_info(data_posterior, result)
    weight_z_autodiffse!(network_info, d_keys_posterior; min_sig=min_sig )

    x_info=network_info["x_info"] 
    x=network_info["x"]

    network_info["f_h"]=create_h(network_info)
    
    
    network_info["f_g"]=create_g(network_info)

    
    # network_info["f_ht"] = vcat(network_info["f_h"], network_info["f_hv"], network_info["f_g"])
    network_info["f_ht"] = vcat(network_info["f_h"], network_info["f_g"])
    network_info["W"]=create_W_autodiff(network_info)

  

    print_info && print_network_info(network_info)


    if "H" in keys(network_info)
        update_gain!(network_info)
    else
        Gain,H=compute_gain(network_info["f_ht"], x, network_info["W"])
        network_info["H"] = H
    end


    data_posterior["posterior"]=Dict{String, Any}()
    m_idx=1
    length_x = length(x_info)
    
    for (i, x_i) in enumerate(network_info["x_info"])
        for (j, x_j) in enumerate(network_info["x_info"])
            if abs(Gain[i, j]) > cutoff
                data_posterior["posterior"]["$m_idx"]=Dict(
                    "a_ij"=>1/(Gain[i, j]), 
                    "cmp_i"=>x_i.cmp,
                    "cmp_j"=>x_j.cmp,
                    "var_i"=>x_i.var,
                    "var_j"=>x_j.var,
                    "i_i"=> i,
                    "i_j"=> j, 
                    "conn_i" => x_i.conn,
                    "conn_j" => x_j.conn,
                    "cmp_id_i" => x_i.cmp_id,
                    "cmp_id_j"=> x_j.cmp_id,
                    "π_i"=>network_info["x"][i],
                    "π_j"=>network_info["x"][j]
                )  
                m_idx += 1          
            end
        end
    end

    return network_info
    
end 


function generate_res_cov_QR!(result,data,d_keys;min_sig::Float64=1e-4,print_info=false,error_res::Float64=1e-12)


    network_info = build_network_info(data, result)
    weight_z_autodiffse!(network_info, d_keys; min_sig=min_sig )

    x_info=network_info["x_info"] 
    x=network_info["x"]

    network_info["f_h"]=create_h(network_info)
    
    
    network_info["f_g"]=create_g(network_info)

    
    # network_info["f_ht"] = vcat(network_info["f_h"], network_info["f_hv"], network_info["f_g"])
    network_info["f_ht"] = vcat(network_info["f_h"], network_info["f_g"])
    network_info["W"]=create_W_autodiff(network_info)

  

    print_info && print_network_info(network_info)


    if "H" in keys(network_info)
        Gain,H=update_gain!(network_info)
    else
        Gain,H=compute_gain(network_info["f_ht"], x, network_info["W"])
        network_info["H"] = H
    end

    W=network_info["W"]
    Whalf=diagm(sqrt.(diag(W)))
    Q,R = qr(Whalf*H )
    Q= Matrix(Q)
    R= Matrix(R)
    Ω= inv(network_info["W"]) - H * inv(R)*inv(R') * H'   


    
    z = [z.z for z in network_info["z"]]
    g = [0 for g in network_info["g"]]
    z_total = vcat(z, g)
    hx = h_vec(x, network_info["f_ht"])


    r = abs.(z_total - hx)
    rn=Vector{Float64}(undef,length(r))

    for i in 1:length(r)
        Ω[i,i] < error_res ? rn[i]=0.0 : rn[i]=r[i]/sqrt(Ω[i,i])
    end
    return Ω,r,W


end




function generate_res_cov!(result,data,d_keys;min_sig::Float64=1e-4,print_info=false,error_res::Float64=1e-12)


    network_info = build_network_info(data, result)
    weight_z_autodiffse!(network_info, d_keys; min_sig=min_sig )

    x_info=network_info["x_info"] 
    x=network_info["x"]

    network_info["f_h"]=create_h(network_info)
    
    
    network_info["f_g"]=create_g(network_info)

    
    # network_info["f_ht"] = vcat(network_info["f_h"], network_info["f_hv"], network_info["f_g"])
    network_info["f_ht"] = vcat(network_info["f_h"], network_info["f_g"])
    network_info["W"]=create_W_autodiff(network_info)

  

    print_info && print_network_info(network_info)


    if "H" in keys(network_info)
        Gain,H=update_gain!(network_info)
        network_info["H"] = H
    else
        Gain,H=compute_gain(network_info["f_ht"], x, network_info["W"])
        network_info["H"] = H
    end

    W=network_info["W"]
    H=network_info["H"]


    Ginv=inv(Gain)
    
    Ω= inv(W) - H * Ginv * H'
    
    z = [z.z for z in network_info["z"]]
    g = [0 for g in network_info["g"]]
    z_total = vcat(z, g)
    hx = h_vec(x, network_info["f_ht"])


    r = abs.(z_total - hx)
    rn=Vector{Float64}(undef,length(r))


    return Ω,r,W,Gain,network_info

end



function generate_res_cov_new!(result,data,d_keys;min_sig::Float64=1e-4,print_info=false,error_res::Float64=1e-12)


    network_info = build_network_info_bad_data(data, result)
    weight_z_autodiffse!(network_info, d_keys; min_sig=min_sig )

    x_info=network_info["x_info"] 
    x=network_info["x"]

    network_info["f_h"]=create_h(network_info)
        
    network_info["f_g"]=create_g(network_info)

    network_info["f_g_Pinj"]=create_g_Pinj(network_info)
    network_info["f_c_Qinj"]=create_g_Qinj(network_info)
    network_info["f_c_Qinj"][3](x)

    network_info["f_ht"] = vcat(network_info["f_h"], network_info["f_g"], network_info["f_g_Pinj"], network_info["f_c_Qinj"])

    network_info["W"]=create_W_autodiff_eg(network_info)

    Gain,H=compute_gain(network_info["f_ht"], x, network_info["W"])

    # Ginv=inv(Gain)
        
    # Ω= inv(network_info["W"]) - H * Ginv * H'

    W=network_info["W"]
    Whalf=diagm(sqrt.(diag(W)))
    Q,R = qr(Whalf*H )
    Q= Matrix(Q)
    R= Matrix(R)
    Ω= inv(network_info["W"]) - H * inv(R)*inv(R') * H'   

    z = [z.z for z in network_info["z"]]
    g = [0 for g in network_info["g"]]
    g_Pinj = [0 for g in network_info["g_Pinj"]]
    g_Qinj = [0 for g in network_info["g_Qinj"]]
    z_total = vcat(z, g, g_Pinj, g_Qinj)
    hx = h_vec(x, network_info["f_ht"])

    hx = h_vec(x, network_info["f_ht"])

    r = abs.(z_total - hx)

    return Ω,r,network_info["W"],Gain,network_info

end




function generate_data_prior_test!(result,data_prior,data_posterior,d_meas_set_prior;d::Float64=1e-2,min_sig::Float64=1e-4,print_info=false)



    network_info = build_network_info(data_prior, result)
    weight_z_autodiffse!(network_info, d_meas_set_prior; min_sig=min_sig )

    x_info=network_info["x_info"] 
    x=network_info["x"]

    network_info["f_h"]=create_h(network_info)
    network_info["f_v"]=create_h_virtual(network_info)

    network_info["f_g"]=create_g(network_info)
    network_info["f_ht"] = vcat(network_info["f_h"], network_info["f_g"])
    network_info["W"]=create_W_autodiff(network_info)

  

    print_info && print_network_info(network_info)

    for (key, bus) in result["solution"]["bus"]
    data_posterior["bus"][key]["vm_start"] = bus["vm"]
    data_posterior["bus"][key]["va_start"] = bus["va"]
    end

    for (key, conv) in result["solution"]["convdc"]
        data_posterior["convdc"][key]["vmf_start"] = conv["vmfilt"]
        data_posterior["convdc"][key]["vaf_start"] = conv["vafilt"]
        data_posterior["convdc"][key]["vmc_start"] = conv["vmconv"]
        data_posterior["convdc"][key]["vac_start"] = conv["vaconv"]
    end

    for (key,budc) in result["solution"]["busdc"]
        data_posterior["busdc"][key]["vm_start"] = [budc["vm"][1],budc["vm"][2],0.0]
    end

   

    return network_info
    
end 

function generate_data_netinfo!(result,data_prior,data_posterior,d_meas_set_prior;d::Float64=1e-2,min_sig::Float64=1e-4,print_info=false)



    network_info = build_network_info(data_prior, result)
    weight_z_autodiffse!(network_info, d_meas_set_prior; min_sig=min_sig )

    x_info=network_info["x_info"] 
    x=network_info["x"]

    network_info["f_h"]=create_h(network_info)
    network_info["f_g"]=create_g(network_info)
    network_info["f_ht"] = vcat(network_info["f_h"], network_info["f_g"])
    network_info["W"]=create_W_autodiff(network_info)

  
    return network_info
    
end 








function generate_data_prior_Gdiag!(result,data_prior,data_posterior,d::Float64=1e-2)



    network_info = build_network_info(data_prior, result)
    weight_z_autodiffse!(network_info, d_meas_set; min_sig=1e-3 )
    h=create_h(network_info)
    g=create_g(network_info)
    W = create_W_autodiff(network_info, h, g)

    x_info=network_info["x_info"] 
    x=network_info["x"]

    open("x_info.txt", "w") do io
        for (i, x_i) in enumerate(x_info)
            println(io, "[$i] cmp=$(x_i.cmp), var=$(x_i.var), conn=$(x_i.conn), cmp_id=$(x_i.cmp_id), label=$(x_i.label)")
        end
    end
    open("z_sigma.txt", "w") do io
        for (i, z_i) in enumerate(network_info["z"])
            println(io, "z[$i] = $(z_i.z_label), sigma = $(z_i.σ)")
        end
    end
    println("Calculating Gain matrix...")

    Gain=compute_gain(vcat(h, g), x, W)


    data_posterior["prior"]=Dict{String, Any}()
    m_idx=1
    for (i, x_i) in enumerate(x_info)
        if abs(Gain[i, i]) > 1e-6
            if Gain[i, i] > 1e8
                Gain[i, i]=1e8
            end
            data_posterior["prior"]["$m_idx"]=Dict(
                "a_ij"=>1/(d*Gain[i, i]), 
                "cmp_i"=>x_i.cmp,
                "cmp_j"=>x_i.cmp,
                "var_i"=>x_i.var,
                "var_j"=>x_i.var,
                "conn_i" => x_i.conn,
                "conn_j" => x_i.conn,
                "cmp_id_i" => x_i.cmp_id,
                "cmp_id_j"=> x_i.cmp_id,
                "π_i"=>x[i],
                "π_j"=>x[i]
                )  
                m_idx += 1          
        end
    end

    
end 