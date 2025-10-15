"""
From power flow results `pf_res`, creates synthetic measurements
and adds them to the network `data`, which is then ready for state estimation
`σ_dict` is a dictionary with the standard deviation of the different measurement types
"""
function powerflow2measurements!(data::Dict, pf_res::Dict, σ_dict::Dict; sample_error::Bool = true, measurements::Vector{String} = ["vm", "vdcm"])
    data["meas"] = Dict{String, Any}() 

    #######################
    ### AC MEASUREMENTS ###
    #######################

    # NODAL measurements (voltages and injections) # creates all measurements with noise
    if "vm" ∈ measurements add_vm!(data, pf_res, σ_dict, sample_error) end     # ac bus voltage magnitude
    if "va" ∈ measurements add_va!(data, pf_res, σ_dict, sample_error) end     # ac bus voltage angle
    if "vr" ∈ measurements add_vr!(data, pf_res, σ_dict, sample_error) end     # ac bus voltage rectangular real
    if "vi" ∈ measurements add_vi!(data, pf_res, σ_dict, sample_error) end     # ac bus voltage rectangular imag
    if "pg" ∈ measurements add_pg!(data, pf_res, σ_dict, sample_error) end     # active power injection from generator
    if "qg" ∈ measurements add_qg!(data, pf_res, σ_dict, sample_error) end     # reactive power injection from generator
    if "pd" ∈ measurements add_pd!(data, pf_res, σ_dict, sample_error) end     # active power injection from load
    if "qd" ∈ measurements add_qd!(data, pf_res, σ_dict, sample_error) end     # reactive power injection from load
    #TODO ADD CURRENT INJECTIONS??? (CGM CGA CDM CDA)
    # NOTE: I can't imagine there are additional injection measurements! (unless we add storage or similar components?)



    # FLOW measurements
    if "p_fr"  ∈ measurements  add_p_fr!(data, pf_res, σ_dict, sample_error) end # ac active power from
    if "q_fr"  ∈ measurements  add_q_fr!(data, pf_res, σ_dict, sample_error) end # ac reactive power from
    if "p_to"  ∈ measurements  add_p_to!(data, pf_res, σ_dict, sample_error) end # ac active power to
    if "q_to"  ∈ measurements  add_q_to!(data, pf_res, σ_dict, sample_error) end # ac reactive power to 
     
    if "cm_fr" ∈ measurements add_cm_fr!(data, pf_res, σ_dict, sample_error) end # ac current magnitude from
    if "cm_to" ∈ measurements add_cm_to!(data, pf_res, σ_dict, sample_error) end # ac current magnitude to


    # ↓ the below to add? (not sure if all need to be added)
    if "ca_fr" ∈ measurements add_ca_fr!(data, pf_res, σ_dict, sample_error) end # ac current angle from
    if "ca_to" ∈ measurements add_ca_to!(data, pf_res, σ_dict, sample_error) end # ac current angle to
    if "cr_fr" ∈ measurements add_cr_fr!(data, pf_res, σ_dict, sample_error) end # ac current rectangular real from
    if "cr_to" ∈ measurements add_cr_to!(data, pf_res, σ_dict, sample_error) end # ac current rectangular real to
    if "ci_fr" ∈ measurements add_ci_fr!(data, pf_res, σ_dict, sample_error) end # ac current rectangular imag from
    if "ci_to" ∈ measurements add_ci_to!(data, pf_res, σ_dict, sample_error) end # ac current rectangular imag to

    #######################
    ### DC MEASUREMENTS ###
    #######################

    # NODAL measurements (voltages and injections)
    if "vdcm" ∈ measurements add_vdcm!(data, pf_res, σ_dict, sample_error) end # dc bus voltage
    #TODO ADD injections (FROM CONVERTERS)

    # FLOW measurements
    if "i_dcgrid_fr" ∈ measurements
        display("adding i_dcgrid_fr")
        add_i_dcgrid_fr!(data, pf_res, σ_dict, sample_error) end # dc branch currents from
    if "i_dcgrid_to" ∈ measurements add_i_dcgrid_to!(data, pf_res, σ_dict, sample_error) end # dc branch currents to

    if "p_dc_fr" ∈ measurements add_p_dc_fr!(data, pf_res, σ_dict, sample_error) end     # dc branch power from
    if "p_dc_to" ∈ measurements add_p_dc_to!(data, pf_res, σ_dict, sample_error) end     # dc branch power to


    # Converter measurements AC

    if "pconv" ∈ measurements add_pconv!(data, pf_res, σ_dict, sample_error) end
    if "qconv" ∈ measurements add_qconv!(data, pf_res, σ_dict, sample_error) end
    if "ppr_fr" ∈ measurements add_ppr_fr!(data, pf_res, σ_dict, sample_error) end
    if "qpr_fr" ∈ measurements add_qpr_fr!(data, pf_res, σ_dict, sample_error) end
    if "pgrid" ∈ measurements add_pgrid!(data, pf_res, σ_dict, sample_error) end
    if "qgrid" ∈ measurements add_qgrid!(data, pf_res, σ_dict, sample_error) end
    if "ptf_to" ∈ measurements add_ptf_to!(data, pf_res, σ_dict, sample_error) end
    if "qtf_to" ∈ measurements add_qtf_to!(data, pf_res, σ_dict, sample_error) end
    if "vmfilt" ∈ measurements add_vmfilt!(data, pf_res, σ_dict, sample_error) end
    if "vmconv" ∈ measurements add_vmconv!(data, pf_res, σ_dict, sample_error) end

    #conveter measurements DC
    
    if "pdc" ∈ measurements add_pdc!(data, pf_res, σ_dict, sample_error) end
    if "pdcg" ∈ measurements add_pdcg!(data, pf_res, σ_dict, sample_error) end
    if "pdcg_shunt" ∈ measurements add_pdcg_shunt!(data, pf_res, σ_dict, sample_error) end

    #converter eletronic MEASUREMENTS

    if "mconv" ∈ measurements add_mconv!(data, pf_res, σ_dict, sample_error) end
    

end



function add_vm!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1 #if it is the first measurement
    for (b, bus) in pf_res["solution"]["bus"]
        σ = maximum([abs(bus["vm"]*σ_dict["vm"]/3), σ_dict["sigma_min"]/3])
        μ = sample_error ? _RAN.rand(_DST.Normal(bus["vm"], σ), ) : bus["vm"]
        dst = [_DST.Normal(μ, σ)]
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :bus,
            "cmp_id"  => parse(Int, b),
            "var"     => :vm,
            "dst"     => dst
        )
        m_idx+=1
    end
end

function add_va!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (b, bus) in pf_res["solution"]["bus"]
        σ = maximum([abs(bus["va"]*σ_dict["va"]/3), σ_dict["sigma_min"]/3])
        μ = sample_error ? _RAN.rand(_DST.Normal(bus["va"], σ), ) : bus["va"]
        dst = [_DST.Normal(μ, σ)]
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :bus,
            "cmp_id"  => parse(Int, b),
            "var"     => :va,
            "dst"     => dst
        )
        m_idx+=1
    end
end

function add_p_fr!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (br, branch) in pf_res["solution"]["branch"]
        σ = maximum([abs(σ_dict["p_ac"]*branch["pf"]/3), σ_dict["sigma_min"]/3])
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["pf"], σ), ) : branch["pf"]
        dst = [_DST.Normal(μ, σ)]
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :branch,
            "cmp_id"  => (data["branch"][br]["index"], data["branch"][br]["f_bus"], data["branch"][br]["t_bus"]),
            "var"     => :p,
            "dst"     => dst,
            "direction" => :from
        )
        m_idx+=1
    end
end

function add_q_fr!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (br, branch) in pf_res["solution"]["branch"]
        σ = maximum([abs(σ_dict["q_ac"]*branch["qf"]/3), σ_dict["sigma_min"]/3])
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["qf"], σ), ) : branch["qf"]
        dst = [_DST.Normal(μ, σ)]
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :branch,
            "cmp_id"  => (data["branch"][br]["index"], data["branch"][br]["f_bus"], data["branch"][br]["t_bus"]),
            "var"     => :q,
            "dst"     => dst,
            "direction" => :from
        )
        m_idx+=1
    end
end

function add_p_to!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (br, branch) in pf_res["solution"]["branch"]
        σ = maximum([abs(σ_dict["p_ac"]*branch["pt"]/3), σ_dict["sigma_min"]/3])
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["pt"], σ), ) : branch["pt"]
        dst = [_DST.Normal(μ, σ)]
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :branch,
            "cmp_id"  => (data["branch"][br]["index"], data["branch"][br]["t_bus"], data["branch"][br]["f_bus"]),
            "var"     => :p,
            "dst"     => dst,
            "direction" => :to
        )
        m_idx+=1
    end
end

function add_q_to!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (br, branch) in pf_res["solution"]["branch"]
        σ = maximum([abs(σ_dict["q_ac"]*branch["qt"]/3), σ_dict["sigma_min"]/3])
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["qt"], σ), ) : branch["qt"]
        dst = [_DST.Normal(μ, σ)]
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :branch,
            "cmp_id"  => (data["branch"][br]["index"], data["branch"][br]["t_bus"], data["branch"][br]["f_bus"]),
            "var"     => :q,
            "dst"     => dst,
            "direction" => :to
        )
        m_idx+=1
    end
end

function add_pg!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (g, gen) in pf_res["solution"]["gen"]
        σ = maximum([abs(σ_dict["pg"]*gen["pg"]/3), σ_dict["sigma_min"]/3])
        μ = sample_error ? _RAN.rand(_DST.Normal(gen["pg"], σ), ) : gen["pg"]
        dst = [_DST.Normal(μ, σ)]
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :gen,
            "cmp_id"  => data["gen"][g]["index"],
            "var"     => :pg,
            "dst"     => dst
        )
        m_idx+=1
    end
end

function add_qg!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (g, gen) in pf_res["solution"]["gen"]
        σ = maximum([abs(σ_dict["qg"]*gen["qg"]/3), σ_dict["sigma_min"]/3])
        μ = sample_error ? _RAN.rand(_DST.Normal(gen["qg"], σ), ) : gen["qg"]
        dst = [_DST.Normal(μ, σ)]
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :gen,
            "cmp_id"  => data["gen"][g]["index"],
            "var"     => :qg,
            "dst"     => dst
        )
        m_idx+=1
    end
end

function add_pd!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (l, load) in data["load"] # contrary to generators, loads are not (o)pf variables, so can't find them in the solution dict
        σ = maximum([abs(σ_dict["pd"]*load["pd"]/3), σ_dict["sigma_min"]/3])
        μ = sample_error ? _RAN.rand(_DST.Normal(load["pd"], σ), ) : load["pd"]
        dst = [_DST.Normal(μ, σ)]
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :load,
            "cmp_id"  => data["load"][l]["index"],
            "var"     => :pd,
            "dst"     => dst
        )
        m_idx+=1
    end
end

function add_qd!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (l, load) in data["load"] # contrary to generators, loads are not (o)pf variables, so can't find them in the solution dict
        σ = maximum([abs(σ_dict["qd"]*load["qd"]/3), σ_dict["sigma_min"]/3])
        μ = sample_error ? _RAN.rand(_DST.Normal(load["qd"], σ), ) : load["qd"]
        dst = [_DST.Normal(μ, σ)]
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :load,
            "cmp_id"  => data["load"][l]["index"],
            "var"     => :qd,
            "dst"     => dst
        )
        m_idx+=1
    end
end

function add_cm_fr!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (br, branch) in pf_res["solution"]["branch"]
        σ = maximum([abs(branch["cmf"]*σ_dict["cm"]/3), σ_dict["sigma_min"]/3])
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["cmf"], σ), ) : branch["cmf"]
        dst = [_DST.Normal(μ, σ)]
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :branch,
            "cmp_id"  => (data["branch"][br]["index"], data["branch"][br]["f_bus"], data["branch"][br]["t_bus"]),
            "var"     => :cm,
            "dst"     => dst,
            "direction" => :from
        )
        m_idx+=1
    end
end

function add_cm_to!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (br, branch) in pf_res["solution"]["branch"]
        σ = maximum([abs(branch["cmt"]*σ_dict["cm"]/3), σ_dict["sigma_min"]/3])
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["cmt"], σ), ) : branch["cmt"]
        dst = [_DST.Normal(μ, σ)]
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :branch,
            "cmp_id"  => (data["branch"][br]["index"], data["branch"][br]["t_bus"], data["branch"][br]["f_bus"]),
            "var"     => :cm,
            "dst"     => dst,
            "direction" => :to
        )
        m_idx+=1
    end
end


function add_ca_fr!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (br, branch) in pf_res["solution"]["branch"]
        σ = maximum([abs(branch["caf"]*σ_dict["ca"]/3), σ_dict["sigma_min"]/3])
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["caf"], σ), ) : branch["caf"]
        dst = [_DST.Normal(μ, σ)]
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :branch,
            "cmp_id"  => (data["branch"][br]["index"], data["branch"][br]["f_bus"], data["branch"][br]["t_bus"]),
            "var"     => :ca,
            "dst"     => dst,
            "direction" => :from
        )
        m_idx+=1
    end
end

function add_ca_to!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (br, branch) in pf_res["solution"]["branch"]
        σ = maximum([abs(branch["cat"]*σ_dict["ca"]/3), σ_dict["sigma_min"]/3])
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["cat"], σ), ) : branch["cat"]
        dst = [_DST.Normal(μ, σ)]
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :branch,
            "cmp_id"  => (data["branch"][br]["index"], data["branch"][br]["t_bus"], data["branch"][br]["f_bus"]),
            "var"     => :ca,
            "dst"     => dst,
            "direction" => :to
        )
        m_idx+=1
    end
end


function add_cr_fr!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (br, branch) in pf_res["solution"]["branch"]
        σ = maximum([abs(branch["crf"]*σ_dict["cr"]/3), σ_dict["sigma_min"]/3])
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["crf"], σ), ) : branch["crf"]
        dst = [_DST.Normal(μ, σ)]
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :branch,
            "cmp_id"  => (data["branch"][br]["index"], data["branch"][br]["f_bus"], data["branch"][br]["t_bus"]),
            "var"     => :cr,
            "dst"     => dst,
            "direction" => :from
        )
        m_idx+=1
    end
end

function add_ci_fr!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (br, branch) in pf_res["solution"]["branch"]
        σ = maximum([abs(branch["cif"]*σ_dict["ci"]/3), σ_dict["sigma_min"]/3])
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["cif"], σ), ) : branch["cif"]
        dst = [_DST.Normal(μ, σ)]
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :branch,
            "cmp_id"  => (data["branch"][br]["index"], data["branch"][br]["f_bus"], data["branch"][br]["t_bus"]),
            "var"     => :ci,
            "dst"     => dst,
            "direction" => :from
        )
        m_idx+=1
    end
end

function add_cr_to!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (br, branch) in pf_res["solution"]["branch"]
        σ = maximum([abs(branch["crt"]*σ_dict["cr"]/3), σ_dict["sigma_min"]/3])
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["crt"], σ), ) : branch["crt"]
        dst = [_DST.Normal(μ, σ)]
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :branch,
            "cmp_id"  => (data["branch"][br]["index"], data["branch"][br]["t_bus"], data["branch"][br]["f_bus"]),
            "var"     => :cr,
            "dst"     => dst,
            "direction" => :to
        )
        m_idx+=1
    end
end



function add_ci_to!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (br, branch) in pf_res["solution"]["branch"]
        σ = maximum([abs(branch["cit"]*σ_dict["ci"]/3), σ_dict["sigma_min"]/3])
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["cit"], σ), ) : branch["cit"]
        dst = [_DST.Normal(μ, σ)]
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :branch,
            "cmp_id"  => (data["branch"][br]["index"], data["branch"][br]["t_bus"], data["branch"][br]["f_bus"]),
            "var"     => :ci,
            "dst"     => dst,
            "direction" => :to
        )
        m_idx+=1
    end
end



# DC NODAL MEASUREMENTS

function add_vdcm!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (b, bus) in pf_res["solution"]["busdc"]
        σ = [maximum([abs(bus["vm"][i]* σ_dict["vdcm"]/3), σ_dict["sigma_min"]/3]) for i in 1:3]
        μ = sample_error ? [_RAN.rand(_DST.Normal(bus["vm"][i], σ[i]), ) for i in 1:3] : bus["vm"]
        dst = [_DST.Normal(μ[i], σ[i]) for i in 1:3]
        data["meas"]["$m_idx"] = Dict(
            "cmp"       => :busdc,
            "cmp_id"    => parse(Int, b),
            "var"       => :vdcm,
            "dst"       => dst
        )
        m_idx+=1
    end
end

# DC FLOW MEASUREMENTS

function add_i_dcgrid_fr!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (b, branch) in pf_res["solution"]["branchdc"]
        c = length(branch["i_from"])
        σ =[maximum([abs(branch["i_from"][i]*σ_dict["i_dcgrid"]/3), σ_dict["sigma_min"]/3]) for i in 1:c]
        μ = sample_error ? [_RAN.rand(_DST.Normal(branch["i_from"][i], σ[i]), ) for i in 1:c] : branch["i_from"][1:c]
        dst = [_DST.Normal(μ[i], σ[i]) for i in 1:c]
        data["meas"]["$m_idx"] = Dict(
            "cmp"       => :branchdc,
            "cmp_id"    => (data["branchdc"][b]["index"], data["branchdc"][b]["fbusdc"], data["branchdc"][b]["tbusdc"]), 
            "var"       => :i_dcgrid,
            "dst"       => dst,
            "direction" => :from
        )
        m_idx+=1
    end
end

function add_i_dcgrid_to!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (b, branch) in pf_res["solution"]["branchdc"]
        c = length(branch["i_to"])
        σ = [maximum([abs(branch["i_to"][i]*σ_dict["i_dcgrid"]/3), σ_dict["sigma_min"]/3]) for i in 1:c]
        μ = sample_error ? [_RAN.rand(_DST.Normal(branch["i_to"][i], σ[i]), ) for i in 1:c] : branch["i_to"][1:c]
        dst = [_DST.Normal(μ[i], σ[i]) for i in 1:c]
        data["meas"]["$m_idx"] = Dict(
            "cmp"       => :branchdc,
            "cmp_id"    => (data["branchdc"][b]["index"], data["branchdc"][b]["tbusdc"], data["branchdc"][b]["fbusdc"]),
            "var"       => :i_dcgrid,
            "dst"       => dst,
            "direction" => :to
        )
        m_idx+=1
    end
end

function add_p_dc_fr!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (b, branch) in pf_res["solution"]["branchdc"]
        c = length(branch["pf"])
        μ = sample_error ? [_RAN.rand(_DST.Normal(branch["pf"][i], abs(branch["pf"][i])*σ_dict["p_dc"]/3), ) for i in 1:c] : branch["pf"][1:c]
        dst = [_DST.Normal(μ[i], σ_dict["p_dc"]) for i in 1:c]
        data["meas"]["$m_idx"] = Dict(
            "cmp"       => :branchdc,
            "cmp_id"    => (data["branchdc"][b]["index"], data["branchdc"][b]["fbusdc"], data["branchdc"][b]["tbusdc"]),
            "var"       => :p_dcgrid,
            "dst"       => dst,
            "direction" => :from
        )
        m_idx+=1
    end
end

function add_p_dc_to!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (b, branch) in pf_res["solution"]["branchdc"]
        c = length(branch["pt"])
        μ = sample_error ? [_RAN.rand(_DST.Normal(branch["pt"][i], abs(branch["pt"][i])*σ_dict["p_dc"]/3), ) for i in 1:c] : branch["pt"][1:c]
        dst = [_DST.Normal(μ[i], σ_dict["p_dc"]) for i in 1:c]
        data["meas"]["$m_idx"] = Dict(
            "cmp"       => :branchdc,
            "cmp_id"    => (data["branchdc"][b]["index"], data["branchdc"][b]["tbusdc"], data["branchdc"][b]["fbusdc"]),
            "var"       => :p_dcgrid,
            "dst"       => dst,
            "direction" => :to
        )
        m_idx+=1
    end
end



# Converter measurements AC part

function add_pconv!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (c, conv) in pf_res["solution"]["convdc"]
        nconvs=length(conv["pconv"])
        σ = [maximum([abs(conv["pconv"][i]* σ_dict["pconv"]/3), σ_dict["sigma_min"]/3]) for i in 1:nconvs]
        μ = sample_error ? [_RAN.rand(_DST.Normal(conv["pconv"][i], σ[i]), ) for i in 1:nconvs] : conv["pconv"]
        dst = [_DST.Normal(μ[i], σ[i]) for i in 1:nconvs]
        data["meas"]["$m_idx"] = Dict(
            "cmp"       => :convdc,
            "cmp_id"    => parse(Int, c),
            "var"       => :pconv_ac,
            "dst"       => dst
        )
        m_idx+=1
    end
end

function add_qconv!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (c, conv) in pf_res["solution"]["convdc"]
        nconvs=length(conv["qconv"])
        σ = [maximum([abs(conv["qconv"][i]* σ_dict["qconv"]/3), σ_dict["sigma_min"]/3]) for i in 1:nconvs]
        μ = sample_error ? [_RAN.rand(_DST.Normal(conv["qconv"][i], σ[i]), ) for i in 1:nconvs] : conv["qconv"]
        dst = [_DST.Normal(μ[i], σ[i]) for i in 1:nconvs]
        data["meas"]["$m_idx"] = Dict(
            "cmp"       => :convdc,
            "cmp_id"    => parse(Int, c),
            "var"       => :qconv_ac,
            "dst"       => dst
        )
        m_idx+=1
    end
end

function add_vmconv!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (c, conv) in pf_res["solution"]["convdc"]
        nconvs=length(conv["vmconv"])
        σ = [maximum([abs(conv["vmconv"][i]* σ_dict["vmconv"]/3), σ_dict["sigma_min"]/3]) for i in 1:nconvs]
        μ = sample_error ? [_RAN.rand(_DST.Normal(conv["vmconv"][i], σ[i]), ) for i in 1:nconvs] : conv["vmconv"]
        dst = [_DST.Normal(μ[i], σ[i]) for i in 1:nconvs]
        data["meas"]["$m_idx"] = Dict(
            "cmp"       => :convdc,
            "cmp_id"    => parse(Int, c),
            "var"       => :vmc,
            "dst"       => dst
        )
        m_idx+=1
    end
end




function add_ppr_fr!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (c, conv) in pf_res["solution"]["convdc"]
        nconvs=length(conv["ppr_fr"])
        σ = [maximum([abs(conv["ppr_fr"][i]* σ_dict["ppr_fr"]/3), σ_dict["sigma_min"]/3]) for i in 1:nconvs]
        μ = sample_error ? [_RAN.rand(_DST.Normal(conv["ppr_fr"][i], σ[i]), ) for i in 1:nconvs] : conv["ppr_fr"]
        dst = [_DST.Normal(μ[i], σ[i]) for i in 1:nconvs]
        data["meas"]["$m_idx"] = Dict(
            "cmp"       => :convdc,
            "cmp_id"    => parse(Int, c),
            "var"       => :pconv_pr_fr,
            "dst"       => dst
        )
        m_idx+=1
    end
end

function add_qpr_fr!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (c, conv) in pf_res["solution"]["convdc"]
        nconvs=length(conv["qpr_fr"])
        σ = [maximum([abs(conv["qpr_fr"][i]* σ_dict["qpr_fr"]/3), σ_dict["sigma_min"]/3]) for i in 1:nconvs]
        μ = sample_error ? [_RAN.rand(_DST.Normal(conv["qpr_fr"][i], σ[i]), ) for i in 1:nconvs] : conv["qpr_fr"]
        dst = [_DST.Normal(μ[i], σ[i]) for i in 1:nconvs]
        data["meas"]["$m_idx"] = Dict(
            "cmp"       => :convdc,
            "cmp_id"    => parse(Int, c),
            "var"       => :qconv_pr_fr,
            "dst"       => dst
        )
        m_idx+=1
    end
end


function add_vmfilt!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (c, conv) in pf_res["solution"]["convdc"]
        nconvs=length(conv["vmfilt"])
        σ = [maximum([abs(conv["vmfilt"][i]* σ_dict["vmfilt"]/3), σ_dict["sigma_min"]/3]) for i in 1:nconvs]
        μ = sample_error ? [_RAN.rand(_DST.Normal(conv["vmfilt"][i], σ[i]), ) for i in 1:nconvs] : conv["vmfilt"]
        dst = [_DST.Normal(μ[i], σ[i]) for i in 1:nconvs]
        data["meas"]["$m_idx"] = Dict(
            "cmp"       => :convdc,
            "cmp_id"    => parse(Int, c),
            "var"       => :vmf,
            "dst"       => dst
        )
        m_idx+=1
    end
end


function add_pgrid!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (c, conv) in pf_res["solution"]["convdc"]
        nconvs=length(conv["pgrid"])
        σ = [maximum([abs(conv["pgrid"][i]* σ_dict["pgrid"]/3), σ_dict["sigma_min"]/3]) for i in 1:nconvs]
        μ = sample_error ? [_RAN.rand(_DST.Normal(conv["pgrid"][i], σ[i]), ) for i in 1:nconvs] : conv["pgrid"]
        dst = [_DST.Normal(μ[i], σ[i]) for i in 1:nconvs]
        data["meas"]["$m_idx"] = Dict(
            "cmp"       => :convdc,
            "cmp_id"    => parse(Int, c),
            "var"       => :pconv_tf_fr,
            "dst"       => dst
        )
        m_idx+=1
    end
end

function add_qgrid!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (c, conv) in pf_res["solution"]["convdc"]
        nconvs=length(conv["qgrid"])
        σ = [maximum([abs(conv["qgrid"][i]* σ_dict["qgrid"]/3), σ_dict["sigma_min"]/3]) for i in 1:nconvs]
        μ = sample_error ? [_RAN.rand(_DST.Normal(conv["qgrid"][i], σ[i]), ) for i in 1:nconvs] : conv["qgrid"]
        dst = [_DST.Normal(μ[i], σ[i]) for i in 1:nconvs]
        data["meas"]["$m_idx"] = Dict(
            "cmp"       => :convdc,
            "cmp_id"    => parse(Int, c),
            "var"       => :qconv_tf_fr,
            "dst"       => dst
        )
        m_idx+=1
    end
end


function add_ptf_to!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (c, conv) in pf_res["solution"]["convdc"]
        nconvs=length(conv["ptf_to"])
        σ = [maximum([abs(conv["ptf_to"][i]* σ_dict["ptf_to"]/3), σ_dict["sigma_min"]/3]) for i in 1:nconvs]
        μ = sample_error ? [_RAN.rand(_DST.Normal(conv["ptf_to"][i], σ[i]), ) for i in 1:nconvs] : conv["ptf_to"]
        dst = [_DST.Normal(μ[i], σ[i]) for i in 1:nconvs]
        data["meas"]["$m_idx"] = Dict(
            "cmp"       => :convdc,
            "cmp_id"    => parse(Int, c),
            "var"       => :pconv_tf_to,
            "dst"       => dst
        )
        m_idx+=1
    end
end

function add_qtf_to!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (c, conv) in pf_res["solution"]["convdc"]
        nconvs=length(conv["qtf_to"])
        σ = [maximum([abs(conv["qtf_to"][i]* σ_dict["qtf_to"]/3), σ_dict["sigma_min"]/3]) for i in 1:nconvs]
        μ = sample_error ? [_RAN.rand(_DST.Normal(conv["qtf_to"][i], σ[i]), ) for i in 1:nconvs] : conv["qtf_to"]
        dst = [_DST.Normal(μ[i], σ[i]) for i in 1:nconvs]
        data["meas"]["$m_idx"] = Dict(
            "cmp"       => :convdc,
            "cmp_id"    => parse(Int, c),
            "var"       => :qconv_tf_to,
            "dst"       => dst
        )
        m_idx+=1
    end
end


# Converter measurements DC part



function add_pdc!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (c, conv) in pf_res["solution"]["convdc"]
        nterm=length(conv["pdc"])
        σ = [maximum([abs(conv["pdc"][i]* σ_dict["pdc"]/3), σ_dict["sigma_min"]/3]) for i in 1:nterm]
        μ = sample_error ? [_RAN.rand(_DST.Normal(conv["pdc"][i], σ[i]), ) for i in 1:nterm] : conv["pdc"]
        dst = [_DST.Normal(μ[i], σ[i]) for i in 1:nterm]
        data["meas"]["$m_idx"] = Dict(
            "cmp"       => :convdc,
            "cmp_id"    => parse(Int, c),
            "var"       => :pconv_dc,
            "dst"       => dst
        )
        m_idx+=1
    end
end


function add_pdcg!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (c, conv) in pf_res["solution"]["convdc"]
        nterm=length(conv["pdcg"])
        σ = [maximum([abs(conv["pdcg"][i]* σ_dict["pdcg"]/3), σ_dict["sigma_min"]/3]) for i in 1:nterm]
        μ = sample_error ? [_RAN.rand(_DST.Normal(conv["pdcg"][i], σ[i]), ) for i in 1:nterm] : conv["pdcg"]
        dst = [_DST.Normal(μ[i], σ[i]) for i in 1:nterm]
        data["meas"]["$m_idx"] = Dict(
            "cmp"       => :convdc,
            "cmp_id"    => parse(Int, c),
            "var"       => :pconv_dcg,
            "dst"       => dst
        )
        m_idx+=1
    end
end

function add_pdcg_shunt!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (c, conv) in pf_res["solution"]["convdc"]
        σ = maximum([abs(σ_dict["pdcg_shunt"]*conv["pdcg_shunt"]/3), σ_dict["sigma_min"]/3])
        μ = sample_error ? _RAN.rand(_DST.Normal(conv["pdcg_shunt"], σ), ) : conv["pdcg_shunt"]
        dst = [_DST.Normal(μ, σ)]
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :convdc,
            "cmp_id"  => parse(Int, c),
            "var"     => :pconv_dcg_shunt,
            "dst"     => dst
        )
        m_idx+=1
    end
end

function add_mconv!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (c, conv) in pf_res["solution"]["convdc"]
        nconvs=length(conv["m"])
        σ = [maximum([abs(conv["m"][i]* σ_dict["mconv"]/3), σ_dict["sigma_min"]/3]) for i in 1:nconvs]
        μ = sample_error ? [_RAN.rand(_DST.Normal(conv["m"][i], σ[i]), ) for i in 1:nconvs] : conv["m"]
        dst = [_DST.Normal(μ[i], σ[i]) for i in 1:nconvs]
        data["meas"]["$m_idx"] = Dict(
            "cmp"       => :convdc,
            "cmp_id"    => parse(Int, c),
            "var"       => :mconv,
            "dst"       => dst
        )
        m_idx+=1
    end
end