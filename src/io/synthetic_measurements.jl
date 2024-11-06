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

    # NODAL measurements (voltages and injections)
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
    # if "ca_fr" ∈ measurements add_ca_fr!(data, pf_res, σ_dict, sample_error) end # ac current angle from
    # if "ca_to" ∈ measurements add_ca_to!(data, pf_res, σ_dict, sample_error) end # ac current angle to
    # if "cr_fr" ∈ measurements add_cr_fr!(data, pf_res, σ_dict, sample_error) end # ac current rectangular real from
    # if "cr_to" ∈ measurements add_cr_to!(data, pf_res, σ_dict, sample_error) end # ac current rectangular real to
    # if "ci_fr" ∈ measurements add_ci_fr!(data, pf_res, σ_dict, sample_error) end # ac current rectangular imag from
    # if "ci_to" ∈ measurements add_ci_to!(data, pf_res, σ_dict, sample_error) end # ac current rectangular imag to

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

end

function add_vm!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))
    for (b, bus) in pf_res["solution"]["bus"]
        μ = sample_error ? _RAN.rand(_DST.Normal(bus["vm"], σ_dict["vm"]), ) : bus["vm"]
        dst = [_DST.Normal(μ, σ_dict["vm"])]
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
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))
    for (b, bus) in pf_res["solution"]["bus"]
        μ = sample_error ? _RAN.rand(_DST.Normal(bus["va"], σ_dict["va"]), ) : bus["va"]
        dst = [_DST.Normal(μ, σ_dict["va"])]
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
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["pf"], σ_dict["p_ac"]), ) : branch["pf"]
        dst = [_DST.Normal(μ, σ_dict["p_ac"])]
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
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["qf"], σ_dict["q_ac"]), ) : branch["qf"]
        dst = [_DST.Normal(μ, σ_dict["q_ac"])]
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
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["pt"], σ_dict["p_ac"]), ) : branch["pt"]
        dst = [_DST.Normal(μ, σ_dict["p_ac"])]
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
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["qt"], σ_dict["q_ac"]), ) : branch["qt"]
        dst = [_DST.Normal(μ, σ_dict["q_ac"])]
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
        μ = sample_error ? _RAN.rand(_DST.Normal(gen["pg"], σ_dict["pg"]), ) : gen["pg"]
        dst = [_DST.Normal(μ, σ_dict["pg"])]
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
        μ = sample_error ? _RAN.rand(_DST.Normal(gen["qg"], σ_dict["qg"]), ) : gen["qg"]
        dst = [_DST.Normal(μ, σ_dict["qg"])]
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
        μ = sample_error ? _RAN.rand(_DST.Normal(load["pd"], σ_dict["pd"]), ) : load["pd"]
        dst = [_DST.Normal(μ, σ_dict["pd"])]
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
        μ = sample_error ? _RAN.rand(_DST.Normal(load["qd"], σ_dict["qd"]), ) : load["qd"]
        dst = [_DST.Normal(μ, σ_dict["qd"])]
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
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["cmf"], σ_dict["cm"]), ) : branch["cmf"]
        dst = [_DST.Normal(μ, σ_dict["cm"])]
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
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["cmf"], σ_dict["cm"]), ) : branch["cmf"]
        dst = [_DST.Normal(μ, σ_dict["cm"])]
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

# DC NODAL MEASUREMENTS

function add_vdcm!(data, pf_res, σ_dict, sample_error)
    m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))+1
    for (b, bus) in pf_res["solution"]["busdc"]
        μ = sample_error ? [_RAN.rand(_DST.Normal(bus["vm"][i], σ_dict["vdcm"]), ) for i in 1:3] : bus["vm"]
        dst = [_DST.Normal(μ[i], σ_dict["vdcm"]) for i in 1:3]
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
        μ = sample_error ? [_RAN.rand(_DST.Normal(branch["i_from"][i], σ_dict["i_dcgrid"]), ) for i in 1:c] : branch["i_from"][1:c]
        dst = [_DST.Normal(μ[i], σ_dict["i_dcgrid"]) for i in 1:c]
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
        μ = sample_error ? [_RAN.rand(_DST.Normal(branch["i_to"][i], σ_dict["i_dcgrid"]), ) for i in 1:c] : branch["i_to"][1:c]
        dst = [_DST.Normal(μ[i], σ_dict["i_dcgrid"]) for i in 1:c]
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
        μ = sample_error ? [_RAN.rand(_DST.Normal(branch["pf"][i], σ_dict["p_dc"]), ) for i in 1:c] : branch["pf"][1:c]
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
        μ = sample_error ? [_RAN.rand(_DST.Normal(branch["pt"][i], σ_dict["p_dc"]), ) for i in 1:c] : branch["pt"][1:c]
        dst = [_DST.Normal(μ[i], σ_dict["i_dcgrid"]) for i in 1:c]
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