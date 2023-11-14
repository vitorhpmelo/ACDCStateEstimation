"""
From power flow results `pf_res`, creates synthetic measurements
and adds them to the network `data`, which is then ready for state estimation
`σ_dict` is a dictionary with the standard deviation of the different measurement types
"""
function powerflow2measurements!(data::Dict, pf_res::Dict, σ_dict::Dict; sample_error::Bool = true)
    data["meas"] = Dict{String, Any}()
    m_idx = 1
    # ac bus voltage magnitude
    for (b, bus) in pf_res["solution"]["bus"]
        μ = sample_error ? _RAN.rand(_DST.Normal(bus["vm"], σ_dict["vm"]), ) : bus["vm"]
        dst = _DST.Normal(μ, σ_dict["vm"])
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :bus,
            "cmp_id"  => parse(Int, b),
            "var"     => :vm,
            "dst"     => dst
        )
        m_idx += 1
    end
    # ac branch: active and reactive power from
    for (br, branch) in pf_res["solution"]["branch"]
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["pf"], σ_dict["p_ac"]), ) : branch["pf"]
        dst = _DST.Normal(μ, σ_dict["p_ac"])
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :branch,
            "cmp_id"  => (data["branch"][br]["index"], data["branch"][br]["f_bus"], data["branch"][br]["t_bus"]),
            "var"     => :p,
            "dst"     => dst
        )
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["qf"], σ_dict["q_ac"]), ) : branch["qf"]
        dst = _DST.Normal(μ, σ_dict["q_ac"])
        data["meas"]["$(m_idx+1)"] = Dict(
            "cmp"     => :branch,
            "cmp_id"  => (data["branch"][br]["index"], data["branch"][br]["f_bus"], data["branch"][br]["t_bus"]),
            "var"     => :q,
            "dst"     => dst
        )
        m_idx += 2
    end
    # ac branch: active and reactive power to
    for (br, branch) in pf_res["solution"]["branch"]
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["pt"], σ_dict["p_ac"]), ) : branch["pt"]
        dst = _DST.Normal(μ, σ_dict["p_ac"])
        data["meas"]["$m_idx"] = Dict(
            "cmp"     => :branch,
            "cmp_id"  => (data["branch"][br]["index"], data["branch"][br]["t_bus"], data["branch"][br]["f_bus"]),
            "var"     => :p,
            "dst"     => dst
        )
        μ = sample_error ? _RAN.rand(_DST.Normal(branch["qt"], σ_dict["q_ac"]), ) : branch["qt"]
        dst = _DST.Normal(μ, σ_dict["q_ac"])
        data["meas"]["$(m_idx+1)"] = Dict(
            "cmp"     => :branch,
            "cmp_id"  => (data["branch"][br]["index"], data["branch"][br]["t_bus"], data["branch"][br]["f_bus"]),
            "var"     => :q,
            "dst"     => dst
        )
        m_idx += 2
    end
    # dc bus voltages (dc buses are always three-terminal in pmmcdc)
    for (b, bus) in pf_res["solution"]["busdc"]
        for conductor in 1:length(bus["vm"]) 
            μ = sample_error ? _RAN.rand(_DST.Normal(bus["vm"][conductor], σ_dict["vdcm"]), ) : bus["vm"][conductor]
            dst = _DST.Normal(μ, σ_dict["vdcm"])
            data["meas"]["$m_idx"] = Dict(
                "cmp"       => :bus,
                "cmp_id"    => parse(Int, b),
                "var"       => :vdcm,
                "dst"       => dst,
                "conductor" => conductor
            )
            m_idx += 1
        end
    end
end

# m_idx = isempty(data["meas"]) ? 1 : maximum(parse.(Int, keys(data["meas"])))