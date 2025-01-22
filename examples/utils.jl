
function generate_data_basic_acdcse(data_pf, data_se, nlp_optimizer; sample_error::Bool = true)
    
    result = _PMMCDC.solve_mcdcopf(data_pf, _PM.ACPPowerModel, nlp_optimizer)
    _ACDCSE.get_dc_power!(result, data_pf)
    
    σ_dict = Dict{String, Float64}(
        "vm"       => 0.005, # (AC) voltage magnitude
        "va"       => 0.005, # (AC) voltage angle
        "pg"       => 0.01, # Active power injection from (AC) generator
        "qg"       => 0.01, # Reactive power injection from (AC) generator
        "pd"       => 0.02, # Active power injection from (AC) load
        "qd"       => 0.02, # Rective power injection from (AC) load
        "p_ac"     => 0.015, # Active power flow on an AC branch   (either injection or flow, really)
        "q_ac"     => 0.015, # Reactive power flow on an AC branch (either injection or flow, really)
        "vdcm"     => 0.01, # DC voltage magnitude
        "p_dc"     => 0.01, # DC power (either injection or flow)
        "i_dcgrid" => 0.01,  # DC current (either injection or flow)
        "cm"       => 0.01
    )

    #_ACDCSE.add_current_flows_to_se_result!(result)

    _ACDCSE.powerflow2measurements!(data_se, result, σ_dict, sample_error = sample_error, measurements = ["vm", "va", "p_to", "q_to", "vdcm", "pg", 
                                "pd", "qg", "qd"])#, "p_dc_to", "p_dc_fr"])

    _ACDCSE.prepare_data_for_se_default!(data_se, exceptions = [1])

    return result, σ_dict, data_se
end

function introduce_error_conv_losses!(data::Dict)
    for (_, conv) in data["convdc"]
        conv["LossA"] = conv["LossA"].*1.1 
        conv["LossB"] = conv["LossB"].*1.1
        conv["LossCinv"] = conv["LossCinv"].*1.1
    end
end