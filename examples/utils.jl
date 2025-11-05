


function generate_data_basic_acdcse(data_pf, data_se, nlp_optimizer; sample_error::Bool = true)
    
    result = _PMMCDC.solve_mcdcopf(data_pf, _PM.ACPPowerModel, nlp_optimizer)
    "branchdc" ∈ keys(data_pf) ? _ACDCSE.get_dc_power!(result, data_pf) : nothing
    "branch" ∈ keys(data_pf) ? _ACDCSE.get_ac_currents!(result, data_pf) : nothing
    "convdc" ∈ keys(data_pf) ? _ACDCSE.get_convacdc_m!(result, data_pf) : nothing

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
        "cm"       => 0.005,  # (AC)
        "ca" => 0.005,
        "cr"=> 0.005,  # (AC) current real
        "ci"=> 0.005,
        "mconv"=> 0.01, # converter Vac/Vdc
        "pconv" => 0.01, # active power into the converters
        "qconv" => 0.01, # reactive power into the converters
        "vmconv" => 0.01,
        "vmfilt" => 0.01,
        "ppr_fr"=> 0.01,
        "qpr_fr" =>0.01,
        "vmfilt" => 0.01,
        "pgrid"=> 0.01,
        "qgrid"=> 0.01,
        "ptf_to"=> 0.01,
        "qtf_to"=> 0.01,
        "pdc"=>0.01,
        "pdcg"=>0.01,
        "pdcg_shunt"=>0.01,
        "sigma_min"=>1e-10

    ) #defines the maximum sgima of measurements ?

    #_ACDCSE.add_current_flows_to_se_result!(result)


    # _ACDCSE.powerflow2measurements!(data_se, result, σ_dict, sample_error = sample_error, measurements = ["vm", "va", "p_to", "q_to", "vdcm", "pg", 
                                # "pd", "qg", "qd","i_dcgrid_to","i_dcgrid_fr","cm_fr"])#, "p_dc_to", "p_dc_fr"])

    # _ACDCSE.powerflow2measurements!(data_se, result, σ_dict, sample_error = sample_error, measurements = ["vm", "va", "p_to", "q_to","p_fr","q_fr",  "pg", 
                                # "pd", "qg", "qd","cr_fr","cr_to","ci_fr","ci_to"])#, "p_dc_to", "p_dc_fr"])
    _ACDCSE.powerflow2measurements!(data_se, result, σ_dict, sample_error = sample_error, measurements = ["vm", "va", "p_to", "q_to","p_fr","q_fr", "vdcm", "pg", 
                                "pd", "qg", "qd","mconv","pconv","qconv","vmconv","ppr_fr","qpr_fr","vmfilt","pgrid","qgrid","ptf_to","qtf_to","pdc","p_dc_fr","p_dc_to","i_dcgrid_to","i_dcgrid_fr","cr_to","cr_fr","ci_to","ci_fr","pdc"])#, "p_dc_to", "p_dc_fr"])
# 
    _ACDCSE.prepare_data_for_se_default!(data_se, exceptions = [1]) #transforms all PV buses in PQ , and remove slack buses

    return result, σ_dict, data_se
end



function generate_data_basic_acdcse(data_pf, data_se, nlp_optimizer,reference::Vector{Int}, measurements::Vector{String}; sample_error::Bool = true)

    result = _PMMCDC.solve_mcdcopf(data_pf, _PM.ACPPowerModel, nlp_optimizer)
    "branchdc" ∈ keys(data_pf) ? _ACDCSE.get_dc_power!(result, data_pf) : nothing
    "branch" ∈ keys(data_pf) ? _ACDCSE.get_ac_currents!(result, data_pf) : nothing
    "convdc" ∈ keys(data_pf) ? _ACDCSE.get_convacdc_m!(result, data_pf) : nothing

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
        "cm"       => 0.005,  # (AC)
        "ca" => 0.005,
        "cr"=> 0.005,  # (AC) current real
        "ci"=> 0.005,
        "mconv"=> 0.01, # converter Vac/Vdc
        "pconv" => 0.01, # active power into the converters
        "qconv" => 0.01, # reactive power into the converters
        "vmconv" => 0.01,
        "vmfilt" => 0.01,
        "ppr_fr"=> 0.01,
        "qpr_fr" =>0.01,
        "vmfilt" => 0.01,
        "pgrid"=> 0.01,
        "qgrid"=> 0.01,
        "ptf_to"=> 0.01,
        "qtf_to"=> 0.01,
        "pdc"=>0.01,
        "pdcg"=>0.01,
        "pdcg_shunt"=>0.01,
        "sigma_min"=>1e-5

    ) #defines the maximum sgima of measurements ?

    #_ACDCSE.add_current_flows_to_se_result!(result)


    _ACDCSE.powerflow2measurements!(data_se, result, σ_dict, sample_error = sample_error, measurements = measurements)#, "p_dc_to", "p_dc_fr"])
# 
    _ACDCSE.prepare_data_for_se_default!(data_se, exceptions = reference) #transforms all PV buses in PQ , and remove slack buses

    return result, σ_dict, data_se
end


function generate_data_basic_acdcse(data_pf, data_se, nlp_optimizer, meas_types; sample_error::Bool = true)
    
    result = _PMMCDC.solve_mcdcopf(data_pf, _PM.ACPPowerModel, nlp_optimizer)
    "branchdc" ∈ keys(data_pf) ? _ACDCSE.get_dc_power!(result, data_pf) : nothing
    "branch" ∈ keys(data_pf) ? _ACDCSE.get_ac_currents!(result, data_pf) : nothing
    "convdc" ∈ keys(data_pf) ? _ACDCSE.get_convacdc_m!(result, data_pf) : nothing

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
        "cm"       => 0.005,  # (AC)
        "ca" => 0.005,
        "cr"=> 0.005,  # (AC) current real
        "ci"=> 0.005,
        "mconv"=> 0.01, # converter Vac/Vdc
        "pconv" => 0.01, # active power into the converters
        "qconv" => 0.01, # reactive power into the converters
        "vmconv" => 0.01,
        "vmfilt" => 0.01,
        "ppr_fr"=> 0.01,
        "qpr_fr" =>0.01,
        "vmfilt" => 0.01,
        "pgrid"=> 0.01,
        "qgrid"=> 0.01,
        "ptf_to"=> 0.01,
        "qtf_to"=> 0.01,
        "pdc"=>0.01,
        "pdcg"=>0.01,
        "pdcg_shunt"=>0.01,
        "sigma_min"=>1e-10

    ) #defines the maximum sgima of measurements ?

    measurements=[]
    if meas_types == "all"
       measurements = ["vm", "va", "p_to", "q_to","p_fr","q_fr", "vdcm", "pg", 
                                "pd", "qg", "qd","mconv","pconv","qconv","vmconv","ppr_fr","qpr_fr","vmfilt","pgrid","qgrid","ptf_to","qtf_to","pdc","p_dc_fr","p_dc_to","i_dcgrid_to","i_dcgrid_fr","cr_to","cr_fr","ci_to","ci_fr","pdc"]
    elseif meas_types == "ac_only"
        measurements = ["vm", "va", "p_to", "q_to","p_fr","q_fr",  "pg", 
                        "pd", "qg", "qd","cr_fr","cr_to","ci_fr","ci_to"]
    elseif meas_types == "no_branch"
        measurements = ["vm", "va", "vdcm", "pg", 
                                "pd", "qg", "qd","mconv","pconv","qconv","vmconv","ppr_fr","qpr_fr","vmfilt","pgrid","qgrid","ptf_to","qtf_to","pdc","p_dc_fr","p_dc_to","i_dcgrid_to","i_dcgrid_fr","pdc"]
    else
        error("Unknown measurement type: $meas_types")
    end

    _ACDCSE.powerflow2measurements!(data_se, result, σ_dict, sample_error = sample_error, measurements = measurements)

    _ACDCSE.prepare_data_for_se_default!(data_se, exceptions = [1,6,11]) #transforms all PV buses in PQ , and remove slack buses

    return result, σ_dict, data_se
end

function generate_data_basic_acdcse(data_pf, data_se, nlp_optimizer, meas_types,reference::Vector{Int64}; sample_error::Bool = true)
    
    result = _PMMCDC.solve_mcdcopf(data_pf, _PM.ACPPowerModel, nlp_optimizer)
    "branchdc" ∈ keys(data_pf) ? _ACDCSE.get_dc_power!(result, data_pf) : nothing
    "branch" ∈ keys(data_pf) ? _ACDCSE.get_ac_currents!(result, data_pf) : nothing
    "convdc" ∈ keys(data_pf) ? _ACDCSE.get_convacdc_m!(result, data_pf) : nothing

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
        "cm"       => 0.005,  # (AC)
        "ca" => 0.005,
        "cr"=> 0.005,  # (AC) current real
        "ci"=> 0.005,
        "mconv"=> 0.01, # converter Vac/Vdc
        "pconv" => 0.01, # active power into the converters
        "qconv" => 0.01, # reactive power into the converters
        "vmconv" => 0.01,
        "vmfilt" => 0.01,
        "ppr_fr"=> 0.01,
        "qpr_fr" =>0.01,
        "vmfilt" => 0.01,
        "pgrid"=> 0.01,
        "qgrid"=> 0.01,
        "ptf_to"=> 0.01,
        "qtf_to"=> 0.01,
        "pdc"=>0.01,
        "pdcg"=>0.01,
        "pdcg_shunt"=>0.01,
        "sigma_min"=>1e-10

    ) #defines the maximum sgima of measurements ?

    measurements=[]
    if meas_types == "all"
       measurements = ["vm", "va", "p_to", "q_to","p_fr","q_fr", "vdcm", "pg", 
                                "pd", "qg", "qd","mconv","pconv","qconv","vmconv","ppr_fr","qpr_fr","vmfilt","pgrid","qgrid","ptf_to","qtf_to","pdc","p_dc_fr","p_dc_to","i_dcgrid_to","i_dcgrid_fr","cr_to","cr_fr","ci_to","ci_fr","pdc"]
    elseif meas_types == "ac_only"
        measurements = ["vm", "va", "p_to", "q_to","p_fr","q_fr",  "pg", 
                        "pd", "qg", "qd","cr_fr","cr_to","ci_fr","ci_to"]
    elseif meas_types == "no_branch"
        measurements = ["vm", "va", "vdcm", "pg", 
                                "pd", "qg", "qd","mconv","pconv","qconv","vmconv","ppr_fr","qpr_fr","vmfilt","pgrid","qgrid","ptf_to","qtf_to","pdc","p_dc_fr","p_dc_to","i_dcgrid_to","i_dcgrid_fr","pdc"]
    else
        error("Unknown measurement type: $meas_types")
    end

    _ACDCSE.powerflow2measurements!(data_se, result, σ_dict, sample_error = sample_error, measurements = measurements)

    _ACDCSE.prepare_data_for_se_default!(data_se, exceptions = reference) #transforms all PV buses in PQ , and remove slack buses

    return result, σ_dict, data_se
end

function introduce_error_conv_losses!(data::Dict)
    for (_, conv) in data["convdc"]
        conv["LossA"] = conv["LossA"].*1.1 
        conv["LossB"] = conv["LossB"].*1.1
        conv["LossCinv"] = conv["LossCinv"].*1.1
    end
end

function introduce_error_conv_losses2!(data::Dict,per=1.1)
    for (_, conv) in data["convdc"]
        conv["LossA"] = conv["LossA"].*per
        conv["LossB"] = conv["LossB"].*per
        conv["LossCinv"] = conv["LossCinv"].*per
    end
end
