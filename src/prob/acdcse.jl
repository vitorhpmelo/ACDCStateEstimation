export solve_acdcse

"""
    solve_acdcse(data, model_type, optimizer; <keyword arguments>)

Build and solve the AC/DC SE problem over a hybrid AC/DC network, using a multi-conductor model for the DC part.

Input can be a Matpower `file` or a `data` dictionary.
The SE problem being built is the one defined in `build_acdcse`.
Keyword arguments, if any, are forwarded to `PowerModels.solve_model`.
"""
function solve_acdcse end

function solve_acdcse(file::String, model_type::Type, optimizer; kwargs...)
    data = _PMMCDC.parse_file(file)
    return solve_acdcse(data, model_type, optimizer; ref_extensions=[_PMMCDC.add_ref_dcgrid!], kwargs...)
end

function solve_acdcse(data::Dict{String,Any}, model_type::Type, optimizer; kwargs...)
    if haskey(data["se_settings"], "criterion")
        assign_unique_individual_criterion!(data)
    end
    if !haskey(data["se_settings"], "rescaler")
        data["se_settings"]["rescaler"] = 1
        @warn "Rescaler set to default value, edit data dictionary if you wish to change it."
    end
    return _PM.solve_model(data, model_type, optimizer, build_acdcse; ref_extensions=[_PMMCDC.add_ref_dcgrid!], kwargs...)
end

function build_acdcse(pm::_PM.AbstractPowerModel)

    ## VARIABLES
    # PowerModels variables
    _PM.variable_bus_voltage(pm, bounded=true)
    _PM.variable_gen_power(pm, bounded = true)
    _PM.variable_branch_power(pm, report = false)

    # PowerModelsMCDC variables 
    # _PMMCDC.variable_mc_active_dcbranch_flow(pm, bounded=true)
    _PMMCDC.variable_mcdcgrid_voltage_magnitude(pm, bounded=true)
    _PMMCDC.variable_mcdc_converter(pm, bounded=true)
    _PMMCDC.variable_mc_dcbranch_current(pm, bounded=true)
    
    # state estimation variables 
    variable_load(pm)
    variable_measurement(pm)
    variable_residual(pm)

    ## CONSTRAINTS
    # PowerModels constraints
    for i in _PM.ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end
    for i in _PM.ids(pm, :branch)
        _PM.constraint_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)
        _PM.constraint_voltage_angle_difference(pm, i) # angle difference across transformer and reactor - useful for LPAC if available?
        _PM.constraint_thermal_limit_from(pm, i)
        _PM.constraint_thermal_limit_to(pm, i)
    end

    # PowerModelsMCDC constraints
    for i in _PM.ids(pm, :busdc)
        _PMMCDC.constraint_kcl_shunt_dcgrid(pm, i) # this will have to be replaced if we want to allow dc loads
    end
    for i in _PM.ids(pm, :branchdc)
        _PMMCDC.constraint_ohms_dc_branch(pm, i)
    end
    for i in _PM.ids(pm, :convdc)
        _PMMCDC.constraint_converter_losses(pm, i)
        _PMMCDC.constraint_converter_current(pm, i)
        _PMMCDC.constraint_converter_dc_current(pm, i)
        _PMMCDC.constraint_conv_transformer(pm, i)
        _PMMCDC.constraint_conv_reactor(pm, i)
        _PMMCDC.constraint_conv_filter(pm, i)
        if pm.ref[:it][_PM.pm_it_sym][:nw][_PM.nw_id_default][:convdc][i]["islcc"] == 1
            _PMMCDC.constraint_conv_firing_angle(pm, i)
        end
    end
    _PMMCDC.constraint_converter_dc_ground_shunt_ohm(pm)

    # State estimation constraints
    for i in _PM.ids(pm, :bus)
        constraint_kcl_shunt_se(pm, i)
    end

    for (i, _) in _PM.ref(pm, :meas)
        constraint_residual(pm, i)
    end

    ## OBJECTIVE
    objective_minimize_residuals(pm)

end