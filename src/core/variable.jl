"""
    variable_residual

This is the residual variable, which is later associated to the residual (in)equality constraint(s), depending on the chosen
state estimation criterion.
If bounded, the lower bound is set to zero, while the upper bound defaults to Inf, unless the user provides
a different value in the measurement dictionary.
"""
function variable_residual( pm::_PM.AbstractPowerModel;
                            nw::Int=_PM.nw_id_default, bounded::Bool = true, report::Bool = true  )

    connections = Dict(i => length(meas["dst"]) for (i,meas) in _PM.ref(pm, nw, :meas) ) #creates a dictionary with all meas connections (vector meas)

    res = _PM.var(pm, nw)[:res] = Dict(
        i => JuMP.@variable(pm.model,
        [c in 1:connections[i]], 
        base_name = "$(nw)_res_$(i)",
        start = _PMMCDC.comp_start_value(_PM.ref(pm, nw, :meas, i), "res_start", 0.0)
        ) for i in _PM.ids(pm, nw, :meas)
    )

    if bounded
        for i in _PM.ids(pm, nw, :meas), c in 1:connections[i]
            JuMP.set_lower_bound(res[i][c], 0.0)
            res_max = haskey(_PM.ref(pm, nw, :meas, i), "res_max") ? meas["res_max"] : 1e8
            JuMP.set_upper_bound(res[i][c], res_max)
        end
    end

    report && _PM.sol_component_value(pm, nw, :meas, :res, _PM.ids(pm, nw, :meas), res)
end

"""
    set_variable_start_value!(pm::_PM.AbstractPowerModel, var_symbol::Symbol, start_values::Dict)

Sets the initial value (start) of existing JuMP variables in the model.
- `var_symbol`: Symbol of the variable (e.g., :res, :prior, :pd, etc.)
- `start_values`: Dict mapping variable indices to start values (can be scalar or vector)
"""
function set_variable_start_value(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default)
    
    for i in _PM.ids(pm, nw, :bus)
        if haskey(_PM.ref(pm, nw, :bus, i), "vm_start") && haskey(_PM.ref(pm, nw, :bus, i), "va_start")
            vm_start = _PM.ref(pm, nw, :bus, i, "vm_start")
            va_start = _PM.ref(pm, nw, :bus, i, "va_start")

            JuMP.set_start_value(_PM.var(pm, nw, :vm, i), vm_start)
            JuMP.set_start_value(_PM.var(pm, nw, :va, i), va_start)
        end
    end
    for i in _PM.ids(pm, nw, :busdc)
        if haskey(_PM.ref(pm, nw, :busdc, i), "vm_start")
            vm_start = _PM.ref(pm, nw, :busdc, i, "vm_start")
            for c in 1:length(_PM.var(pm, nw, :vdcm, i))
                JuMP.set_start_value(_PM.var(pm, nw, :vdcm, i)[c], vm_start[c])
            end
        end
    end

    for i in _PM.ids(pm, nw, :convdc)
        if haskey(_PM.ref(pm, nw, :convdc, i), "vmf_start")
            vmf_start = _PM.ref(pm, nw, :convdc, i, "vmf_start")
            for c in 1:length(_PM.var(pm, nw, :vmf, i))
                JuMP.set_start_value(_PM.var(pm, nw, :vmf, i)[c], vmf_start[c])
            end
        end
        if haskey(_PM.ref(pm, nw, :convdc, i), "vaf_start")
            vaf_start = _PM.ref(pm, nw, :convdc, i, "vaf_start")
            for c in 1:length(_PM.var(pm, nw, :vaf, i))
                JuMP.set_start_value(_PM.var(pm, nw, :vaf, i)[c], vaf_start[c])
            end
        end
        if haskey(_PM.ref(pm, nw, :convdc, i), "vmc_start")
            vmc_start = _PM.ref(pm, nw, :convdc, i, "vmc_start")
            for c in 1:length(_PM.var(pm, nw, :vmc, i))
                JuMP.set_start_value(_PM.var(pm, nw, :vmc, i)[c], vmc_start[c])
            end
        end
        if haskey(_PM.ref(pm, nw, :convdc, i), "vac_start")
            vac_start = _PM.ref(pm, nw, :convdc, i, "vac_start")
            for c in 1:length(_PM.var(pm, nw, :vac, i))
                JuMP.set_start_value(_PM.var(pm, nw, :vac, i)[c], vac_start[c])
            end
        end

    end
end

function variable_prior( pm::_PM.AbstractPowerModel;
                            nw::Int=_PM.nw_id_default, bounded::Bool = true, report::Bool = true  ) #TODO

    connections = Dict(i => 1 for (i,meas) in _PM.ref(pm, nw, :prior) ) 

    prior = _PM.var(pm, nw)[:prior] = Dict(
        i => JuMP.@variable(pm.model, [c in 1:connections[i]],
        base_name = "$(nw)_prior_$(i)",
        start = _PMMCDC.comp_start_value(_PM.ref(pm, nw, :prior, i), "prior_start", 0.0) # set initial value
        ) for i in _PM.ids(pm, nw, :prior)
    )
    if bounded
        for i in _PM.ids(pm, nw, :prior) #might make sense to put an inferior bound
            JuMP.set_lower_bound(prior[i][1], -1.e1)
            JuMP.set_upper_bound(prior[i][1], 1.e1)
        end
    end

    report && _PM.sol_component_value(pm, nw, :prior, :prior, _PM.ids(pm, nw, :prior), prior) #TODO do it
end

"""
    variable_measurement
checks for every measurement if the measured
quantity belongs to the formulation's variable space. If not, the function
`create_conversion_constraint' is called, that adds a constraint that
associates the measured quantity to the formulation's variable space.
"""
function variable_measurement(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=false)
    for i in _PM.ids(pm, nw, :meas)
        msr_var = _PM.ref(pm, nw, :meas, i, "var")
        cmp_id = _PM.ref(pm, nw, :meas, i, "cmp_id")
        cmp_type = _PM.ref(pm, nw, :meas, i, "cmp")
        direction = haskey(_PM.ref(pm, nw, :meas, i), "direction") ? _PM.ref(pm, nw, :meas, i, "direction") : :none
        connections = get_active_connections(pm, nw, cmp_type, cmp_id)

        if no_conversion_needed(pm, msr_var, cmp_type)
            #no additional variable is created, it is already by default in the formulation
        else
            if haskey(_PM.var(pm, nw), msr_var)
                push!(_PM.var(pm, nw)[msr_var], cmp_id => JuMP.@variable(pm.model,
                    [c in connections], base_name="$(nw)_$(String(msr_var))_$(cmp_id)"))
            else
                _PM.var(pm, nw)[msr_var] = Dict(cmp_id => JuMP.@variable(pm.model,
                    [c in connections], base_name="$(nw)_$(String(msr_var))_$(cmp_id)"))
            end
            msr_type = assign_conversion_type_to_msr(pm, i, msr_var, cmp_type, direction; nw=nw)
            create_conversion_constraint(pm, _PM.var(pm, nw)[msr_var], msr_type; nw=nw)
        end
    end
end
"""
    variable_load in terms of power, for ACR and ACP
"""
function variable_load(pm::_PM.AbstractPowerModel; kwargs...)
    # NB: currently, only (single-phase eq.) loads on the ac side are supported
    variable_load_active(pm; kwargs...)
    variable_load_reactive(pm; kwargs...)
end

function variable_load_active(pm::_PM.AbstractPowerModel;
                                 nw::Int=_PM.nw_id_default, bounded::Bool=true, report::Bool=true)

    pd = _PM.var(pm, nw)[:pd] = Dict(i => JuMP.@variable(pm.model,
            base_name="$(nw)_pd_$(i)"
            #,start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "pd_start", c, 0.0) #findall(idx -> idx == c, connections[i])[1]
        ) for i in _PM.ids(pm, nw, :load)
    )

    if bounded
        for (i,load) in _PM.ref(pm, nw, :load)
            if haskey(load, "pmin")
                for (idx, c) in enumerate(connections[i])
                    JuMP.set_lower_bound(pd[i][c], load["pmin"][idx])
                end
            end
            if haskey(load, "pmax")
                for (idx, c) in enumerate(connections[i])
                    JuMP.set_upper_bound(pd[i][c], load["pmax"][idx])
                end
            end
        end
    end
    report && _PM.sol_component_value(pm, nw, :load, :pd, _PM.ids(pm, nw, :load), pd)
end

function variable_load_reactive(pm::_PM.AbstractPowerModel;
                                   nw::Int=_PM.nw_id_default, bounded::Bool=true, report::Bool=true)

    qd = _PM.var(pm, nw)[:qd] = Dict(i => JuMP.@variable(pm.model,
            base_name="$(nw)_qd_$(i)"
            #,start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "qd_start", c, 0.0) #findall(idx -> idx == c, connections[i])[1]
        ) for i in _PM.ids(pm, nw, :load)
    )

    if bounded
        for (i,load) in _PM.ref(pm, nw, :load)
            if haskey(load, "qmin")
                for (idx, c) in enumerate(connections[i])
                    JuMP.set_lower_bound(qd[i][c], load["qmin"][idx])
                end
            end
            if haskey(load, "qmax")
                for (idx, c) in enumerate(connections[i])
                    JuMP.set_upper_bound(qd[i][c], load["qmax"][idx])
                end
            end
        end
    end
    report && _PM.sol_component_value(pm, nw, :load, :qd, _PM.ids(pm, nw, :load), qd)
end
"""
    variable_load_current, IVR current equivalent of variable_load
"""
function variable_load_current(pm::_PM.AbstractIVRModel; kwargs...)
    # NB: currently, only (single-phase eq.) loads on the ac side are supported
    variable_load_current_real(pm; kwargs...)
    variable_load_current_imag(pm; kwargs...)
end

function variable_load_current_real(pm::_PM.AbstractIVRModel;
                                 nw::Int=_PM.nw_id_default, bounded::Bool=true, report::Bool=true)

    crd = _PM.var(pm, nw)[:crd] = Dict(i => JuMP.@variable(pm.model,
            base_name="$(nw)_crd_$(i)"
            #,start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "crd_start", c, 0.0)
        ) for i in _PM.ids(pm, nw, :load)
    )

    report && _PM.sol_component_value(pm, nw, :load, :crd, _PM.ids(pm, nw, :load), crd)
end

function variable_load_current_imag(pm::_PM.AbstractIVRModel; nw::Int=_PM.nw_id_default, bounded::Bool=true, report::Bool=true, meas_start::Bool=false)

    cid = _PM.var(pm, nw)[:cid] = Dict(i => JuMP.@variable(pm.model,
           base_name="$(nw)_cid_$(i)"
            #,start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "cid_start",c, 0.0)
        ) for i in _PM.ids(pm, nw, :load)
    )

    report && _PM.sol_component_value(pm, nw, :load, :cid, _PM.ids(pm, nw, :load), cid)
end
