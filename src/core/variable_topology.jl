function variable_dc_branch_state(pm::_PM.AbstractPowerModel;
                                 nw::Int=_PM.nw_id_default, relax::Bool=false, report::Bool=true)

    candidate_branches = isempty(_PM.ref(pm, nw, :dc_branches_to_consider)) ? _PM.ids(pm, nw, :load) : _PM.ref(pm, nw, :dc_branches_to_consider)

    s = _PM.var(pm, nw)[:s] = Dict(i => JuMP.@variable(pm.model,
            base_name="$(nw)_s_$(i)",
            start = _PM.comp_start_value(_PM.ref(pm, nw, :branch, i), "state", 1),
            upper_bound = 1,
            lower_bound = 0,
            binary = !relax
        ) for i in candidate_branches
    )

    report && _PM.sol_component_value(pm, nw, :branch, :s, candidate_branches, s)
end

function variable_ac_branch_state(pm::_PM.AbstractPowerModel;
                                 nw::Int=_PM.nw_id_default, relax::Bool=false, report::Bool=true)

    candidate_branches = isempty(_PM.ref(pm, nw, :ac_branches_to_consider)) ? _PM.ids(pm, nw, :branchdc) : _PM.ref(pm, nw, :ac_branches_to_consider)

    s = _PM.var(pm, nw)[:s] = Dict(i => JuMP.@variable(pm.model,
            base_name="$(nw)_s_$(i)",
            start = _PM.comp_start_value(_PM.ref(pm, nw, :branchdc, i), "state", 1),
            upper_bound = 1,
            lower_bound = 0,
            binary = !relax
        ) for i in candidate_branches
    )

    report && _PM.sol_component_value(pm, nw, :branchdc, :s, candidate_branches, s)
end