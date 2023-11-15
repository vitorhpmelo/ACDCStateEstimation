"""
    constraint_residual

Equality constraint that describes the residual definition, which depends on the
criterion assigned to each individual measurement in data["meas"]["m"]["crit"].
"""
function constraint_residual(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)

    cmp_id = get_cmp_id(pm, nw, i)
    res = _PM.var(pm, nw, :res, i)
    var = _PM.var(pm, nw, _PM.ref(pm, nw, :meas, i, "var"), cmp_id)
    dst = _PM.ref(pm, nw, :meas, i, "dst")
    rsc = _PM.ref(pm, nw, :se_settings)["rescaler"]
    crit = _PM.ref(pm, nw, :meas, i, "crit")
    
    conns = get_active_connections(pm, nw, _PM.ref(pm, nw, :meas, i, "cmp"), cmp_id)

    for (idx, c) in enumerate(conns)
        if (occursin("ls", crit) || occursin("lav", crit)) && isa(dst[idx], _DST.Normal)
            μ, σ = occursin("w", crit) ? (_DST.mean(dst[idx]), _DST.std(dst[idx])) : (_DST.mean(dst[idx]), 1.0)
        end
        if isa(dst[idx], Float64)
            JuMP.@constraint(pm.model, var[c] == dst[idx])
            JuMP.@constraint(pm.model, res[idx] == 0.0)         
        elseif crit ∈ ["wls", "ls"] && isa(dst[idx], _DST.Normal)
            JuMP.@constraint(pm.model,
                res[idx] * rsc^2 * σ^2 == (var[c] - μ)^2 
            )
        elseif crit == "rwls" && isa(dst[idx], _DST.Normal)
            JuMP.@constraint(pm.model,
                res[idx] * rsc^2 * σ^2 >= (var[c] - μ)^2
            )
        elseif crit ∈ ["wlav", "lav"] && isa(dst[idx], _DST.Normal)
            JuMP.@NLconstraint(pm.model,
                res[idx] * rsc * σ == abs(var[c] - μ)
            )
        elseif crit == "rwlav" && isa(dst[idx], _DST.Normal)
            JuMP.@constraint(pm.model,
                res[idx] * rsc * σ >= (var[c] - μ) 
            )
            JuMP.@constraint(pm.model,
                res[idx] * rsc * σ >= - (var[c] - μ)
            )
        else
            error("SE criterion of measurement $(i) not recognized")
        end
    end
end

# adaptation for se of _PMMCDC.constraint_kcl_shunt
# the difference here is that load demand is a variable, not a parameter read from the dictionary
function constraint_kcl_shunt_se(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    bus = _PM.ref(pm, nw, :bus, i)
    bus_arcs = _PM.ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = _PM.ref(pm, nw, :bus_arcs_dc, i)
    bus_gens = _PM.ref(pm, nw, :bus_gens, i)
    bus_convs_ac = _PM.ref(pm, nw, :bus_convs_ac, i)
    bus_loads = _PM.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PM.ref(pm, nw, :bus_shunts, i)

    gs = Dict(k => _PM.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bs = Dict(k => _PM.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    vm = _PM.var(pm, nw, :vm, i)
    p = _PM.var(pm, nw, :p)
    q = _PM.var(pm, nw, :q)
    pg = _PM.var(pm, nw, :pg)
    qg = _PM.var(pm, nw, :qg)
    pd = _PM.var(pm, nw, :pd)
    qd = _PM.var(pm, nw, :qd)
    pconv_grid_ac = _PM.var(pm, nw, :pconv_tf_fr)
    qconv_grid_ac = _PM.var(pm, nw, :qconv_tf_fr)

    JuMP.@NLconstraint(pm.model, sum(p[a] for a in bus_arcs) + sum(sum(pconv_grid_ac[c][d] for d in 1:length(_PM.var(pm, nw, :pconv_tf_fr, c))) for c in bus_convs_ac) == sum(pg[g] for g in bus_gens) - sum(pd[d] for d in bus_loads) - sum(gs[s] for s in bus_shunts) * vm^2)
    JuMP.@NLconstraint(pm.model, sum(q[a] for a in bus_arcs) + sum(sum(qconv_grid_ac[c][d] for d in 1:length(_PM.var(pm, nw, :qconv_tf_fr, c))) for c in bus_convs_ac) == sum(qg[g] for g in bus_gens) - sum(qd[d] for d in bus_loads) + sum(bs[s] for s in bus_shunts) * vm^2)
end