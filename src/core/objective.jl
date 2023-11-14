function objective_minimize_residuals(pm::_PM.AbstractPowerModel)
    return JuMP.@objective(pm.model, Min,
    sum(
        sum(
            sum(_PM.var(pm, nw, :res, i)[idx] for idx in 1:length(_PM.var(pm, nw, :res, i)) )
        for i in _PM.ids(pm, nw, :meas))
    for (nw, nw_ref) in _PM.nws(pm) )
    )
end