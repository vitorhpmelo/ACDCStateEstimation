"""
    get_cmp_id(pm, nw, i)
Retrieves the id of the component in measurement i. This is a tuple if the component is a branch. Otherwise, it is a scalar.
"""
function get_cmp_id(pm::_PM.AbstractPowerModel, nw::Int, i::Int)
    if  _PM.ref(pm, nw, :meas, i, "cmp") == :branch
        branch_id = _PM.ref(pm, nw, :meas, i, "cmp_id")
        cmp_id = (branch_id, _PM.ref(pm,nw,:branch, branch_id)["f_bus"], _PM.ref(pm,nw,:branch,branch_id)["t_bus"])
    else
        cmp_id = _PM.ref(pm, nw, :meas, i, "cmp_id")
    end
    return cmp_id
end
"""
    function get_active_connections(pm::_PM.AbstractPowerModel, nw::Int, cmp_type::Symbol, cmp_id::Int)
Returns the list of terminals, connections or t_ and f_connections, depending on the type of the component.
"""
function get_active_connections(pm::_PM.AbstractPowerModel, nw::Int, cmp_type::Symbol, cmp_id::Int)
    if cmp_type == :busdc
       active_conn = [1, 2, 3]
   elseif cmp_type ∈ [:gen, :load]
       active_conn = _PM.ref(pm, nw, cmp_type, cmp_id)["connections"]
   elseif cmp_type == :branch
       active_conn = intersect(_PM.ref(pm, nw, :branch, cmp_id)["f_connections"], _PM.ref(pm, nw, :branch, cmp_id)["t_connections"])
   else
       error("Measurements for component of type $cmp_type are not supported")
   end
   return active_conn
end
"""
    assign_unique_individual_criterion!(data::Dict)
    - data: data model of the network (dictionary)
    Assigns the criterion in data["se_settings"]["criterion"] to all individual measurements.
"""
function assign_unique_individual_criterion!(data::Dict)
    for (_, meas) in data["meas"]
        meas["crit"] = data["se_settings"]["criterion"]
    end
end