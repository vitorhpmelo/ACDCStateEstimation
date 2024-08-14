"""
    get_cmp_id(pm, nw, i)
Retrieves the id of the component in measurement i. This is a tuple if the component is a branch. Otherwise, it is a scalar.
"""
function get_cmp_id(pm::_PM.AbstractPowerModel, nw::Int, i::Int)
    # TODO add clause for converters, if we end up considering them
    cmp_id = _PM.ref(pm, nw, :meas, i, "cmp_id")
    if _PM.ref(pm, nw, :meas, i, "cmp") == :branchdc 
        cmp_id = cmp_id[1]
        bus1, bus2 = _PM.ref(pm, nw, :meas, i, "direction") == :from ? (_PM.ref(pm,nw,:branchdc, cmp_id)["fbusdc"], _PM.ref(pm,nw,:branchdc, cmp_id)["tbusdc"]) : (_PM.ref(pm,nw,:branchdc, cmp_id)["tbusdc"], _PM.ref(pm,nw,:branchdc, cmp_id)["fbusdc"])
        cmp_id = (cmp_id, bus1, bus2)
    elseif _PM.ref(pm, nw, :meas, i, "cmp") == :branch
        cmp_id = cmp_id[1] 
        bus1, bus2 = _PM.ref(pm, nw, :meas, i, "direction") == :from ? (_PM.ref(pm,nw,:branch, cmp_id)["f_bus"], _PM.ref(pm,nw,:branch, cmp_id)["t_bus"]) : (_PM.ref(pm,nw,:branch, cmp_id)["t_bus"], _PM.ref(pm,nw,:branch, cmp_id)["f_bus"])
        cmp_id = (cmp_id, bus1, bus2)
    end
    return cmp_id
end
"""
    function get_active_connections(pm::_PM.AbstractPowerModel, nw::Int, cmp_type::Symbol, cmp_id)
Returns the list of terminals, connections or t_ and f_connections, depending on the type of the component.
"""
function get_active_connections(pm::_PM.AbstractPowerModel, nw::Int, cmp_type::Symbol, cmp_id)
    if cmp_type == :busdc
       return [1, 2, 3]
    elseif cmp_type == :branchdc
       return _PM.ref(pm, nw, cmp_type, cmp_id[1])["line_confi"] == 2 ? [1,2,3] : [1]
#    elseif cmp_type == :branchdc
#         active_conn =  _PM.ref(pm, nw, cmp_type, cmp_id)["line_confi"] == 2 ? [1,2,3] : get_monopolar_conn(_PM.ref(pm, nw, cmp_type, cmp_id)["connect_at"])
#     elseif cmp_type == :convdc
#         active_conn =  _PM.ref(pm, nw, cmp_type, cmp_id)["conv_confi"] == 2 ? [1,2,3] : get_monopolar_conn(_PM.ref(pm, nw, cmp_type, cmp_id)["connect_at"])
   else
       return [1]
   end
end
"""
    get_monopolar_conn(connect_at::Int)
"""
function get_monopolar_conn(connect_at::Int)::Vector{Int}
    if connect_at == 0
        return [1,2]
    elseif connect_at == 1
        return [1,3]
    elseif connect_at == 2
        return [2,3]
    end
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
