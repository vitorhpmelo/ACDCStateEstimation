




function objective_minimize_residuals(pm::_PM.AbstractPowerModel)
    return JuMP.@objective(pm.model, Min,
    sum(
        sum(
            sum(_PM.var(pm, nw, :res, i)[idx] for idx in 1:length(_PM.var(pm, nw, :res, i)) )
        for i in _PM.ids(pm, nw, :meas))
    for (nw, nw_ref) in _PM.nws(pm) )
    )
end




# function objective_minimize_fase(pm::_PM.AbstractPowerModel)

    
#     meas_exp=sum(
#         sum(
#             sum(_PM.var(pm, nw, :res, i)[idx] for idx in 1:length(_PM.var(pm, nw, :res, i)) )
#         for i in _PM.ids(pm, nw, :meas))
#     for (nw, nw_ref) in _PM.nws(pm) )

#     pirori_exp=sum(
#         sum(
#             sum(_PM.var(pm, nw, :prior, i)[idx] for idx in 1:length(_PM.var(pm, nw, :prior, i)) )
#         for i in _PM.ids(pm, nw, :prior))
#     for (nw, nw_ref) in _PM.nws(pm) )
    
#     return  JuMP.@objective(pm.model, Min, meas_exp+pirori_exp)
# end

function objective_minimize_fase(pm::_PM.AbstractPowerModel)

    
    meas_exp=sum(
        sum(
            sum(_PM.var(pm, nw, :res, i)[idx] for idx in 1:length(_PM.var(pm, nw, :res, i)) )
        for i in _PM.ids(pm, nw, :meas))
    for (nw, nw_ref) in _PM.nws(pm) )

    pirori_exp=sum(
        sum(
            sum(_PM.var(pm, nw, :prior, i)[idx] for idx in 1:length(_PM.var(pm, nw, :prior, i)) )
        for i in _PM.ids(pm, nw, :prior))
    for (nw, nw_ref) in _PM.nws(pm) )
    

    return  JuMP.@objective(pm.model, Min, meas_exp+pirori_exp)
end






function calculate_objective(pm::_PM.AbstractPowerModel)

    
    JuMP.@expression(pm.model,meas_exp,
        sum(
            sum(
                sum(_PM.var(pm, nw, :res, i)[idx] for idx in 1:length(_PM.var(pm, nw, :res, i)) )
            for i in _PM.ids(pm, nw, :meas))
        for (nw, nw_ref) in _PM.nws(pm) )
    )

    JuMP.@expression(pm.model,pirori_exp,sum(
        sum(
            sum(_PM.var(pm, nw, :prior, i)[idx] for idx in 1:length(_PM.var(pm, nw, :prior, i)) )
        for i in _PM.ids(pm, nw, :prior))
    for (nw, nw_ref) in _PM.nws(pm) ))
        
    open("objective_values.txt", "w") do io
        println(io, "meas_exp: ", JuMP.value(meas_exp))
        println(io, "pirori_exp: ", JuMP.value(pirori_exp))
    end


end




function objective_minimize_residuals_alt(pm::_PM.AbstractPowerModel)
    
    residual_exp=JuMP.QuadExpr()
    for (nw, nw_ref) in _PM.nws(pm) 
        nw_id=nw
        for (i, _) in _PM.ref(pm, nw_id, :meas)
            cmp_id = _PM.ref(pm, nw_id, :meas, i, "cmp_id")
            var = _PM.var(pm, nw_id, _PM.ref(pm, nw_id, :meas, i, "var"), cmp_id)
            dst = _PM.ref(pm, nw_id, :meas, i, "dst")
            rsc = _PM.ref(pm, nw_id, :se_settings)["rescaler"]
            crit = _PM.ref(pm, nw_id, :meas, i, "crit")
            varmeas=_PM.ref(pm, nw_id, :meas, i,"var")

            conns = get_active_connections(pm, nw_id, _PM.ref(pm, nw_id, :meas, i, "cmp"), cmp_id, var=varmeas)
           
            for (idx, c) in enumerate(conns) #check via this
                μ,σ=(_DST.mean(dst[idx]), _DST.std(dst[idx])) 
                residual=JuMP.@expression(pm.model,(var[c] - μ)^2/(σ^2))
                JuMP.add_to_expression!(residual_exp,
                        residual
                    )
            end
        
        end
    end

    return  JuMP.@objective(pm.model, Min, residual_exp)
end

function objective_minimize_fase_alt(pm::_PM.AbstractPowerModel)

    
    meas_exp=sum(
        sum(
            sum(_PM.var(pm, nw, :res, i)[idx] for idx in 1:length(_PM.var(pm, nw, :res, i)) )
        for i in _PM.ids(pm, nw, :meas))
    for (nw, nw_ref) in _PM.nws(pm) )
    
    priori_exp=JuMP.QuadExpr()
    
    for (nw, nw_ref) in _PM.nws(pm) 
        nw_id=nw
        for (i,_) in _PM.ref(pm,nw_id,:prior)
            cmp_id_i = _PM.ref(pm, nw_id, :prior, i, "cmp_id_i") 
            cmp_id_j = _PM.ref(pm, nw_id, :prior, i, "cmp_id_j")
            conn_i= _PM.ref(pm, nw_id, :prior, i, "conn_i")
            conn_j= _PM.ref(pm, nw_id, :prior, i, "conn_j")

            var_i = _PM.var(pm, nw_id, _PM.ref(pm, nw_id, :prior, i, "var_i"), cmp_id_i)[conn_i] #gets the var_i connections might be a problem
            var_j = _PM.var(pm, nw_id, _PM.ref(pm, nw_id, :prior, i, "var_j"), cmp_id_j)[conn_j] #gets the var_j connections might be a problem

            π_i = _PM.ref(pm, nw_id, :prior, i, "π_i") # previous value of the variable i
            π_j = _PM.ref(pm, nw_id, :prior, i, "π_j") # previous value of the variable j
            a_ij= _PM.ref(pm, nw_id, :prior, i, "a_ij") 

            rsc = _PM.ref(pm, nw_id, :se_settings)["rescaler"]
            priori=JuMP.@expression(pm.model,(π_i-var_i)*(π_j-var_j)/(a_ij*rsc^2))
            JuMP.add_to_expression!(priori_exp, priori)
        end
    end


    return  JuMP.@objective(pm.model, Min, meas_exp+priori_exp)
end



function objective_minimize_fase_alt2(pm::_PM.AbstractPowerModel)
    
    residual_exp=JuMP.QuadExpr()
    for (nw, nw_ref) in _PM.nws(pm) 
        nw_id=nw
        for (i, _) in _PM.ref(pm, nw_id, :meas)
            cmp_id = _PM.ref(pm, nw_id, :meas, i, "cmp_id")
            var = _PM.var(pm, nw_id, _PM.ref(pm, nw_id, :meas, i, "var"), cmp_id)
            dst = _PM.ref(pm, nw_id, :meas, i, "dst")
            rsc = _PM.ref(pm, nw_id, :se_settings)["rescaler"]
            crit = _PM.ref(pm, nw_id, :meas, i, "crit")
            varmeas=_PM.ref(pm, nw_id, :meas, i,"var")

            conns = get_active_connections(pm, nw_id, _PM.ref(pm, nw_id, :meas, i, "cmp"), cmp_id, var=varmeas)
           
            for (idx, c) in enumerate(conns) #check via this
                μ,σ=(_DST.mean(dst[idx]), _DST.std(dst[idx]))
                residual=JuMP.@expression(pm.model,(var[c] - μ)^2/(σ^2 *rsc^2))
                JuMP.add_to_expression!(residual_exp,
                        residual)
            end
        
        end
    end

    
    priori_exp=JuMP.QuadExpr()
    
    for (nw, nw_ref) in _PM.nws(pm) 
        nw_id=nw
        for (i,_) in _PM.ref(pm,nw_id,:prior)
            cmp_id_i = _PM.ref(pm, nw_id, :prior, i, "cmp_id_i") 
            cmp_id_j = _PM.ref(pm, nw_id, :prior, i, "cmp_id_j")
            conn_i= _PM.ref(pm, nw_id, :prior, i, "conn_i")
            conn_j= _PM.ref(pm, nw_id, :prior, i, "conn_j")

            var_i = _PM.var(pm, nw_id, _PM.ref(pm, nw_id, :prior, i, "var_i"), cmp_id_i)[conn_i] #gets the var_i connections might be a problem
            var_j = _PM.var(pm, nw_id, _PM.ref(pm, nw_id, :prior, i, "var_j"), cmp_id_j)[conn_j] #gets the var_j connections might be a problem

            π_i = _PM.ref(pm, nw_id, :prior, i, "π_i") # previous value of the variable i
            π_j = _PM.ref(pm, nw_id, :prior, i, "π_j") # previous value of the variable j
            a_ij= _PM.ref(pm, nw_id, :prior, i, "a_ij") 

            rsc = _PM.ref(pm, nw_id, :se_settings)["rescaler"]
            priori=JuMP.@expression(pm.model,(π_i-var_i)*(π_j-var_j)/(a_ij*rsc^2))
            JuMP.add_to_expression!(priori_exp, priori)
        end
    end
    return  JuMP.@objective(pm.model, Min, residual_exp+priori_exp)
end
