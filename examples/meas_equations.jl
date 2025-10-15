#powers 

"""
Creates a branch power flow function, with the state variables mapped in the x vector
"""
function create_branch_p(vm_ruler::Dict{String, Int}, va_ruler::Dict{String, Int}, direction::Symbol, branch_i::Dict{String, Any})

    if direction == :from
        i = branch_i["f_bus"]
        j = branch_i["t_bus"]
        gii = branch_i["gii"]
        gij = branch_i["gij"]
        bij = branch_i["bij"]
        bii = branch_i["bii"]
    else
        i = branch_i["t_bus"]
        j = branch_i["f_bus"]
        gii = branch_i["gjj"]
        gij = branch_i["gji"]
        bij = branch_i["bji"]
        bii = branch_i["bjj"]
    end

    vmi = vm_ruler[i]
    vmj = vm_ruler[j]
    vai = va_ruler[i]
    vaj = va_ruler[j]

    return function(x)
        return (x[vmi]^2)*gii + x[vmi]* x[vmj]*(gij*cos(x[vai]-x[vaj])+bij*sin(x[vai]-x[vaj]))
    end
end





"""
Creates a branch reactive power flow function, with the state variables mapped in the x vector
"""
function create_branch_q(vm_ruler::Dict{String, Int}, va_ruler::Dict{String, Int}, direction::Symbol,  branch_i::Dict{String, Any})

    if direction == :from
        i = branch_i["f_bus"]
        j = branch_i["t_bus"]
        gii = branch_i["gii"]
        gij = branch_i["gij"]
        bij = branch_i["bij"]
        bii = branch_i["bii"]
    else
        i = branch_i["t_bus"]
        j = branch_i["f_bus"]
        gii = branch_i["gjj"]
        gij = branch_i["gji"]
        bij = branch_i["bji"]
        bii = branch_i["bjj"]
    end

    vmi = vm_ruler[i]
    vmj = vm_ruler[j]
    vai = va_ruler[i]
    vaj = va_ruler[j]
    
    return function(x)
        return -(x[vmi]^2)*bii - x[vmi]* x[vmj]*(bij*cos(x[vai]-x[vaj]) -gij*sin(x[vai]-x[vaj]))
    end
end

function create_bus_p_sh(vm_ruler::Dict{String, Int}, bus_i::Dict{String, Any})

    vmi = vm_ruler["$(bus_i["bus_i"])"]
    gs = bus_i["gs"]
    return function(x)
        return (x[vmi]^2)*gs 
    end
end



function create_bus_q_sh(vm_ruler::Dict{String, Int}, bus_i::Dict{String, Any})
    vmi = vm_ruler["$(bus_i["bus_i"])"]
    bs = bus_i["bs"]
    return function(x)
        return -(x[vmi]^2)*bs
    end
end

function create_conv_q_sh(vm_ruler::Dict{String, Int}, conv_i::Dict{String, Any},c::Int)
    bus_i=conv_i["tf"][c]["t_bus"]
    vmi = vm_ruler[bus_i]
    bs = conv_i["f"][c]["bf"]

    return function(x)
        return -(x[vmi]^2)*bs
    end
end


function create_bus_p(vm_ruler, va_ruler, bus_calc, branch_calc, conv_calc, bus)
    h_int = Function[]
    for branch in bus_calc["$bus"]["branch_fr"]
        push!(h_int, create_branch_p(vm_ruler, va_ruler, :from, branch_calc["$branch"]))
    end
    for branch in bus_calc["$bus"]["branch_to"]
        push!(h_int, create_branch_p(vm_ruler, va_ruler, :to, branch_calc["$branch"]))
    end
    push!(h_int, create_bus_p_sh(vm_ruler, bus_calc["$bus"]))
    for conv in bus_calc["$bus"]["convdc"]
        haskey(conv_calc["$conv"], "tf") ? br = "tf" : br = "pr" # Do not work for converters connected directly to the bus
        for c in 1:conv_calc["$conv"]["conn"]
            push!(h_int, create_branch_p(vm_ruler, va_ruler, :from, conv_calc["$conv"][br][c]))
        end
    end
    return sum_functions(h_int)
end



function create_bus_q(vm_ruler, va_ruler, bus_calc, branch_calc, conv_calc, bus)
    
    h_int = Function[]
    for branch in bus_calc["$bus"]["branch_fr"]
        push!(h_int, create_branch_q(vm_ruler, va_ruler, :from, branch_calc["$branch"]))
    end
    for branch in bus_calc["$bus"]["branch_to"]
        push!(h_int, create_branch_q(vm_ruler, va_ruler, :to, branch_calc["$branch"]))
    end
    push!(h_int, create_bus_q_sh(vm_ruler, bus_calc["$bus"]))
    for conv in bus_calc["$bus"]["convdc"]
        haskey(conv_calc["$conv"],"tf") ? br="tf" : br="pr" # Do not work for converters connected directly to the bus
        for c in 1:conv_calc["$conv"]["conn"]
        push!(h_int, create_branch_q(vm_ruler, va_ruler, :from, conv_calc["$conv"][br][c]))
        end
    end
    return sum_functions(h_int)
end

#currents 


"""
Creates a branch current (real) flow function, with the state variables mapped in the x vector
"""
function create_branch_cr(vm_ruler::Dict{String, Int}, va_ruler::Dict{String, Int}, direction::Symbol, branch_i::Dict{String, Any})

    if direction == :from
        i = branch_i["f_bus"]
        j = branch_i["t_bus"]
        gii = branch_i["gii"]
        gij = branch_i["gij"]
        bij = branch_i["bij"]
        bii = branch_i["bii"]
    else
        i = branch_i["t_bus"]
        j = branch_i["f_bus"]
        gii = branch_i["gjj"]
        gij = branch_i["gji"]
        bij = branch_i["bji"]
        bii = branch_i["bjj"]
    end

    vmi = vm_ruler[i]
    vmj = vm_ruler[j]
    vai = va_ruler[i]
    vaj = va_ruler[j]

    return function(x)
        return -x[vmi]*bii*sin(x[vai])+x[vmi]*gii*cos(x[vai])-x[vmj]*bij*sin(x[vaj])+x[vmj]*gij*cos(x[vaj])
    end
end

"""
Creates a branch branch current (imaginary) flow function, with the state variables mapped in the x vector
"""
function create_branch_ci(vm_ruler::Dict{String, Int}, va_ruler::Dict{String, Int}, direction::Symbol,  branch_i::Dict{String, Any})

    if direction == :from
        i = branch_i["f_bus"]
        j = branch_i["t_bus"]
        gii = branch_i["gii"]
        gij = branch_i["gij"]
        bij = branch_i["bij"]
        bii = branch_i["bii"]
    else
        i = branch_i["t_bus"]
        j = branch_i["f_bus"]
        gii = branch_i["gjj"]
        gij = branch_i["gji"]
        bij = branch_i["bji"]
        bii = branch_i["bjj"]
    end

    vmi = vm_ruler[i]
    vmj = vm_ruler[j]
    vai = va_ruler[i]
    vaj = va_ruler[j]
    
    return function(x)
        return x[vmi]*bii*cos(x[vai]) + x[vmi]*gii*sin(x[vai]) + x[vmj]*bij*cos(x[vaj]) + x[vmj]*gij*sin(x[vaj])
    end
end





"""
Creates a branch branch current (module) flow function, with the state variables mapped in the x vector
"""

function create_branch_cm(vm_ruler::Dict{String, Int}, va_ruler::Dict{String, Int}, direction::Symbol, branch::Dict{String, Any})

    if direction == :from
        i = branch["f_bus"]

    else
        i = branch["t_bus"]
    end

    p = create_branch_p(vm_ruler, va_ruler, direction, branch)
    q = create_branch_q(vm_ruler, va_ruler, direction, branch)
    v = create_bus_vm(vm_ruler[i])

    return function(x)
        return sqrt(p(x)^2 + q(x)^2)/v(x)
    end
end


"""Creates a branch current (real) shunt injection function, with the state variables mapped in the x vector
"""
function create_bus_cr_sh(vm_ruler::Dict{String, Int}, bus_i::Dict{String, Any})

    vmi = vm_ruler["$(bus_i["bus_i"])"]
    vai = va_ruler[i]
    gs = bus_i["gs"]
    bs = bus_i["bs"]
    return function(x)
        return x[vmi]*(-bs*sin(x[vai])+gs*cos(x[vai]))
    end
end

"""Creates a branch current (imaginary) shunt injection function, with the state variables mapped in the x vector
"""
function create_bus_ci_sh(vm_ruler::Dict{String, Int}, bus_i::Dict{String, Any})
    vmi = vm_ruler["$(bus_i["bus_i"])"]
    vai = va_ruler[i]
    gs = bus_i["gs"]
    bs = bus_i["bs"]
    return function(x)
        return x[vmi]*(-bs*sin(x[vai])+gs*cos(x[vai]))
    end
end


"""Creates a branch current (imaginary) injection function, with the state variables mapped in the x vector
"""
function create_bus_cr(vm_ruler, va_ruler, bus_calc, branch_calc, conv_calc, bus)
    h_int = Function[]
    for branch in bus_calc["$bus"]["branch_fr"]
        push!(h_int, create_branch_cr(vm_ruler, va_ruler, :from, branch_calc["$branch"]))
    end
    for branch in bus_calc["$bus"]["branch_to"]
        push!(h_int, create_branch_cr(vm_ruler, va_ruler, :to, branch_calc["$branch"]))
    end
    push!(h_int, create_bus_cr_sh(vm_ruler, bus_calc["$bus"]))
    for conv in bus_calc["$bus"]["convdc"]
        haskey(conv_calc["$conv"], "tf") ? br = "tf" : br = "pr" # Do not work for converters connected directly to the bus
        for c in 1:conv_calc["$conv"]["conn"]
            push!(h_int, create_branch_cr(vm_ruler, va_ruler, :from, conv_calc["$conv"][br][c]))
        end
    end
    return sum_functions(h_int)
end


"""
Creates a branch current (imaginary) injection function, with the state variables mapped in the x vector
"""
function create_bus_ci(vm_ruler, va_ruler, bus_calc, branch_calc, conv_calc, bus)
    
    h_int = Function[]
    for branch in bus_calc["$bus"]["branch_fr"]
        push!(h_int, create_branch_ci(vm_ruler, va_ruler, :from, branch_calc["$branch"]))
    end
    for branch in bus_calc["$bus"]["branch_to"]
        push!(h_int, create_branch_ci(vm_ruler, va_ruler, :to, branch_calc["$branch"]))
    end
    push!(h_int, create_bus_ci_sh(vm_ruler, bus_calc["$bus"]))
    for conv in bus_calc["$bus"]["convdc"]
        haskey(conv_calc["$conv"],"tf") ? br="tf" : br="pr" # Do not work for converters connected directly to the bus
        for c in 1:conv_calc["$conv"]["conn"]
        push!(h_int, create_branch_ci(vm_ruler, va_ruler, :from, conv_calc["$conv"][br][c]))
        end
    end
    return sum_functions(h_int)
end




function create_conv_ploss(vm_ruler::Dict{String, Int}, va_ruler::Dict{String, Int}, conv_i::Dict{String, Any},c::Int)


    haskey(conv_i,"pr") ? br="pr" : br="tf"
    Iabs= create_branch_cm(vm_ruler, va_ruler, :to, conv_i[br][c])


    return function(x)
        return conv_i["LossA"][c] +  conv_i["LossB"][c] * Iabs(x) +  conv_i["LossC"][c] * Iabs(x)^2
    end
end





function create_bus_vm(vmi::Int)
    return function(x)
        return x[vmi]
    end
end



function create_bus_va(vai::Int)
    return function(x)
        return x[vai]
    end
end


function create_bus_vm(vm_ruler::Dict{String, Int}, bus_i::Dict{String, Any})
    return function(x)
        return x[vm_ruler["$(bus_i["bus_i"])"]]
    end
end

function create_bus_va(vm_ruler::Dict{String, Int}, bus_i::Dict{String, Any})
    return function(x)
        return x[vm_ruler["$(bus_i["bus_i"])"]]
    end
end



function create_branch_idc(vdc_ruler::Dict{String, Int}, direction::Symbol, branchdc_i::Dict{String, Any})


    (i,j) = direction == :from ? (branchdc_i["f_busdc"], branchdc_i["t_busdc"]) : (branchdc_i["t_busdc"], branchdc_i["f_busdc"])

    vdcmi=vdc_ruler[i]
    vdcmj=vdc_ruler[j]

    return function(x)
        return (x[vdcmi]- x[vdcmj])/branchdc_i["r"]
    end
end

function create_branch_pdc(vdc_ruler::Dict{String, Int}, direction::Symbol, branchdc_i::Dict{String, Any})

    
    (i,j) = direction == :from ? (branchdc_i["f_busdc"], branchdc_i["t_busdc"]) : (branchdc_i["t_busdc"], branchdc_i["f_busdc"])

    vdcmi=vdc_ruler[i]
    vdcmj=vdc_ruler[j]

    return function(x)
        return x[vdcmi]*(x[vdcmi]- x[vdcmj])/branchdc_i["r"]
    end
end

function create_conv_m(vm_ruler::Dict{String, Int}, vdc_ruler::Dict{String, Int}, conv_i::Dict{String, Any}, c::Int)

    vmac = vm_ruler["$(conv_i["conv_i"])-c-$c"]
    vmdcp = vdc_ruler["$(conv_i["busdc_i"])-dc-$(conv_i["vdcp"][c])"]
    vmdcn = vdc_ruler["$(conv_i["busdc_i"])-dc-$(conv_i["vdcn"][c])"]

    return function(x)
        return x[vmac]/(x[vmdcp]-x[vmdcn])
    end
end


function create_bus_pdc(vmdc_ruler::Dict{String, Int}, busdc_calc::Dict{String, Any}, branchdc_calc::Dict{String, Any},busdc_i::Int, c::Int)
    h_int = Function[]
    for branch in busdc_calc["$busdc_i"]["branchdc_fr"]
        push!(h_int, create_branch_pdc(vmdc_ruler, :from, branchdc_calc["$branch"][c]))
    end
    for branch in busdc_calc["$busdc_i"]["branchdc_to"]
        push!(h_int, create_branch_pdc(vmdc_ruler, :to, branchdc_calc["$branch"][c]))
    end
    
    return sum_functions(h_int)
end




function create_conv_idc(vmdc_ruler::Dict{String, Int}, busdc_calc::Dict{String, Any}, branchdc_calc::Dict{String, Any}, conv_calc_i::Dict{String, Any}, c::Int)
    
    h_int = Function[]
    vdcp=conv_calc_i["vdcp"][c]
    vdcn=conv_calc_i["vdcn"][c]
    busdc_i = conv_calc_i["busdc_i"]
    pol = vdcp == 3 ? vdcn : vdcp

    for branch in busdc_calc["$busdc_i"]["branchdc_fr"]
        for cond in branchdc_calc["$branch"]["conductor"]
            if cond["f_busdc"] == "$busdc_i-dc-$pol"
                push!(h_int, create_branch_idc(vmdc_ruler, :from, cond))
            end        
        end
    end 
    
    for branch in busdc_calc["$busdc_i"]["branchdc_to"]
        for cond in branchdc_calc["$branch"]["conductor"]
            if cond["t_busdc"] == "$busdc_i-dc-$pol"
                push!(h_int, create_branch_idc(vmdc_ruler, :to, cond))
            end
        end
    end
    
    return sum_functions(h_int)
end



#TODO improve this function, it has to call the dictorionaries bus dc and branch dc for the whole network, but the conv for the specific converter, this is not consistent with the remaining implementations
function create_conv_pdc(vmdc_ruler::Dict{String, Int}, busdc_calc::Dict{String, Any}, branchdc_calc::Dict{String, Any}, conv_calc_i::Dict{String, Any}, c::Int)

    vdcp=conv_calc_i["vdcp"][c]
    vdcn=conv_calc_i["vdcn"][c]
    busdc_i = conv_calc_i["busdc_i"]

    key_vdcp = "$busdc_i-dc-$vdcp"
    key_vdcn = "$busdc_i-dc-$vdcn"
    (i,j) = vdcp == 3 ? (vmdc_ruler[key_vdcp], vmdc_ruler[key_vdcn]) : (vmdc_ruler[key_vdcn], vmdc_ruler[key_vdcp])
    v_drop(x) = x[i] - x[j]
    idc = create_conv_idc(vmdc_ruler, busdc_calc, branchdc_calc, conv_calc_i, c)

    return prod_functions([idc, v_drop])
    
end

function prod_functions(flist::Vector{Function})
    return function(x)
        p = 1.0
        for f in flist
            p *= f(x)
        end
        return p
    end
end

function sum_functions(flist::Vector{Function})
    return function(x)
        s = 0.0
        for f in flist
            s += f(x)
        end
        return s
    end
end

function create_c_conv_loss(vm_ruler::Dict{String, Int}, va_ruler::Dict{String, Int}, vdc_ruler::Dict{String, Int}, conv_calc_i::Dict{String, Any}, busdc_calc::Dict{String, Any}, branchdc_calc::Dict{String, Any}, c::Int)

    p_loss=create_conv_ploss(vm_ruler, va_ruler, conv_calc_i ,c)
    haskey(conv_calc_i,"pr") ? br="pr" : br="tf"
    pac_conv=create_branch_p(vm_ruler, va_ruler, :to, conv_calc_i[br][c])
    pdc_conv=create_conv_pdc(vdc_ruler, busdc_calc, branchdc_calc, conv_calc_i, c)
    f(x) = p_loss(x) + pac_conv(x) - pdc_conv(x)
    return f
    
end


function create_c_conv_Pinj(vm_ruler::Dict{String, Int}, va_ruler::Dict{String, Int}, conv_calc_i::Dict{String, Any}, c::Int)

    p_pr =create_branch_p(vm_ruler, va_ruler, :from, conv_calc_i["pr"][c])
    p_tf =create_branch_p(vm_ruler, va_ruler, :to, conv_calc_i["tf"][c])
    f(x) = p_pr(x) + p_tf(x)
    return f
    
end

function create_c_conv_Qinj(vm_ruler::Dict{String, Int}, va_ruler::Dict{String, Int}, conv_calc_i::Dict{String, Any}, c::Int)

    q_pr =create_branch_q(vm_ruler, va_ruler, :from, conv_calc_i["pr"][c])
    q_tf =create_branch_q(vm_ruler, va_ruler, :to, conv_calc_i["tf"][c])
    q_sh = create_conv_q_sh(vm_ruler, conv_calc_i, c)
    f(x) = q_pr(x) + q_tf(x) + q_sh(x)  
    return f
    
end

function create_h(data::Dict{String, Any}, network_info)
    

    (vm_ruler, va_ruler, vdc_ruler) = (network_info["vm_ruler"], network_info["va_ruler"], network_info["vdc_ruler"])

    (branch_calc, bus_calc, conv_calc, branchdc_calc, busdc_calc) = (network_info["branch_calc"], network_info["bus_calc"], network_info["conv_calc"], network_info["branchdc_calc"], network_info["busdc_calc"])

    h = Vector{Function}()
    bus_P_injc_added=Vector{Int64}()
    bus_Q_injc_added=Vector{Int64}()

    for (key, meas) in data["meas"]
        if meas["var"] == :p
            cmp_id = meas["cmp_id"][1]
            push!(h, create_branch_p(vm_ruler, va_ruler, meas["direction"], branch_calc["$cmp_id"]))
        elseif  meas["var"] == :q #branch power flows
            cmp_id = meas["cmp_id"][1]
            push!(h, create_branch_q(vm_ruler, va_ruler, meas["direction"], branch_calc["$cmp_id"]))
        elseif meas["var"] == :vm
            cmp_id = meas["cmp_id"]
            push!(h, create_bus_vm(vm_ruler,bus_calc["$cmp_id"]))
        elseif meas["var"] == :va
            cmp_id = meas["cmp_id"]
            push!(h, create_bus_va(va_ruler,bus_calc["$cmp_id"]))
        elseif  meas["var"] in [:pd, :pg]  # add injection measurements
            cmp_id = meas["cmp_id"] # gets the load or gen id
            bus = get_bus_d_o_g(data, meas["var"], cmp_id) #gets the bus
            (bus in bus_P_injc_added) ? continue : nothing 
            #loads and gens meas need to be condensed into a single measurement
            push!(bus_P_injc_added, bus)
            push!(h, create_bus_p(vm_ruler, va_ruler, bus_calc, branch_calc, conv_calc, bus)) 
        elseif meas["var"] in [:qd, :qg] # reactive power injections
            cmp_id = meas["cmp_id"]
            bus = get_bus_d_o_g(data, meas["var"], cmp_id) #gets the bus
            (bus in bus_Q_injc_added) ? continue : nothing
            push!(bus_Q_injc_added, bus)
            push!(h, create_bus_q(vm_ruler, va_ruler, bus_calc, branch_calc, conv_calc, bus))
        elseif meas["var"] == :pconv_ac #inverted sing
            conv=meas["cmp_id"]
            haskey(conv_calc["$conv"],"pr") ? br="pr" : br="tf"
            for c in 1:conv_calc["$conv"]["conn"]
                push!(h, create_branch_p(vm_ruler, va_ruler, :to, conv_calc["$conv"][br][c])) # attention
            end
        elseif meas["var"] == :qconv_ac #inverted sing
            conv=meas["cmp_id"]
            haskey(conv_calc["$conv"],"pr") ? br="pr" : br="tf"
            for c in 1:conv_calc["$conv"]["conn"]
                push!(h, create_branch_q(vm_ruler, va_ruler, :to, conv_calc["$conv"][br][c])) # attention
            end
        elseif meas["var"] == :pconv_tf_to
            conv=meas["cmp_id"]
            haskey(conv_calc["$conv"],"tf") ? br="tf" : continue
            for c in 1:conv_calc["$conv"]["conn"]
                push!(h, create_branch_p(vm_ruler, va_ruler, :to, conv_calc["$conv"][br][c])) # attention
            end
        elseif meas["var"] == :qconv_tf_to
            conv=meas["cmp_id"]
            haskey(conv_calc["$conv"],"tf") ? br="tf" : continue
            for c in 1:conv_calc["$conv"]["conn"]
                push!(h, create_branch_q(vm_ruler, va_ruler, :to, conv_calc["$conv"][br][c])) # attention
            end
        elseif meas["var"] == :pconv_tf_fr
            conv=meas["cmp_id"]
            haskey(conv_calc["$conv"],"tf") ? br="tf" : continue
            for c in 1:conv_calc["$conv"]["conn"]
                push!(h, create_branch_p(vm_ruler, va_ruler, :from, conv_calc["$conv"][br][c])) # attention
            end
        elseif meas["var"] == :qconv_tf_fr
            conv=meas["cmp_id"]
            haskey(conv_calc["$conv"],"tf") ? br="tf" : continue
            for c in 1:conv_calc["$conv"]["conn"]
                push!(h, create_branch_q(vm_ruler, va_ruler, :from, conv_calc["$conv"][br][c])) # attention
            end
        elseif meas["var"] == :pconv_pr_fr
            conv=meas["cmp_id"]
            haskey(conv_calc["$conv"],"pr") ? br="pr" : continue
            for c in 1:conv_calc["$conv"]["conn"]
                push!(h, create_branch_p(vm_ruler, va_ruler, :from, conv_calc["$conv"][br][c])) # attention
            end
        elseif meas["var"] == :qconv_pr_fr
            conv=meas["cmp_id"]
            haskey(conv_calc["$conv"],"pr") ? br="pr" : continue
            for c in 1:conv_calc["$conv"]["conn"]
                push!(h, create_branch_q(vm_ruler, va_ruler, :from, conv_calc["$conv"][br][c])) # attention
            end
        elseif meas["var"] == :vmc
            conv=meas["cmp_id"]
            for c in 1:conv_calc["$conv"]["conn"]
                push!(h, create_bus_vm(vm_ruler["$conv-c-$c"]))
            end
        elseif meas["var"] == :vac
            conv=meas["cmp_id"]
            for c in 1:conv_calc["$conv"]["conn"]
                push!(h, create_bus_va(va_ruler["$conv-c-$c"]))
            end 
        elseif meas["var"] == :vmf
            conv=meas["cmp_id"]
            for c in 1:conv_calc["$conv"]["conn"]
                push!(h, create_bus_vm(vm_ruler["$conv-f-$c"]))
            end
        elseif meas["var"] == :vaf
            conv=meas["cmp_id"]
            for c in 1:conv_calc["$conv"]["conn"]
                push!(h, create_bus_va(va_ruler["$conv-f-$c"]))
            end
        elseif meas["var"] == :i_dcgrid
            branchdc = meas["cmp_id"][1]
            direction = meas["direction"]
            for c in 1:branchdc_calc["$branchdc"]["conn"]
                push!(h, create_branch_idc(vdc_ruler, direction, branchdc_calc["$branchdc"]["conductor"][c]))
            end
        elseif meas["var"] == :p_dcgrid
            branchdc = meas["cmp_id"][1]
            direction = meas["direction"]
            for c in 1:branchdc_calc["$branchdc"]["conn"]
                push!(h, create_branch_pdc(vdc_ruler, direction, branchdc_calc["$branchdc"]["conductor"][c]))
            end
        elseif meas["var"] == :vdcm
            busdc = meas["cmp_id"]
            for c in 1:busdc_calc["$busdc"]["conn"]
                push!(h, create_bus_vm(vdc_ruler["$busdc-dc-$c"]))
            end
        elseif meas["var"] == :mconv
            conv = meas["cmp_id"]
            for c in 1:conv_calc["$conv"]["conn"]
                push!(h, create_conv_m(vm_ruler, vdc_ruler, conv_calc["$(conv)"], c))
            end
        elseif meas["var"] == :cr
            cmp_id = meas["cmp_id"][1]
            direction = meas["direction"]
            push!(h, create_branch_cr(vm_ruler, va_ruler, direction, branch_calc["$(cmp_id)"]))
        elseif meas["var"] == :ci
            cmp_id = meas["cmp_id"][1]
            direction = meas["direction"]
            push!(h, create_branch_ci(vm_ruler, va_ruler, direction, branch_calc["$(cmp_id)"]))
        end
    end
    return h
end


function create_h(network_info::Dict{String, Any})
    

    (vm_ruler, va_ruler, vdc_ruler) = (network_info["vm_ruler"], network_info["va_ruler"], network_info["vdc_ruler"])

    (branch_calc, bus_calc, conv_calc, branchdc_calc, busdc_calc) = (network_info["branch_calc"], network_info["bus_calc"], network_info["conv_calc"], network_info["branchdc_calc"], network_info["busdc_calc"])

    h = Vector{Function}()
    

    for meas in network_info["z"]
        if meas.var == :p
            push!(h, create_branch_p(vm_ruler, va_ruler, meas.direction, branch_calc["$(meas.cmp_id)"]));
        elseif  meas.var == :q #branch power flows 
            push!(h, create_branch_q(vm_ruler, va_ruler, meas.direction, branch_calc["$(meas.cmp_id)"]));
        elseif meas.var == :vm
            push!(h, create_bus_vm(vm_ruler,bus_calc["$(meas.cmp_id)"]));
        elseif meas.var == :va
            push!(h, create_bus_va(va_ruler,bus_calc["$(meas.cmp_id)"]));
        elseif  meas.var == :pinj  # add injection measurements
            push!(h, create_bus_p(vm_ruler, va_ruler, bus_calc, branch_calc, conv_calc, meas.cmp_id)); #loads and gens meas need to be condensed into a single measurement
        elseif meas.var == :qinj # reactive power injections 
            push!(h, create_bus_q(vm_ruler, va_ruler, bus_calc, branch_calc, conv_calc, meas.cmp_id)); #loads and gens meas need to be condensed into a single measurement
        elseif meas.var == :pconv_ac #inverted sing
            push!(h, create_branch_p(vm_ruler, va_ruler, :to, conv_calc["$(meas.cmp_id)"]["pr"][meas.conn])); # attention
        elseif meas.var == :qconv_ac #inverted sing
            push!(h, create_branch_q(vm_ruler, va_ruler, :to, conv_calc["$(meas.cmp_id)"]["pr"][meas.conn])); # attention
        elseif meas.var == :pconv_tf_to
            push!(h, create_branch_p(vm_ruler, va_ruler, :to, conv_calc["$(meas.cmp_id)"]["tf"][meas.conn])); # attention
        elseif meas.var == :qconv_tf_to
            push!(h, create_branch_q(vm_ruler, va_ruler, :to, conv_calc["$(meas.cmp_id)"]["tf"][meas.conn])); # attention
        elseif meas.var == :pconv_tf_fr
            push!(h, create_branch_p(vm_ruler, va_ruler, :from, conv_calc["$(meas.cmp_id)"]["tf"][meas.conn])); # attention
        elseif meas.var == :qconv_tf_fr
            push!(h, create_branch_q(vm_ruler, va_ruler, :from, conv_calc["$(meas.cmp_id)"]["tf"][meas.conn])); # attention
        elseif meas.var == :pconv_pr_fr
            push!(h, create_branch_p(vm_ruler, va_ruler, :from, conv_calc["$(meas.cmp_id)"]["pr"][meas.conn])); # attention
        elseif meas.var == :qconv_pr_fr
            push!(h, create_branch_q(vm_ruler, va_ruler, :from, conv_calc["$(meas.cmp_id)"]["pr"][meas.conn])); # attention
        elseif meas.var == :vmc
            push!(h, create_bus_vm(vm_ruler["$(meas.cmp_id)-c-$(meas.conn)"]));
        elseif meas.var == :vac
            push!(h, create_bus_va(va_ruler["$(meas.cmp_id)-c-$(meas.conn)"]));
        elseif meas.var == :vmf
            push!(h, create_bus_vm(vm_ruler["$(meas.cmp_id)-f-$(meas.conn)"]));
        elseif meas.var == :vaf
            push!(h, create_bus_va(va_ruler["$(meas.cmp_id)-f-$(meas.conn)"]));
        elseif meas.var == :i_dcgrid
            push!(h, create_branch_idc(vdc_ruler, meas.direction, branchdc_calc["$(meas.cmp_id)"]["conductor"][meas.conn]));
        elseif meas.var == :p_dcgrid
            push!(h, create_branch_pdc(vdc_ruler, meas.direction, branchdc_calc["$(meas.cmp_id)"]["conductor"][meas.conn]));
        elseif meas.var == :vdcm
            push!(h, create_bus_vm(vdc_ruler["$(meas.cmp_id)-dc-$(meas.conn)"]));
        elseif meas.var == :mconv
            push!(h, create_conv_m(vm_ruler, vdc_ruler, conv_calc["$(meas.cmp_id)"], meas.conn));
        elseif meas.var == :cr
            push!(h, create_branch_cr(vm_ruler, va_ruler, meas.direction, branch_calc["$(meas.cmp_id)"]));
        elseif meas.var == :ci
            push!(h, create_branch_ci(vm_ruler, va_ruler, meas.direction, branch_calc["$(meas.cmp_id)"]));
        else
            @warn "Measurement type $(meas.var) not recognized, skipping."
        end
    end
    return h
end



function create_h_virtual(network_info::Dict{String, Any})
    

    (vm_ruler, va_ruler, vdc_ruler) = (network_info["vm_ruler"], network_info["va_ruler"], network_info["vdc_ruler"])

    (branch_calc, bus_calc, conv_calc, branchdc_calc, busdc_calc) = (network_info["branch_calc"], network_info["bus_calc"], network_info["conv_calc"], network_info["branchdc_calc"], network_info["busdc_calc"])

    h = Vector{Function}()
    

    for meas in network_info["z_virtual"]
        if  meas.var == :pinj  # add injection measurements
            push!(h, create_bus_p(vm_ruler, va_ruler, bus_calc, branch_calc, conv_calc, meas.cmp_id)); #loads and gens meas need to be condensed into a single measurement
        elseif meas.var == :qinj # reactive power injections 
            push!(h, create_bus_q(vm_ruler, va_ruler, bus_calc, branch_calc, conv_calc, meas.cmp_id)); #loads and gens meas need to be condensed into a single measurement
        else
            @warn "Virtual Measurement type $(meas.var) not recognized, skipping."
        end
    end
    return h
end



function create_g(network_info)
    (vm_ruler, va_ruler, vdc_ruler) = (network_info["vm_ruler"], network_info["va_ruler"], network_info["vdc_ruler"])
    (busdc_calc, branchdc_calc,conv_calc) = (network_info["busdc_calc"], network_info["branchdc_calc"], network_info["conv_calc"])
    g=network_info["g"]
    g_out=Vector{Function}()
    for g_i in g
        conv_i= g_i.cmp_id   
        push!(g_out, create_c_conv_loss(vm_ruler,va_ruler, vdc_ruler,conv_calc[conv_i],busdc_calc, branchdc_calc, g_i.conn))
    end
    return g_out
end


function create_g_Pinj(network_info)
    (vm_ruler, va_ruler, vdc_ruler) = (network_info["vm_ruler"], network_info["va_ruler"], network_info["vdc_ruler"])
    (busdc_calc, branchdc_calc,conv_calc) = (network_info["busdc_calc"], network_info["branchdc_calc"], network_info["conv_calc"])
    g=network_info["g_Pinj"]
    g_out=Vector{Function}()
    for g_i in g
        conv_i= g_i.cmp_id   
        push!(g_out, create_c_conv_Pinj(vm_ruler,va_ruler,conv_calc[conv_i], g_i.conn))
    end
    return g_out
end


function create_g_Qinj(network_info)
    (vm_ruler, va_ruler, vdc_ruler) = (network_info["vm_ruler"], network_info["va_ruler"], network_info["vdc_ruler"])
    (busdc_calc, branchdc_calc,conv_calc) = (network_info["busdc_calc"], network_info["branchdc_calc"], network_info["conv_calc"])
    g=network_info["g_Qinj"]
    g_out=Vector{Function}()
    for g_i in g
        conv_i= g_i.cmp_id
        push!(g_out, create_c_conv_Qinj(vm_ruler,va_ruler,conv_calc[conv_i], g_i.conn))
    end
    return g_out
end