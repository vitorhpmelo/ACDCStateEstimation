abstract type ConversionType end

struct SquareFraction<:ConversionType
    msr_id::Int64
    cmp_type::Symbol
    cmp_id::Int64
    bus_ind::Int64
    numerator::Array
    denominator::Array
end

struct Square<:ConversionType
    msr_id::Int64
    cmp_type::Symbol
    cmp_id::Int64
    bus_ind::Int64
    elements::Array
end

struct Multiplication<:ConversionType
    msr_sym::Symbol
    msr_id::Int64
    cmp_type::Symbol
    cmp_id::Int64
    bus_ind::Int64
    mult1::Array
    mult2::Array
end

struct Tangent<:ConversionType
    msr_id::Int64
    cmp_type::Symbol
    cmp_id::Int64
    numerator::Symbol
    denominator::Symbol
end

struct Fraction<:ConversionType
    msr_type::Symbol
    msr_id::Int64
    cmp_type::Symbol
    cmp_id::Int64
    bus_ind::Int64
    numerator::Array
    denominator::Symbol
end

struct MultiplicationFraction<:ConversionType
    msr_type::Symbol
    msr_id::Int64
    cmp_type::Symbol
    cmp_id::Int64
    bus_ind::Int64
    power::Array
    voltage::Array
end

function assign_conversion_type_to_msr(pm::_PM.AbstractACPModel,i,msr::Symbol, cmp_type;nw=nw)
    cmp_id = _PM.ref(pm, nw, :meas, i, "cmp_id")
    if cmp_type ∈ [:branchdc, :convdc, :busdc]
      msr_type = assign_conversion_type_to_msr_dc(pm, i, msr, cmp_id; nw=nw)
    else
      msr_type = assign_conversion_type_to_msr_ac(pm, i, msr, cmp_id; nw=nw)
    end
    return msr_type
end

function assign_conversion_type_to_msr_dc(pm::_PM.AbstractPowerModel, i, msr::Symbol, cmp_id; nw=nw)
  if msr == :p_dcgrid
      return Multiplication(msr, i,:branchdc, cmp_id, _PM.ref(pm,nw,:branchdc,cmp_id)["fbusdc"], [:i_dcgrid], [:vdcm])
  else
      error("the chosen measurement $(msr) at $(_PM.ref(pm, nw, :meas, i, "cmp")) $(_PM.ref(pm, nw, :meas, i, "cmp_id")) is not supported and should be removed")
  end
end

function assign_conversion_type_to_msr_ac(pm::_PM.AbstractACPModel, i, msr::Symbol, cmp_id; nw=nw)
  if msr == :cm
      msr_type = SquareFraction(i,:branch, cmp_id, _PM.ref(pm,nw,:branch,cmp_id)["f_bus"], [:p, :q], [:vm])
  elseif msr == :cmg
      msr_type = SquareFraction(i,:gen, cmp_id, _PM.ref(pm,nw,:gen,cmp_id)["gen_bus"], [:pg, :qg], [:vm])
  elseif msr == :cmd
      msr_type = SquareFraction(i,:load, cmp_id, _PM.ref(pm,nw,:load,cmp_id)["load_bus"], [:pd, :qd], [:vm])
  elseif msr == :cr
      msr_type = Fraction(msr, i,:branch, cmp_id, _PM.ref(pm,nw,:branch,cmp_id)["f_bus"], [:p, :q, :va], :vm)
  elseif msr == :crg
      msr_type = Fraction(msr, i,:gen, cmp_id, _PM.ref(pm,nw,:gen,cmp_id)["gen_bus"], [:pg, :qg, :va], :vm)
  elseif msr == :crd
      msr_type = Fraction(msr, i,:load, cmp_id, _PM.ref(pm,nw,:load,cmp_id)["load_bus"], [:pd, :qd, :va], :vm)
  elseif msr == :ci
      msr_type = Fraction(msr, i,:branch, cmp_id, _PM.ref(pm,nw,:branch,cmp_id)["f_bus"], [:p, :q, :va], :vm)
  elseif msr == :cig
      msr_type = Fraction(msr, i,:gen, cmp_id, _PM.ref(pm,nw,:gen,cmp_id)["gen_bus"], [:pg, :qg, :va], :vm)
  elseif msr == :cid
      msr_type = Fraction(msr, i,:load, cmp_id, _PM.ref(pm,nw,:load,cmp_id)["load_bus"], [:pd, :qd, :va], :vm)
  else
      error("the chosen measurement $(msr) at $(_PM.ref(pm, nw, :meas, i, "cmp")) $(_PM.ref(pm, nw, :meas, i, "cmp_id")) is not supported and should be removed")
  end
end

function assign_conversion_type_to_msr_ac(pm::_PM.AbstractACRModel,i,msr::Symbol, cmp_type;nw=nw)
    cmp_id = _PM.ref(pm, nw, :meas, i, "cmp_id")
    if msr == :vm
        msr_type = Square(i,:bus, cmp_id, _PM.ref(pm,nw,:bus,cmp_id)["index"], [:vi, :vr])
    elseif msr == :va
        msr_type = Tangent(i, :bus, cmp_id, :vi, :vr)
    elseif msr == :cm
        msr_type = SquareFraction(i,:branch, cmp_id, _PM.ref(pm,nw,:branch,cmp_id)["f_bus"], [:p, :q], [:vi, :vr])
    elseif msr == :cmg
        msr_type = SquareFraction(i,:gen, cmp_id, _PM.ref(pm,nw,:gen,cmp_id)["gen_bus"], [:pg, :qg], [:vi, :vr])
    elseif msr == :cmd
        msr_type = SquareFraction(i,:load, cmp_id, _PM.ref(pm,nw,:load,cmp_id)["load_bus"], [:pd, :qd], [:vi, :vr])
    elseif msr == :cr
        msr_type = MultiplicationFraction(msr, i,:branch, cmp_id, _PM.ref(pm,nw,:branch,cmp_id)["f_bus"], [:p, :q], [:vr, :vi])
    elseif msr == :crg
        msr_type = MultiplicationFraction(msr, i,:gen, cmp_id, _PM.ref(pm,nw,:gen,cmp_id)["gen_bus"], [:pg, :qg], [:vr, :vi])
    elseif msr == :crd
        msr_type = MultiplicationFraction(msr, i,:load, cmp_id, _PM.ref(pm,nw,:load,cmp_id)["load_bus"], [:pd, :qd], [:vr, :vi])
    elseif msr == :ci
        msr_type = MultiplicationFraction(msr, i,:branch, cmp_id, _PM.ref(pm,nw,:branch,cmp_id)["f_bus"], [:p, :q], [:vr, :vi])
    elseif msr == :cig
        msr_type = MultiplicationFraction(msr, i,:gen, cmp_id, _PM.ref(pm,nw,:gen,cmp_id)["gen_bus"], [:pg, :qg], [:vr, :vi])
    elseif msr == :cid
        msr_type = MultiplicationFraction(msr, i,:load, cmp_id, _PM.ref(pm,nw,:load,cmp_id)["load_bus"], [:pd, :qd], [:vr, :vi])
    else
       error("the chosen measurement $(msr) at $(_PM.ref(pm, nw, :meas, i, "cmp")) $(_PM.ref(pm, nw, :meas, i, "cmp_id")) is not supported and will be ignored")
    end
    return msr_type
end

function assign_conversion_type_to_msr_ac(pm::_PM.AbstractIVRModel,i,msr::Symbol, cmp_type;nw=nw)
    cmp_id = _PM.ref(pm, nw, :meas, i, "cmp_id")
    if msr == :vm
        msr_type = Square(i,:bus, cmp_id, _PM.ref(pm,nw,:bus,cmp_id)["index"], [:vi, :vr])
    elseif msr == :va
        msr_type = Tangent(i, :bus, cmp_id, :vi, :vr)
    elseif msr == :cm
        msr_type = Square(i,:branch, cmp_id, _PM.ref(pm,nw,:branch,cmp_id)["f_bus"], [:cr, :ci])
    elseif msr == :cmg
        msr_type = Square(i,:gen, cmp_id, _PM.ref(pm,nw,:gen,cmp_id)["gen_bus"], [:crg, :cig])
    elseif msr == :cmd
        msr_type = Square(i,:load, cmp_id, _PM.ref(pm,nw,:load,cmp_id)["load_bus"], [:crd, :cid])
    elseif msr == :ca
        msr_type = Tangent(i, :branch, cmp_id, :ci, :cr)
    elseif msr == :cag
        msr_type = Tangent(i, :gen, cmp_id, :cig, :crg)
    elseif msr == :cad
        msr_type = Tangent(i, :load, cmp_id, :cid, :crd)
    elseif msr == :p
        msr_type = Multiplication(msr, i,:branch, cmp_id, _PM.ref(pm,nw,:branch,cmp_id)["f_bus"], [:cr, :ci], [:vr, :vi])
    elseif msr == :pg
        msr_type = Multiplication(msr, i,:gen, cmp_id, _PM.ref(pm,nw,:gen,cmp_id)["gen_bus"], [:crg, :cig], [:vr, :vi])
    elseif msr == :pd
        msr_type = Multiplication(msr, i,:load, cmp_id, _PM.ref(pm,nw,:load,cmp_id)["load_bus"], [:crd, :cid], [:vr, :vi])
    elseif msr == :q
        msr_type = Multiplication(msr, i,:branch, cmp_id, _PM.ref(pm,nw,:branch,cmp_id)["f_bus"], [:cr, :ci], [:vr, :vi])
    elseif msr == :qg
        msr_type = Multiplication(msr, i,:gen, cmp_id, _PM.ref(pm,nw,:gen,cmp_id)["gen_bus"], [:crg, :cig], [:vr, :vi])
    elseif msr == :qd
        msr_type = Multiplication(msr, i,:load, cmp_id, _PM.ref(pm,nw,:load,cmp_id)["load_bus"], [:crd, :cid], [:vr, :vi])
    else
       error("the chosen measurement $(msr) at $(_PM.ref(pm, nw, :meas, i, "cmp")) $(_PM.ref(pm, nw, :meas, i, "cmp_id")) is not supported and should be removed")
    end
    return msr_type
end

function no_conversion_needed(pm::_PM.AbstractACPModel, msr_var::Symbol, cmp_type::Symbol)
  if cmp_type == :branchdc
    return msr_var ∈ [:i_dcgrid]
  elseif cmp_type == :convdc
    error("dc converter support not added yet!")
  elseif cmp_type == :busdc
    return msr_var ∈ [:vdcm]
  else
    return msr_var ∈ [:vm, :va, :pd, :qd, :pg, :qg, :p, :q]
  end
end
  
function no_conversion_needed(pm::_PM.AbstractACRModel, msr_var::Symbol, cmp_type::Symbol)
  return msr_var ∈ [:vr, :vi, :pd, :qd, :pg, :qg, :p, :q]
end

function no_conversion_needed(pm::_PM.AbstractIVRModel, msr_var::Symbol, cmp_type::Symbol)
  return msr_var ∈ [:vr, :vi, :cr, :ci, :crg, :cig, :crd, :cid]
end

function create_conversion_constraint(pm::_PM.AbstractPowerModel, original_var, msr::SquareFraction; nw=nw)

  new_var_num = []
  for nvn in msr.numerator
      if occursin("v", String(nvn)) && msr.cmp_type != :bus
          push!(new_var_num, _PM.var(pm, nw, nvn, msr.bus_ind))
      elseif msr.cmp_type == :branch
          push!(new_var_num, _PM.var(pm, nw, nvn, (msr.cmp_id, msr.bus_ind, _PM.ref(pm,nw,:branch,msr.cmp_id)["t_bus"])))
      else
          push!(new_var_num, _PM.var(pm, nw, nvn, msr.cmp_id))
      end
  end

  msr.cmp_type == :branch ? id = (msr.cmp_id,  msr.bus_ind, _PM.ref(pm,nw,:branch,msr.cmp_id)["t_bus"]) : id = msr.cmp_id
  conn = get_active_connections(pm, nw, msr.cmp_type, msr.cmp_id)

  new_var_den = []
  for nvd in msr.denominator
      if !isa(nvd, Symbol) #only case is when I have an array of ones
          push!(new_var_den, nvd)
      elseif occursin("v", String(nvd)) && msr.cmp_type != :bus
          push!(new_var_den, _PM.var(pm, nw, nvd, msr.bus_ind))
      else
          push!(new_var_den, _PM.var(pm, nw, nvd, msr.cmp_id))
      end
  end

  JuMP.@NLconstraint(pm.model, [c in conn],
      original_var[id][c]^2 == (sum( n[c]^2 for n in new_var_num ))/
                  (sum( d[c]^2 for d in new_var_den))
      )
end

function create_conversion_constraint(pm::_PM.AbstractPowerModel, original_var, msr::Square; nw=nw)

    new_var = []
    for nvn in msr.elements
        if msr.cmp_type == :branch
            push!(new_var, _PM.var(pm, nw, nvn, (msr.cmp_id, _PM.ref(pm,nw,:branch,msr.cmp_id)["f_bus"], _PM.ref(pm,nw,:branch,msr.cmp_id)["t_bus"])))
        else
            push!(new_var, _PM.var(pm, nw, nvn, msr.cmp_id))
        end
    end

    msr.cmp_type == :branch ? id = (msr.cmp_id,  _PM.ref(pm,nw,:branch,msr.cmp_id)["f_bus"], _PM.ref(pm,nw,:branch,msr.cmp_id)["t_bus"]) : id = msr.cmp_id
    conn = get_active_connections(pm, nw, msr.cmp_type, msr.cmp_id)

    JuMP.@constraint(pm.model, [c in conn],
        original_var[id][c]^2 == (sum( n[c]^2 for n in new_var ))
        )
end

function create_conversion_constraint(pm::_PM.AbstractPowerModel, original_var, msr::Multiplication; nw=nw)

    m1 = []
    m2 = []

    for m in msr.mult1
        if occursin("v", String(m)) && msr.cmp_type != :bus
            push!(m1, _PM.var(pm, nw, m, msr.bus_ind))
        elseif msr.cmp_type == :branch
            push!(m1, _PM.var(pm, nw, m, (msr.cmp_id, msr.bus_ind, _PM.ref(pm, nw, :branch,msr.cmp_id)["t_bus"])))
        elseif msr.cmp_type == :branchdc
          push!(m1, _PM.var(pm, nw, m, (msr.cmp_id, msr.bus_ind, _PM.ref(pm, nw, :branchdc,msr.cmp_id)["tbusdc"])))
        else
            push!(m1, _PM.var(pm, nw, m, msr.cmp_id))
        end
    end

    for mm in msr.mult2
        if occursin("v", String(mm)) && msr.cmp_type != :bus
            if msr.cmp_type == :branchdc && _PM.ref(pm, nw, :branchdc, msr.cmp_id)["line_confi"] == 1 # we need to align to the fact that all +-0 voltages exist but only two currents
                if _PM.ref(pm, nw, :branchdc, msr.cmp_id)["connect_at"] == 0
                    push!(m2, _PM.var(pm, nw, mm, msr.bus_ind)[[1,2]])    
                elseif _PM.ref(pm, nw, :branchdc, msr.cmp_id)["connect_at"] == 1
                    push!(m2, _PM.var(pm, nw, mm, msr.bus_ind)[[1,3]])
                else
                    push!(m2, _PM.var(pm, nw, mm, msr.bus_ind)[[2,3]])
                end
            else
                push!(m2, _PM.var(pm, nw, mm, msr.bus_ind))
            end
        elseif msr.cmp_type == :branch
            push!(m2, _PM.var(pm, nw, mm, (msr.cmp_id, msr.bus_ind, _PM.ref(pm, nw, :branch, msr.cmp_id)["t_bus"])))
        elseif msr.cmp_type == :branchdc
            push!(m2, _PM.var(pm, nw, m, (msr.cmp_id, msr.bus_ind, _PM.ref(pm, nw, :branchdc,msr.cmp_id)["tbusdc"])))
        else
            push!(m2, _PM.var(pm, nw, mm, msr.cmp_id))
        end
    end

    id = if msr.cmp_type == :branch 
      (msr.cmp_id,  msr.bus_ind, _PM.ref(pm,nw,:branch,msr.cmp_id)["t_bus"])
    elseif msr.cmp_type == :branchdc 
      (msr.cmp_id,  msr.bus_ind, _PM.ref(pm,nw,:branchdc,msr.cmp_id)["tbusdc"])
    else
      msr.cmp_id
    end
    
    conn = get_active_connections(pm, nw, msr.cmp_type, msr.cmp_id)
    
    if msr.cmp_type == :branchdc
      JuMP.@constraint(pm.model, [c in conn],
              original_var[id][c] == m1[1][c]*m2[1][c]
              )
    else
      if occursin("p", String(msr.msr_sym))
          JuMP.@constraint(pm.model, [c in conn],
              original_var[id][c] == m1[1][c]*m2[1][c]+m1[2][c]*m2[2][c]
              )
      elseif occursin("q", String(msr.msr_sym))
          JuMP.@constraint(pm.model, [c in conn],
              original_var[id][c] == -m1[2][c]*m2[1][c]+m1[1][c]*m2[2][c]
              )
      end
    end
end

function create_conversion_constraint(pm::_PM.AbstractPowerModel, original_var, msr::Tangent; nw=nw)
    #TODO for v0.2.0 this needs to be general to every distribution or we need to provide an exception
    @warn "Performing a Tangent conversion only makes sense for Normal distributions and is in general not advised"
    conn = get_active_connections(pm, nw, msr.cmp_type, msr.cmp_id)
    for c in conn
        if _PM.ref(pm, nw, :meas, msr.msr_id, "dst")[c] != 0.0

            μ_tan = tan(_DST.mean(_PM.ref(pm, nw, :meas, msr.msr_id, "dst")[c]))
            σ_tan = abs(sec(μ_tan)*_DST.std(_PM.ref(pm, nw, :meas, msr.msr_id, "dst")[c]))
            _PM.ref(pm, nw, :meas, msr.msr_id, "dst")[c] = _DST.Normal( μ_tan, σ_tan )

            if msr.cmp_type == :branch
                num = _PM.var(pm, nw, msr.numerator, (msr.cmp_id, _PM.ref(pm, nw, :branch,msr.cmp_id)["f_bus"], _PM.ref(pm, nw, :branch,msr.cmp_id)["t_bus"]))
                den = _PM.var(pm, nw, msr.denominator, (msr.cmp_id, _PM.ref(pm, nw, :branch,msr.cmp_id)["f_bus"], _PM.ref(pm, nw, :branch,msr.cmp_id)["t_bus"]))
            else
                num = _PM.var(pm, nw, msr.numerator, msr.cmp_id)
                den = _PM.var(pm, nw, msr.numerator, msr.cmp_id)
            end
            msr.cmp_type == :branch ? id = (msr.cmp_id, _PM.ref(pm,nw,:branch,msr.cmp_id)["f_bus"], _PM.ref(pm,nw,:branch,msr.cmp_id)["t_bus"]) : id = msr.cmp_id
            JuMP.@NLconstraint(pm.model,
                original_var[id][c]*den[c] == num[c]
                )
        end
    end
end

function create_conversion_constraint(pm::_PM.AbstractPowerModel, original_var, msr::Fraction; nw=nw)
    num = []
    for n in msr.numerator
        if occursin("v", String(n))
            push!(num, _PM.var(pm, nw, n, msr.bus_ind))
        elseif !occursin("v", String(n)) && msr.cmp_type != :branch
            push!(num, _PM.var(pm, nw, n, msr.cmp_id))
        else
            push!(num, _PM.var(pm, nw, n, (msr.cmp_id, _PM.ref(pm, nw, :branch,msr.cmp_id)["f_bus"], _PM.ref(pm, nw, :branch,msr.cmp_id)["t_bus"])))
        end
    end

    den = _PM.var(pm, nw, msr.denominator, msr.bus_ind)
    msr.cmp_type == :branch ? id = (msr.cmp_id,  msr.bus_ind, _PM.ref(pm,nw,:branch,msr.cmp_id)["t_bus"]) : id = msr.cmp_id
    conn = get_active_connections(pm, nw, msr.cmp_type, msr.cmp_id)

    if occursin("r", String(msr.msr_type))
        JuMP.@NLconstraint(pm.model, [c in conn],
            original_var[id][c]*den[c] == num[1][c]*cos(num[3][c])+num[2][c]*sin(num[3][c])
            )
    elseif occursin("i", String(msr.msr_type))
        JuMP.@NLconstraint(pm.model, [c in conn],
            original_var[id][c]*den[c] == -num[2][c]*cos(num[3][c])+num[1][c]*sin(num[3][c])
            )
    else
        error("wrong measurement association")
    end
end

function create_conversion_constraint(pm::_PM.AbstractPowerModel, original_var, msr::MultiplicationFraction; nw=nw)

    p = []
    v = []
    for pw in msr.power
        if msr.cmp_type == :branch
            push!(p, _PM.var(pm, nw, pw, (msr.cmp_id, _PM.ref(pm, nw, :branch,msr.cmp_id)["f_bus"], _PM.ref(pm, nw, :branch,msr.cmp_id)["t_bus"])))
        else
            push!(p, _PM.var(pm, nw, pw, msr.cmp_id))
        end
    end
    for vl in msr.voltage
        push!(v, _PM.var(pm, nw, vl, msr.bus_ind))
    end

    msr.cmp_type == :branch ? id = (msr.cmp_id,  _PM.ref(pm, nw, :branch,msr.cmp_id)["f_bus"], _PM.ref(pm,nw,:branch,msr.cmp_id)["t_bus"]) : id = msr.cmp_id
    conn = get_active_connections(pm, nw, msr.cmp_type, msr.cmp_id)

    if occursin("cr", string(msr.msr_type))
        JuMP.@NLconstraint(pm.model, [c in conn],
            original_var[id][c] == (p[1][c]*v[1][c]+p[2][c]*v[2][c])/(v[1][c]^2+v[2][c]^2)
            )
    elseif occursin("ci", string(msr.msr_type))
        JuMP.@NLconstraint(pm.model, [c in conn],
            original_var[id][c] == (-p[2][c]*v[1][c]+p[1][c]*v[2][c])/(v[1][c]^2+v[2][c]^2)
            )
    end
end