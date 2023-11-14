function no_conversion_needed(pm::_PM.AbstractACPModel, msr_var::Symbol)
    return msr_var ∈ [:vm, :va, :pd, :qd, :pg, :qg, :p, :q]
  end
  
  function no_conversion_needed(pm::_PM.AbstractACRModel, msr_var::Symbol)
    return msr_var ∈ [:vr, :vi, :pd, :qd, :pg, :qg, :p, :q]
  end
  
  function no_conversion_needed(pm::_PM.AbstractIVRModel, msr_var::Symbol)
    return msr_var ∈ [:vr, :vi, :cr, :ci, :crg, :cig, :crd, :cid]
  end