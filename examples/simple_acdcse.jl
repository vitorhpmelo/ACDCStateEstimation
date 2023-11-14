import ACDCStateEstimation as _ACDCSE

import Ipopt
import PowerModels as _PM
import PowerModelsMCDC as _PMMCDC

nlp_optimizer = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 0, "sb" => "yes"
)

_PMMCDC_dir = dirname(dirname(pathof(_PMMCDC)))

data = _PMMCDC.parse_file(joinpath(_PMMCDC_dir, "test/data/matacdc_scripts/case5_2grids_MC.m"))
result = _PMMCDC.solve_mcdcopf(data, _PM.ACPPowerModel, nlp_optimizer)
get_dc_power!(result, data)

σ_dict = Dict{String, Float64}(
    "vm"   => 0.01,
    "p_ac" => 0.01,
    "q_ac" => 0.01,
    "vdcm" => 0.01
)

_ACDCSE.powerflow2measurements!(data, result, σ_dict, sample_error=false)

data["se_settings"] = Dict(
    "rescaler" => 1,
    "criterion" => "rwlav"
)

se_res = _ACDCSE.solve_acdcse(data, _PM.ACPPowerModel, nlp_optimizer)