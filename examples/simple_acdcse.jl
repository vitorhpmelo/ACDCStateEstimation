import ACDCStateEstimation as _ACDCSE

import Ipopt
import PowerModels as _PM
import PowerModelsMCDC as _PMMCDC

nlp_optimizer = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 0, "sb" => "yes"
)

_PMMCDC_dir = dirname(dirname(pathof(_PMMCDC)))


data_se = _PMMCDC.parse_file(joinpath(_PMMCDC_dir, "test/data/matacdc_scripts/case5_2grids_MC.m"))
data_pf = _PMMCDC.parse_file(joinpath(_PMMCDC_dir, "test/data/matacdc_scripts/case5_2grids_MC.m"))

# if data_se is the same as data_pf (which is the "real system"), no problem
# below, some room to play around and see what happens if it is not the case (e.g., branch has been open but se people don't know)

deenergize_dc_branches!(data_pf, 1)

result = _PMMCDC.solve_mcdcopf(data_pf, _PM.ACPPowerModel, nlp_optimizer)
_ACDCSE.get_dc_power!(result, data_pf)

σ_dict = Dict{String, Float64}(
    "vm"       => 0.01,
    "va"       => 0.01,
    "pg"       => 0.01,
    "qg"       => 0.01,
    "pd"       => 0.01,
    "qd"       => 0.01,
    "p_ac"     => 0.01,
    "q_ac"     => 0.01,
    "vdcm"     => 0.01,
    "p_dc"     => 0.01,
    "i_dcgrid" => 0.01
)

_ACDCSE.powerflow2measurements!(data_se, result, σ_dict, sample_error=false, measurements = ["vm", "va", "p_g","p_fr", "q_fr", "vdcm", "p_dc_fr", "pg", 
                                "pd", "qg", "qd", "i_dcgrid_fr"])

_ACDCSE.prepare_data_for_se_default!(data_se)

se_res = _ACDCSE.solve_acdcse(data_se, _PM.ACPPowerModel, nlp_optimizer)
_ACDCSE.get_dc_power!(se_res, data_se)