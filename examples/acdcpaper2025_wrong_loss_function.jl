import ACDCStateEstimation as _ACDCSE

import Ipopt
import PowerModels as _PM
import PowerModelsMCDC as _PMMCDC

include("utils.jl")

nlp_optimizer = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 0, "sb" => "yes"
)

data_se = _ACDCSE.quickget_case5()
data_pf = _ACDCSE.quickget_case5()

result, σ_dict, data_se = generate_data_basic_acdcse(data_pf, data_se, nlp_optimizer, sample_error = true);

introduce_error_conv_losses!(data_se)

se_res = _ACDCSE.solve_acdcse(data_se, _PM.ACPPowerModel, nlp_optimizer)

Plots.plot()

