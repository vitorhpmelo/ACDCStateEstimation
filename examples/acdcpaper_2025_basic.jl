import ACDCStateEstimation as _ACDCSE

import Ipopt
import PowerModels as _PM
import PowerModelsMCDC as _PMMCDC
import Plots

using LaTeXStrings

include("utils.jl")

nlp_optimizer = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 5, "sb" => "yes"
)

data_se_noiseless = _ACDCSE.quickget_case5()
data_se_noisy = _ACDCSE.quickget_case5()
data_pf = _ACDCSE.quickget_case5()

result, σ_dict, data_se_noiseless = generate_data_basic_acdcse(data_pf, data_se_noiseless, nlp_optimizer, sample_error = false);
result, σ_dict, data_se_noisy = generate_data_basic_acdcse(data_pf, data_se_noisy, nlp_optimizer, sample_error = true);

se_res_noiseless = _ACDCSE.solve_acdcse(data_se_noiseless, _PM.ACPPowerModel, nlp_optimizer)
se_res_noisy = _ACDCSE.solve_acdcse(data_se_noisy, _PM.ACPPowerModel, nlp_optimizer)

ground_truth_1 = [busdc["vm"][1] for (b,busdc) in result["solution"]["busdc"]]
ground_truth_2 = [busdc["vm"][2] for (b,busdc) in result["solution"]["busdc"]]
ground_truth_0 = [busdc["vm"][3] for (b,busdc) in result["solution"]["busdc"]]

se_noiseless_1 = [busdc["vm"][1] for (b,busdc) in se_res_noiseless["solution"]["busdc"]]
se_noiseless_2 = [busdc["vm"][2] for (b,busdc) in se_res_noiseless["solution"]["busdc"]]
se_noiseless_0 = [busdc["vm"][3] for (b,busdc) in se_res_noiseless["solution"]["busdc"]]

se_noisy_1 = [busdc["vm"][1] for (b,busdc) in se_res_noisy["solution"]["busdc"]]
se_noisy_2 = [busdc["vm"][2] for (b,busdc) in se_res_noisy["solution"]["busdc"]]
se_noisy_0 = [busdc["vm"][3] for (b,busdc) in se_res_noisy["solution"]["busdc"]]


Plots.scatter(vcat([ground_truth_1, ground_truth_2, ground_truth_0]), marker = :X)

Plots.scatter!([se_noiseless_1, se_noiseless_2, se_noiseless_0], marker = :circle)

Plots.scatter!([se_noisy_1, se_noisy_2, se_noisy_0], marker = :diamond)

Plots.scatter([bus["vm"] for (b,bus) in result["solution"]["bus"]], marker = :circle, label = L"\textrm{True}", ylabel = L"\textrm{U}^{\textrm{mag}} ~~ \textrm{[p.u.]}", color="white", ms=9)
Plots.scatter!([bus["vm"] for (b,bus) in se_res_noiseless["solution"]["bus"]], marker = :X, label = L"\textrm{Noiseless~~SE~~results}", color="black", ms=8)
Plots.scatter!([bus["vm"] for (b,bus) in se_res_noisy["solution"]["bus"]], marker = :diamond, label = L"\textrm{Noisy~~SE~~results}", color="grey", ms=5)
Plots.plot!(xlabel=L"\textrm{AC~~bus~~id~~[-]}", xticks=(1:11, [L"%$i" for i in 1:11]), yticks = (1.06:0.01:1.1, [L"%$i" for i in 1.06:0.01:1.1]),
            legendfontsize=12, xlabelfontsize=13, ylabelfontsize=13, tickfontsize=12)

Plots.savefig("AC-bus-results.pdf")

data_se_noiseless_39 = _ACDCSE.quickget_case39()
data_se_noisy_39 = _ACDCSE.quickget_case39()
data_pf_39 = _ACDCSE.quickget_case39()

result, σ_dict, data_se_noiseless_39 = generate_data_basic_acdcse(data_pf_39, data_se_noiseless_39, nlp_optimizer, sample_error = false);
result, σ_dict, data_se_noisy_39 = generate_data_basic_acdcse(data_pf_39, data_se_noisy_39, nlp_optimizer, sample_error = true);

se_res_noiseless_39 = _ACDCSE.solve_acdcse(data_se_noiseless_39, _PM.ACPPowerModel, nlp_optimizer)
se_res_noisy_39 = _ACDCSE.solve_acdcse(data_se_noisy_39, _PM.ACPPowerModel, nlp_optimizer)

data_se_noiseless_67 = _ACDCSE.quickget_case67()
data_se_noisy_67 = _ACDCSE.quickget_case67()
data_pf_67 = _ACDCSE.quickget_case67()

result, σ_dict, data_se_noiseless_67 = generate_data_basic_acdcse(data_pf_67, data_se_noiseless_67, nlp_optimizer, sample_error = false);
result, σ_dict, data_se_noisy_67 = generate_data_basic_acdcse(data_pf_67, data_se_noisy_67, nlp_optimizer, sample_error = true);

se_res_noiseless_67 = _ACDCSE.solve_acdcse(data_se_noiseless_67, _PM.ACPPowerModel, nlp_optimizer)
se_res_noisy_67 = _ACDCSE.solve_acdcse(data_se_noisy_67, _PM.ACPPowerModel, nlp_optimizer)