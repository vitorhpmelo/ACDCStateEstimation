import ACDCStateEstimation as _ACDCSE
import PowerModels as _PM
import PowerModelsMCDC as _PMMCDC
import Ipopt
using Test

nlp_optimizer = _PMMCDC.optimizer_with_attributes(
    Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 0, "sb" => "yes"
)

σ_dict = Dict{String, Float64}(
    "vm"       => 0.01, # (AC) voltage magnitude
    "va"       => 0.01, # (AC) voltage angle
    "pg"       => 0.01, # Active power injection from (AC) generator
    "qg"       => 0.01, # Reactive power injection from (AC) generator
    "pd"       => 0.01, # Active power injection from (AC) load
    "qd"       => 0.01, # Rective power injection from (AC) load
    "p_ac"     => 0.01, # Active power flow on an AC branch   (either injection or flow, really)
    "q_ac"     => 0.01, # Reactive power flow on an AC branch (either injection or flow, really)
    "vdcm"     => 0.01, # DC voltage magnitude
    "p_dc"     => 0.01, # DC power (either injection or flow)
    "i_dcgrid" => 0.01  # DC current (either injection or flow)
)

# TODO ↑ add ac current magnitude information


@testset "ACDCStateEstimation.jl" begin
    include("acdcse.jl")
end
