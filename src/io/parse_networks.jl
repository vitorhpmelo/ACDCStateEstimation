"""
    prepare_data_for_se_default!(data::Dict)

    Simple function that modifies the network data in the `data` dictionary
        adding standard assumptions and settings for SE calculations.
"""
function prepare_data_for_se_default!(data::Dict; exceptions::Vector{Int} = Int[])
    data["se_settings"] = Dict(
        "rescaler" => 1,
        "criterion" => "rwlav"
    )
    PV2PQbuses!(data)
    remove_slack_buses!(data, exceptions=exceptions)
end
"""
All PV buses replaced as PQ buses.
"""
function PV2PQbuses!(data::Dict)
    for (b, bus) in data["bus"]
        if bus["bus_type"] == 2
            bus["bus_type"] = 1
        end
    end
end
"""
Removes all buses of type 3 except those is the `exception` 
Voltage angles of the former slack buses are estimated, not fixed).
"""
function remove_slack_buses!(data::Dict; exceptions::Vector{Int} = Int[])
    slackbuses = [bus["index"] for (_, bus) in data["bus"] if bus["bus_type"] == 3]
    exceptions = isempty(exceptions) ? slackbuses[1] : exceptions
    for b in setdiff(Set(slackbuses), Set(exceptions)) 
        data["bus"]["$b"]["bus_type"] = 1
    end
end
"""
Shortcut to get the network data of case5_2grids_MC.m from local PowerModelsMCDC.jl
"""
quickget_case5() = _PMMCDC.parse_file(joinpath(_PMMCDC_dir, "test/data/matacdc_scripts/case5_2grids_MC.m"))
"""
Shortcut to get the network data of case67mcdc_scopf4.m from local PowerModelsMCDC.jl
"""
quickget_case67() = _PMMCDC.parse_file(joinpath(_PMMCDC_dir, "test/data/matacdc_scripts/case67mcdc_scopf4.m"))
"""
Shortcut to get the network data of case39_mcdc.m from local PowerModelsMCDC.jl
"""
quickget_case39() = _PMMCDC.parse_file(joinpath(_PMMCDC_dir, "test/data/matacdc_scripts/case39_mcdc.m"))
"""
Shortcut to get the network data of case3120sp_mcdc.m from local PowerModelsMCDC.jl
"""
quickget_case3120() = _PMMCDC.parse_file(joinpath(_PMMCDC_dir, "test/data/matacdc_scripts/case3120sp_mcdc.m"))