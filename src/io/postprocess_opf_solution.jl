"""
Adds the dc power values to a PMMCDC (o)pf result dictionary
`sol`. The original network dictionary `data` is used to look up the
bus connection of the dc branch
"""
function get_dc_power!(sol::Dict, data::Dict)
    for (br, dcbranch) in sol["solution"]["branchdc"]
        bf = data["branchdc"][br]["fbusdc"]
        bt = data["branchdc"][br]["tbusdc"]
        dcbranch["pf"] = dcbranch["i_from"].*sol["solution"]["bus"]["$bf"]["vm"]
        dcbranch["pt"] = dcbranch["i_to"].*sol["solution"]["bus"]["$bt"]["vm"]
    end
end