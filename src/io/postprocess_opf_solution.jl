"""
Adds the dc power values to a PMMCDC (o)pf result dictionary
`sol`. The original network dictionary `data` is used to look up the
bus connection of the dc branch
"""
function get_dc_power!(sol::Dict, data::Dict)
    for (br, dcbranch) in sol["solution"]["branchdc"]
        bf = data["branchdc"][br]["fbusdc"]
        bt = data["branchdc"][br]["tbusdc"]
        if data["branchdc"][br]["line_confi"] == 2
            dcbranch["pf"] = dcbranch["i_from"].*sol["solution"]["busdc"]["$bf"]["vm"]
            dcbranch["pt"] = dcbranch["i_to"].*sol["solution"]["busdc"]["$bt"]["vm"]
        else
            if data["branchdc"][br]["connect_at"] == 0
                dcbranch["pf"] = dcbranch["i_from"].*sol["solution"]["busdc"]["$bf"]["vm"][[1,2]]
                dcbranch["pt"] = dcbranch["i_to"].*sol["solution"]["busdc"]["$bt"]["vm"][[1,2]]
            elseif data["branchdc"][br]["connect_at"] == 1
                dcbranch["pf"] = dcbranch["i_from"].*sol["solution"]["busdc"]["$bf"]["vm"][[1,3]]
                dcbranch["pt"] = dcbranch["i_to"].*sol["solution"]["busdc"]["$bt"]["vm"][[1,3]]
            else
                dcbranch["pf"] = dcbranch["i_from"].*sol["solution"]["busdc"]["$bf"]["vm"][[2,3]]
                dcbranch["pt"] = dcbranch["i_to"].*sol["solution"]["busdc"]["$bt"]["vm"][[2,3]]
            end
        end
    end
end