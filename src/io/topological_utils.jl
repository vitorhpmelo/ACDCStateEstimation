function deenergize_dc_branches!(data::Dict, id::Int)
    delete!(data["branchdc"], "$id")
end

function deenergize_dc_branches!(data::Dict, id::Vector{Int})
    for i in id
        delete!(data["branchdc"], "$i")
    end
end

function deenergize_ac_branches!(data::Dict, id::Int)
    delete!(data["branch"], "$id")
end

function deenergize_ac_branches!(data::Dict, id::Vector{Int})
    for i in id
        delete!(data["branch"], "$i")
    end
end

# TODO: open pole/contingency?