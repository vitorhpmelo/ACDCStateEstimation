"""
Adds the dc power values to a PMMCDC (o)pf result dictionary
`sol`. The original network dictionary `data` is used to look up the
bus connection of the dc branch
"""
function get_dc_power!(sol::Dict, data::Dict)
    if length(data["branchdc"])==0
        @warn "DC branches dictionary exists but it is empty."
        return
    end
    for (br, dcbranch) in sol["solution"]["branchdc"]
        bf = data["branchdc"][br]["fbusdc"]
        bt = data["branchdc"][br]["tbusdc"]
        if data["branchdc"][br]["line_confi"] == 2
            dcbranch["pf"] = dcbranch["i_from"].*sol["solution"]["busdc"]["$bf"]["vm"]
            dcbranch["pt"] = dcbranch["i_to"].*sol["solution"]["busdc"]["$bt"]["vm"]
        else
            if data["branchdc"][br]["connect_at"] == 0 #what does the connect at field means
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

function get_ac_currents!(sol::Dict, data::Dict)

    j=complex(0,1)
    if length(data["branch"])==0
        @warn "AC branches dictionary exists but it is empty."
        return
    end
    
    for (br, branch) in sol["solution"]["branch"]
        bf = data["branch"][br]["f_bus"]
        bt = data["branch"][br]["t_bus"]
        

        
        If=conj((branch["pf"] + j*branch["qf"])/(sol["solution"]["bus"]["$bf"]["vm"]*exp(j*sol["solution"]["bus"]["$bf"]["va"])))
        
        It=conj((branch["pt"] + j*branch["qt"])/(sol["solution"]["bus"]["$bt"]["vm"]*exp(j*sol["solution"]["bus"]["$bt"]["va"])))

        branch["cmf"] = abs(If)
        branch["caf"] = angle(If)
        branch["cmt"] = abs(It)
        branch["cat"] = angle(It)


        branch["crf"] = real(If)
        branch["cif"] = imag(If)
        branch["crt"] = real(It)
        branch["cit"] = imag(It)
    

    end

    
end



function get_convacdc_m!(sol::Dict, data::Dict)

    if length(data["convdc"])==0
        @warn "converters dictionary exists but it is empty."
        return
    end
    for (key, conv) in sol["solution"]["convdc"]
        bus_dc = data["convdc"][key]["busdc_i"]
        conv_confi = data["convdc"][key]["conv_confi"]
        connect_at =data["convdc"][key]["connect_at"]
        # Find the key in data_pf["bus_dc"] where "busdc_i" equals bus_dc
        key_busdc_match = nothing
        for (key_busdc, busdc) in data["busdc"]
            if busdc["busdc_i"] == bus_dc
                key_busdc_match = key_busdc
                break
            end
        end
        if conv_confi == 2
            m1 = conv["vmconv"][1]/(sol["solution"]["busdc"][key_busdc_match]["vm"][1]-sol["solution"]["busdc"][key_busdc_match]["vm"][3])
            m2 = conv["vmconv"][2]/(sol["solution"]["busdc"][key_busdc_match]["vm"][3]-sol["solution"]["busdc"][key_busdc_match]["vm"][2])
            conv["m"]=Vector([m1,m2])
        end
        if conv_confi == 1
            if connect_at ==0
                m=conv["vmconv"][1]/(sol["solution"]["busdc"][key_busdc_match]["vm"][1]-sol["solution"]["busdc"][key_busdc_match]["vm"][2])
            elseif connect_at ==1
                m=conv["vmconv"][1]/(sol["solution"]["busdc"][key_busdc_match]["vm"][1]-sol["solution"]["busdc"][key_busdc_match]["vm"][3])
            else
                m=conv["vmconv"][1]/(sol["solution"]["busdc"][key_busdc_match]["vm"][3]-sol["solution"]["busdc"][key_busdc_match]["vm"][2])
            end
            conv["m"]=Vector([m])
        end
    end
end

