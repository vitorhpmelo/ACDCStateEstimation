

function flat_start(x, network_info)
    for (key, value) in network_info["vm_ruler"]
        if !occursin("-c", key)
            x[value] = 1.0
        else
            x[value] = 1.0  # perturb voltage magnitudes
        end
    end
    for (key, value) in network_info["vdc_ruler"]
        if occursin("dc-1", key)
            x[value] = 1.0
        elseif occursin("dc-2", key)
            x[value] = -1.0
        else
            x[value] = 0.0
        end
    end
end


function determine_reference_buses(reference_buses, network_info)
    cols_to_remove = Vector{Int}()
    for bus in reference_buses
        push!(cols_to_remove, network_info["va_ruler"][string(bus)])
    end
    return cols_to_remove
end

function determine_reference_buses(network_info)
    cols_to_remove = Vector{Int}()
    for bus in network_info["reference_buses"]
        push!(cols_to_remove, network_info["va_ruler"][string(bus)])
    end
    return cols_to_remove
end


function update_x!(x, dx, cols_to_remove)
    j = 1
    for i in eachindex(x)
        if i in cols_to_remove
            continue
        else
            x[i] = x[i] + dx[j]
            j += 1
        end
    end
end


function SE_WLS(network_info, d_meas_set; min_sig=1e-7, itmax=10,tol=1e-6)
    weight_z_autodiffse!(network_info, d_meas_set; min_sig=min_sig)

    h = create_h(network_info)
    g = create_g(network_info)
    W = create_W_autodiff(network_info, h, g)

    h_vec(x) = [f(x) for f in vcat(h, g)]

    x = zeros(length(network_info["x"]))

    flat_start(x, network_info)

    cols_to_remove = determine_reference_buses(network_info["reference_buses"], network_info)
    cols_to_keep = setdiff(1:length(network_info["x"]), cols_to_remove)

    z1_val = [m.z for m in network_info["z"]]
    z2_val = [0.0 for _ in network_info["g"]]

    z_val = vcat(z1_val, z2_val)
    m = length(network_info["z"]) + length(network_info["g"])
    n = length(network_info["x"])
    it = 0
    H=zeros(m,n)
    hx=vcat(h, g)
    while it < itmax
        it += 1

        compute_Jacobian!(hx, x, W, H, m)

        # Hx = H[:, cols_to_keep]
        Gain = H[:, cols_to_keep]' * W * H[:, cols_to_keep]
        dz = z_val - h_vec(x)
        b = H[:, cols_to_keep]' * W * dz
        dx = Gain \ b
        update_x!(x, dx, cols_to_remove)
        println("Norm of dx: ", norm(dx))
        if norm(dx) < tol
            break
        end
    end

    return x
end


function SE_WLS_timebench(network_info, d_meas_set; min_sig=1e-7, itmax=10,tol=1e-6)
    weight_z_autodiffse!(network_info, d_meas_set; min_sig=min_sig)

    h = create_h(network_info)
    g = create_g(network_info)
    W = create_W_autodiff(network_info)

    h_vec(x) = [f(x) for f in vcat(h, g)]
    network_info["h_vec"] = h_vec
    x = zeros(length(network_info["x"]))

    flat_start(x, network_info)

    cols_to_remove = determine_reference_buses(network_info["reference_buses"], network_info)
    cols_to_keep = setdiff(1:length(network_info["x"]), cols_to_remove)
    println(cols_to_remove)

    z1_val = [m.z for m in network_info["z"]]
    z2_val = [0.0 for _ in network_info["g"]]

    z_val = vcat(z1_val, z2_val)
    m = length(network_info["z"]) + length(network_info["g"])
    n = length(network_info["x"])
    it = 0
    H=zeros(m,n)
    hx=vcat(h, g)

    t_start = time_ns()
    dx_norms = zeros(Float64, itmax)
    while it < itmax
        it += 1

        compute_Jacobian!(hx, x, H, m)

        Gain = H[:, cols_to_keep]' * W * H[:, cols_to_keep]
        dz = z_val - h_vec(x)
        b = H[:, cols_to_keep]' * W * dz
        println("Rank of Gain Matrix: ", rank(Gain))
        dx = Gain \ b
        update_x!(x, dx, cols_to_remove)
        
        dx_norms[it] = norm(dx)
        println("Norm of dx: ", dx_norms[it])
        if dx_norms[it] < tol
            break
        end
    end
    t_end = time_ns()
    println((t_end - t_start) / 1e6, "ms")


    return x,dx_norms[1:it]
end


function build_prior_autodiff!(network_info_prior,network_info_posterior; min_sig=1e-7, itmax=10,tol=1e-6)
    
    weight_z_autodiffse!(network_info_prior, d_meas_set; min_sig=min_sig)

    h = create_h(network_info_prior)
    g = create_g(network_info_prior)
    W = create_W_autodiff(network_info_prior, h, g)

    h_vec(x) = [f(x) for f in vcat(h, g)]

    x = network_info_prior["x"]


    m = length(network_info_prior["z"]) + length(network_info_prior["g"])
    n = length(network_info_prior["x"])

    H=zeros(m,n)
    hx=vcat(h, g)
    compute_Jacobian!(hx, x, W, H, m)

    Pinv = H' * W * H
    network_info_posterior["prior"] = Dict("Pinv" => Pinv,
                                    "x" => x,
    )

end