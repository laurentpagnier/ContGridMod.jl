export hm_plot, time_plot, hm_movie, disc_plot, cmp_plot

using Plots


function hm_plot(
    contmodel::ContModel,
    field::Symbol;
    border::Array{Float64,2} = zeros(0,2),
    clim::Tuple{Float64, Float64} = (0.0,0.0),
    cs = :inferno,
    cb_title::String = "",
    colorbar = true,
)
    
    cval = copy(getfield(contmodel,field))
    if clim == (0.0, 0.0)
        clim = (minimum(cval), maximum(cval))
    end
    
    cval .-= clim[1]
    cval ./= (clim[2] - clim[1])

    plot(legend=false)
    for (k, v) in enumerate(contmodel.mesh.cell_vertices)
        x = 1:length(v) .|> i -> v[i][1]
        y = 1:length(v) .|> i -> v[i][2]
        plot!(Shape(x,y), linecolor=:transparent, color=cgrad(cs)[cval[k]])
    end

    current()
end


function time_plot(
    contmod::ContModel,
    t::Vector{Float64},
    values::Matrix{Float64}, # timeseries from the continous model
    gps_coord::Vector{Tuple{Float64,Float64}}; 
    xlabel::String = String("\$t\\;[s]\$"),
    ylabel::String = String("\$\\omega \$"),
    tlim::Tuple{Float64,Float64} = (0.0, 0.0),
)
    tlim[1] != 0.0 ? idstart = findall(tlim[1] .< t)[1] : idstart = 1
    tlim[2] != 0.0 ? idend = findall(t .< tlim[2])[end] : idend = length(t) 
    t_id = idstart:idend
    cp = palette(:tab10)
    plot_id = find_node_from_gps(contmod, gps_coord)

    plot()
    for k in 1:length(plot_id)
        plot!(t[t_id], values[t_id, plot_id[k]], color=cp[k])
    end

    p1 = plot!(legend = false, xlabel = xlabel, ylabel = ylabel,
        grid = false, linewidth = 1)
    
    plot()
    for (k, id) in enumerate(plot_id)
        scatter!([contmod.mesh.node_coord[id][1]], [contmod.mesh.node_coord[id][2]],
            color = cp[k], markerstrokecolor = cp[k], markersize = 5)
    end
    temp = to_mat(contmod.mesh.border)
    p2 = plot!(temp[:, 1], temp[:, 2], color=:black, grid=false, showaxis=:hide,
        xaxis=nothing, yaxis=nothing, linewidth=3.0, aspect_ratio=:equal, legend = false)
    plot(p1, p2, layout=(1, 2), size=(800,300))
end
