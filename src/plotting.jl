using Plots
using FerriteViz

export hm_plot, time_plot, hm_movie, disc_plot, cmp_plot, nodal_plot


function nodal_plot(model::ContModelFer, fieldname::Symbol; colormap=:inferno, colorbar=true)::Figure
    f = Figure()
    Axis(f[1, 1])
    sp = solutionplot!(model.dh₁, getproperty(model, fieldname), colormap=colormap)
    if colorbar
        Colorbar(f[1, 2], sp)
    end
    return f
end

function nodal_plot(model::ContModelFer, val::Vector{T}; colormap=:inferno, colorbar=true)::Figure where {T}
    f = Figure()
    Axis(f[1, 1])
    sp = solutionplot!(model.dh₁, val, colormap=colormap)
    if colorbar
        Colorbar(f[1, 2], sp)
    end
    f
end

function hm_plot(
    contmodel::ContModel,
    field::String;
    border::Array{Float64,2} = zeros(0,2),
    clim::Tuple{Float64, Float64} = (0.0,0.0),
    c = :inferno,
    cb_title::String = "",
    colorbar = true,
)
    global cm = contmodel
    eval(Meta.parse("temp = cm.$(field)"))
    value = zeros(contmodel.mesh.Ny, contmodel.mesh.Nx)
    value[contmodel.mesh.is_grid] = temp
    value[.!contmodel.mesh.is_grid] .= NaN
    if(clim == (0.0, 0.0))
        Plots.heatmap(contmodel.mesh.xrange, contmodel.mesh.yrange,
            reshape(value, contmodel.mesh.Ny, contmodel.mesh.Nx),
            #c = c, colorbar_title=cb_title)
            c = c, colorbar = colorbar)
    else
        Plots.heatmap(contmodel.mesh.xrange, contmodel.mesh.yrange,
        reshape(value, contmodel.mesh.Ny, contmodel.mesh.Nx),
        #c = c, clim = clim, colorbar_title=cb_title)
        c = c, clim = clim, colorbar = colorbar)
    end
    if(size(border,1) > 0)
        x1 = 0.42; x2 = -.62; y1 = 0.37; y2 = -0.4
        plot!(Shape([border[:, 1]; border[end,1];x1;x1;x2;x2;border[end,1]],
        [border[:, 2];y1; y1;y2;y2;y1;y1]); color = :white, linecolor = :white)
    end
    plot!(legend = false, grid = false, showaxis = :hide, xaxis = nothing,
        yaxis = nothing, aspect_ratio = :equal)
end


function time_plot(
    contmod::ContModel,
    time::Vector{Float64},
    values::Matrix{Float64}, # timeseries from the continous model
    gps_coord::Matrix{Float64}; 
    borders::Vector{Matrix{Float64}} = Matrix{Float64}[],
    xlabel::String = String("\$t\\;[s]\$"),
    ylabel::String = String("\$\\omega \$"),
    tlim::Tuple{Float64,Float64} = (0.0, 0.0),
)
    tlim[1] != 0.0 ? idstart = findall(tlim[1] .< time)[1] : id1= 1
    tlim[2] != 0.0 ? idend = findall(time .< tlim[2])[end] : id2 = length(time) 
    t_id = id1:id2
    
    cp = palette(:tab10)
    grid_coord = contmod.mesh.coord
    plot_id = find_node_from_gps(contmod, gps_coord)
    
    p1 = plot(time[t_id], values[plot_id[1], t_id], color=cp[1])
    for k in 2:length(plot_id)
        p1 = plot!(time[t_id], values[plot_id[k], t_id], color=cp[k])
    end
    
    plot!(legend = false, xlabel = xlabel, ylabel = ylabel,
        grid = false, linewidth = 1)
    
    p2 = scatter([grid_coord[plot_id[1], 1]], [grid_coord[plot_id[1], 2]],
        color = cp[1], markerstrokecolor = cp[1], markersize = 5)
    for k in 1:length(plot_id)
        p2 = scatter!([grid_coord[plot_id[k], 1]], [grid_coord[plot_id[k], 2]],
            color = cp[k], markerstrokecolor = cp[k], markersize = 5)
    end
    for k in 1:length(borders)
        p2 = plot!(borders[k][:, 1], borders[k][:, 2], color=:black,
        grid=false, showaxis=:hide, xaxis=nothing, yaxis=nothing, linewidth=3.0, aspect_ratio=:equal)
    end
    plot!(legend = false)
    plot(p1, p2, layout=(1, 2), size=(800,300))
end
function disc_plot(
    coord::Array{Float64, 2},
    values::Array{Float64, 1};
    borders::Array{Array{Float64,2},1} = Array{Float64,2}[],
    clim::Tuple{Float64, Float64} = (0.0,0.0),
    c::Symbol = :inferno,
    cbar::Bool = true,
    cb_title::String = "",
    markersize = 6.0
)
    Plots.plot()
    g = :inferno
    for k in 1:length(borders)
        p2 = Plots.plot!(borders[k][:, 1], borders[k][:, 2], color=:black, lw=3.0)
    end
    if (clim==(0.0,0.0))
        Plots.scatter!(coord[:,2], coord[:,1], zcolor=values, legend=false, grid=false,
        msw=0, showaxis=:hide, xaxis=nothing, yaxis=nothing, markersize=markersize, c=c, cbar=cbar,
        cb_title=cb_title)
    else
        Plots.scatter!(coord[:,2], coord[:,1], zcolor=values, legend=false, clim=clim, grid=false,
        msw=0, showaxis=:hide, xaxis=nothing, yaxis=nothing, markersize=markersize, c=c, cbar=:bottom,
        cb_title=cb_title)
    end

end
