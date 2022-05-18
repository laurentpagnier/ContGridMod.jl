using Plots

export hm_plot, time_plot, hm_movie, disc_plot, cmp_plot


function hm_plot(
    contmodel::ContModel,
    field::String;
    borders::Array{Array{Float64,2},1} = Array{Float64,2}[],
    clim::Tuple{Float64, Float64} = (0.0,0.0),
    c = :inferno,
    cb_title::String = ""
)
    global cm = contmodel
    eval(Meta.parse("temp = cm.$(field)"))
    value = zeros(contmodel.mesh.Ny, contmodel.mesh.Nx)
    value[contmodel.mesh.isgrid] = temp
    value[.!contmodel.mesh.isgrid] .= NaN
    if(clim == (0.0, 0.0))
        Plots.heatmap(contmodel.mesh.xrange, contmodel.mesh.yrange,
            reshape(value, contmodel.mesh.Ny, contmodel.mesh.Nx),
            #c = c, colorbar_title=cb_title)
            c = c, colorbar=false)
    else
        Plots.heatmap(contmodel.mesh.xrange, contmodel.mesh.yrange,
        reshape(value, contmodel.mesh.Ny, contmodel.mesh.Nx),
        #c = c, clim = clim, colorbar_title=cb_title)
        c = c, clim = clim, colorbar=false)
    end
    for k in 1:length(borders)
        p2 = plot!(borders[k][:, 1], borders[k][:, 2], color=:black,linewidth=3.0)
    end
    plot!(legend=false, grid=false, showaxis=:hide, xaxis=nothing,
    yaxis=nothing, aspect_ratio=:equal)
end

#=
function hm_plot(
    contmodel::ContModel,
    values::Array{Float64, 1}; # timeseries from the continous model
    borders::Array{Array{Float64,2},1} = Array{Float64,2}[],
    clim::Tuple{Float64, Float64} = (0.0,0.0),
    c = :inferno,
    cb_title::String = ""
)
    temp = copy(values)
    temp[.!contmodel.isinside] .= NaN
    if(clim == (0.0, 0.0))
        Plots.heatmap(contmodel.yrange, contmodel.xrange,
            reshape(temp, contmodel.Ny, contmodel.Nx),
            #c = c, colorbar_title=cb_title)
            c = c, colorbar=false)
    else
        Plots.heatmap(contmodel.yrange, contmodel.xrange,
        reshape(temp, contmodel.Ny, contmodel.Nx),
        #c = c, clim = clim, colorbar_title=cb_title)
        c = c, clim = clim, colorbar=false)
    end
    for k in 1:length(borders)
        p2 = plot!(borders[k][:, 1], borders[k][:, 2], color=:black,linewidth=3.0)
    end
    plot!(legend=false, grid=false, showaxis=:hide, xaxis=nothing,
    yaxis=nothing, aspect_ratio=:equal)
end
=#


function time_plot(
    contmod::ContModel,
    time::Array{Float64, 1},
    values::Array{Float64, 2}, # timeseries from the continous model
    coord::Array{Float64, 2}; 
    borders::Array{Array{Float64,2},1} = Array{Float64,2}[],
    xlabel::String = String("\$t\\;[s]\$"),
    ylabel::String = String("\$\\omega \$"),
    tstart::Float64 = 0.0,
    tend::Float64 = 0.0
)
    if(tstart != 0.0)
        idstart = findall(tstart .< ts)[1]
    else
        idstart = 1 
    end
    if(tend != 0.0)
        idend = findall(time .< tend)[end]
    else
        idend = length(time) 
    end
    cp = palette(:tab10)
    grid_coord = contmod.coord[contmod.isgrid,:]

    p1 = Plots.Plot()
    for k in 1:size(coord, 1)
        dx = grid_coord[:, 1] .- coord[k, 2]
        dy = grid_coord[:, 2] .- coord[k, 1]
        id = argmin(dx.^2 + dy.^2)
        if(k == 1)
            p1 = plot(time[idstart:idend], values[id, idstart:idend], color=cp[k])
        else
            p1 = plot!(time[idstart:idend], values[id, idstart:idend], color=cp[k])
        end
    end
    plot!(legend = false, xlabel = xlabel, ylabel = ylabel, grid=false, linewidth=1)
     
    p2 = Plots.Plot()
    for k in 1:size(coord, 1)
        if(k == 1)
            p2 = scatter([coord[k, 1]], [coord[k, 2]], color=cp[k], markerstrokecolor=cp[k], markersize=5)
        else
            p2 = scatter!([coord[k, 1]], [coord[k, 2]], color=cp[k], markerstrokecolor=cp[k], markersize=5)
        end
    end
    for k in 1:length(borders)
        p2 = plot!(borders[k][:, 1], borders[k][:, 2], color=:black,
        grid=false, showaxis=:hide, xaxis=nothing, yaxis=nothing, linewidth=3.0, aspect_ratio=:equal)
    end
    plot!(legend = false)
    plot(p1, p2, layout=(1, 2), size=(800,300))
end


function hm_movie(
    contmod::ContModel,
    ts::Array{Float64, 1},
    values::Array{Float64, 2}; # timeseries from the continous model
    tstart::Float64 = 0.0,
    tend::Float64 = 0.0,
    interval::Int64 = 1
)
    if(tstart != 0.0)
        idstart = findall(tstart .< ts)[1]
    else
        idstart = 1 
    end
    if(tend != 0.0)
        idend = findall(ts .< tend)[end]
    else
        idend = length(ts) 
    end

    @gif for t in idstart:interval:idend
        temp = zeros(size(contmod.isgrid))
        temp[contmod.isgrid] = values[:,t]
        temp[.!contmod.isgrid] .= NaN
        clim = (minimum(values), maximum(values))
        heatmap(reshape(temp, contmod.Ny, contmod.Nx), fill=true, clim=clim, aspect_ratio=:equal)
    end
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
    plot()
    g = :inferno
    for k in 1:length(borders)
        p2 = plot!(borders[k][:, 1], borders[k][:, 2], color=:black, lw=3.0)
    end
    if (clim==(0.0,0.0))
        scatter!(coord[:,2], coord[:,1], zcolor=values, legend=false, grid=false,
        msw=0, showaxis=:hide, xaxis=nothing, yaxis=nothing, markersize=markersize, c=c, cbar=cbar,
        cb_title=cb_title)
    else
        scatter!(coord[:,2], coord[:,1], zcolor=values, legend=false, clim=clim, grid=false,
        msw=0, showaxis=:hide, xaxis=nothing, yaxis=nothing, markersize=markersize, c=c, cbar=:bottom,
        cb_title=cb_title)
    end

end


function cmp_plot(
    grid_coord::Array{Float64, 2},
    cont_values::Array{Float64, 2},
    cont_time::Array{Float64, 1},
    disc_coord::Array{Float64, 2},
    disc_values::Array{Float64, 2},
    disc_time::Array{Float64, 1},
    plot_coord::Array{Float64, 2};
    tstart::Float64 = 0.0,
    tend::Float64 = 0.0
)
    if(tstart != 0.0)
        ids1 = findall(tstart .< cont_time)[1]
        ids2 = findall(tstart .< disc_time)[1]
    else
        ids1 = 1 
        ids2 = 1 
    end
    if(tend != 0.0)
        ide1 = findall(cont_time .< tend)[end]
        ide2 = findall(disc_time .< tend)[end]
    else
        ide1 = length(cont_time) 
        ide2 = length(disc_time) 
    end

    for k in 1:size(plot_coord, 1)
        id2 = argmin((disc_coord[:, 1] .- plot_coord[k, 1]).^2 +
            (disc_coord[:, 2] .- plot_coord[k, 2]).^2)
        id1 = argmin((grid_coord[:, 1] .- disc_coord[id2, 1]).^2 +
            (grid_coord[:, 2] .- disc_coord[id2, 2]).^2)
        if(k == 1)
            plot(cont_time[ids1:ide1], cont_values[id1, ids1:ide1])
            plot!(disc_time[ids2:ide2], disc_values[id2, ids2:ide2], linestyle=:dash)
        else
            plot!(cont_time[ids1:ide1], cont_values[id1, ids1:ide1])
            plot!(disc_time[ids2:ide2], disc_values[id2, ids2:ide2], linestyle=:dash)
        end
    end
end
