using Plots


function hm_plot(
    isgrid::BitArray,
    Ny::Int64,
    Nx::Int64,
    values::Array{Float64, 1}; # timeseries from the continous model
    clim::Tuple = (0.0,0.0),
    c = :inferno,
    cb_title::String = ""
)
    temp = copy(values)
    temp[.!isgrid] .= NaN
    if(clim == (0.0, 0.0))
        return Plots.heatmap(reshape(temp, Ny, Nx), c = c, colorbar_title=cb_title)
    else
        return Plots.heatmap(reshape(temp, Ny, Nx), c = c, clim = clim, colorbar_title=cb_title)
    end
end


function time_plot(
    time::Array{Float64, 1},
    cont_value::Array{Float64, 2}, # timeseries from the continous model
    grid_coord::Array{Float64, 2}, # t
    coord::Array{Float64, 2}; # locations from where we wantg to fetch data
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

    p1 = Plots.Plot()
    for k in 1:size(coord, 1)
        dx = grid_coord[:, 1] .- coord[k, 1]
        dy = grid_coord[:, 2] .- coord[k, 2]
        id = argmin(dx.^2 + dy.^2)
        if(k == 1)
            p1 = plot(time[idstart:idend], cont_value[id, idstart:idend])
        else
            p1 = plot!(time[idstart:idend], cont_value[id, idstart:idend])
        end
    end
    plot!(legend = false, xlabel = xlabel, ylabel = ylabel)
     
    p2 = Plots.Plot()
    for k in 1:size(coord, 1)
        if(k == 1)
            p2 = scatter([coord[k, 1]], [coord[k, 2]])
        else
            p2 = scatter!([coord[k, 1]], [coord[k, 2]])
        end
    end
    for k in 1:length(borders)
        p2 = plot!(borders[k][:, 1], borders[k][:, 2], color=:black,)
    end
    plot!(legend = false)
    plot(p1, p2, layout=(1, 2), size=(800,300))
end


function hm_movie(
    isgrid::BitArray,
    Ny::Int64,
    Nx::Int64,
    ts::Array{Float64, 1},
    cont_value::Array{Float64, 2}; # timeseries from the continous model
    tstart::Float64 = 0.0,
    tend::Float64 = 0.0,
    interval::Int64 = 1
)
    # !!!!!!!!!!!!! STILL SOME WORK TO DO
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
        temp = zeros(size(isgrid))
        temp[isgrid] = cont_value[:,t]
        temp[.!isgrid] .= NaN
        clim = (minimum(cont_value), maximum(cont_value))
        heatmap(reshape(temp, Ny, Nx), fill=true, clim=clim)
    end
end

function disc_plot(
    coord::Array{Float64, 2},
    values::Array{Float64, 1}
)
    temp = copy(values)
    temp .-= minimum(temp)
    temp ./= maximum(temp)
    C(g::ColorGradient) = RGB[g[z] for z=temp]
    g = :inferno
    scatter(coord[:,2], coord[:,1], color=(cgrad(g) |> C), legend=false)
    #scatter([0,0], [0,1], zcolor=[0,3], clims=(0.0,1.0),
    #             xlims=(1,1.1), xshowaxis=false, yshowaxis=false, label="", grid=false)
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
    
    #println(size(grid_coord))
    #println(size(cont_values))
    #println(size(cont_time))
    #println(size(disc_coord))
    #println(size(disc_values))
    #println(size(disc_time))
    #println(size(plot_coord))

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
