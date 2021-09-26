using Plots


function hm_plot(
    isgrid::BitArray,
    Ny::Int64,
    Nx::Int64,
    ts::Array{Float64, 1},
    cont_value::Array{Float64, 1}; # timeseries from the continous model
    clim::Tuple = (0.0,0.0),
    c = :viridis
)
    temp = copy(values)
    temp[.!isgrid] .= NaN
    if(clim == (0.0, 0.0))
        return Plots.heatmap(temp, c = c)
    else
        return Plots.heatmap(temp, c = c, clim = clim)
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
