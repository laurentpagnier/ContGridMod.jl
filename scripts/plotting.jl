function ctr_plot(isin::BitMatrix,
    values::Array{Float64,2};
    xlim::Tuple = (0.0, 0.0),
    ylim::Tuple = (0.0, 0.0),
    fill::Bool = true
)
    temp = copy(values)
    temp[.!isin] .= NaN
    if(ylim == (0.0, 0.0) || xlim == (0.0, 0.0))
        return contour(temp, fill=fill)
    else
        return contour(temp, fill=fill, xlim=xlim, ylim=ylim)
    end
end


function hm_plot(
    isin::BitMatrix,
    xrange::Array{Float64,1},
    yrange::Array{Float64,1},
    values::Array{Float64,2};
    clim::Tuple = (0.0,0.0),
    c = :viridis
)
    temp = copy(values)
    temp[.!isin] .= NaN
    if(clim == (0.0, 0.0))
        return heatmap(xrange, yrange, temp, c = c)
    else
        return heatmap(xrange, yrange, temp, c = c, clim = clim)
    end
end


function time_plot(
    time::Array{Float64, 1},
    cont_value::Array{Float64, 3}, # timeseries from the continous model
    coord::Array{Float64, 2}; # locations from where we wantg to fetch data
    borders::Array{Array{Float64,2},1} = Array{Float64,2}[],
    xlab::String = String("\$t[s]\$"),
    tstart::Float64 = 0.0,
    tend::Float64 = 0.0
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
    idin = findall(isinside)
    cont_coord = zeros(length(idin), 2)
    for i = 1:length(idin)
        cont_coord[i,:] = [xrange[idin[i][2]], yrange[idin[i][1]]]
    end
    p1 = Plots.Plot()
    for k in 1:size(coord, 1)
        dx = cont_coord[:, 1] .- coord[k, 1]
        dy = cont_coord[:, 2] .- coord[k, 2]
        id = argmin(dx.^2 + dy.^2)
        if(k == 1)
            p1 = plot(time[idstart:idend], cont_value[idin[id][1], idin[id][2], idstart:idend])
        else
            p1 = plot!(time[idstart:idend], cont_value[idin[id][1], idin[id][2], idstart:idend])
        end
    end
    plot!(legend = false)
     
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


function ctr_movie(
    ts::Array{Float64, 1},
    cont_value::Array{Float64, 3}; # timeseries from the continous model
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
        #ctr_plot(isin, values)
        contour(cont_value[:,:,t],fill=true)
    end
end
