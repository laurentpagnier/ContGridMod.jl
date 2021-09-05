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
    cont_value::Array{Float64, 1}, # values from the continous model
    coord::Array{Float64, 2}; # locations from where we wantg to fetch data
    xlab::String = String("\$t[s]\$")
)
    # !!!!!!!!!!!!! STILL SOME WORK TO DO
    
    idin = findall(isinside)
    cont_coord = zeros(length(idin), 2)
    for i = 1:length(idin)
        cont_coord[i,:] = [xrange[idin[i][2]], yrange[idin[i][1]]]
    end

    for k in 1:size(coord, 1)
        dx = cont_coord[:, 1] .- coord[k, 1]
        dy = cont_coord[:, 2] .- coord[k, 2]
        dist = min.(
            sqrt.(dx.^2 + dy.^2),
            10*minimum(sqrt.(dx.^2 + dy.^2))
            ) # min is here to prevent NaN
        factor = exp.(-dist) ./ sum(exp.(-dist))
        cont_value[is,:]
        plot!(time, cont_value[i,j,:])
    end
end
