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






using Downloads

function get_dataset(
	name::String;
	location::String = ENV["HOME"]
)
	if name == "panta"
		download("url", location * "panta")
    elseif name == "smthesle"
	
	else
		println(name * " is not a registered data set.")
	end
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



function compute_stable_sol_old!(
    cm::ContModel;
    interval::Int64 = 1000,
    Niter::Int64 = 15000,
    tol::Float64 = 1e-7
)
    p2 = cm.mesh.dx^2 * cm.p
    b = -vec(diag(cm.xi))
    xi  = cm.xi + sparse(1:length(b),1:length(b),b)
    @time begin
        for t in 1:Niter
            if(mod(t, interval) == 0)
                temp = copy(cm.th)
            end
            cm.th = (xi * cm.th + p2) ./ b            
            if(mod(t, interval) == 0)
                println( [t maximum(abs.(cm.th - temp))] )
                if( maximum(abs.(cm.th - temp)) < tol )
                    break
                end
            end
        end
    end
end


function get_mesh(
    border::Array{Float64, 2},
    dx::Float64
)
    xlim = [minimum(border[:, 1])-2*dx maximum(border[:, 1])+2*dx]
    ylim = [minimum(border[:, 2])-2*dx maximum(border[:, 2])+2*dx]
    xrange = collect(xlim[1]:dx:xlim[2])
    yrange = collect(ylim[1]:dx:ylim[2])

    Nx = length(xrange)
    Ny = length(yrange)

    isborder = falses(Ny, Nx)
    isinside = falses(Ny, Nx)
    isgrid = falses(Ny, Nx)

    # check if inside the shape
    for j = 1:Nx
        for i = 1:Ny
            p = [xrange[j] yrange[i]]
            isinside[i, j] = inPolygon(p, border)[1]
        end
    end
    
    # check if alone 
    for j=2:Nx-1
        for i=2:Ny-1
            if(isinside[i, j])
                if(isinside[i-1, j] + isinside[i+1, j] + isinside[i, j-1] + isinside[i, j+1] == 0)
                    isinside[i, j] = false
                end
            end
        end
    end
    
    # add a boundary layer encompassing the inside
    for j=2:Nx-1
       for i=2:Ny-1
            nx = isinside[i, j-1] - isinside[i, j+1]
            ny = isinside[i-1, j] - isinside[i+1, j]
            if((nx^2 + ny^2) > 0 && isinside[i,j] == false)
                isgrid[i, j] = true
            end
       end
    end
    
    isgrid = isgrid .| isinside
    # check if inside and considered as outside. This check isnt able
    # to find hole that are bigger than a single point. Som it should
    # be improved in the future
    for j=2:Nx-1
        for i=2:Ny-1
            if(.!isgrid[i, j])
                if(isgrid[i-1, j] + isgrid[i+1, j] + isgrid[i, j-1] + isgrid[i, j+1] == 4)
                    isgrid[i, j] = true
                    isinside[i, j] = true
                end
            end
        end
    end
    
    for j=2:Nx-1
        for i=2:Ny-1
            if(isgrid[i, j])
                if(!isgrid[i-1, j] + !isgrid[i+1, j] + !isgrid[i, j-1] + !isgrid[i, j+1] > 2)
                    isgrid[i, j] = false
                end
            end
        end
    end
    
    
    n = Array{Float64, 2}( reshape([], 0, 4) ) # defined as coordy coordx ny nx
    for j=2:Nx-1
        for i=2:Ny-1
            if(isgrid[i, j])
                if(!isgrid[i,j+1] & isgrid[i,j-1])
                    nx = 1
                elseif(!isgrid[i,j-1] & isgrid[i,j+1])
                    nx = -1
                else
                    nx = 0
                end
                if(!isgrid[i+1,j] & isgrid[i-1,j])
                    ny = 1
                elseif(!isgrid[i-1,j] & isgrid[i+1,j])
                    ny = -1
                else
                    ny = 0
                end
                if((nx^2 + ny^2) > 0)
                    isborder[i, j] = true
                    n = vcat(n, Array{Int64,2}([i j ny nx]))
                end
            end
        end
    end

    coord = reshape([],0,2)
    for j in 1:Nx
        for i in 1:Ny
            coord = [coord; reshape([yrange[i] xrange[j]], 1, 2)]
        end
    end

    isinside = isgrid .& .!isborder
    
    # create incidence matrix
    bus_id = Int64.(zeros(0,2))
    bus_coord = zeros(0,2)
    for j=1:Nx-1
        for i=1:Ny-1
            if(isgrid[i, j])
                bus_id = [bus_id; j i]
                bus_coord = [bus_coord; yrange[i] xrange[j]]
            end
        end
    end

    incidence_mat = Int64.(zeros(0,2))
    line_coord = zeros(0,4)
    for j=1:Nx-1
        for i=1:Ny-1
            if(isgrid[i, j] && isgrid[i, j+1])
                id1 = findfirst(all(bus_id .== [j i],dims=2))[1]
                id2 = findfirst(all(bus_id .== [j+1 i],dims=2))[1]
                incidence_mat = [incidence_mat; id1 id2]
                line_coord = [line_coord; [yrange[i] xrange[j] yrange[i] xrange[j+1]]]
            end
            if(isgrid[i, j] && isgrid[i+1, j])
                id1 = findfirst(all(bus_id .== [j i],dims=2))[1]
                id2 = findfirst(all(bus_id .== [j i+1],dims=2))[1]
                incidence_mat = [incidence_mat; id1 id2]
                line_coord = [line_coord; [yrange[i] xrange[j] yrange[i+1] xrange[j]]]
            end
        end
    end
    
return mesh = Mesh(
        Nx,
        Ny,
        Float64.(bus_coord[:,2:-1:1]),
        Float64.(line_coord),
        incidence_mat,
        vec(isgrid),
        yrange,
        xrange,
        dx,
    )
end


function update_model!(
    contmod::ContModel,
    dm::DiscModel,
)
    # Update the dynamical parameters of the model from a discrete model
    N = size(contmod.mesh.coord, 1)

    m = zeros(N)
    d = zeros(N)
    pl = zeros(N)
    pg = zeros(N)

    Threads.@threads for g in 1:dm.Ngen
        if(dm.p_gen[g] > 0.0)
            k = argmin((contmod.mesh.coord[:, 1] .- dm.coord[dm.id_gen[g],1]).^2 +
                (contmod.mesh.coord[:, 2] .- dm.coord[dm.id_gen[g],2]).^2)
            m[k] += dm.m_gen[g]
            d[k] += dm.d_gen[g]
            pg[k] += dm.p_gen[g]
        end
    end

    Threads.@threads for l in 1:dm.Nbus
        k = argmin((contmod.mesh.coord[:, 1] .- dm.coord[l,1]).^2 +
            (contmod.mesh.coord[:, 2] .- dm.coord[l,2]).^2)
        d[k] += dm.d_load[l]
        pl[k] += dm.p_load[l]
    end
    @time begin
        m = heat_diff(contmod.mesh, m, Niter = contmod.Niter, tau = contmod.tau)
        d = heat_diff(contmod.mesh, d, Niter = contmod.Niter, tau = contmod.tau)
        pl = heat_diff(contmod.mesh, pl, Niter = contmod.Niter, tau = contmod.tau)
        pg = heat_diff(contmod.mesh, pg, Niter = contmod.Niter, tau = contmod.tau)
    end

    # asign minimal values to the quantities
    m .= max.(m, contmod.min_factor * sum(m) / length(m))
    d .= max.(d, contmod.min_factor * sum(d) / length(m))
    
    # ensure that the integrals of the different quantities over
    # the medium is equivalent to their sum over the discrete model
    m = (sum(dm.m_gen) / sum(m) / contmod.mesh.dx^2) .* m 
    d = ((sum(dm.d_gen) + sum(dm.d_load)) / sum(d) / contmod.mesh.dx^2) .* d
    pl = (sum(dm.p_load) / sum(pl) / contmod.mesh.dx^2) .* pl
    pg = (sum(dm.p_gen) / sum(pg) / contmod.mesh.dx^2) .* pg

    # ensure that the integral of the power injection is 0, i.e
    # that generation matches the demand.
    p = pg - pl
    p = p .- sum(p) / length(p) 

    minv = m.^(-1)    
    # Update the paramters
    contmod.minv = minv
    contmod.gamma = d.* minv
    contmod.m = m
    contmod.d = d
    contmod.p = p
end


function update_susceptance!(cm::ContModel, b::Vector{Float64})
    n = size(cm.mesh.inc_mat, 1)
    incMat = sparse([cm.mesh.inc_mat[:,1]; cm.mesh.inc_mat[:,2]], [1:n; 1:n], [-ones(n); ones(n)])
    cm.b = b
    cm.xi = -incMat * (b .* incMat')
end


function update_params!(
    contmod::ContModel
)
    # update xi, minv and gamma
    minv = contmod.m.^(-1)
    gamma = contmod.d.* minv
    contmod.minv = minv
    contmod.gamma = contmod.d .* minv
    n = size(contmod.mesh.inc_mat, 1)
    temp = sparse([contmod.mesh.inc_mat[:,1]; contmod.mesh.inc_mat[:,2]], [1:n; 1:n], [-ones(n); ones(n)])
    contmod.xi = - temp * (contmod.b .* temp')
end
