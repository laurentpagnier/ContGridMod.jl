function inPolygon_new(
    point::Array{Tuple{Float64, Float64}, 1},
    poly::Array{Tuple{Float64, Float64}, 1},
)
    N = length(poly)
    Np = length(point)
    b = falses(Np)
    Threads.@threads for k in 1:Np
        j = N
        for i = 1:N
            if (
                ((poly[i][2] < point[k][2]) & (poly[j][2] >= point[k][2])) |
                ((poly[j][2] < point[k][2]) & (poly[i][2] >= point[k][2]))
            )
                if (
                    poly[i][1] +
                    (point[k][2] - poly[i][2]) / (poly[j][2] - poly[i][2]) *
                    (poly[j][1] - poly[i][1]) < point[k][1]
                )
                    b[k] = !b[k]
                end
            end
            j = i
        end
    end
    return b
end



function vectorize(
    isinside::BitMatrix,
    isborder::BitMatrix,
    n::Array{Float64, 2},
    bx::Array{Float64, 2},
    by::Array{Float64, 2},
    p::Array{Float64, 2},
    m::Array{Float64, 2},
    d::Array{Float64, 2},
)
    # flatten the m, d, p, bx and by matrices
    isinsideflat = vec(isinside)
    isgridflat = vec(isinside .| isborder)
    bxflat = vec(bx)
    byflat = vec(by)
    pflat = vec(p)
    mflat = vec(m)
    dflat = vec(d)
    Ny = size(bx,1)
    Nx = size(bx,2)

    id1 = Array{Int64,1}()
    id2 = Array{Int64,1}()
    v = Array{Float64,1}()
    for k in 1:Nx*Ny
        if(isinsideflat[k])
            append!(id1, k)
            append!(id2, k)
            append!(v, - (bxflat[k] +
                bxflat[k-Ny] +
                byflat[k] +
                byflat[k-1])
                )
            append!(id1, k)
            append!(id2, k-Ny)
            append!(v, bxflat[k-Ny])
            append!(id1, k)
            append!(id2, k+Ny)
            append!(v, bxflat[k])
            append!(id1, k)
            append!(id2, k-1)
            append!(v, byflat[k-1])
            append!(id1, k)
            append!(id2, k+1)
            append!(v, byflat[k])            
        end
    end
    
    for id in 1:size(n, 1)
        k = (Int64(n[id, 2]) - 1) * Ny + Int64(n[id, 1])
        ny = n[id, 3]
        nx = n[id, 4] 
        etamx = 1 - nx/2-nx^2/2
        etapx = 1 + nx/2-nx^2/2
        etamy = 1 - ny/2-ny^2/2
        etapy = 1 + ny/2-ny^2/2
        append!(id1, k)
        append!(id2, k)
        append!(v, - (etamx * bxflat[k] +
            etapx * bxflat[k-Ny] +
            etamy * byflat[k] +
            etapy * byflat[k-1]))
        append!(id1, k)
        append!(id2, k-Ny)
        append!(v, etapx * bxflat[k-Ny])
        append!(id1, k)
        append!(id2, k+Ny)
        append!(v, etamx * bxflat[k])
        append!(id1, k)
        append!(id2, k-1)
        append!(v, etapy * byflat[k-1])   
        append!(id1, k)
        append!(id2, k+1)
        append!(v, etamy * byflat[k])
    end
    
    
    xi = sparse(id1, id2, v, Ny * Nx, Ny * Nx)
    minvflat = mflat.^(-1)
    minvflat = minvflat[isgridflat]
    gammaflat = dflat[isgridflat] .* minvflat
    xi = SparseMatrixCSC{Float64, Int64}(xi[isgridflat,isgridflat])
    return isinsideflat, pflat, minvflat, gammaflat, xi
end



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


function local_disturbance_old(
    isinside::BitMatrix,
    xrange::Array{Float64,1},
    yrange::Array{Float64,1},
    location::Array{Float64,1},
    dP::Float64,
    sigma::Float64
)
    dp = zeros(Ny, Nx)
    for i = 1:Nx
        for j = 1:Ny
            if(isinside[j,i])
                dp[j, i] = exp(-((xrange[i] - location[1])^2 +(yrange[j] - location[2])^2)/2/sigma^2)
            end
        end
    end
    return dP .* dp ./ trapz((yrange, xrange), dp)
end



function ctr_movie_old(
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
        heatmap(cont_value[:,:,t], clim = (-10., 10.))
    end
end


function time_plot_old(
    time::Array{Float64, 1},
    cont_value::Array{Float64, 3}, # timeseries from the continous model
    coord::Array{Float64, 2}; # locations from where we wantg to fetch data
    borders::Array{Array{Float64,2},1} = Array{Float64,2}[],
    xlabel::String = String("\$t\\;[s]\$"),
    ylabel::String = String("\$\\omega \$"),
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
        idend = findall(time .< tend)[end]
    else
        idend = length(time) 
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





function get_params(
    isinside::BitMatrix,
    isborder::BitMatrix,
    n_vector::Array{Float64, 2},
    dx::Float64,
    yrange::Array{Float64, 1},
    xrange::Array{Float64, 1},
    scaling_factor::Float64,
    dataname::String,
    savename::String;
    sigma::Float64 = 100.0,
    min_factor::Float64 = 0.1,
    bmin::Float64 = 1.0,
    patch::Float64 = 1.0 # a temporary patch so that the discrete and continuous
                  #stables solutions more or less coincide
)
    gen, dem, bline, idb, idgen, coord, mg, dg, dl, th = load_discrete_model(dataname::String, scaling_factor)
    
    isgrid = BitArray(isinside + isborder)
    
    Nx = length(xrange)
    Ny = length(yrange)
    
    #bx = bmin * ones(Ny, Nx)
    #by = bmin * ones(Ny, Nx)
    bx = zeros(Ny, Nx)
    by = zeros(Ny, Nx)    
    m = zeros(Ny, Nx)
    d = zeros(Ny, Nx)
    pl = zeros(Ny, Nx)
    pg = zeros(Ny, Nx)
    
    for g in 1:length(idgen)
        Threads.@threads for i=2:Ny-1
            Threads.@threads for j=2:Nx-1
                idg1 = coord[idgen[g],1]
                idg2 = coord[idgen[g],2]
                n = norm_dist([idg1, idg2], [xrange[j], yrange[i]], sigma)
                m[i,j] += mg[g] * n
                d[i,j] += dg[g] * n
                pg[i,j] += gen[g] * n
            end
        end
    end

    for l in 1:length(dl)
        Threads.@threads for i=2:Ny-1
            Threads.@threads for j=2:Nx-1
                n = norm_dist(vec(coord[l,:]), [xrange[j], yrange[i]], sigma)
                d[i,j] += dl[l] * n
                pl[i,j] += dem[l] * n
            end
        end
    end
    
    # a point x on line defined by two points a=(x1,y1) and b=(x2,y2)
    # can be written as x = alpha*(dx,dy) + a, with dx-x2-x1 and dy=y2-y1
    # the line prependicular to the first one and passing by c=(x3,y3) is
    # given by y = beta*(-dy,dx) + c. The point where the two line cross
    # each other is given by y = x. Solving for beta and alpha, the
    # distance separating the point c and the line ab is beta*sqrt(dy^2+dx^2) 
    Threads.@threads for l in 1:size(idb, 1)
        x2 = coord[idb[l, 2], 1]
        x1 = coord[idb[l, 1], 1]
        y2 = coord[idb[l, 2], 2]
        y1 = coord[idb[l, 1], 2] 
        dx_l = x2 - x1
        dy_l = y2 - y1
        phi = atan(dy_l, dx_l)
        ds2 = dx_l^2 + dy_l^2
        if (dx_l != 0 && dy_l != 0) # if not a transformer
            Threads.@threads for i = 1:Ny-1
                Threads.@threads for j = 1:Nx-1
                    if (isgrid[i,j]) # if the line bx(i,j) is the grid                     
                        # we set the coordinates of the lines to be the ones of their middle points
                        bx_x = (xrange[j] + xrange[j+1]) / 2
                        bx_y = yrange[i]
                        beta = (dx_l * (y1 - bx_y) + dy_l * (bx_x - x1)) / ds2
                        alpha = (dx_l * (bx_x - x1) + dy_l * (bx_y - y1)) / ds2
                        # if the projection is in the segment ab
                        if ((0 < alpha) & (alpha < 1))
                            dist = abs(beta) * sqrt(dx_l^2 + dy_l^2)
                            # if close enough to the line 
                            f = exp(-dist^2 / 2 / sigma^2) / sigma^2 / (2 * pi)
                            bx[i, j] += bline[l] * abs(cos(phi)) * f * dx^2 * patch
                            # else close enough to an end of the line 
                        else
                            dist = min(
                                (bx_y - y1)^2 + (bx_x - x1)^2,
                                (bx_y - y2)^2 + (bx_x - x2)^2,
                            )
                            f = exp(-dist / 2 / sigma^2) / sigma^2 / (2 * pi)
                            bx[i, j] += bline[l] * abs(cos(phi)) * f * dx^2 * patch
                        end
                    end
                    
                    if(isgrid[i, j]) # if the line by(i,j) is the grid 
                        # we set the coordinates of the lines to be the ones of their middle points
                        by_x = xrange[j]
                        by_y = (yrange[i] + yrange[i+1]) / 2      
                        beta = (dx_l*(y1 - by_y) + dy_l*(by_x - x1)) / ds2
                        alpha = (dx_l*(by_x - x1) + dy_l*(by_y - y1)) / ds2
                        # if the projection is in the segment ab
                        if((0 < alpha) & (alpha < 1))
                            dist = abs(beta) * sqrt(ds2) 
                            # if close enough to the line 
                            f = exp(-dist^2 / 2 / sigma^2) / sigma^2 / (2*pi)
                            by[i,j] += bline[l] * abs(sin(phi)) * f * dx^2 * patch 
                        # else close enough to an end of the line 
                        else
                            dist = min((by_y - y1)^2 + (by_x - x1)^2,
                                (by_y - y2)^2 + (by_x - x2)^2)
                            f = exp(-dist / 2 / sigma^2) / sigma^2 / (2*pi)
                            by[i,j] += bline[l] * abs(sin(phi)) * f * dx^2 * patch
                        end
                    end

                end
            end
        end
    end

    #  assign minimal values to ensure numerical stability
    #=
    Threads.@threads for i=2:Ny-1
        Threads.@threads for j=2:Nx-1
            if(isgrid[i, j-1] & isgrid[i, j] & !(isborder[i, j-1] & isborder[i, j])) # if the line bx(i,j) is the grid
                bx[i, j] = max(bx[i, j], bmin)
            end
            if(isgrid[i, j] & isgrid[i, j] & !(isborder[i, j-1] & isborder[i, j])) # if the line by(i,j) is the grid 
                by[i, j] = max(by[i, j], bmin)
            end
        end
    end
    =#
    
    #=
    Threads.@threads for i=2:Ny-1
        Threads.@threads for j=2:Nx-1
            if(isgrid[i, j-1] & isgrid[i, j]) # if the line bx(i,j) is the grid
                bx[i, j] = max(bx[i, j], bmin)
            end
            if(isgrid[i, j] & isgrid[i, j]) # if the line by(i,j) is the grid 
                by[i, j] = max(by[i, j], bmin)
            end
        end
    end
    =#
    
    bx[isgrid] .= max.(bx[isgrid], bmin)
    by[isgrid] .= max.(by[isgrid], bmin)
    
    # Threads.@threads for k in 1:size(n_vector,1)
    #     i = Int64(n_vector[k,1])
    #     j = Int64(n_vector[k,2])
                
    #     nx = n_vector[k,4]
    #     ny = n_vector[k,3]
        
    #     by[i-1, j] = ny * (1 + ny) / 2  * by[i-2, j]
    #     by[i, j] = ny * (ny - 1) / 2 * by[i+1, j]
    #     bx[i, j-1] = nx * (1 + nx) / 2 * by[i, j-2]
    #     bx[i, j] = nx * (nx - 1) / 2 * by[i, j+1]
    # end
    # due to how the boundary is treated in the code, interia, damping or
    # power injection on boundary won't be taken into account
    m[.!isgrid] .= 0
    d[.!isgrid] .= 0
    #pl[.!isgrid] .= 0
    #pg[.!isgrid] .= 0
    pl[.!isinside] .= 0
    pg[.!isinside] .= 0
    
    # asign minimal values to the quantities
    m[isgrid] .= max.(m[isgrid], min_factor * maximum(m))
    d[isgrid] .= max.(d[isgrid], min_factor * maximum(d))
    #pl[isgrid] .= max.(pl[isgrid], min_factor * maximum(pl))
    #pg[isgrid] .= max.(pg[isgrid], min_factor * maximum(pg))
    pl[isinside] .= max.(pl[isinside], min_factor * maximum(pl))
    pg[isinside] .= max.(pg[isinside], min_factor * maximum(pg))
    
    # ensure that the integrals of the different quantities over
    # the medium is equivalent to their sum over the discrete model

    m = sum(mg) .* m ./ trapz((yrange, xrange), m)
    d = (sum(dg) + sum(dl)) .* d ./ trapz((yrange, xrange), d)
    pl = sum(dem) .* pl ./ trapz((yrange, xrange), pl)
    pg = sum(gen) .* pg ./ trapz((yrange, xrange), pg)
    

    # save the quantities
    fid = h5open(savename, "w")
    write(fid, "bx", bx)
    write(fid, "by", by)
    write(fid, "m", m)
    write(fid, "d", d)
    write(fid, "pl", pl)
    write(fid, "pg", pg)
    write(fid, "xrange", xrange)
    write(fid, "yrange", yrange)
    close(fid)

    # ensure that the integral of the power injection is 0, i.e
    # that generation match the demand.
    p = pg - pl 
    #p[isgrid] = p[isgrid] .- sum(p[isgrid]) / sum(isgrid)
    p[isinside] = p[isinside] .- sum(p[isinside]) / sum(isinside)
    
    return bx, by, p, m, d
end








function compute_stable_sol(
    isinside::BitMatrix,
    n::Array{Float64, 2},
    bx::Array{Float64, 2},
    by::Array{Float64, 2},
    p::Array{Float64, 2};
    interval::Int64 = 1000,
    Niter::Int64 = 15000,
    tol::Float64 = 1e-7
)

    temp = findall(isinside)
    idin = Int64.(zeros(size(temp, 1), 2))
    for i in 1:size(temp,1)
        idin[i,:] = [temp[i][1], temp[i][2]]
    end

    th = zeros(Ny, Nx)
    ths = zeros(Ny, Nx, Int64(floor(Niter/interval)))
    @time begin
        for t in 1:Niter
            if(mod(t, interval) == 0)
                temp = copy(th)
            end
            
            Threads.@threads for k in 1:size(idin, 1)
                i = idin[k, 1]
                j = idin[k, 2]
                bij = by[i-1, j] + by[i, j] + bx[i, j-1] + bx[i, j]
                th[i,j] = (
                    by[i, j] * th[i+1, j] +
                    by[i-1, j] * th[i-1, j] + 
                    bx[i, j] * th[i, j+1] +
                    bx[i, j-1] * th[i, j-1] +
                    dx^2 * p[i, j]
                    ) / bij
            end
            
            Threads.@threads for k in 1:size(n,1)
                i = Int64(n[k,1])
                j = Int64(n[k,2])
                nx = n[k, 4]
                ny = n[k, 3]
                
                bij = (1.0 + ny) * by[i-1, j] +
                    (1.0 - ny) * by[i, j] +
                    (1.0 + nx) * bx[i, j-1] +
                    (1.0 - nx) * bx[i, j]

                th[i, j] = (
                    (1.0 + ny) * by[i-1, j] * th[i-1, j] +
                    (1.0 - ny) * by[i, j] * th[i+1, j] + 
                    (1.0 + nx) * bx[i, j-1] * th[i, j-1] +
                    (1.0 - nx) * bx[i, j] * th[i, j+1] +
                    dx^2 * p[i, j]
                    ) / bij
            end
            
            if(mod(t, interval) == 0)
                println( [t maximum(abs.(th - temp))] )
                ths[:, :, Int64(t/interval)] = th
                if( maximum(abs.(th - temp)) < tol )
                    break
                end
            end
        end
    end
    return th, ths
end



function get_grid_fast(
    border::Array{Float64, 2},
    dx::Float64
)
    xlim = [minimum(border[:, 1])-2*dx maximum(border[:, 1])+2*dx]
    ylim = [minimum(border[:, 2])-2*dx maximum(border[:, 2])+2*dx]
    println(xlim)
    println(ylim)
    xrange = collect(xlim[1]:dx:xlim[2])
    yrange = collect(ylim[1]:dx:ylim[2])
    Nx = length(xrange)
    Ny = length(yrange)
    
    #p = [(yrange[i], xrange[j])  for j in 1:Nx for i in 1:Ny]

    p = Float64.(reshape([],0,2))
    for j in 1:Nx
        for i in 1:Ny
            #if(isgrid[i,j])
            #    coord = [coord; reshape([yrange[i] xrange[j]], 1, 2)]
            #end
            p = [p; reshape([yrange[i] xrange[j]], 1, 2)]
        end
    end
    p = Float64.(p)
    
    isborder = falses(Ny * Nx)
    #isinside = falses(Ny * Nx)
    isgrid = falses(Ny * Nx)
    #border2 = Array{Tuple{Float64, Float64}, 1}()
    #for i in 1:size(border,1)
    #    push!(border2, (border[i,2], border[i,1]))
    #end
    
    # check if inside the shape
    
    #rintln(minimum(border))
    #println(minimum(p))
    #println(maximum(border))
    #println(maximum(p))
    isinside = inPolygon(p, border)
    println(size(isinside))
    scatter(p[isinside,1],p[isinside,2])
    # check if alone 
    for k in Ny+1:Ny*(Nx-1)-1
        if(isinside[k])
            if(isinside[k-1] + isinside[k+1] + isinside[k-Ny] + isinside[k+Ny] == 0)
                isinside[k] = false
            end
        end
    end
    
    
    # add a boundary layer encompassing the inside
    #id = Array{Int64, 2}( reshape([], 0, 2) ) # defined as coordy coordx ny nx
    for k in Ny+1:Ny*(Nx-1)-1
        if(isinside[k] == false)
            nx = isinside[k-Ny] - isinside[k+Ny]
            ny = isinside[k-1] - isinside[k+1]
            if((nx^2 + ny^2) > 0)
                isgrid[k] = true
            end
        end
    end

    isgrid = isgrid .| isinside
    # check if inside and considered as outside. This check isnt able
    # to find hole that are bigger than a single point. Som it should
    # be improved in the future
    for k in Ny+1:Ny*(Nx-1)-1
        if(!isgrid[k])
            if(isgrid[k-1] + isgrid[k+1] + isgrid[k-Ny] + isgrid[k+Ny] == 4)
                isgrid[k] = true
                isinside[k] = true
            end
        end
    end
    
    for k in Ny+1:Ny*(Nx-1)-1
        if(isgrid[k])
            if(!isgrid[k-1] + !isgrid[k+1] + !isgrid[k+Ny] + !isgrid[k+Ny] > 2)
                isgrid[k] = false
            end
        end
    end
  
    n = Array{Float64, 2}( reshape([], 0, 4) ) # defined as coordy coordx ny nx
    for k in Ny+1:Ny*(Nx-1)-1
        if(isgrid[k])
            if(!isgrid[k+Ny] & isgrid[k-Ny])
                nx = 1
            elseif(!isgrid[k-Ny] & isgrid[k+Ny])
                nx = -1
            else
                nx = 0
            end
            if(!isgrid[k+1] & isgrid[k-1])
                ny = 1
            elseif(!isgrid[k-1] & isgrid[k+1])
                ny = -1
            else
                ny = 0
            end
            if((nx^2 + ny^2) > 0)
                isborder[k] = true
                i = Int64(floor(k/Nx)) + 1
                j = Ny*(i-1) 
                n = vcat(n, Array{Int64,2}([i j ny nx]))
            end
        end

    end
    

    isinside = isgrid .& .!isborder
    return Nx, Ny, p, isinside, isborder, n, isgrid
end






function get_grid(border::Array{Float64,2}, dx::Float64)

    xlim = [minimum(border[:,1])-2*dx maximum(border[:,1])+2*dx]
    ylim = [minimum(border[:,2])-2*dx maximum(border[:,2])+2*dx]
    xrange = collect(xlim[1]:dx:xlim[2])
    yrange = collect(ylim[1]:dx:ylim[2])

    Nx = length(xrange)
    Ny = length(yrange)

    isborder = Bool.(zeros(Ny,Nx))
    isinside = Bool.(zeros(Ny,Nx))
    isgrid = Bool.(zeros(Ny,Nx))

    for j=1:Nx
        for i=1:Ny
            p = [xrange[j] yrange[i]]
            tmp = inPolygon(p,border)[1]
            isinside[i,j] = tmp
            isgrid[i, j] = tmp
        end
    end
    # check if alone 
    for j=2:Nx-1
        for i=2:Ny-1
            if(isinside[i, j])
                if(isinside[i-1, j] + isinside[i+1, j] + isinside[i, j-1] + isinside[i, j+1] == 0)
                    isinside[i, j] = false
                    isgrid[i, j] = false
                end
            end
        end
    end
    # Check if nodes are "chained"
    while true
        rem = Array{Float64, 2}(reshape([], 0, 2))
        for j=2:Nx-1
            for i=2:Ny-1
                if(isinside[i, j] && isinside[i+1, j] + isinside[i-1, j] + isinside[i, j+1] + isinside[i, j-1] < 2)
                    rem = vcat(rem, [i j])
                end
            end
        end
        if size(rem, 1) > 0
            for i=1:size(rem, 1)
                isinside[Int64(rem[i, 1]), Int64(rem[i, 2])] = false
                isgrid[Int64(rem[i, 1]), Int64(rem[i, 2])] = false
            end
        else
            break
        end
    end
    
    n = Array{Float64,2}(reshape([],0,4)) # defined as coordy coordx ny nx
    for j=2:Nx-1
       for i=2:Ny-1
            nx = isinside[i, j-1] - isinside[i, j+1]
            ny = isinside[i-1, j] - isinside[i+1, j]
            if((nx^2 + ny^2) > 0 && isinside[i,j] == false)
                isborder[i, j] = true
                isgrid[i, j] = true
                n = vcat(n,Array{Int64,2}([i j ny nx]))
            end
       end
    end
    
    return Nx, Ny, xrange, yrange, isinside, isborder, n, isgrid
end






function vectorize(
    isinside::BitMatrix,
    isborder::BitMatrix,
    n::Array{Float64, 2},
    bx::Array{Float64, 2},
    by::Array{Float64, 2},
    p::Array{Float64, 2},
    m::Array{Float64, 2},
    d::Array{Float64, 2},
)
    # flatten the m, d, p, bx and by matrices
    isinsideflat = vec(isinside)
    isgridflat = vec(isinside .| isborder)
    bxflat = vec(bx)
    byflat = vec(by)
    pflat = vec(p)
    mflat = vec(m)
    dflat = vec(d)
    Ny = size(bx,1)
    Nx = size(bx,2)
    xi = sparse([], [], Float64.([]), Ny * Nx, Ny * Nx)

    #Threads.@threads for i in 1:Nx*Ny
    for k in 1:Nx*Ny
        if(isinsideflat[k])
            xi[k, k] = - (bxflat[k] +
                bxflat[k-Ny] +
                byflat[k] +
                byflat[k-1])
            xi[k, k-Ny] = bxflat[k-Ny]
            xi[k, k+Ny] = bxflat[k]
            xi[k, k-1] = byflat[k-1]
            xi[k, k+1] = byflat[k]
        end
    end
    
    #Threads.@threads for k in 1:size(n, 1)
    for id in 1:size(n, 1)
        k = (Int64(n[id, 2]) - 1) * Ny + Int64(n[id, 1])
        ny = n[id, 3]
        nx = n[id, 4] 
        
        xi[k, k] = - ((1.0 - nx) * bxflat[k] +
            (1.0 + nx) * bxflat[k-Ny] +
            (1.0 - ny) * byflat[k] +
            (1.0 + ny) * byflat[k-1])
        xi[k, k-Ny] = (1.0 + nx) * bxflat[k-Ny]
        xi[k, k+Ny] = (1.0 - nx) * bxflat[k]
        xi[k, k-1] = (1.0 + ny) * byflat[k-1]
        xi[k, k+1] = (1.0 - ny) * byflat[k]
    end  
    
    minvflat = mflat.^(-1)
    minvflat = minvflat[isgridflat]
    gammaflat = dflat[isgridflat] .* minvflat
    xi = SparseMatrixCSC{Float64, Int64}(xi[isgridflat,isgridflat])
    return isinsideflat, pflat, minvflat, gammaflat, xi
end




function find_time_step(isin::BitMatrix, m::Array{Float64,2}, d::Array{Float64,2}, p::Array{Float64,2},
    bx::Array{Float64,2}, by::Array{Float64,2}, dx::Float64, alpha=0.1)
    # !!!!!!! NOT FUNCTIONAL STILL SOME WORK TO DO
    Nx = size(m,2)
    Ny = size(m,1)
    bij = zeros(size(m))
    for i = 2:Ny-2
        for j = 2:Nx-2
            bij[i,j] = bx[i,j] + bx[i,j+1] + by[i-1,j] + by[i,j]
        end
    end
    gamma = d ./ m
    println(alpha*dx^2*minimum(m[isin] ./ bij[isin]))
    println(alpha*minimum(d[isin] ./ abs.(p[isin])))
end


function vectorize_ref(
    isinside::BitMatrix,
    isborder::BitMatrix,
    n::Array{Float64, 2},
    bx::Array{Float64, 2},
    by::Array{Float64, 2},
    p::Array{Float64, 2},
    m::Array{Float64, 2},
    d::Array{Float64, 2},
)
    # flatten the m, d, p, bx and by matrices
    isinsideflat = vec(isinside)
    isgridflat = vec(isinside .| isborder)
    bxflat = vec(bx)
    byflat = vec(by)
    pflat = vec(p)
    mflat = vec(m)
    dflat = vec(d)
    Ny = size(bx,1)
    Nx = size(bx,2)
    xi = sparse([], [], Float64.([]), Ny * Nx, Ny * Nx)

    #Threads.@threads for i in 1:Nx*Ny
    for k in 1:Nx*Ny
        if(isinsideflat[k])
            xi[k, k] = - bxflat[k] - bxflat[k-Ny] -
                byflat[k] - byflat[k-1]
            xi[k, k-Ny] = bxflat[k-Ny]
            xi[k, k+Ny] = bxflat[k]
            xi[k, k-1] = byflat[k-1]
            xi[k, k+1] = byflat[k]
        end
    end
    
    #Threads.@threads for k in 1:size(n, 1)
    for id in 1:size(n, 1)
        k = (Int64(n[id, 2]) - 1) * Ny + Int64(n[id, 1])
        ny = n[id, 3]
        nx = n[id, 4] 
        
        xi[k, k] = - (1.0 - nx) * bxflat[k] - (1.0 + nx) * bxflat[k-Ny] -
            (1.0 - ny) * byflat[k] - (1.0 + ny) * byflat[k-1]
        xi[k, k-Ny] = (1.0 + nx) * bxflat[k-Ny]
        xi[k, k+Ny] = (1.0 - nx) * bxflat[k]
        xi[k, k-1] = (1.0 + ny) * byflat[k-1]
        xi[k, k+1] = (1.0 - ny) * byflat[k]
    end  
    
    minvflat = mflat.^(-1)
    minvflat = minvflat[isgridflat]
    gammaflat = dflat[isgridflat] .* minvflat
    xi = SparseMatrixCSC{Float64, Int64}(xi[isgridflat,isgridflat])
    return isinsideflat, pflat, minvflat, gammaflat, xi
end



function get_cont_values(
    isinside::BitMatrix,
    yrange::Array{Float64, 1},
    xrange::Array{Float64, 1},
    disc_coord::Array{Float64, 2},
    disc_values::Array{Float64, 1};
    Niter = 1000
)
    
    idin = findall(isinside)
    cont_coord = zeros(length(idin), 2)
    for i = 1:length(idin)
        cont_coord[i,:] = [xrange[idin[i][2]], yrange[idin[i][1]]]
    end
    
    cont_values_vec = zeros(length(idin))
    for i = 1:size(idin, 1)
        dx = disc_coord[:,1] .- cont_coord[i,1] 
        dy = disc_coord[:,2] .- cont_coord[i,2] 
        id = argmin(dx.^2 + dy.^2)
        cont_values_vec[i] = disc_values[id]
    end
    
    cont_values = zeros(Ny, Nx)
    cont_values[idin] = cont_values_vec
    
    new_cont_values = zeros(Ny, Nx)
    tau = 0.001
    interval = 100

    @time begin
        for k in 1:Niter
            Threads.@threads for i in 2:Ny-1
                Threads.@threads for j in 2:Nx-1
                    if(isinside[i,j])
                        new_cont_values[i,j] = (1.0 - 4.0 * tau) * cont_values[i,j] + tau * (cont_values[i+1,j] +
                            cont_values[i-1,j] + cont_values[i,j+1] + cont_values[i,j-1])
                    end
                end
            end
            
            Threads.@threads for k in 1:size(n,1)
                i = Int64(n[k,1])
                j = Int64(n[k,2])
                nx = n[k,4]
                ny = n[k,3]
                if(nx == 1)
                    new_cont_values[i,j] = new_cont_values[i,j-2]
                elseif(nx == -1)
                    new_cont_values[i,j] = new_cont_values[i,j+2]
                end
                if(ny == 1)
                    new_cont_values[i,j] = new_cont_values[i-2,j]
                elseif(ny == -1)
                    new_cont_values[i,j] = new_cont_values[i+2,j]
                end
            end
            
            cont_values = copy(new_cont_values)
            
            if(mod(k,interval) == 0)
                println(k)
            end
        end
    end
    return cont_values
end
