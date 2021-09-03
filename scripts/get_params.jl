function get_params(
    isinside::BitMatrix,
    isborder::BitMatrix,
    dx::Float64,
    yrange::Array{Float64,1},
    xrange::Array{Float64,1},
    dataname::String,
    savename::String;
    sigma::Float64 = 50.0,
    min_factor::Float64 = 0.1,
    bmin::Float64 = 1.0,
    patch::Float64 = 1.0 # a temporary patch so that the discrete and continuous
                  #stables solutions more or less coincide
)
    gen, dem, bline, idb, idgen, coord, mg, dg, dl, th = load_discrete_model(dataname::String)
    
    isgrid = isinside .| isborder
    
    Nx = length(xrange)
    Ny = length(yrange)
    
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
    isnew = true
    for l in 1:size(idb,1)
        x2 = coord[idb[l,2], 1]
        x1 = coord[idb[l,1], 1]
        y2 = coord[idb[l,2], 2]
        y1 = coord[idb[l,1], 2] 
        dx_l = x2 - x1
        dy_l = y2 - y1
        ds2 = (dy_l^2 + dx_l^2)
        phi = atan(dy_l, dx_l)
        if(dx_l != 0 && dy_l != 0) # if not a transformer
            Threads.@threads for i = 2:Ny-1
                Threads.@threads for j = 2:Nx-1
                    if((isgrid[i,j-1] & isgrid[i,j]) & !(isborder[i,j-1] & isborder[i,j])) # if the line bx(i,j) is the grid                     
                        # we set the coordinates of the lines to be the ones of their middle points
                        bx_x = (xrange[j-1] + xrange[j]) / 2
                        bx_y = yrange[i]
                        beta = (dx_l*(y1 - bx_y) + dy_l*(bx_x- x1)) / ds2
                        alpha = (dx_l*(bx_x - x1) + dy_l*(bx_y - y1)) / ds2
                        # if the projection is in the segment ab
                        if((0 < alpha) & (alpha < 1))
                            dist = abs(beta) * sqrt(ds2) 
                            # if close enough to the line 
                            f = exp(-dist^2/2/sigma^2)/sigma^2 / (2*pi)
                            bx[i,j] += bline[l] * abs(cos(phi)) * f * dx^2 * patch 
                        # else close enough to an end of the line 
                        else
                            dist = min((bx_y - y1)^2 + (bx_x - x1)^2,
                            (bx_y - y2)^2 + (bx_x - x2)^2)
                            f = exp(-dist / 2 / sigma^2) / sigma^2 / (2*pi)
                            bx[i,j] += bline[l] * abs(cos(phi)) * f * dx^2 * patch
                        end
                    end
                    
                    if((isgrid[i,j] & isgrid[i+1,j]) & !(isborder[i,j] & isborder[i+1,j])) # if the line by(i,j) is the grid 
                        # we set the coordinates of the lines to be the ones of their middle points
                        by_x = xrange[j]
                        by_y = (yrange[i] + yrange[i+1]) / 2      
                        beta = (dx_l*(y1 - by_y) + dy_l*(by_x - x1)) / ds2
                        alpha = (dx_l*(by_x - x1) + dy_l*(by_y - y1)) / ds2
                        # if the projection is in the segment ab
                        if((0 < alpha) & (alpha < 1))
                            dist = abs(beta) * sqrt(ds2) 
                            # if close enough to the line 
                            f = exp(-dist^2/2/sigma^2) / sigma^2 / (2*pi)
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
    Threads.@threads for i=2:Ny-1
        Threads.@threads for j=2:Nx-1
            if(isgrid[i,j-1] & isgrid[i,j] & !(isborder[i,j-1] & isborder[i,j])) # if the line bx(i,j) is the grid
                bx[i,j] = max(bx[i,j], bmin)
            end
            if(isgrid[i,j] & isgrid[i+1,j] & !(isborder[i,j] & isborder[i+1,j])) # if the line by(i,j) is the grid 
                by[i,j] = max(by[i,j], bmin)
            end
        end
    end
    
    println(minimum(m[isinside]))
    m[isinside] .= max.(m[isinside], min_factor * sum(m[isinside]) / sum(isinside))
    d[isinside] .= max.(d[isinside], min_factor * sum(d[isinside]) / sum(isinside))
    println(minimum(m[isinside]))
    
    # due to how the boundary is treated in the code, interia, damping or
    # power injection on boundary won't be taken into account
    m[.!isinside] .= 0
    d[.!isinside] .= 0
    pl[.!isinside] .= 0
    pg[.!isinside] .= 0
    
    # asign minimal values to the quantities
    m[isinside] .= max.(m[isinside], min_factor * maximum(m))
    d[isinside] .= max.(d[isinside], min_factor * maximum(d))
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
    p[isinside] = p[isinside] .- sum(p[isinside]) / sum(isinside)

    return bx, by, p, m, d
end


function norm_dist(a::Array{Float64,1}, b::Array{Float64,1}, sigma::Float64)
    return exp(-((a[1] - b[1])^2 + (a[2] - b[2])^2) / 2 / sigma^2) / (2*pi) / sigma^2
end


function get_params_diff(
    isinside::BitMatrix,
    isborder::BitMatrix,
    dx::Float64,
    yrange::Array{Float64,1},
    xrange::Array{Float64,1},
    dataname::String,
    savename::String;
    def_val::Float64 = 1E-2,
    dmax::Float64 = 100.0,
    Niter::Int64 = 5000,
    tau::Float64 = 0.001,
    patch::Float64 = 1.0,
    bmin::Float64 = 1.0)
    # uses heat equation to "diffuse" the discrete parameters over the lattice
    # for heat equation tau=kappa*dt/dx^2
    
    isgrid = isinside .| isborder
    println(patch)
    gen, dem, bline, idb, idgen, coord, mg, dg, dl, th = load_discrete_model(dataname)

    Nx = length(xrange)
    Ny = length(yrange)
    
    idin = findall(isinside)
    Nbus = length(idin)
    lat_coord = zeros(Nbus, 2)
    for i = 1:Nbus
        lat_coord[i, :] = [yrange[idin[i][1]], xrange[idin[i][2]]]
    end
    
    m_vec = def_val * ones(Nbus)
    d_vec = def_val * ones(Nbus)
    pl_vec = def_val * ones(Nbus)
    pg_vec = def_val * ones(Nbus)
    bx_vec = def_val * ones(Nbus)
    by_vec = def_val * ones(Nbus)
    
    for g in 1:length(idgen)
        k = argmin((lat_coord[:, 2] .- coord[idgen[g],1]).^2 +
            (lat_coord[:, 1] .- coord[idgen[g],2]).^2)
        m_vec[k] += mg[g]
        d_vec[k] += dg[g]
        pg_vec[k] += gen[g]
    end
    
    for l in 1:length(dl)
        k = argmin((lat_coord[:, 2] .- coord[l,1]).^2 +
            (lat_coord[:, 1] .- coord[l,2]).^2)
        d_vec[k] += dl[l]
        pl_vec[k] += dem[l]
    end

    for l in 1:size(idb,1)
        x2 = coord[idb[l,2], 1]
        x1 = coord[idb[l,1], 1]
        y2 = coord[idb[l,2], 2]
        y1 = coord[idb[l,1], 2] 
        dx_l = x2 - x1
        dy_l = y2 - y1
        ds2 = (dy_l^2 + dx_l^2)
        phi = atan(dy_l, dx_l)
        if(dx_l != 0 && dy_l != 0) # if not a transformer
            x = lat_coord[:, 2] # using node locations instead of centers of lines, 
            y = lat_coord[:, 1] # it's "less precise" but way simpler to implement 
        
            beta = (dx_l .* (y1 .- y) + dy_l .* (x .- x1)) ./ ds2
            alpha = (dx_l .* (x .- x1) + dy_l .* (y .- y1)) ./ ds2
            in_seg = (0 .< alpha) .& (alpha .< 1)
            dist = abs.(beta) .* in_seg .* sqrt(ds2) + # if close to the segment
                .!in_seg .* min.(sqrt.((x .- x1).^2 + (y .- y1).^2), # if close to the ends
                sqrt.((x .- x2).^2 + (y .- y2).^2)) 
            bx_vec[dist .< dmax] .+= bline[l] * abs(cos(phi)) * dx^2 * patch
            by_vec[dist .< dmax] .+= bline[l] * abs(sin(phi)) * dx^2 * patch
        end
    end
    
    m = zeros(Ny, Nx)
    d = zeros(Ny, Nx)
    pg = zeros(Ny, Nx)
    pl = zeros(Ny, Nx)
    bx = zeros(Ny, Nx)
    by = zeros(Ny, Nx)
  
    m_new = zeros(Ny, Nx)
    d_new = zeros(Ny, Nx)
    pg_new = zeros(Ny, Nx)
    pl_new = zeros(Ny, Nx)
    bx_new = zeros(Ny, Nx)
    by_new = zeros(Ny, Nx)
    
    m[isinside] = m_vec
    d[isinside] = d_vec
    pg[isinside] = pg_vec
    pl[isinside] = pl_vec
    bx[isinside] = bx_vec
    by[isinside] = by_vec
    
    interval = 100
    @time begin
        for k in 1:Niter
            Threads.@threads for i in 2:Ny-1
                Threads.@threads for j in 2:Nx-1
                    if(isinside[i,j])
                        m_new[i,j] = (1.0 - 4.0 * tau) * m[i,j] + tau * (m[i+1,j] +
                            m[i-1,j] + m[i,j+1] + m[i,j-1])
                        d_new[i,j] = (1.0 - 4.0 * tau) * d[i,j] + tau * (d[i+1,j] +
                            d[i-1,j] + d[i,j+1] + d[i,j-1])
                        pg_new[i,j] = (1.0 - 4.0 * tau) * pg[i,j] + tau * (pg[i+1,j] +
                            pg[i-1,j] + pg[i,j+1] + pg[i,j-1])
                        pl_new[i,j] = (1.0 - 4.0 * tau) * pl[i,j] + tau * (pl[i+1,j] +
                            pl[i-1,j] + pl[i,j+1] + pl[i,j-1])
                        bx_new[i,j] = (1.0 - 4.0 * tau) * bx[i,j] + tau * (bx[i+1,j] +
                            bx[i-1,j] + bx[i,j+1] + bx[i,j-1])
                        by_new[i,j] = (1.0 - 4.0 * tau) * by[i,j] + tau * (by[i+1,j] +
                            by[i-1,j] + by[i,j+1] + by[i,j-1])
                    end
                end
            end
            
            Threads.@threads for k in 1:size(n,1)
                i = Int64(n[k,1])
                j = Int64(n[k,2])
                nx = n[k,4]
                ny = n[k,3]
                if(nx == 1)
                    m_new[i,j] = m_new[i,j-2]
                    d_new[i,j] = d_new[i,j-2]
                    pg_new[i,j] = pg_new[i,j-2]
                    pl_new[i,j] = pl_new[i,j-2]
                    bx_new[i,j] = bx_new[i,j-2]
                    by_new[i,j] = by_new[i,j-2]
                elseif(nx == -1)
                    m_new[i,j] = m_new[i,j+2]
                    d_new[i,j] = d_new[i,j+2]
                    pg_new[i,j] = pg_new[i,j+2]
                    pl_new[i,j] = pl_new[i,j+2]
                    bx_new[i,j] = bx_new[i,j+2]
                    by_new[i,j] = by_new[i,j+2]
                end
                if(ny == 1)
                    m_new[i,j] = m_new[i-2,j]
                    d_new[i,j] = d_new[i-2,j]
                    pg_new[i,j] = pg_new[i-2,j]
                    pl_new[i,j] = pl_new[i-2,j]
                    bx_new[i,j] = bx_new[i-2,j]
                    by_new[i,j] = by_new[i-2,j]
                elseif(ny == -1)
                    m_new[i,j] = m_new[i+2,j]
                    d_new[i,j] = d_new[i+2,j]
                    pg_new[i,j] = pg_new[i+2,j]
                    pl_new[i,j] = pl_new[i+2,j]
                    bx_new[i,j] = bx_new[i+2,j]
                    by_new[i,j] = by_new[i+2,j]
                end
            end
            
            m = copy(m_new)
            d = copy(d_new)
            pg = copy(pg_new)
            pl = copy(pl_new)
            bx = copy(bx_new)
            by = copy(by_new)
            
            if(mod(k,interval) == 0)
                println(k)
            end
        end
    end
    
    # due to how the boundary is treated in the code, interia, damping or
    # power injection on boundary won't be taken into account
    m[.!isinside] .= 0
    d[.!isinside] .= 0
    pl[.!isinside] .= 0
    pg[.!isinside] .= 0
    
    
    #  assign minimal values to ensure numerical stability
    Threads.@threads for i=2:Ny-1
        Threads.@threads for j=2:Nx-1
            if(isgrid[i,j-1] & isgrid[i,j] & !(isborder[i,j-1] & isborder[i,j])) # if the line bx(i,j) is the grid
                bx[i,j] = max(bx[i,j], bmin)
            end
            if(isgrid[i,j] & isgrid[i+1,j] & !(isborder[i,j] & isborder[i+1,j])) # if the line by(i,j) is the grid 
                by[i,j] = max(by[i,j], bmin)
            end
        end
    end
    
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
    p[isinside] = p[isinside] .- sum(p[isinside]) / sum(isinside)

    return bx, by, p, m, d
end


function load_discrete_model(dataname::String)
    data = h5read(dataname, "/")
    mg = vec(data["gen_inertia"])
    dg = vec(data["gen_prim_ctrl"])
    idgen = Int64.(vec(data["gen"][:, 1]))
    coord = alberts_projection(
        data["bus_coord"] ./ (180 / pi),
        13.37616 / 180 * pi,
        46.94653 / 180 * pi,
        10 / 180 * pi,
        50 / 180 * pi,
    )
    dl = vec(data["load_freq_coef"])
    idb = Int64.(data["branch"][:, 1:2])
    bline = 1 ./ data["branch"][:, 4]

    dem = vec(data["bus"][:, 3]) / 100.0
    th = vec(data["bus"][:, 9]) / 180.0 * pi
    gen = vec(data["gen"][:, 2]) / 100.0
    return gen, dem, bline, idb, idgen, coord, mg, dg, dl, th
end


function get_params(isinside::BitMatrix, filename::String)
    data = h5read(filename, "/")
    m = data["m"]
    d = data["d"]
    bx = data["bx"]
    by = data["by"]
    pl = data["pl"]
    pg = data["pg"]

    p = pg - pl
    p[isinside] = p[isinside] .- sum(p[isinside]) / sum(isinside)
    return bx, by, p, m, d
end
