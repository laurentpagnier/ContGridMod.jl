using HDF5




function get_params_diff_fast(
    isinside::BitMatrix,
    n::Array{Float64, 2},
    dx::Float64,
    lat_coord::Array{Float64, 2},
    scale_factor::Float64,
    dataname::String,
    savename::String;
    def_val::Float64 = 1E-2,
    dmax::Float64 = 100.0,
    Niter::Int64 = 5000,
    tau::Float64 = 0.001,
    patch::Float64 = 0.001,
    min_factor::Float64 = 0.1,
    bmin::Float64 = 1.0
)
    # this function uses the heat equation to "diffuse" the discrete 
    # parameters over the lattice (tau=kappa*dt/dx^2)

    # load the discrete model
    gen, dem, bline, idb, idgen, coord, mg, dg, dl, th = load_discrete_model(dataname, scale_factor)
    
    isgrid = isinside .| isborder
    isgridflat = vec(isgrid)
    
    # get the coordinates of every point in the grid
    #temp = findall(isgrid)
    #Nbus = length(temp)
    #lat_coord = zeros(Nbus, 2)
    #idgrid = Int64.(zeros(size(temp, 1), 2))
    #for i = 1:Nbus
    #    idgrid[i,:] = [temp[i][1], temp[i][2]]
    #    lat_coord[i, :] = [yrange[temp[i][1]], xrange[temp[i][2]]]
    #end
    grid_coord = lat_coord[vec(isgrid),:]
    idgrid = findall(vec(isgrid))
    idin = findall(vec(isinside))
    N = size(lat_coord, 1)
    #temp = findall(isinside)
    #idin = Int64.(zeros(size(temp, 1), 2))
    #for i = 1:size(temp,1)
    #    idin[i,:] = [temp[i][1], temp[i][2]]
    #end
    m = zeros(N)
    d = zeros(N)
    pl = zeros(N)
    pg = zeros(N)
    bx = zeros(N)
    by = zeros(N)
    println(size(m))
    println("mmm")
    println(tau)
    println(size(grid_coord))
    
    println(minimum(grid_coord[:,1]))
    println(minimum(coord[:,1]))
    println(maximum(grid_coord[:,1]))
    println(maximum(coord[:,1]))
    println(minimum(grid_coord[:,2]))
    println(minimum(coord[:,2]))
    println(maximum(grid_coord[:,2]))
    println(maximum(coord[:,2])) 
    for g in 1:length(idgen)
        k = argmin((grid_coord[:, 1] .- coord[idgen[g],1]).^2 +
            (grid_coord[:, 2] .- coord[idgen[g],2]).^2)
        m[idgrid[k]] += mg[g]
        d[idgrid[k]] += dg[g]
        pg[idgrid[k]] += gen[g]
    end
    
    for l in 1:length(dl)
        k = argmin((grid_coord[:, 1] .- coord[l,1]).^2 +
            (grid_coord[:, 2] .- coord[l,2]).^2)
        d[idgrid[k]] += dl[l]
        pl[idgrid[k]] += dem[l]
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
            y = grid_coord[:, 1] # it's "less precise" but way simpler to implement 
            x = grid_coord[:, 2] # using node locations instead of centers of lines, 
            beta = (dx_l .* (y1 .- y) + dy_l .* (x .- x1)) ./ ds2
            alpha = (dx_l .* (x .- x1) + dy_l .* (y .- y1)) ./ ds2
            in_seg = (0 .< alpha) .& (alpha .< 1)
            dist = abs.(beta) .* in_seg .* sqrt(ds2) + # if close to the segment
                .!in_seg .* min.(sqrt.((x .- x1).^2 + (y .- y1).^2), # if close to the ends
                sqrt.((x .- x2).^2 + (y .- y2).^2))
            bx[idgrid[dist .< dmax]] .+= bline[l] * abs(cos(phi)) * dx^2 * patch
            by[idgrid[dist .< dmax]] .+= bline[l] * abs(sin(phi)) * dx^2 * patch
        end
    end
    m0 = copy(m)
    @time begin
        m = heat_diff_fast(vec(isinside), n, m, Nx, Ny, Niter = Niter, tau = tau)
        println("Done with inertia.")
        #d = heat_diff(idin, n, d, Niter = Niter, tau = tau)
        #println("Done with damping.")
        #pg = heat_diff(idin, n, pg, Niter = Niter, tau = tau)
        #println("Done with with generation.")
        #pl = heat_diff(idin, n, pl, Niter = Niter, tau = tau)
        #println("Done with load.")
        #bx = heat_diff(idin, n, bx, Niter = Niter, tau = tau)
        #by = heat_diff(idin, n, by, Niter = Niter, tau = tau)
    end
    println(sum(abs.(m-m0)))
    # due to how the boundary is treated in the code, interia, damping or
    # power injection on boundary won't be taken into account
    # this is just in case, but should'nt be need as parameters are not
    # allowed to diffuse ouside of the grid
    #println(size(m))
    #println(size(vec(isgrid)))
    #println(Nx*Ny)
    #println(sum(m[.!vec(isgrid)]))
    #println("till here it's fine")
    #m[.!isgrid] .= 0
    #d[.!isgrid] .= 0
    #pl[.!isgrid] .= 0
    #pg[.!isgrid] .= 0
    
    # asign minimal values to the quantities
    m[isgridflat] .= max.(m[isgridflat], min_factor * sum(m) / sum(isgridflat))
    d[isgridflat] .= max.(d[isgridflat], min_factor * sum(d) / sum(isgridflat))
    bx[isgridflat] .= max.(bx[isgridflat], min_factor * sum(bx) / sum(isgridflat))
    by[isgridflat] .= max.(by[isgridflat], min_factor * sum(by) / sum(isgridflat))
    println("till here it's fine")
    # ensure that the integrals of the different quantities over
    # the medium is equivalent to their sum over the discrete model

    #m = sum(mg) .* m ./ trapz((yrange, xrange), m)
    #d = (sum(dg) + sum(dl)) .* d ./ trapz((yrange, xrange), d)
    #pl = sum(dem) .* pl ./ trapz((yrange, xrange), pl)
    #pg = sum(gen) .* pg ./ trapz((yrange, xrange), pg)

    m = (sum(mg) / sum(m) / dx^2) .* m 
    d = ((sum(dg) + sum(dl)) / sum(d) / dx^2) .* d
    pl = (sum(dem) / sum(pl) / dx^2) .* pl
    pg = (sum(gen) / sum(pg) / dx^2) .* pg
    
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
    # that generation matches the demand.
    p = pg - pl 
    println(sum(p))
    p[isgridflat] = p[isgridflat] .- sum(p[isgridflat]) / sum(isgridflat)

    return bx, by, p, m, m0
end





function heat_diff_fast(
    isinsideflat::BitVector,
    n::Array{Float64,2},
    v0::Array{Float64,1},
    Ny::Int64,
    Nx::Int64;
    Niter::Int64 = 5000,
    tau::Float64 = 0.001,
)
    #v_new = copy(v0)
    x = copy(v0)
    
    id1 = Array{Int64,1}()
    id2 = Array{Int64,1}()
    v = Array{Float64,1}()
    for k in 1:length(v0)
        if(isinsideflat[k])
            append!(id1, k)
            append!(id2, k)
            append!(v, -4.0)
            append!(id1, k)
            append!(id2, k-Ny)
            append!(v, 1.0)
            append!(id1, k)
            append!(id2, k+Ny)
            append!(v, 1.0)
            append!(id1, k)
            append!(id2, k-1)
            append!(v, 1.0)
            append!(id1, k)
            append!(id2, k+1)
            append!(v, 1.0)            
        end
    end
    
    #= 
    for id in 1:size(n, 1)
        k = (Int64(n[id, 2]) - 1) * Ny + Int64(n[id, 1])
        ny = n[id, 3]
        nx = n[id, 4] 
        append!(id1, k)
        append!(id2, k)
        append!(v, - ((1.0 - nx) +
            (1.0 + nx) +
            (1.0 - ny) +
            (1.0 + ny)))
        append!(id1, k)
        append!(id2, k-Ny)
        append!(v, (1.0 + nx) )
        append!(id1, k)
        append!(id2, k+Ny)
        append!(v, (1.0 - nx) )
        append!(id1, k)
        append!(id2, k-1)
        append!(v, (1.0 + ny) )   
        append!(id1, k)
        append!(id2, k+1)
        append!(v, (1.0 - ny) )
    end
    =#
    
    for id in 1:size(n, 1)
        k = (Int64(n[id, 2]) - 1) * Ny + Int64(n[id, 1])
        ny = n[id, 3]
        nx = n[id, 4] 
        #append!(id1, k)
        #append!(id2, k)
        #append!(v, - ((1.0 - nx) +
        #    (1.0 + nx) +
        #    (1.0 - ny) +
        #    (1.0 + ny)))
        #a = nx^2 + ny^2
        
        if(nx == 1 && ny == 1)
            append!(id1, k)
            append!(id2, k-Ny)
            append!(v, -2/3)
            append!(id1, k)
            append!(id2, k-2*Ny)
            append!(v, 2/3)
            append!(id1, k)
            append!(id2, k-1)
            append!(v, -2/3)
            append!(id1, k)
            append!(id2, k-2)
            append!(v, 2/3)
        elseif(nx == 1 && ny == -1)
            append!(id1, k)
            append!(id2, k-Ny)
            append!(v, -2/3)
            append!(id1, k)
            append!(id2, k-2*Ny)
            append!(v, 2/3)
            append!(id1, k)
            append!(id2, k+1)
            append!(v, -2/3)
            append!(id1, k)
            append!(id2, k+2)
            append!(v, 2/3)
        elseif(nx == -1 && ny == 1)
            append!(id1, k)
            append!(id2, k+Ny)
            append!(v, -2/3)
            append!(id1, k)
            append!(id2, k+2*Ny)
            append!(v, 2/3)
            append!(id1, k)
            append!(id2, k-1)
            append!(v, -2/3)
            append!(id1, k)
            append!(id2, k-2)
            append!(v, 2/3)
        elseif(nx == -1 && ny == -1)
            append!(id1, k)
            append!(id2, k+Ny)
            append!(v, -2/3)
            append!(id1, k)
            append!(id2, k+2*Ny)
            append!(v, 2/3)
            append!(id1, k)
            append!(id2, k+1)
            append!(v, -2/3)
            append!(id1, k)
            append!(id2, k+2)
            append!(v, 2/3)
        elseif(nx == 1 && ny == 0)
            append!(id1, k)
            append!(id2, k-Ny)
            append!(v, -10/3)
            append!(id1, k)
            append!(id2, k-2*Ny)
            append!(v, 4/3)
            append!(id1, k)
            append!(id2, k-1)
            append!(v, 1.0)
            append!(id1, k)
            append!(id2, k+1)
            append!(v, 1.0)
        elseif(nx == -1 && ny == 0)
            append!(id1, k)
            append!(id2, k+Ny)
            append!(v, -10/3)
            append!(id1, k)
            append!(id2, k+2*Ny)
            append!(v, 4/3)
            append!(id1, k)
            append!(id2, k-1)
            append!(v, 1.0)
            append!(id1, k)
            append!(id2, k+1)
            append!(v, 1.0)
        elseif(nx == 0 && ny == 1)
            append!(id1, k)
            append!(id2, k-Ny)
            append!(v, 1.0)
            append!(id1, k)
            append!(id2, k+Ny)
            append!(v, 1.0)
            append!(id1, k)
            append!(id2, k-1)
            append!(v, -10/3)
            append!(id1, k)
            append!(id2, k-2)
            append!(v, 4/3)
        elseif(nx == 0 && ny == -1)
            append!(id1, k)
            append!(id2, k-Ny)
            append!(v, 1.0)
            append!(id1, k)
            append!(id2, k+Ny)
            append!(v, 1.0)
            append!(id1, k)
            append!(id2, k+1)
            append!(v, -10/3)
            append!(id1, k)
            append!(id2, k+2)
            append!(v, 4/3)
        else
            println("problem")
        end
    end
    println("machin")
    xi = sparse(id1, id2, v, Ny * Nx, Ny * Nx)
    N = Ny * Nx
    I = sparse(1:N, 1:N, ones(N))
    #A = 2.0 .* I - tau .* xi
    #B = 2.0 .* I + tau .* xi
    A = I - tau .* xi
    
    println("machin")
    #println(size(I))
    #println(size(xi))
    #println(size(A))
    #println(A)
    #println(Ny * Nx)
    #println(v)
    for t in 1:Niter
        #x = A \ (B * x ) #way slower when dx -> 0
        x = A \ x
        #gmres!(x, A , B * x)
    end
    println("done")
    return x 
end




function get_params_diff(
    isinside::BitMatrix,
    n::Array{Float64, 2},
    dx::Float64,
    yrange::Array{Float64, 1},
    xrange::Array{Float64, 1},
    scaling_factor::Float64,
    dataname::String,
    savename::String;
    def_val::Float64 = 1E-2,
    dmax::Float64 = 100.0,
    Niter::Int64 = 5000,
    tau::Float64 = 0.001,
    patch::Float64 = 0.001,
    min_factor::Float64 = 0.1,
    bmin::Float64 = 1.0
)
    # this function uses the heat equation to "diffuse" the discrete 
    # parameters over the lattice (tau=kappa*dt/dx^2)

    # load the discrete model
    gen, dem, bline, idb, idgen, coord, mg, dg, dl, th = load_discrete_model(dataname, scaling_factor)
    
    isgrid = isinside .| isborder

    Nx = length(xrange)
    Ny = length(yrange)
    
    # get the coordinates of every point in the grid
    temp = findall(isgrid)
    Nbus = length(temp)
    lat_coord = zeros(Nbus, 2)
    idgrid = Int64.(zeros(size(temp, 1), 2))
    for i = 1:Nbus
        idgrid[i,:] = [temp[i][1], temp[i][2]]
        lat_coord[i, :] = [yrange[temp[i][1]], xrange[temp[i][2]]]
    end
    
    temp = findall(isinside)
    idin = Int64.(zeros(size(temp, 1), 2))
    for i = 1:size(temp,1)
        idin[i,:] = [temp[i][1], temp[i][2]]
    end

    m = zeros(Ny, Nx)
    d = zeros(Ny, Nx)
    pl = zeros(Ny, Nx)
    pg = zeros(Ny, Nx)
    bx = zeros(Ny, Nx)
    by = zeros(Ny, Nx)
    
    for g in 1:length(idgen)
        k = argmin((lat_coord[:, 2] .- coord[idgen[g],1]).^2 +
            (lat_coord[:, 1] .- coord[idgen[g],2]).^2)
        m[idgrid[k,1], idgrid[k,2]] += mg[g]
        d[idgrid[k,1], idgrid[k,2]] += dg[g]
        pg[idgrid[k,1], idgrid[k,2]] += gen[g]
    end
    
    for l in 1:length(dl)
        k = argmin((lat_coord[:, 2] .- coord[l,1]).^2 +
            (lat_coord[:, 1] .- coord[l,2]).^2)
        d[idgrid[k,1], idgrid[k,2]] += dl[l]
        pl[idgrid[k,1], idgrid[k,2]] += dem[l]
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
            bx[idgrid[dist .< dmax,1], idgrid[dist .< dmax,2]] .+= bline[l] * abs(cos(phi)) * dx^2 * patch
            by[idgrid[dist .< dmax,1], idgrid[dist .< dmax,2]] .+= bline[l] * abs(sin(phi)) * dx^2 * patch
        end
    end

    
    @time begin
        m = heat_diff(idin, n, m, Niter = Niter, tau = tau)
        println("Done with inertia.")
        d = heat_diff(idin, n, d, Niter = Niter, tau = tau)
        println("Done with damping.")
        pg = heat_diff(idin, n, pg, Niter = Niter, tau = tau)
        println("Done with with generation.")
        pl = heat_diff(idin, n, pl, Niter = Niter, tau = tau)
        println("Done with load.")
        bx = heat_diff(idin, n, bx, Niter = Niter, tau = tau)
        by = heat_diff(idin, n, by, Niter = Niter, tau = tau)
    end

    # due to how the boundary is treated in the code, interia, damping or
    # power injection on boundary won't be taken into account
    # this is just in case, but should'nt be need as parameters are not
    # allowed to diffuse ouside of the grid
    m[.!isgrid] .= 0
    d[.!isgrid] .= 0
    pl[.!isgrid] .= 0
    pg[.!isgrid] .= 0
    
    # asign minimal values to the quantities
    m[isgrid] .= max.(m[isgrid], min_factor * sum(m) / sum(isgrid))
    d[isgrid] .= max.(d[isgrid], min_factor * sum(d) / sum(isgrid))
    bx[isgrid] .= max.(bx[isgrid], min_factor * sum(bx) / sum(isgrid))
    by[isgrid] .= max.(by[isgrid], min_factor * sum(by) / sum(isgrid))
    
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
    # that generation matches the demand.
    p = pg - pl 
    p[isgrid] = p[isgrid] .- sum(p[isgrid]) / sum(isgrid)

    return bx, by, p, m, d
end




function heat_diff(
    idin::Array{Int64, 2},
    n::Array{Float64,2},
    v0::Array{Float64,2};
    Niter::Int64 = 5000,
    tau::Float64 = 0.001,
)
    v_new = copy(v0)
    v = copy(v0)
    for t in 1:Niter
        Threads.@threads for k in 1:size(idin, 1)
            i = idin[k, 1]
            j = idin[k, 2]
            v_new[i, j] = (1.0 - 4.0 * tau) * v[i, j] + tau * (v[i+1, j] +
                v[i-1, j] + v[i, j+1] + v[i, j-1])

        end
            
        Threads.@threads for k in 1:size(n,1)
            i = Int64(n[k,1])
            j = Int64(n[k,2])
            nx = n[k,4]
            ny = n[k,3]
            v_new[i, j] = (1.0 - 4.0 * tau) * v[i, j] + tau * (
                (1 - ny) * v[i+1, j] +
                (1 - ny) * v[i-1, j] +
                (1 - nx) * v[i, j+1] +
                (1 + nx) *  v[i, j-1]
            )
        end
        v = copy(v_new)
    end
    return v 
end


function load_discrete_model(
    dataname::String,
    scaling_factor::Float64
)
    data = h5read(dataname, "/")
    mg = vec(data["gen_inertia"])
    dg = vec(data["gen_prim_ctrl"])
    idgen = Int64.(vec(data["gen"][:, 1]))
    coord = alberts_projection( data["bus_coord"] ./ (180 / pi) )
    dl = vec(data["load_freq_coef"])
    idb = Int64.(data["branch"][:, 1:2])
    bline = 1 ./ data["branch"][:, 4]
    dem = vec(data["bus"][:, 3]) / 100.0
    th = vec(data["bus"][:, 9]) / 180.0 * pi
    gen = vec(data["gen"][:, 2]) / 100.0
    return gen, dem, bline, idb, idgen, coord[:,[2,1]] / scaling_factor, mg, dg, dl, th
end


function get_params(
    isinside::BitMatrix,
    filename::String
)
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


function norm_dist(
    a::Array{Float64,1},
    b::Array{Float64,1},
    sigma::Float64
)
    return exp(-((a[1] - b[1])^2 + (a[2] - b[2])^2) / 2 / sigma^2) / (2*pi) / sigma^2
end




