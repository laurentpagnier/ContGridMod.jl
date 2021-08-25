function get_params(
    isinside::BitMatrix,
    sigma::Float64,
    dataname::String,
    savename::String,
)
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

    dem = vec(data["bus"][:, 3]) / 100
    gen = vec(data["gen"][:, 2]) / 100

    bx = zeros(Ny, Nx)
    by = zeros(Ny, Nx)
    m = zeros(Ny, Nx)
    d = zeros(Ny, Nx)
    pl = zeros(Ny, Nx)
    pg = zeros(Ny, Nx)

    Threads.@threads for j=2:Nx-1
        Threads.@threads for i=2:Ny-1
            Threads.@threads for g in 1:length(idgen)
                idg1 = coord[idgen[g],1]
                idg2 = coord[idgen[g],2]
                n = norm_dist([idg1, idg2], [xrange[j], yrange[i]], sigma)
                m[i,j] += mg[g] * n
                d[i,j] += dg[g] * n
                pg[i,j] += gen[g] * n
            end
            Threads.@threads for l in 1:length(dl)
                n = norm_dist(vec(coord[l,:]), [xrange[j], yrange[i]], sigma)
                d[i,j] += dl[l] * n
                pl[i,j] += dem[l] * n
            end

            Threads.@threads for l in 1:size(idb,1)
                # a point x on line defined by two points a=(x1,y1) and b=(x2,y2)
                # can be written as x = alpha*(dx,dy) + a, with dx-x2-x1 and dy=y2-y1
                # the line prependicular to the first one and passing by c=(x3,y3) is
                # given by y = beta*(-dy,dx) + c. The point where the two line cross
                # each other is given by y = x. Solving for beta and alpha, the
                # distance separating the point c and the line ab is beta*sqrt(dy^2+dx^2) 
                x2 = coord[idb[l,2],1]
                x1 = coord[idb[l,1],1]
                y2 = coord[idb[l,2],2]
                y1 = coord[idb[l,1],2] 
                dx = x2 - x1
                dy = y2 - y1
                phi = atan(dy,dx)
                # we set the coordinates of the lines to be the ones of the middle points
                bx_x = (xrange[j-1] + xrange[j]) / 2
                bx_y = yrange[i]
                by_x = xrange[j]
                by_y = (yrange[i] + yrange[i+1]) / 2           
                beta_x = (dx*(y1 - bx_y) + dy*(bx_x- x1)) / (dy^2 + dx^2)
                alpha_x = (dx*(bx_x - x1) + dy*(bx_y - y1)) / (dy^2 + dx^2)
                dist1_x = sqrt( (bx_y - y1)^2 + (bx_x - x1)^2 )
                dist2_x = sqrt( (bx_y - y2)^2 + (bx_x - x2)^2 )
                beta_y = ( dx*(y1 - by_y) + dy*(by_x - x1) ) / (dy^2 + dx^2)
                alpha_y = ( dx*(by_x - x1) + dy*(by_y - y1) ) / (dy^2 + dx^2)
                dist1_y = sqrt( (by_y - y1)^2 + (by_x - x1)^2 )
                dist2_y = sqrt( (by_y - y2)^2 + (by_x - x2)^2 )
                # if the line bx(i,j) is the grid
                if(inPolygon([xrange[j-1] yrange[i]],border)[1] & inPolygon([xrange[j] yrange[i]],border)[1])
                
                    # if the projection is in the segment ab
                    if((0 < alpha_x) & (alpha_x < 1))
                        dist = abs(beta_x) * sqrt(dx^2 + dy^2) 
                        # if close enough to the line 
                        if(dist < 4*sigma)
                            f = exp(-dist^2/2/sigma^2)/sigma/sqrt(2*pi)
                            bx[i,j] += bline[l] * abs(cos(phi)) * f * dx^2 / 1000 # 1000 is a temporary patch
                        end
                    # else close enough to an end of the line 
                    elseif(dist1_x < 4*sigma)
                        f = exp(-dist1_x/2/sigma^2)/sigma/sqrt(2*pi)
                        bx[i,j] += bline[l] * abs(cos(phi)) * f * dx^2 / 1000 # 1000 is a temporary patch
                    elseif(dist2_x < 4*sigma)
                        f = exp(-dist2_x/2/sigma^2)/sigma/sqrt(2*pi)
                        bx[i,j] += bline[l] * abs(cos(phi)) * f * dx^2 / 1000 # 1000 is a temporary patch
                    end
                end
                # if the line by(i,j) is the grid
                if(inPolygon([xrange[j] yrange[i]],border)[1] & inPolygon([xrange[j] yrange[i+1]],border)[1])
                    if((0 < alpha_y) & (alpha_y < 1))
                        dist = abs(beta_y) * sqrt(dx^2 + dy^2)
                        if(dist < 4*sigma)
                            f = exp(-dist^2/2/sigma^2)/sigma/sqrt(2*pi)
                            by[i,j] += bline[l] * abs(sin(phi)) * f * dx^2 / 1000 # 1000 is a temporary patch
                        end
                    elseif(dist1_y < 4*sigma)
                        f = exp(-dist1_y/2/sigma^2)/sigma/sqrt(2*pi)
                        by[i,j] += bline[l] * abs(sin(phi)) * f * dx^2 / 1000 # 1000 is a temporary patch
                    elseif(dist2_y < 4*sigma)
                        f = exp(-dist2_y/2/sigma^2)/sigma/sqrt(2*pi)
                        by[i,j] += bline[l] * abs(sin(phi)) * f * dx^2 / 1000 # 1000 is a temporary patch
                    end
                end
            end
        end
    end


    # due to how the boundary is treated in the code, interia, damping or
    # power injection on boundary won't be taken into account
    m[.!isinside] .= 0
    d[.!isinside] .= 0
    pl[.!isinside] .= 0
    pg[.!isinside] .= 0
    
    
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
