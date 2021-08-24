function get_params(
    isborder::BitMatrix,
    sigmal_l::Float64,
    dataname::String,
    savename::String,
)
    sigma_d = 100# Sigma for the population dynamics
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
    println(size(xrange))
    println(size(yrange))
    println(size(m))
    for j = 2:Nx-1
        println(j)
        for i = 2:Ny-1
            for g = 1:length(idgen)
                idg1 = coord[idgen[g], 1]
                idg2 = coord[idgen[g], 2]
                m[i, j] +=
                    mg[g] / 2 / pi / sigma_l^2 *
                    exp(-((idg1 - xrange[j])^2 + (idg2 - yrange[i])^2) / 2 / sigma_l^2)
                d[i, j] +=
                    dg[g] / 2 / pi / sigma_l^2 *
                    exp(-((idg1 - xrange[j])^2 + (idg2 - yrange[i])^2) / 2 / sigma_l^2)
                pg[i, j] +=
                    gen[g] / 2 / pi / sigma_l^2 *
                    exp(-((idg1 - xrange[j])^2 + (idg2 - yrange[i])^2) / 2 / sigma_l^2)
            end
            for l = 1:length(dl)
                d[i, j] +=
                    dl[l] * exp(
                        -((coord[l, 1] - xrange[j])^2 + (coord[l, 2] - yrange[i])^2) / 2 /
                        sigma_d^2,
                    ) / 2 / pi / sigma_d^2
                # This might need to be discussed 
                pl[i, j] +=
                    dem[l] / 2 / pi / sigma_l^2 * exp(
                        -((coord[l, 1] - xrange[j])^2 + (coord[l, 2] - yrange[i])^2) / 2 /
                        sigma_l^2,
                    )
            end

            for l = 1:size(idb, 1)
                dx = coord[idb[l, 2], 1] - coord[idb[l, 1], 1]
                dy = coord[idb[l, 2], 2] - coord[idb[l, 1], 2]
                phi = atan(dy, dx)
                cxx = (xrange[j-1] + xrange[j]) / 2
                cxy = yrange[i]
                cyx = xrange[j]
                cyy = (yrange[i] + yrange[i+1]) / 2
                beta_x =
                    (dx * (coord[idb[l, 1], 2] - cxy) + dy * (cxx - coord[idb[l, 1], 1])) /
                    (dy^2 + dx^2)
                alpha_x =
                    (dx * (cxx - coord[idb[l, 1], 1]) + dy * (cxy - coord[idb[l, 1], 2])) /
                    (dy^2 + dx^2)
                dist1_x = (cxy - coord[idb[l, 1], 2])^2 + (cxx - coord[idb[l, 1], 1])^2
                dist2_x = (cxy - coord[idb[l, 2], 2])^2 + (cxx - coord[idb[l, 2], 1])^2
                beta_y =
                    (dx * (coord[idb[l, 1], 2] - cyy) + dy * (cyx - coord[idb[l, 1], 1])) /
                    (dy^2 + dx^2)
                alpha_y =
                    (dx * (cyx - coord[idb[l, 1], 1]) + dy * (cyy - coord[idb[l, 1], 2])) /
                    (dy^2 + dx^2)
                dist1_y = (cyy - coord[idb[l, 1], 2])^2 + (cyx - coord[idb[l, 1], 1])^2
                dist2_y = (cyy - coord[idb[l, 2], 2])^2 + (cyx - coord[idb[l, 2], 1])^2
                if (
                    inPolygon([xrange[j-1] yrange[i]], border)[1] &
                    inPolygon([xrange[j] yrange[i]], border)[1]
                )
                    if ((0 < alpha_x) & (alpha_x < 1))
                        dist = abs(beta_x) * sqrt(dx^2 + dy^2)
                        if (dist < 4 * sigma_l)
                            f = exp(-dist^2 / 2 / sigma_l^2) / sigma_l^2 / 2 / pi
                            bx[i, j] += bline[l] * abs(cos(phi)) * f
                        end
                    elseif (dist1_x < 16 * sigma_l^2)
                        f = exp(-dist1_x / 2 / sigma_l^2) / sigma_l^2 / 2 / pi
                        bx[i, j] += bline[l] * abs(cos(phi)) * f
                    elseif (dist2_x < 16 * sigma_l^2)
                        f = exp(-dist2_x / 2 / sigma_l^2) / sigma_l^2 / 2 / pi
                        bx[i, j] += bline[l] * abs(cos(phi)) * f
                    end
                end
                if (
                    inPolygon([xrange[j] yrange[i]], border)[1] &
                    inPolygon([xrange[j] yrange[i+1]], border)[1]
                )
                    if ((0 < alpha_y) & (alpha_y < 1))
                        dist = abs(beta_y) * sqrt(dx^2 + dy^2)
                        if (dist < 4 * sigma_l)
                            f = exp(-dist^2 / 2 / sigma_l^2) / sigma_l^2 / 2 / pi
                            by[i, j] += bline[l] * abs(sin(phi)) * f
                        end
                    elseif (dist1_y < 16 * sigma_l^2)
                        f = exp(-dist1_y / 2 / sigma_l^2) / sigma_l^2 / 2 / pi
                        by[i, j] += bline[l] * abs(sin(phi)) * f
                    elseif (dist2_y < 16 * sigma_l^2)
                        f = exp(-dist2_y / 2 / sigma_l^2) / sigma_l^2 / 2 / pi
                        by[i, j] += bline[l] * abs(sin(phi)) * f
                    end
                end
            end
        end
    end

    m[.!isgrid] .= 0
    d[.!isgrid] .= 0
    pl[.!isgrid] .= 0
    pg[.!isgrid] .= 0
    mm = trapz((yrange, xrange), m)
    dd = trapz((yrange, xrange), d)
    ppl = trapz((yrange, xrange), pl)
    ppg = trapz((yrange, xrange), pg)
    m = sum(mg) .* m ./ mm
    d = (sum(dg) + sum(dl)) .* d ./ dd
    pl = sum(dem) .* pl ./ ppl
    pg = sum(gen) .* pg ./ ppg
    mm = trapz((yrange, xrange), m)
    dd = trapz((yrange, xrange), d)
    ppl = trapz((yrange, xrange), pl)
    ppg = trapz((yrange, xrange), pg)
    println("Numeric inertia: ", mm)
    println("Actual inertia: ", sum(mg))
    println("Numeric damping: ", dd)
    println("Actual damping: ", sum(dg) + sum(dl))
    println("Numeric power draw: ", ppl)
    println("Actual power draw: ", sum(dem))
    println("Numeric power gen: ", ppg)
    println("Actual power gen: ", sum(gen))

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

    isinside = isgrid .& .!isborder
    pl[isinside] = sum(pl) .* pl[isinside] ./ sum(pl[isinside])
    pg[isinside] = sum(pg) .* pg[isinside] ./ sum(pg[isinside])
    pl[.!isinside] .= 0
    pg[.!isinside] .= 0
    p = pl - pg
    p[isinside] = p[isinside] .- sum(p[isinside]) / sum(isinside)

    return bx, by, p, m, d
end

function get_params(filename::String)
    data = h5read(filename, "/")
    m = data["m"]
    d = data["d"]
    bx = data["bx"]
    by = data["by"]
    pl = data["pl"]
    pg = data["pg"]
    isinside = isgrid .& .!isborder
    pl[isinside] = sum(pl) .* pl[isinside] ./ sum(pl[isinside])
    pg[isinside] = sum(pg) .* pg[isinside] ./ sum(pg[isinside])
    pl[.!isinside] .= 0
    pg[.!isinside] .= 0
    p = pl - pg
    p[isinside] = p[isinside] .- sum(p[isinside]) / sum(isinside)
    return bx, by, p, m, d
end
