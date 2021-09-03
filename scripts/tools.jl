function inPolygon(p::Array{Float64,2}, poly::Array{Float64,2})
    N = size(poly, 1)
    b = Bool.(zeros(size(p,1), 1))
    for k in 1:size(p,1)
        j = N
        for i in 1:N
            if(( (poly[i,2] < p[k,2]) & (poly[j,2]>=p[k,2]) ) | 
                    ( (poly[j,2] < p[k,2]) & (poly[i,2] >=p[k,2]) ))
                if(poly[i,1]+(p[k,2]-poly[i,2])/(poly[j,2]-poly[i,2])*(poly[j,1]-poly[i,1])<p[k,1])
                    b[k] =  ! b[k]
                end
            end
            j = i
        end
    end
    return b
end

function alberts_projection(coord::Array{Float64,2}, lon0::Float64, lat0::Float64, lat1::Float64, lat2::Float64)
    R = 6371 #Earth radius in km
    n = 1/2*(sin(lat1) + sin(lat2))
    theta = n * (coord[:,1] .- lon0)
    c = cos(lat1)^2 + 2*n*sin(lat1)
    rho = R/n*sqrt.(c .- 2*n*sin.(coord[:,2]))
    rho0 = R/n*sqrt.(c - 2*n*sin(lat0)) 
    x = rho .* sin.(theta)
    y = rho0 .- rho .* cos.(theta)
    return [x y]
end
    
function import_border(filename::String)
    data = JSON.parsefile(filename)
    N = size(data["border"],1)
    
    b = zeros(N,2)
    for i = 1:N
        b[i,1] = data["border"][i][1]/180*pi
        b[i,2] = data["border"][i][2]/180*pi
    end
    
    if(b[1,:] != b[end,:])
        b = vcat(b,reshape(b[1,:],1,2))
    end
    
    lon0 = 13.37616/180*pi # German parliament
    lat0 = 46.94653/180*pi # Swiss parliament
    lat1 = 10/180*pi
    lat2 = 50/180*pi
    
    b = alberts_projection(b, lon0, lat0, lat1, lat2)
    return b
end

function do_plot(isin::BitMatrix, values::Array{Float64,2}, xlim::Tuple{Int64, Int64}, ylim::Tuple{Int64, Int64})
    temp = copy(values)
    temp[.!isin] .= NaN
    return contour(temp,fill=true, xlim=xlim, ylim=ylim)
end

function do_plot(isin::BitMatrix, values::Array{Float64,2})
    temp = copy(values)
    temp[.!isin] .= NaN
    return contour(temp,fill=true)
end


function set_ref_phase(isinside::BitMatrix, yrange::Array{Float64,1},
    xrange::Array{Float64,1}, theta::Array{Float64,2}, coord::Array{Float64, 1},
    th_ref::Float64 = 0.0)
    # set the phases according to some reference phase at at a given location 
    idin = findall(isinside)
    cont_coord = zeros(length(idin), 2)
    for i = 1:length(idin)
        cont_coord[i,:] = [xrange[idin[i][2]], yrange[idin[i][1]]]
    end
    dx = cont_coord[:,1] .- coord[1]
    dy = cont_coord[:,2] .- coord[2]
    id = argmin(dx.^2 + dy.^2)
    println(cont_coord[id,:], theta[idin[id]])
    theta[isinside] .= theta[isinside] .- theta[idin[id]] .+ th_ref
    return theta
end



function get_discrete_values(isinside::BitMatrix, yrange::Array{Float64,1},
    xrange::Array{Float64,1}, cont_values::Array{Float64,2}, disc_coord::Array{Float64, 2})
    
    idin = findall(isinside)
    cont_coord = zeros(length(idin), 2)
    for i = 1:length(idin)
        cont_coord[i,:] = [xrange[idin[i][2]], yrange[idin[i][1]]]
    end
    
    disc_values = zeros(size(coord, 1))
    for i = 1:size(disc_coord, 1)
        dx = cont_coord[:,1] .- disc_coord[i,1]
        dy = cont_coord[:,2] .- disc_coord[i,2]
        dist = min.(sqrt.(dx.^2 + dy.^2), 10*minimum(sqrt.(dx.^2 + dy.^2))) # max is here to prevent NaN
        factor = exp.(-dist) ./ sum(exp.(-dist)) # here we use softmax for no particular reason
        if(isnan(sum(factor)))
            println(minimum(dist))
            println(maximum(dist))
        end
        disc_values[i] = sum(factor .* cont_values[idin])
    end
    
    return disc_values
end


function get_cont_values(isinside::BitMatrix, yrange::Array{Float64,1},
    xrange::Array{Float64,1}, disc_coord::Array{Float64, 2}, disc_values::Array{Float64, 1})
    
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
    Niter = 1000
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


function power_flow(isinside::BitMatrix, theta::Array{Float64,2}, bx::Array{Float64,2},
    by::Array{Float64,2}, dx::Float64)
    dx_th = zeros(size(theta))
    dy_th = zeros(size(theta))
    for i in 2:size(theta, 1)-1
        for j in 2:size(theta, 2)-1
            if(isinside[i,j])
                dx_th[i,j] = (theta[i,j+1] - theta[i,j-1]) / 2.0 / dx
                dy_th[i,j] = (theta[i+1,j] - theta[i-1,j]) / 2.0 / dx
            end
        end
    end
    fx = bx .* dx_th
    fy = by .* dy_th
    return fx, fy
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
