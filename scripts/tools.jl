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

function get_discrete_values(rangeyrange::Array{Float64,1}, xrange::Array{Float64,1},
    cont_values::Array{Float64,2}, coord::Array{Float64, 2})
    x = repeat(reshape(xrange,1,Nx),Ny,1)
    y = repeat(reshape(yrange,Ny,1),1,Nx)
    disc_values = zeros(size(coord, 1), 1)
    for i = 1:size(coord, 1)
        # TODO
        dx = xrange .- coord[i, 1]
        dy = yrange .- coord[i, 2]
        dist = sqrt.(dx.^2 + dy.^2)
        factor = exp.(-dist) / sum(exp.(-dist))
        v
    end
end
    
function find_time_step(isin::BitMatrix, m::Array{Float64,2}, d::Array{Float64,2}, p::Array{Float64,2},
    bx::Array{Float64,2}, by::Array{Float64,2}, dx::Float64, alpha=0.1)
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
