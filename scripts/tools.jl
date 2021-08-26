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
    
function import_border(filename)
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

function do_plot(isin::BitMatrix, values::Array{Float64,2})
    temp = copy(values)
    temp[.!isin] .= NaN
    contour(temp,fill=true)
end
