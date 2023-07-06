export nodal_plot


function nodal_plot(model::ContModel, fieldname::Symbol; colormap=:inferno, colorbar=true)::Figure
    f = Figure()
    Axis(f[1, 1])
    sp = solutionplot!(model.dh₁, getproperty(model, fieldname), colormap=colormap)
    if colorbar
        Colorbar(f[1, 2], sp)
    end
    return f
end

function nodal_plot(model::ContModel, val::Vector{T}; colormap=:inferno, colorbar=true)::Figure where {T}
    f = Figure()
    Axis(f[1, 1])
    sp = solutionplot!(model.dh₁, val, colormap=colormap)
    if colorbar
        Colorbar(f[1, 2], sp)
    end
    f
end
