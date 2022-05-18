export compute_stable_sol!

function compute_stable_sol!(
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
