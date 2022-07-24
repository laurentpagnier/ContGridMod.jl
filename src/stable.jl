export compute_stable_sol!, ompute_stable_sol

function compute_stable_sol!(
    cm::ContModel,
)
    cm.th = compute_stable_sol(cm)
end


function compute_stable_sol(
    cm::ContModel,
)
    B = sparse([cm.mesh.id_edge[:,1]; cm.mesh.id_edge[:,2]],
        [1:cm.mesh.Nedge; 1:cm.mesh.Nedge], [-ones(cm.mesh.Nedge); ones(cm.mesh.Nedge)])
    id = setdiff(1:cm.mesh.Nnode, cm.id_slack)
    I = sparse(id, 1:cm.mesh.Nnode-1, ones(cm.mesh.Nnode-1),
        cm.mesh.Nnode, cm.mesh.Nnode-1)
    Bns = B[id,:]
    Bnst = SparseMatrixCSC{Float64, Int64}(Bns')
    M = -Bns * (cm.b .* Bnst)
    return I * (M \ - cm.p[id] * cm.mesh.dx^2)
end

