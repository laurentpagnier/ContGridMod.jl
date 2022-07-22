export compute_stable_sol!

function compute_stable_sol!(
    cm::ContModel
)
    # Find slack bus and remove from ids, create map from grounded to full model
    n = size(cm.mesh.coord, 1)
    idSlack = 1
    idsWoSlack = setdiff(1:n, idSlack)
    unground = sparse(idsWoSlack, 1:n-1, ones(n-1), n, n-1)
    pGround = cm.p[idsWoSlack]
    cm.th = -unground * (cm.xi[idsWoSlack, idsWoSlack] \ pGround * cm.mesh.dx^2)
end
