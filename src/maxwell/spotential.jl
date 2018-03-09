

mutable struct MWSingleLayerPotential3D{K}
    wavenumber::K
end

#quaddata(op::MWSingleLayerPotential3D,rs,els) = quaddata(rs,els,(2,))
quaddata(op::MWSingleLayerPotential3D,rs,els) = quadpoints(rs,els,(3,))
quadrule(op::MWSingleLayerPotential3D,refspace,p,y,q,el,qdata) = qdata[1,q]

function kernelvals(op::MWSingleLayerPotential3D,y,p)
    k = op.wavenumber
    r = y - cartesian(p)
    R = norm(r)

    expn = exp(im * k * R)
    green = expn / (4pi*R)
    gradgreen = -(im*k + 1/R) * green / R * r

    KernelValsMaxwell3D(im*k, r, R, green, gradgreen)
end

function integrand(op::MWSingleLayerPotential3D, krn, y, fp, p)
    divJ = fp[2]
    γ = krn.gamma
    G = krn.green

    -G*divJ/γ
end
