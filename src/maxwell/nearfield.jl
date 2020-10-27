

mutable struct MWSingleLayerField3D{K}
    wavenumber::K
end

mutable struct MWDoubleLayerField3D{K}
    wavenumber::K
end

"""
    MWSingleLayerField3D(wavenumber = error())

Create the single layer near field operator, for use with `potential`.
"""
MWSingleLayerField3D(;wavenumber = error("Missing arg: `wavenumber`")) = MWSingleLayerField3D(wavenumber)
MWDoubleLayerField3D(;wavenumber = error("Missing arg: `wavenumber`")) = MWDoubleLayerField3D(wavenumber)


const MWField3D = Union{MWSingleLayerField3D,MWDoubleLayerField3D}
quaddata(op::MWField3D,rs,els) = quadpoints(rs,els,(2,))
quadrule(op::MWField3D,refspace,p,y,q,el,qdata) = qdata[1,q]
function kernelvals(op::MWField3D,y,p)

        k = op.wavenumber
        r = y - cartesian(p)
        R = norm(r)

        kr = k * R
        expn = exp(-im * k * R)
        green = expn / (4pi*R)
        gradgreen = -(im*k + 1/R) * green / R * r

        krn = (trans=r, dist=R, green=green, gradgreen=gradgreen)
end

function integrand(op::MWSingleLayerField3D, krn, y, fp, p)

    γ = im * op.wavenumber

    j = fp.value
    ρ = fp.divergence

    G  = krn.green
    ∇G = krn.gradgreen

    -γ*G*j + ∇G*ρ/γ
end

function integrand(op::MWDoubleLayerField3D, krn, y, fp, p)

    j = fp.value
    ∇G = krn.gradgreen

    ∇G × j
end
