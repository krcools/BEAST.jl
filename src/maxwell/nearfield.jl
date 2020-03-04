

mutable struct MWSingleLayerField3D{K}
    wavenumber::K
end

"""
    MWSingleLayerField3D(wavenumber = error())

Create the single layer near field operator, for use with `potential`.
"""
MWSingleLayerField3D(;wavenumber = error("Missing arg: `wavenumber`")) = MWSingleLayerField3D(wavenumber)


quaddata(op::MWSingleLayerField3D,rs,els) = quadpoints(rs,els,(2,))
quadrule(op::MWSingleLayerField3D,refspace,p,y,q,el,qdata) = qdata[1,q]

function kernelvals(op::MWSingleLayerField3D,y,p)

        k = op.wavenumber
        r = y - cartesian(p)
        R = norm(r)

        kr = k * R
        expn = exp(-im * k * R)
        green = expn / (4pi*R)
        gradgreen = -(im*k + 1/R) * green / R * r

        KernelValsMaxwell3D(im*k, r, R, green, gradgreen)
end

function integrand(op::MWSingleLayerField3D, krn, y, fp, p)

    j = fp[1]
    ρ = fp[2]

    #κ = krn.
    γ = krn.gamma
    G = krn.green
    ∇G = krn.gradgreen

    -γ*G*j + ∇G.*ρ/γ
end

mutable struct MWDoubleLayerField3D{K}
    wavenumber::K
end

"""
    MWDoubleLayerField3D(wavenumber = error())

Create the double layer near field operator, for use with `potential`.
"""
MWDoubleLayerField3D(;wavenumber = error("Missing arg: `wavenumber`")) = MWDoubleLayerField3D(wavenumber)


quaddata(op::MWDoubleLayerField3D,rs,els) = quadpoints(rs,els,(2,))
quadrule(op::MWDoubleLayerField3D,refspace,p,y,q,el,qdata) = qdata[1,q]

function kernelvals(op::MWDoubleLayerField3D,y,p)

        k = op.wavenumber
        r = y - cartesian(p)
        R = norm(r)

        kr = k * R
        expn = exp(-im * k * R)
        green = expn / (4pi*R)
        gradgreen = -(im*k + 1/R) * green / R * r

        KernelValsMaxwell3D(im*k, r, R, green, gradgreen)
end

function integrand(op::MWDoubleLayerField3D, krn, y, fp, p)

    j = fp[1]
    ∇G = krn.gradgreen

    j×∇G
end
