

mutable struct MWSingleLayerField3D{T, U}
    gamma::T
    α::U
    β::U
end

mutable struct MWDoubleLayerField3D{T}
    gamma::T
end

"""
    MWSingleLayerField3D(;gamma, wavenumber, alpha, beta)

Create the single layer near field operator, for use with `potential`.
"""
function MWSingleLayerField3D(;
    gamma=nothing,
    wavenumber=nothing,
    alpha=nothing,
    beta=nothing
)
    gamma, _ = gamma_wavenumber_handler(gamma, wavenumber)

    @assert !isstatic(gamma)

    alpha === nothing && (alpha = -gamma)
    beta  === nothing && (beta  = -1/gamma)

    MWSingleLayerField3D(gamma, alpha, beta)
end

MWSingleLayerField3D(op::MWSingleLayer3D{T,U}) where {T,U} = MWSingleLayerField3D(op.gamma, op.α, op.β)

"""
    MWDoubleLayerField3D(; gamma, wavenumber)

Create the double layer near field operator, for use with `potential`.
"""
function MWDoubleLayerField3D(;
    gamma=nothing,
    wavenumber=nothing
)
    gamma, _ = gamma_wavenumber_handler(gamma, wavenumber)
    @assert !isstatic(gamma)

    MWDoubleLayerField3D(gamma)
end

MWDoubleLayerField3D(op::MWDoubleLayer3D) = MWDoubleLayerField3D(op.gamma)

const MWField3D = Union{MWSingleLayerField3D,MWDoubleLayerField3D}
defaultquadstrat(op::MWField3D, basis) = SingleNumQStrat(2)
quaddata(op::MWField3D,rs,els,qs::SingleNumQStrat) = quadpoints(rs,els,(qs.quad_rule,))
quadrule(op::MWField3D,refspace,p,y,q,el,qdata,qs::SingleNumQStrat) = qdata[1,q]

function kernelvals(op::MWField3D,y,p)

        γ = op.gamma
        r = y - cartesian(p)
        R = norm(r)

        γR = γ*R
        expn = exp(-γR)
        green = expn / (4pi*R)
        gradgreen = -(γ + 1/R) * green / R * r

        krn = (trans=r, dist=R, green=green, gradgreen=gradgreen)
end

function integrand(op::MWSingleLayerField3D, krn, y, fp, p)

    γ = op.gamma

    j = fp.value
    ρ = fp.divergence

    G  = krn.green
    ∇G = krn.gradgreen

    op.α*G*j - op.β*∇G*ρ
end

function integrand(op::MWDoubleLayerField3D, krn, y, fp, p)

    j = fp.value
    ∇G = krn.gradgreen

    ∇G × j
end
