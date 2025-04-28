# -------- used packages
using FastGaussQuadrature

using LinearAlgebra
using StaticArrays

# -------- included files
include("gqlog.jl")


abstract type SauterSchwabStrategy1D end

struct CommonEdge{A} <: SauterSchwabStrategy1D
    qps::A
end
struct CommonVertex{A} <: SauterSchwabStrategy1D
    qps::A
end



"""
	(::CommonEdge)(f, ξ, η)

Regularizing coordinate transform for parametrization on the unit line: [0,1] ↦ Γ.
based on Boundary Element Methods by Sauter and Schwab, example, 5.2.3, p.308

Common face case.
"""
function (::CommonEdge)(f, w, z)

    return (1-z) *
            (
            f( (1 - w) * (1 - z), (1 - w) * (1 - z) + z )  +
            f( 1 - (1 - w) * (1 - z),  1 - (1 - w) * (1 - z) - z)
            )
end


#=
# Belgian version: for comparizon purpose
# MSc thesis:"2D electromagnetic field MoM calculations
# using well conditioned higher order polynimials" by Denturck,
# eq. 3.38, p.26
function (::CommonEdge)(f, u, v)

    return (1 - v) *
            (
            f( v + (1 - v) * u, (1 - v) * u )  +  
            f( (1 - v) * u, v + (1 - v) * u )
            )
end
=#

"""
	(::CommonVertex)(f, ξ, η)

Regularizing coordinate transform for parametrization on the unit line: [0,1] ↦ Γ.
based on Boundary Element Methods by Sauter and Schwab, example, 5.2.3, p.308

Common vertex case.
"""
function (::CommonVertex)(f, w, z)
    return z * (
        f((1 - w)*z, 1 - w*z) +
        f(1 - (1 - w)*z, w*z)
    )
end

#=
# Belgian version: for comparizon purpose
# MSc thesis:"2D electromagnetic field MoM calculations
# using well conditioned higher order polynimials" by Denturck,
# eq. 3.39, p.26
function (::CommonVertex)(f, u, v)

    return v * (
        f( 1 - ( 1 - u) * v, u * v) +
        f( u * v, 1 - (1 - u) * v)
        )
end
=#

function _legendre(n, a, b)
    x, w = FastGaussQuadrature.gausslegendre(n)
    w .*= (b - a) / 2
    x = (x .+ 1) / 2 * (b - a) .+ a
    collect(zip(x, w))
end

# NRWRules abbreviation of Ma-Rocklin-Wandzura rule
# Ma, J., Rocklin, D., & Wandzura, S. (1996).
#"Generalized Gaussian Quadrature Rules for Systems of Arbitrary Functions."
function _NRWrules(n,a,b)

    x, w = generalizedquadrature(n)
    return collect(zip(x,w))
end

function sauterschwab_parameterized1D(integrand, strategy::SauterSchwabStrategy1D)

    qps_z = strategy.qps # MRW outer

    # TODO: consider allowing different order for inner integral
    n = length(qps_z)
    qps_w = _legendre(n, 0, 1) # gl_jg outer

    return sum(w1 * w2 * strategy(integrand, ξ, η) for (ξ, w1) in qps_w, (η, w2) in qps_z)
end