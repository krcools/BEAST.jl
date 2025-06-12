# -------- used packages
using FastGaussQuadrature

using LinearAlgebra
using StaticArrays

# -------- included files
include("gqlog.jl")


abstract type SauterSchwabStrategy1D end

struct CommonEdge{A} <: SauterSchwabStrategy1D
    qpso::A # MRW quadrature for the outer integral
    qpsi::A # GL quadrature for the inner integral
end
struct CommonVertex{A} <: SauterSchwabStrategy1D
    qpso::A # MRW quadrature for the outer integral
    qpsi::A # GL quadrature for the inner integral
end



"""
	(::CommonEdge)(f, η, ξ)

Regularizing coordinate transform for parametrization on the unit line: [0,1] ↦ Γ.
based on Boundary Element Methods by Sauter and Schwab, example, 5.2.3, p.308

Common face case.
"""
function (::CommonEdge)(f, η, ξ)

    return (1-ξ) *
            (
            f( (1 - η) * (1 - ξ), (1 - η) * (1 - ξ) + ξ )  +
            f( 1 - (1 - η) * (1 - ξ),  1 - (1 - η) * (1 - ξ) - ξ)
            )
end



"""
	(::CommonVertex)(f, η, ξ)

Regularizing coordinate transform for parametrization on the unit line: [0,1] ↦ Γ.
based on Boundary Element Methods by Sauter and Schwab, example, 5.2.3, p.308

Common vertex case.
"""


# We apply the Duffy trick at vertex (0,0) of the two triangles
# created by splitting along the diagonal (0,0)-(1,1)
function (::CommonVertex)(f, η, ξ)

    return ξ * (
        f( ξ, η * ξ) +
        f( η * ξ, ξ)
        )
end


function _legendre(n, a, b)
    x, w = FastGaussQuadrature.gausslegendre(n)
    w .*= (b - a) / 2
    x = (x .+ 1) / 2 * (b - a) .+ a
    collect(zip(x, w))
end

# MRWRules abbreviation of Ma-Rocklin-Wandzura rule
# Ma, J., Rocklin, D., & Wandzura, S. (1996).
#"Generalized Gaussian Quadrature Rules for Systems of Arbitrary Functions."
function _MRWrules(n,a,b)

    x, w = mrwquadrature(n)
    return collect(zip(x,w))
end

function sauterschwab_parameterized1D(integrand, strategy::SauterSchwabStrategy1D)
    return sum(w1 * w2 * strategy(integrand, η, ξ) for (η, w1) in strategy.qpsi, (ξ, w2) in strategy.qpso)
end

# In the reference domain [0, 1] x [0,1], we assume
# that the singularity is at [0, 0] and then apply the
# Duffy trick at the triangles created by splitting along
# the diagonal (0,0)-(1,1).
# ==> Therefore, the vertices must be mapped that the segments
# with vertices [vt1, vt2] and [vs1, vs2] meet at vertices
# vt2 and vs2 (since the barycentric coordinate ξ = 0 is at vt2/vs2)
function reorder(t, s, strat::CommonVertex)

    T = eltype(t[1])
    tol = 1e3 * eps(T)
    # tol = 1e5 * eps(T)
    # tol = sqrt(eps(T))

    # Find the permutation P of t and s that make
    # Pt = [P, A1, A2]
    # Ps = [P, B1, B2]
    I = zeros(Int, 1)
    J = zeros(Int, 1)
    e = 1
    for i in 1:2
        v = t[i]
        for j in 1:2
            w = s[j]
            if norm(w - v) < tol
                I[e] = i
                J[e] = j
                e += 1
                break
            end
        end
        e == 2 && break
    end

    prepend!(I, setdiff([1, 2], I))
    prepend!(J, setdiff([1, 2], J))

    K = zeros(Int, 2)
    for i in 1:2
        for j in 1:2
            if I[j] == i
                K[i] = j
                break
            end
        end
    end

    L = zeros(Int, 2)
    for i in 1:2
        for j in 1:2
            if J[j] == i
                L[i] = j
                break
            end
        end
    end

    return I, J, K, L
end

function reorder(t, s, strat::CommonEdge)

    T = eltype(t[1])
    tol = 1e3 * eps(T)
    # tol = 1e5 * eps(T)
    # tol = sqrt(eps(T))

    I = [1, 2]
    J = zeros(Int, 2)
    v = t[1]
    w = s[1]
    if norm(w - v) < tol
        J[:] = I[:]
    else # If first vertices do not coincide -> swap
        J[1] = 2
        J[2] = 1
    end

    K = zeros(Int, 2)
    for i in 1:2
        for j in 1:2
            if I[j] == i
                K[i] = j
                break
            end
        end
    end

    L = zeros(Int, 2)
    for i in 1:2
        for j in 1:2
            if J[j] == i
                L[i] = j
                break
            end
        end
    end

    return I, J, K, L
end