import SpecialFunctions: hankelh2
import SpecialFunctions: besselk

abstract type HelmholtzOperator2D{T, K} <: IntegralOperator end

# Some unfortunate code duplication from mwops.jl
# TODO: Maybe we can unify by introducing another layer of abstract types
# or using Holy traits?
scalartype(op::HelmholtzOperator2D{T, K}) where {T, K <: Val{0}} = T
scalartype(op::HelmholtzOperator2D{T, K}) where {T, K} = promote_type(T, K)


# TODO: Another example of code duplication from mwops.jl
# Maybe any integral operator has a gamma field?
# Do we have counter examples?
function isstatic(op::HelmholtzOperator2D)
    return typeof(op.gamma) == Val{0}
end

"""
    hh2d_makegammacomplexifneeded(gamma)

Returns a complexified gamma. Unlike the 3D-case, the handling of the Green's
function is more complicated as, for example, the static kernel is not the limit
of the dynamic kernel for k → 0

First, recall that throughout BEAST, we assume a dependency of exp(+iωt).
The wavenumber k and gamma are related via  γ = ik (and thus k = -iγ).
Accordingly, the homogeneous 2D Helmholtz equation reads
    - Δu - k² u = - Δu + γ² u = 0

For physically meaning full real gamma (i.e., γ >= 0), the solutions of
the BIEs are real-valued. For this reason, we will use gamma and not k to 
deduce the underlying scalartype of the operator.

Note that if γ < 0 -- eventhough it is real valued -- the fundamental solution
is no longer the modified Bessel function K, and instead we have to resort to
the general solution, the Hankel function (i.e., we could use the Hankel function for
γ > 0, but it is slower than then modified Bessel function)

The Green's functions are:

2D Laplace:
    G(x, y) = -1/(2π)* ln|x - y|

2D modified Helmholtz equation:
    G(x, y) = 1/(2π)* K₀(γ |x - y|)

2D Helmholtz equation
    G(x, y) = -i/4 * H₀⁽²⁾(-iγ |x - y|)
where γ = ik (and thus k = -iγ)
"""
function hh2d_makegammacomplexifneeded(gamma::Real)
    if gamma < 0.0
        return Complex(gamma)
    else
        return gamma
    end
end

function hh2d_makegammacomplexifneeded(gamma)
    return gamma
end

struct HH2DSingleLayerFDBIO{T, K} <: HelmholtzOperator2D{T, K}
    alpha::T
    gamma::K

    function HH2DSingleLayerFDBIO(alpha, gamma)
        gamma = hh2d_makegammacomplexifneeded(gamma)
        new{typeof(alpha),typeof(gamma)}(alpha, gamma)
    end
end

struct HH2DDoubleLayerFDBIO{T, K} <: HelmholtzOperator2D{T, K}
    alpha::T
    gamma::K

    function HH2DDoubleLayerFDBIO(alpha, gamma)
        gamma = hh2d_makegammacomplexifneeded(gamma)
        new{typeof(alpha),typeof(gamma)}(alpha, gamma)
    end
end

struct HH2DDoubleLayerTransposedFDBIO{T, K} <: HelmholtzOperator2D{T, K}
    alpha::T
    gamma::K

    function HH2DDoubleLayerTransposedFDBIO(alpha, gamma)
        gamma = hh2d_makegammacomplexifneeded(gamma)
        new{typeof(alpha),typeof(gamma)}(alpha, gamma)
    end
end

struct HH2DHyperSingularFDBIO{T, K} <: HelmholtzOperator2D{T, K}
    alpha::T # Weakly singular term
    beta::T # Hypersingular term
    gamma::K

    function HH2DHyperSingularFDBIO(alpha, beta, gamma)
        gamma = hh2d_makegammacomplexifneeded(gamma)
        new{typeof(alpha),typeof(gamma)}(alpha, beta, gamma)
    end
end

mutable struct KernelValsHelmholtz2D{K1,K2,F}
    gamma::K1
    vect::SVector{2,F}
    dist::F
    green::K2
    gradgreen::SVector{2,K2}
    n::SVector{2,F}
    t::SVector{2,F}
end

# Kernel values for 2D Laplace BIE operators
function kernelvals(biop::HelmholtzOperator2D{T, K}, y, p) where {T, K <: Val{0}}

    xc = cartesian(p)
    yc = cartesian(y)
    r = yc - xc
    R = norm(r)

    #Treatment: static limit (HelmholtzOperator2D ~> "LaplaceOperator2D") Δu = 0 (k, γ -> 0)
    #Logarithmic singularity of Hankel function: Gk = -1/(2π) * log(kR) (at static case, lnk = 0)
    green = -1/(2π) * log(R)
    gradgreen = -1/(2π) * (1/R) * (r/R)

    #txty = dot(normal(tgeo), normal(bgeo))
    n = normal(p)
    t = tangents(p, 1) / norm(tangents(p, 1))

    return KernelValsHelmholtz2D(nothing, r, R, green, gradgreen, n, t)
end

# Kernel values for 2D modified Helmholtz equation BIE operators
function kernelvals(biop::HelmholtzOperator2D{T, K}, y, p) where {T, K <: Real}

    γ = biop.gamma
    xc = cartesian(p)
    yc = cartesian(y)
    r = yc - xc
    R = norm(r)

    # Modified Helmholtz equation assumes: Δu - k² u = Δu + γ² u = f
    # Green's function: Gγ = 1/(2π)* K0(γ |x - y|)
    green = 1/(2π) * besselk(0, γ * R)
    # Jin (pg. 713, C.5.9): dK0/dn = -dK1/dn
    gradgreen = -γ/(2π) * besselk(1, γ * R) * (r/R)

    n = normal(p)
    t = tangents(p, 1) / norm(tangents(p, 1))

    return KernelValsHelmholtz2D(γ, r, R, green, gradgreen, n, t)
end

# Kernel values for 2D Helmholtz equation BIE operators
function kernelvals(biop::HelmholtzOperator2D{T, K}, y, p) where {T, K <: Complex}

    # Even though the evaluation of the Hankel function delivers
    # in general a complex output (sole exception if wavenumber is purely imaginary)
    # the Hankel function is evaluated much faster if the wavenumber, and thus
    # the argument of the Hankelfunction, is purely real
    γ = biop.gamma
    if iszero(real(γ))
        k = imag(γ)
    else
        k = -im*γ
    end
    xc = cartesian(p)
    yc = cartesian(y)
    r = yc - xc
    R = norm(r)

    hankels = hankelh2.([0 1], k * R)
    green = - im / 4 * hankels[1]

    # Gradient with respect to observation point
    gradgreen = k * im / 4 * hankels[2] * (r/R)

    n = normal(p)
    t = tangents(p, 1) / norm(tangents(p, 1))

    return KernelValsHelmholtz2D(γ, r, R, green, gradgreen, n, t)
end

function integrand(op::HH2DSingleLayerFDBIO, krn, f, x, g, y)

    α = op.alpha
    G = krn.green
    fx = f.value
    gy = g.value

    return α * G * fx * gy
end

function integrand(op::HH2DDoubleLayerFDBIO, krn, f, x, g, y)

    α = op.alpha
    ∇G = -krn.gradgreen
    ny = krn.n
    fx = f.value
    gy = g.value
    ∂G∂n = dot(ny, ∇G)

    return α * ∂G∂n * fx * gy
end

function integrand(op::HH2DDoubleLayerTransposedFDBIO, krn, f, x, g, y)

    α = op.alpha
    ∇G = krn.gradgreen
    nx = krn.n
    fx = f.value
    gy = g.value
    ∂G∂n = dot(nx, ∇G)

    return -α * ∂G∂n * fx * gy
end

function integrand(op::HH2DHyperSingularFDBIO, krn, f, x, g, y)

    α = op.alpha
    β = op.beta
    G = krn.green
    fx = f.value
    dfx = f.derivative
    gy = g.value
    dgy = g.derivative
    nxny = dot(normal(x), normal(y))

    # Based on (2.86) in Kumagai et al, “Integral equation methods for electromagnetics”
    # The formula returns the electric field, but we like to match the behavior of the
    # dyadic hypersingular operator and this requires an additional rotation
    return (α * fx * gy * nxny + β * dfx * dgy) * G
end

function integrand(op::HH2DHyperSingularFDBIO{T, K}, krn, f, x, g, y) where {T, K <: Val{0}}

    β = op.beta
    G = krn.green
    dfx = f.derivative
    dgy = g.derivative

    return (β * dfx * dgy) * G
end

function quaddata(op::HelmholtzOperator2D,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_charts, trial_charts, qs::DoubleNumSauterQstrat)

    T = coordtype(test_charts[1])

    tqd = quadpoints(test_local_space,  test_charts,  (qs.outer_rule,))
    bqd = quadpoints(trial_local_space, trial_charts, (qs.inner_rule,))

    leg = (
      convert.(NTuple{2,T},_legendre(qs.sauter_schwab_common_vert,0,1)),
      convert.(NTuple{2,T},_legendre(qs.sauter_schwab_common_edge,0,1)),
      convert.(NTuple{2,T},_legendre(qs.sauter_schwab_common_face,0,1)),
    )

    mrw = (
     convert.(NTuple{2,T},BEAST.SauterSchwabQuadrature1D._MRWrules(qs.sauter_schwab_common_vert,0,1)),
     convert.(NTuple{2,T},BEAST.SauterSchwabQuadrature1D._MRWrules(qs.sauter_schwab_common_edge,0,1)),
     convert.(NTuple{2,T},BEAST.SauterSchwabQuadrature1D._MRWrules(qs.sauter_schwab_common_face,0,1)),
    )

    return (tpoints=tqd, bpoints=bqd, gausslegendre=leg, marokhlinwandura=mrw)
end

function quadrule(op::HelmholtzOperator2D, g::LagrangeRefSpace, f::LagrangeRefSpace,
    i, τ::CompScienceMeshes.Simplex{<:Any, 1},
    j, σ::CompScienceMeshes.Simplex{<:Any, 1},
    qd, qs::DoubleNumSauterQstrat)

    hits = _numhits(τ, σ)
    @assert hits <= 2

    hits == 2 && return BEAST.SauterSchwabQuadrature1D.CommonEdge(qd.marokhlinwandura[2], qd.gausslegendre[2])
    hits == 1 && return BEAST.SauterSchwabQuadrature1D.CommonVertex(qd.marokhlinwandura[1], qd.gausslegendre[1])

    return DoubleQuadRule(
        qd.tpoints[1,i],
        qd.bpoints[1,j],
    )
end