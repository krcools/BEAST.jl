import SpecialFunctions: hankelh2

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
    txty::F
end

# Kernel values for 2D Laplace BIE operators
function kernelvals(biop::HelmholtzOperator2D{T, K}, tgeo, bgeo) where {T, K <: Val{0}}

    error("We currently do not support 2D Laplace equation")
end

# Kernel values for 2D modified Helmholtz equation BIE operators
function kernelvals(biop::HelmholtzOperator2D{T, K}, tgeo, bgeo) where {T, K <: Real}

    error("We currently do not support 2D modified Helmholtz equation")
end

# Kernel values for 2D Helmholtz equation BIE operators
function kernelvals(biop::HelmholtzOperator2D{T, K}, tgeo, bgeo) where {T, K <: Complex}

    # Even though the evaluation of the Hankel function delivers
    # in general a complex output (sole exception if wavenumber is purely imaginary)
    # the Hankel function is evaluated much faster if the wavenumber, and thus
    # the argument of the Hankelfunction, is purely real
    if iszero(real(biop.gamma))
        k = imag(biop.gamma)
    else
        k = -im*gamma
    end
    r = tgeo.cart - bgeo.cart
    R = norm(r)

    kr = k * R
    hankels = hankelh2.([0 1], kr)
    green = - im / 4 * hankels[1]

    # Gradient with respect to observation point
    gradgreen = k * im / 4 * hankels[2] * r / R

    txty = dot(normal(tgeo), normal(bgeo))

    KernelValsHelmholtz2D(gamma, r, R, green, gradgreen, txty)
end

shapevals(op::HelmholtzOperator2D, ϕ, ts) = shapevals(ValDiff(), ϕ, ts)

function integrand(biop::HH2DSingleLayerFDBIO, kerneldata, tvals,
    tgeo, bvals, bgeo)

    α = biop.alpha

    gx = tvals[1]
    fy = bvals[1]

    return α * gx * kerneldata.green * fy
end

function integrand(biop::HH2DDoubleLayerFDBIO, kernel, fp, mp, fq, mq)
    nq = normal(mq)

    α = biop.alpha

    return α * fp[1] * dot(nq, -kernel.gradgreen) * fq[1]
end

function integrand(biop::HH2DDoubleLayerTransposedFDBIO, kernel, fp, mp, fq, mq)
    np = normal(mp)

    α = biop.alpha

    return α * fp[1] * dot(np, kernel.gradgreen) * fq[1]
end

function integrand(biop::HH2DHyperSingularFDBIO, kernel, tvals, tgeo,
    bvals, bgeo)

    α = biop.alpha
    β = biop.beta

    gx = tvals[1]
    fy = bvals[1]

    dgx = tvals[2]
    dfy = bvals[2]

    G    = kernel.green
    txty = kernel.txty

    return (α * txty * gx * fy + β * dgx * dfy) * G
end

function cellcellinteractions!(biop::HelmholtzOperator2D, tshs, bshs, tcell, bcell, z)

    regularcellcellinteractions!(biop, tshs, bshs, tcell, bcell, z)

end

defaultquadstrat(op::HelmholtzOperator2D, tfs, bfs) = DoubleNumSauterQstrat(30,30,0,4,25,25)

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