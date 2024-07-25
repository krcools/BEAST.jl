import SpecialFunctions: hankelh2

abstract type HelmholtzOperator2D <: IntegralOperator end
scalartype(::HelmholtzOperator2D) = ComplexF64

struct SingleLayer{T} <: HelmholtzOperator2D
    wavenumber::T
end

struct HyperSingular{T} <: HelmholtzOperator2D
    wavenumber::T
end

struct DoubleLayer{T} <: HelmholtzOperator2D
    wavenumber::T
end

struct DoubleLayerTransposed{T} <: HelmholtzOperator2D
    wavenumber::T
end



mutable struct KernelValsHelmholtz2D
    wavenumber
    vect
    dist
    green
    gradgreen
    txty
end


function kernelvals(biop::HelmholtzOperator2D, tgeo, bgeo)

    k = biop.wavenumber
    r = tgeo.cart - bgeo.cart
    R = norm(r)

    kr = k * R
    hankels = hankelh2.([0 1], kr)
    green = -im / 4 * hankels[1]
    gradgreen = k * im / 4 * hankels[2] * r / R

    txty = dot(normal(tgeo), normal(bgeo))

    KernelValsHelmholtz2D(k, r, R, green, gradgreen, txty)
end


shapevals(op::HelmholtzOperator2D, ϕ, ts) = shapevals(ValDiff(), ϕ, ts)


function integrand(biop::SingleLayer, kerneldata, tvals,
    tgeo, bvals, bgeo)

    gx = tvals[1]
    fy = bvals[1]

    gx * kerneldata.green * fy
end

function integrand(biop::HyperSingular, kernel, tvals, tgeo,
    bvals, bgeo)

    gx = tvals[1]
    fy = bvals[1]

    dgx = tvals[2]
    dfy = bvals[2]

    k    = kernel.wavenumber
    G    = kernel.green
    txty = kernel.txty

    (dgx * dfy - k*k * txty * gx * fy) * G
end

function integrand(biop::DoubleLayer, kernel, fp, mp, fq, mq)
    nq = normal(mq)
    fp[1] * dot(nq, -kernel.gradgreen) * fq[1]
end

function integrand(biop::DoubleLayerTransposed, kernel, fp, mp, fq, mq)
    np = normal(mp)
    fp[1] * dot(np, kernel.gradgreen) * fq[1]
end

function cellcellinteractions!(biop::HelmholtzOperator2D, tshs, bshs, tcell, bcell, z)

    regularcellcellinteractions!(biop, tshs, bshs, tcell, bcell, z)

end

defaultquadstrat(op::HelmholtzOperator2D, tfs, bfs) = DoubleNumSauterQstrat(3,3,0,4,10,10)

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
     convert.(NTuple{2,T},BEAST.SauterSchwabQuadrature1D._NRWrules(qs.sauter_schwab_common_vert,0,1)),
     convert.(NTuple{2,T},BEAST.SauterSchwabQuadrature1D._NRWrules(qs.sauter_schwab_common_edge,0,1)),
     convert.(NTuple{2,T},BEAST.SauterSchwabQuadrature1D._NRWrules(qs.sauter_schwab_common_face,0,1)),
    )

    return (tpoints=tqd, bpoints=bqd, gausslegendre=leg, marokhlinwandura=mrw)
end

function quadrule(op::HelmholtzOperator2D, g::LagrangeRefSpace, f::LagrangeRefSpace,
    i, τ::CompScienceMeshes.Simplex{<:Any, 1},
    j, σ::CompScienceMeshes.Simplex{<:Any, 1},
    qd, qs::DoubleNumSauterQstrat)

    hits = _numhits(τ, σ)
    @assert hits <= 2

    #hits == 2 && return BEAST.SauterSchwabQuadrature1D.CommonEdge(qd.marokhlinwandura[2])
    #hits == 1 && return BEAST.SauterSchwabQuadrature1D.CommonVertex(qd.marokhlinwandura[1])
    hits == 2 && return BEAST.SauterSchwabQuadrature1D.CommonEdge(qd.marokhlinwandura[2])
    hits == 1 && return BEAST.SauterSchwabQuadrature1D.CommonVertex(qd.marokhlinwandura[1])

    return DoubleQuadRule(
        qd.tpoints[1,i],
        qd.bpoints[1,j],
    )
end