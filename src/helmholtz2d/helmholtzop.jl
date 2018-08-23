using SpecialFunctions

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

function cellcellinteractions!(biop::HelmholtzOperator2D, tshs, bshs, tcell, bcell, z)

    regularcellcellinteractions!(biop, tshs, bshs, tcell, bcell, z)

end

function testfunc1()
    print("test function!")
end

function quaddata(op::HelmholtzOperator2D, g::LagrangeRefSpace, f::LagrangeRefSpace, tels, bels)

  tqd = quadpoints(g, tels, (4,))
  bqd = quadpoints(f, bels, (3,))

  return (tpoints=tqd, bpoints=bqd)
end


function quadrule(op::HelmholtzOperator2D, g::LagrangeRefSpace, f::LagrangeRefSpace, i, τ, j, σ, qd)

    DoubleQuadStrategy(
        qd.tpoints[1,i],
        qd.bpoints[1,j]
    )

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

mutable struct PlaneWaveDirichlet{T,P} <: Functional
    wavenumber::T
    direction::P
end

mutable struct PlaneWaveNeumann{T,P} <: Functional
    wavenumber::T
    direction::P
end

mutable struct ScalarTrace{F} <: Functional
    f::F
end

strace(f, mesh::Mesh) = ScalarTrace(f)

(s::ScalarTrace)(x) = s.f(cartesian(x))
integrand(s::ScalarTrace, tx, fx) = dot(tx.value, fx)

shapevals(f::Functional, ϕ, ts) = shapevals(ValOnly(), ϕ, ts)

function (field::PlaneWaveDirichlet)(mp)

    wavenumber = field.wavenumber
    direction  = field.direction

    cart = cartesian(mp)
    exp(-im * wavenumber * dot(direction, cart))
end


function(field::PlaneWaveNeumann)(mp)

    wavenumber = field.wavenumber
    direction  = field.direction

    cart = cartesian(mp)
    norm = normal(mp)

    wave = exp(-im * wavenumber * dot(direction, cart))
    grad = -im * wavenumber * direction * wave

    d = norm[1] * grad[1]
    for i in 2:length(norm)  d += norm[i] * grad[i]  end
    return d
end

function integrand(pw::PlaneWaveDirichlet, sv, fx)
    tx = sv[1]
    return dot(tx, fx)
end

function integrand(pw::PlaneWaveNeumann, sv, fx)
    tx = sv[1]
    d = tx[1] * fx[1]
    for i in 2:length(tx) d += tx[i]*fx[i] end
    return d
end
