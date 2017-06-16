## Global initialisation
using CompScienceMeshes
using BoundaryElements
import JLD

## Auxiliary scripts
include(Pkg.dir("CompScienceMeshes","examples","mesh_box_width_lid.jl"))
include(Pkg.dir("CompScienceMeshes","examples","matlab_patches.jl"))
include(Pkg.dir("BoundaryElements","examples","matlab_plotting.jl"))
@matlab cd("C:\\Temp")


## operators and excitations
const λ = 3.0
const κ = 2π/λ
const γ = im*κ
const μ = 1.0
const ϵ = 1.0
const η = sqrt(μ / ϵ)
const ω = κ / sqrt(μ * ϵ)
const α = 1 / γ

const I = BoundaryElements.Identity()
const N = BoundaryElements.NCross()
const T = BoundaryElements.MWSingleLayer3D(γ)
const K = BoundaryElements.MWDoubleLayer3D(γ)
const P = BoundaryElements.SingleLayerTrace(γ)
const nxK = BoundaryElements.DoubleLayerRotatedMW3D(γ)

const F = BoundaryElements.MWFarField3D(γ)
const Θ = linspace(0,π,100)
const Φ = [0.0]
const ffpts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ for ϕ in Φ]

const d = point(1,0,1); d = normalize(d)
const p = point(0,1,0); p = cross(p,d); p = normalize(p)
const E = BoundaryElements.planewavemw3d(direction = d, polarization = p, wavenumber = κ)
const H = 1 / (im*ω*μ) * curl(E)

const et = (n × E) × n
const hx = n × H
const ht = (n × H) × n


## define the solvers
function efie(m)
    X = raviartthomas(m)
    A = assemble(T, X, X)
    b = assemble(et, X)
    u = A \ b
    return u, X
end

function mxmfie(m)
    X = raviartthomas(m)
    Y = buffachristiansen(m)
    A = assemble(K + 0.5N, Y, X)
    b = assemble(ht, Y)
    u = A \ b
    return u, X
end

function clmfie(m)
    X = raviartthomas(m)
    A = assemble(nxK - 0.5I, X, X)
    b = assemble(hx, X)
    u = A \ b
    return u, X
end

##shorthand
farfield(u,X) = potential(F, ffpts, u, X)

## reference solution
m0 = meshcube(1.0, 0.04)
u0, X0 = efie(m0)
ff0 = farfield(u0, X0)
JLD.@save "refsolution.jld" u0

## load the reference solution
m0 = meshcube(1.0, 0.04)
X0 = raviartthomas(m0)
JLD.@load "refsolution.jld"
ff0 = farfield(u0, X0)


## the mfie solutions
hs = [0.2 0.15 0.14 0.13 0.12 0.11 0.1 0.08]
#hs = hs[1:2]

ff1, us1, Xs1 = [], [], []
ff2, us2, Xs2 = [], [], []
ff3, us3, Xs3 = [], [], []

for (i,h) in enumerate(hs)
    h = hs[i]

    m = meshcube(1.0, h)
    @show i, numcells(m)

    u1, X1 = efie(m)
    u2, X2 = clmfie(m)
    u3, X3 = mxmfie(m)

    push!(us1, u1); push!(Xs1, X1)
    push!(us2, u2); push!(Xs2, X2)
    push!(us3, u3); push!(Xs3, X3)

    push!(ff1, farfield(u1, X1))
    push!(ff2, farfield(u2, X2))
    push!(ff3, farfield(u3, X3))
end
