using CompScienceMeshes, BEAST

# Γ = readmesh(joinpath(dirname(pathof(BEAST)),"../examples/sphere2.in"))
Γ = meshsphere(radius=1.0, h=0.15)
X, Y = raviartthomas(Γ), buffachristiansen(Γ)

ϵ, μ, ω = 1.0, 1.0, 1.0; κ, η = ω * √(ϵ*μ), √(μ/ϵ)
# K, N = Maxwell3D.doublelayer(wavenumber=κ), NCross()
nxK = BEAST.DoubleLayerRotatedMW3D(1.0, κ*im)
Id = BEAST.Identity()
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
# E = -η/(im*κ)*BEAST.CurlCurlGreen(κ, ẑ, point(2,0,0))
H = -1/(im*μ*ω)*curl(E)

h = (n × H) × n
nxh = n × H
h = nxh × n

@hilbertspace j
@hilbertspace m
mfie = @discretise (nxK-0.5Id)[m,j] == nxh[m]  j∈X m∈X
u = solve(mfie)

include("utils/postproc.jl")
include("utils/plotresults.jl")

# freeze, store = BEAST.allocatestorage(nxK, X, X, Val{:bandedstorage}, BEAST.LongDelays{:compress})
# @enter BEAST.assemble!(nxK, X, X, store, BEAST.Threading{:single})
# Z1 = freeze()
# Z2 = assemble(Id, X, X)
# Z = Z1 - 0.5*Z2
# b = assemble(nxh, X)
# u = Z \ b
