using CompScienceMeshes, BEAST

# Γ = readmesh(joinpath(dirname(pathof(BEAST)),"../examples/sphere2.in"))
Γ = meshsphere(radius=1.0, h=0.2)
X = BEAST.gwpdiv(Γ; order=3)
Y = BEAST.gwpdiv(Γ; order=3)

ϵ, μ, ω = 1.0, 1.0, 1.0; κ = ω * √(ϵ*μ)
# κ = 3.0
NK, Id = BEAST.DoubleLayerRotatedMW3D(1.0, im*κ), Identity()
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
H = -1/(im*μ*ω)*curl(E)
h = n × H

BEAST.@defaultquadstrat (NK,X,X) BEAST.DoubleNumSauterQstrat(7,8,6,6,6,6)
BEAST.@defaultquadstrat (Id,X,X) BEAST.SingleNumQStrat(7)
BEAST.@defaultquadstrat (h,X) BEAST.SingleNumQStrat(12)

@hilbertspace j
@hilbertspace m
mfie = @discretise (NK-0.5Id)[m,j] == h[m]  j∈X m∈Y
u = gmres(mfie)

include("utils/postproc.jl")
include("utils/plotresults.jl")
