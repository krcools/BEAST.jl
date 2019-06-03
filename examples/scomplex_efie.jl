using CompScienceMeshes
using BEAST

Γ1 = readmesh(joinpath(@__DIR__,"sphere2.in"))
Γ2 = CompScienceMeshes.SComplex2D(Γ1)
# Y = raviartthomas(Γ1)
Nd1 = BEAST.nedelec(Γ1)
Nd2 = BEAST.nedelec(Γ2)
RT1 = n × Nd1
RT2 = n × Nd2

N1 = skeleton(Γ1,0)
N2 = skeleton(Γ2,0)

for (n1,n2) in zip(N1.faces, N2.nodes) @assert n1[1] == n2[1] end

E1 = skeleton(Γ1,1)
E2 = skeleton(Γ2,1)

using LinearAlgebra
for (e1,e2) in zip(E1.faces, E2.edges)
    q = getindex.(Ref(E2.nodes), e2)
    @assert q[1][1] == e1[1]
    @assert q[2][1] == e1[2]
    ctr1 = cartesian(center(chart(E1, e1)))
    ctr2 = cartesian(center(chart(E2, e2)))
    if norm(ctr1-ctr2) > 1e-8
        @show ctr1
        @show ctr2
        error("stop")
    end
end

κ = 1.0
t = Maxwell3D.singlelayer(wavenumber=κ)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n

Z1 = assemble(t,RT1,RT1)
Z2 = assemble(t,RT2,RT2)

@hilbertspace j
@hilbertspace k
efie = @discretise t[k,j]==e[k]  j∈RT2 k∈RT2
u = gmres(efie)

Γ = Γ2
X = RT2

include("utils/postproc.jl")
include("utils/plotresults.jl")
