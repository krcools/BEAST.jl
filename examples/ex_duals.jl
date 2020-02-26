using LinearAlgebra

using CompScienceMeshes
using BEAST

include("utils/showfn.jl")

Tetrs = CompScienceMeshes.tetmeshsphere(1.0, 0.35)
Faces = CompScienceMeshes.skeleton_fast(Tetrs, 2)
Edges = CompScienceMeshes.skeleton_fast(Tetrs, 1)

bnd_Faces = boundary(Tetrs)
bnd_Edges = CompScienceMeshes.skeleton_fast(bnd_Faces, 1)

int_Faces = submesh(!in(bnd_Faces), Faces)
int_Edges = submesh(!in(bnd_Edges), Edges)

primal1 = BEAST.nedelecc3d(Tetrs, int_Edges)
primal2 = BEAST.nedelecd3d(Tetrs, int_Faces)

Dir = bnd_Faces
Neu = submesh(!in(Dir), bnd_Faces)
dual1 = BEAST.dual1forms(Tetrs, int_Faces, Dir)
dual2 = BEAST.dual2forms(Tetrs, int_Edges, Neu)

TF = connectivity(int_Faces, Tetrs)
FE = connectivity(int_Edges, int_Faces)

Id = BEAST.Identity()
G = assemble(Id, dual1, primal2)
A = assemble(Id, primal2, primal2)
B = assemble(Id, curl(primal1), curl(primal1))
@show norm(G)
@show norm(TF*G*FE)

B2 = FE'*A*FE;
@show norm(B - B2)

EV = eigen(B)
idx = findfirst(abs.(EV.values) .> 1e-6)
u = real(EV.vectors[:,idx])

hemi = submesh(Tetrs) do Tetr
    cartesian(center(chart(Tetrs, Tetr)))[3] < 0
end
