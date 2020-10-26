using LinearAlgebra

using CompScienceMeshes
using BEAST

include("utils/showfn.jl")

Tetrs = CompScienceMeshes.tetmeshsphere(1.0, 0.40)
Faces = CompScienceMeshes.skeleton_fast(Tetrs, 2)
Edges = CompScienceMeshes.skeleton_fast(Tetrs, 1)
Nodes = CompScienceMeshes.skeleton_fast(Tetrs, 0)

hemi = submesh(Tetrs) do Tetr
    cartesian(CompScienceMeshes.center(chart(Tetrs, Tetr)))[3] < 0
end

bnd_Faces = boundary(Tetrs)
bnd_Edges = CompScienceMeshes.skeleton_fast(bnd_Faces, 1)
bnd_Nodes = CompScienceMeshes.skeleton_fast(bnd_Faces, 0)

int_Faces = submesh(!in(bnd_Faces), Faces)
int_Edges = submesh(!in(bnd_Edges), Edges)
int_Nodes = submesh(!in(bnd_Nodes), Nodes)

@show length(int_Edges)
@show length(int_Faces)

primal1 = BEAST.nedelecc3d(Tetrs, int_Edges)
primal2 = BEAST.nedelecd3d(Tetrs, Faces)

Dir = bnd_Faces
Neu = submesh(!in(Dir), bnd_Faces)

@time dual1 = BEAST.dual1forms(Tetrs, Faces, Dir)
error("stop here")
dual2 = BEAST.dual2forms(Tetrs, int_Edges, Dir)

@assert numfunctions(primal1) == numfunctions(dual2)
@assert numfunctions(primal2) == numfunctions(dual1)

Id = BEAST.Identity()
G = assemble(Id, dual1, primal2)
H = assemble(Id, dual2, primal1)

TF = connectivity(Faces, Tetrs)
FE = connectivity(int_Edges, Faces)
EN = connectivity(int_Nodes, int_Edges)

@show norm(G)
@show norm(H)
@show norm(TF*G*FE)
@show norm(FE*H*EN)

# A = assemble(Id, primal2, primal2)
# B = assemble(Id, curl(primal1), curl(primal1))
# B2 = FE'*A*FE;
# @show norm(B - B2)

# EV = eigen(B)
# idx = findfirst(abs.(EV.values) .> 1e-6)
# u = real(EV.vectors[:,idx+1])

# primal1_hemi = restrict(primal1, hemi)
# trace_primal1_hemi = BEAST.ttrace(primal1_hemi, boundary(hemi))
# fcr, geo = facecurrents(u, trace_primal1_hemi)
# Plotly.plot(patch(geo, norm.(fcr)))


# tetrs, bnd, dir, v2t, v2n = BEAST.dual1forms_init(Tetrs, Dir)
# duals11 = BEAST.dual1forms_body(int_Faces[collect(1:10)], tetrs, bnd, dir, v2t, v2n)

# Ggf = assemble(Id, dual1, primal2);
Cgg = assemble(Id, dual1, curl(primal1))

# Cfg = Ggf \ Cgg;
# Conn = Matrix(connectivity(int_Edges, Faces))
# @show norm(Cfg - Conn)

Zpp = zeros(numfunctions(primal1), numfunctions(primal1))
Zdd = zeros(numfunctions(dual1), numfunctions(dual1))

ϵ = BEAST.Multiplicative(p -> 1.0)
μ = BEAST.Multiplicative(p -> 1.0)

Ipp = assemble(ϵ, primal1, primal1)
Idd = assemble(μ, dual1, dual1)

Zpd = zeros(numfunctions(primal1), numfunctions(dual1))

C = [Zpp Cgg'; -Cgg Zdd];
J = [Ipp Zpd; Zpd' Idd];
γ = 1.0 * im

# P = inv(2J - C) * (2J + C);
# Q = (2J - C) * inv(2J + C);

jg = assemble(BEAST.ScalarTrace(p -> point(1,0,0)), primal1)
mG = zeros(numfunctions(dual1))

b = [jg; mG]
u = (C - γ*J) \ b

eg = u[1:numfunctions(primal1)]
hG = u[numfunctions(primal1)+1:end]

primal1_hemi = restrict(primal1, hemi)
trace_primal1_hemi = BEAST.ttrace(primal1_hemi, boundary(hemi))
fcr, geo = facecurrents(eg, trace_primal1_hemi)
Plotly.plot(patch(geo, norm.(fcr)))

# tetrs, bnd, dir, v2t, v2n = BEAST.dualforms_init(Tetrs, Dir)
# @profview BEAST.dual1forms_body(Faces[collect(1:10)], tetrs, bnd, dir, v2t, v2n)