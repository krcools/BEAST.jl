using Base.Test

using CompScienceMeshes
using BEAST
using StaticArrays

T = Float64
P = SVector{3,T}

omega = 1.0;

id = Identity()
nc = NCross()

fn = Pkg.dir("BEAST","test","sphere2.in")
m = readmesh(fn)
rt = raviartthomas(m)
bc = buffachristiansen(m)

G1 = zeros(Complex128, numfunctions(rt), numfunctions(rt))
BEAST.assemble_local_matched!(id, rt, rt, (v,m,n)->(G1[m,n]+=v))

G2 = zeros(Complex128, numfunctions(rt), numfunctions(rt))
BEAST.assemble_local_mixed!(id, rt, rt, (v,m,n)->(G2[m,n]+=v))

# G1 = BEAST.assemble_local_matched(id, rt, rt)
# G2 = BEAST.assemble_local_mixed(id, rt, rt)

@test norm(G1-G2,  Inf) ≈ 0 atol = sqrt(eps())
@test norm(G1-G1', Inf) ≈ 0 atol = sqrt(eps())
@test norm(G2-G2', Inf) ≈ 0 atol = sqrt(eps())

G = zeros(Complex128, numfunctions(rt), numfunctions(bc))
BEAST.assemble_local_mixed!(nc, rt, bc, (v,m,n)->(G[m,n]+=v))
@test cond(G) ≈ 2.919026332637947

G = zeros(Complex128, numfunctions(rt), numfunctions(rt))
BEAST.assemble_local_mixed!(nc, rt, rt, (v,m,n)->(G[m,n]+=v))
@test norm(G+G',  Inf) ≈ 0 atol = sqrt(eps())


v = [
    point(0.0, 0.0, 0.0),
    point(1.0, 0.0, 0.0),
    point(1.0, 1.0, 0.0),
    point(0.0, 1.0, 0.0),]
f = [
    index(1,2,3),
    index(1,3,4)
]
Γ = Mesh(v,f)
X = raviartthomas(Γ)
I = Identity()
G = assemble(I,X,X)
@test size(G) == (1,1)
@test norm(G[1,1] - 1/3) < sqrt(eps())
