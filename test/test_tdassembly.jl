using Test

using CompScienceMeshes
using BEAST
using LinearAlgebra

sol = 1.0
G = BEAST.HH3DSingleLayerTDBIO(sol)

S = readmesh(joinpath(@__DIR__, "assets", "sphere2.in"))
S1 = meshrectangle(1.0, 1.0, 0.25, 3)
S2 = CompScienceMeshes.translate(S1, point(1.0,1.0,1.0))

@show numcells(S1)
@show numcells(S2)

X1 = lagrangecxd0(S1)
x1 = refspace(X1)

X2 = lagrangecxd0(S2)
x2 = refspace(X2)

Δt = 0.2
Nt = 300
ΔR = sol*Δt
P = timebasiscxd0(Δt, Nt)
H = timebasisc0d1(Δt, Nt)

Q = convolve(P, H)
q = refspace(Q)

# Overwrite the default quadrature strategy for
# this operator/space combination
τ1 = chart(S1, first(Iterators.drop(cells(S1),0)))
τ2 = chart(S2, first(Iterators.drop(cells(S2),0)))
ι = BEAST.ring(first(BEAST.rings(τ1, τ2, ΔR)), ΔR)

struct DoubleQuadTimeDomainRule end

function momintegrals!(z, op::typeof(G),
    g::typeof(x1), f::typeof(x2), T::typeof(q),
    τ::typeof(τ1), σ::typeof(τ2), ι::typeof(ι),
    qr::DoubleQuadTimeDomainRule)

    qp_τ = quadpoints(g, [τ], (3,))[1]
    qp_σ = quadpoints(f, [σ], (4,))[1]

    for wpv_τ in qp_τ
        x = wpv_τ.point
        gx = wpv_τ.value[1]
        for wpv_σ in qp_σ
            y = wpv_σ.point
            fy = wpv_σ.value[1]

            R = norm(cartesian(x)-cartesian(y))
            first(ι) <= R <= last(ι) || continue
            maxdegree = numfunctions(T)-1
            TR = [(-R)^d for d in 0:maxdegree]
            dxdy = wpv_τ.weight * wpv_σ.weight

            for i in 1 : numfunctions(g)
                for j in 1 : numfunctions(f)
                    for k in 1 : numfunctions(T)
                        z[i,j,k] += dxdy * gx[i] * fy[j] * TR[k] / (4π*R)
end end end end end end

# Compare results for a single monomial
z1 = zeros(numfunctions(x1), numfunctions(x2), numfunctions(q))
for r in BEAST.rings(τ1, τ2, ΔR)
    ι = BEAST.ring(r, ΔR)
    momintegrals!(z1, G, x1, x2, q, τ1, τ2, ι, DoubleQuadTimeDomainRule())
end

qd = quaddata(G, x1, x2, q, [τ1], [τ2], nothing)
z2 = zeros(numfunctions(x1), numfunctions(x2), numfunctions(q))
for r in BEAST.rings(τ1, τ2, ΔR)
    ι = BEAST.ring(r, ΔR)
    quad_rule = quadrule(G, x1, x2, q, 1, τ1, 1, τ2, r, ι, qd)
    BEAST.momintegrals!(z2, G, x1, x2, q, τ1, τ2, ι, quad_rule)
end

@test z1≈z2 rtol=1e-6

# Compare results for a single basis function
CSM = CompScienceMeshes
distance = norm(cartesian(CSM.center(τ1)) - cartesian(CSM.center(τ2)))
k = floor(Int, distance/Δt/sol)

timead = BEAST.temporalassemblydata(Q, kmax=k)
