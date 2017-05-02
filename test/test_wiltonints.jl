import BEAST; BE = BEAST;
using CompScienceMeshes
using Base.Test

T = Float64
P = SVector{3,T}
tol = eps(T)*10^3

# s1 = int d/R^3
# s2 = int 1/R
# s3 = int 1
# s4 = int R
# s5 = int R^2
# s6 = int R^3
# v1 = int (rho'-rho)/R^3
# v2 = int (rho'-rho)/R
# v3 = int (rho'-rho)
# v4 = int (rho'-rho)*R
# v5 = int (rho'-rho)*R
# g1 = int grad' (1/R)
# g2 = int grad' R
# g3 = int grad' R^2
# g4 = int grad' R^3
function withquadrules(triangle, r, n)

	u, w = BE.trgauss(n)

    n = (triangle[1] - triangle[3]) × (triangle[2] - triangle[3])
    n /= norm(n)
    d = (r - triangle[1]) ⋅ n
    ρ = r - d*n

    scalar = zeros(T,6)
    vector = zeros(P,6)
    gradgr = zeros(P,4)

	m = simplex(triangle[1], triangle[2], triangle[3])

    for i in 1 : length(w)
        q = neighborhood(m, u[:,i])
        s = cartesian(q)
        dq = w[i] * jacobian(q)

        R = norm(s - r)
        scalar[1] += d / R^3 * dq
        scalar[2] += 1 / R * dq
        scalar[3] += 1 * dq
        scalar[4] += R * dq
        scalar[5] += R^2 * dq
        scalar[6] += R^3 * dq
        vector[1] += (s - ρ) / R^3 * dq
        vector[2] += (s - ρ) / R * dq
        vector[3] += (s - ρ) * dq
        vector[4] += (s - ρ) * R * dq
        vector[5] += (s - ρ) * R^2 * dq
    end

    gradgr[1] = -1( vector[1] - scalar[1] * n )
    gradgr[2] = +1( vector[2] - d*scalar[2]*n )
    gradgr[3] = +2( vector[3] - d*scalar[3]*n )
    gradgr[4] = +3( vector[4] - d*scalar[4]*n )

    return scalar, vector, gradgr
end

tr = [
	point(0.0, 0.0, 0.0),
	point(2.0, 0.0, 0.0),
	point(0.0, 3.0, 0.0)]

pts = [
    (tr[1]+2*tr[2]-1*tr[3]),
    (tr[1]+tr[2])/2 + point(0.0, 0.0, 0.3),
    point(3.0, 0.0, 0.0),
    (tr[1]+tr[2]+tr[3])/3 + point(0.0, 0.0, 0.3),
    tr[2] + point(0.0, 0.0, 0.3),
    (tr[1]+tr[2]+tr[3])/3 + point(0.0, 0.0, 100.0)
]

for p in pts
    s1, v1, g1 = BE.wiltonints(tr[1],tr[2],tr[3],p)
    s2, v2, g2 = withquadrules(tr,p,13)

    for i in length(s1)
        @test norm(s1[i]-s2[i]) < 1.0e-5
    end

    for i in length(v1)
        @test norm(v1[i]-v2[i]) < 1.0e-5
    end

    for i in 1:length(g1)
        for j in 1:3
            @test norm(g1[i][j]-g2[i][j]) < 1.0e-5
        end
    end
end


Γ = readmesh(joinpath(dirname(@__FILE__),"assets","sphere2.in"))
nc = numcells(Γ)
t = chart(Γ, first(cells(Γ)))
s = chart(Γ, last(cells(Γ)))

X = BEAST.raviartthomas(Γ)
x = BEAST.refspace(X)

κ = 1.0
op = BEAST.MWSingleLayer3D(κ)

n = BE.numfunctions(x)
z1 = zeros(Complex128, n, n)
z2 = zeros(z1)

BE = BEAST

tqd = BE.quadpoints(x, [t], (12,))
bqd = BE.quadpoints(x, [s], (13,))

DQ_strategy = BE.DoubleQuadStrategy(tqd[1,1], bqd[1,1])
BEAST.momintegrals!(op, x, x, t, s, z1, DQ_strategy)

SE_strategy = BE.WiltonSEStrategy(
  tqd[1,1],
  BE.DoubleQuadStrategy(
	tqd[1,1],
	bqd[1,1],
  ),
)
BEAST.momintegrals!(op, x, x, t, s, z2, SE_strategy)

@test norm(z1-z2) < 1.0e-7
