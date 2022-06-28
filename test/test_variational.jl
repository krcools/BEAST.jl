using BEAST
using CompScienceMeshes
using Test

@hilbertspace j[1:3]
@hilbertspace k[1:3]

n = BEAST.n

SL = BEAST.diag(BEAST.Maxwell3D.singlelayer(wavenumber=1.0))
EF = BEAST.Maxwell3D.planewave(polarization=point(1,0,0), direction=point(0,0,1), wavenumber=1.0)
ef = (n × EF) × n

bf = SL[j,k]
lf = ef[j]

@test length(bf.terms) == 3
for term in bf.terms
    @test term.test_id == term.trial_id
end

mesh1 = Mesh(
    [
        point(0,0,0),
        point(1,0,0),
        point(1,1,0),
        point(0,1,0),
    ],
    [
        index(1,2,3),
        index(1,3,4)
    ]
)

mesh2 = CompScienceMeshes.translate(mesh1, point(0,0,1))
mesh3 = CompScienceMeshes.translate(mesh1, point(0,0,2))

X1 = raviartthomas(mesh1)
X2 = raviartthomas(mesh2)
X3 = raviartthomas(mesh3)

X = X1 × X2 × X3

dbf = BEAST.discretise(bf, j=>X, k=>X)
dlf = BEAST.discretise(lf, j=>X)

space_mappings = (j=>X, k=>X,)
for sm in space_mappings
    @show sm
end

M = assemble(dbf)
A = Matrix(M)

@test A[1,2] == A[1,3] == A[2,1] == A[2,3] == A[3,1] == A[3,2] == 0
@test A[1,1] ≈ A[2,2] ≈ A[3,3]

N = assemble(SL, X, X)
B = Matrix(N)
@test A == B