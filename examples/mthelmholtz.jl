using CompScienceMeshes, BEAST

function blkdiagm(blks...)
    m = sum(size(b,1) for b in blks)
    n = sum(size(b,2) for b in blks)
    A = zeros(eltype(first(blks)), m, n)
    o1 = 1
    o2 = 1
    for b in blks
        m1 = size(b,1)
        n1 = size(b,2)
        A[o1:o1+m1-1,o2:o2+n1-1] .= b
        o1 += m1
        o2 += n1
    end
    A
end

width, height, h = 1.0, 0.5, 0.05
G1 = meshrectangle(width, height, h)
G2 = CompScienceMeshes.rotate(G1, 0.5π * x̂)
G3 = CompScienceMeshes.rotate(G1, 1.0π * x̂)

G12 = weld(G1,-G2)
G23 = weld(G2,-G3)
G31 = weld(G3,-G1)

X12 = lagrangec0d1(G12)
X23 = lagrangec0d1(G23)
X31 = lagrangec0d1(G31)

Y12 = duallagrangecxd0(G12)
Y23 = duallagrangecxd0(G23)
Y31 = duallagrangecxd0(G31)

κ = 1.0
SL = Helmholtz3D.singlelayer(wavenumber=κ)
HS = Helmholtz3D.hypersingular(gamma=κ*im)
N = Identity()

X = X12 × X23 × X31
Y = Y12 × Y23 × Y31

E = Helmholtz3D.planewave(direction=ẑ, wavenumber=κ)
e = strace(E, G12)
ex = assemble(e, X)

Sxx = assemble(SL, X, X)
Syy = assemble(HS, Y, Y)

N1 = Matrix(assemble(Identity(), X12, Y12))
N2 = Matrix(assemble(Identity(), X23, Y23))
N3 = Matrix(assemble(Identity(), X31, Y31))
Dyx = blkdiagm(inv(N1),inv(N2),inv(N3))

Q = transpose(Dyx) * Syy * Dyx * Sxx
R = transpose(Dyx) * Syy * Dyx * ex

u1, ch1 = solve(BEAST.GMRESSolver(Sxx),ex)
u2, ch2 = solve(BEAST.GMRESSolver(Q),R)

@show ch1.iters
@show ch2.iters

rad, ϕ = 10, π/2
thetas = range(0,π,length=200)
pts = [rad*point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in thetas]

DLp = BEAST.HH3DDoubleLayerNear(wavenumber=κ)
near1 = BEAST.potential(DLp, pts, u1, X, type=ComplexF64)
near2 = BEAST.potential(DLp, pts, u2, X, type=ComplexF64)
inc = E.(pts)

using Plots
plot(thetas, abs.(vec(near1)))
scatter!(thetas, abs.(vec(near2)))

using LinearAlgebra
cond(Sxx)
cond(Q)

w1 = eigvals(Sxx)
w2 = eigvals(Q)

using Plots
plot()
scatter!(w1)
scatter!(w2)

