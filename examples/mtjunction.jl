using CompScienceMeshes, BEAST

width, height, h = 1.0, 0.5, 0.05
G1 = meshrectangle(width, height, h)
G2 = CompScienceMeshes.rotate(G1, 1.0π * x̂)
G3 = CompScienceMeshes.rotate(G1, 1.5π * x̂)

G12 = weld(G1,-G2)
G23 = weld(G2,-G3)
G31 = weld(G3,-G1)

G21 = weld(G2,-G1)

X12 = raviartthomas(G12, sort=:none)
X23 = raviartthomas(G23, sort=:none)
X31 = raviartthomas(G31, sort=:none)

X21 = raviartthomas(G21, sort=:none)

Y12 = buffachristiansen(G12, sort=:none)
Y23 = buffachristiansen(G23, sort=:none)
Y31 = buffachristiansen(G31, sort=:none)

Y21 = buffachristiansen(G21)

κ = 1.0
SL = Maxwell3D.singlelayer(wavenumber=κ)
N = NCross()

X = X12 × X23 × X31
Y = Y12 × Y23 × Y31

# X = X12 × X21
# Y = Y12 × Y21

E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n

@hilbertspace j
@hilbertspace k
mtefie = @discretise(
    SL[k,j] == e[k],
    j ∈ X, k ∈ X)

# u = gmres(mtefie)
#
# offset = 1; u12 = u[offset:offset+numfunctions(X12)-1]
# offset += numfunctions(X12); u23 = u[offset:offset+numfunctions(X23)-1]
# offset += numfunctions(X23); u31 = u[offset:offset+numfunctions(X31)-1]
#
# fcr12, _ = facecurrents(u12, X12)
# fcr23, _ = facecurrents(u23, X23)
# fcr31, _ = facecurrents(u31, X31)

import Plotly
using LinearAlgebra

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

Sxx = BEAST.sysmatrix(mtefie)
N1 = assemble(N, X12, Y12)
N2 = assemble(N, X23, Y23)
N3 = assemble(N, X31, Y31)
Nxy = blkdiagm(N1,N2,N3)
Syy = assemble(SL, Y, Y)
iNxy = blkdiagm(inv(Matrix(N1)), inv(Matrix(N2)), inv(Matrix(N3)))

# N21 = assemble(N, X21, Y21)
# Nxy = blkdiagm(N1, N21)
# iNxy = blkdiagm(inv(Matrix(N1)), inv(Matrix(N21)))

ex = assemble(e,X)
P = transpose(iNxy) * Syy * iNxy
Q = P * Sxx;
R = P * ex;


u1, ch1 = solve(BEAST.GMRESSolver(Sxx), ex)
u2, ch2 = solve(BEAST.GMRESSolver(Q), R) 

@show ch1.iters
@show ch2.iters

G123 = weld(G1,G2,G3)
X123 = raviartthomas(G123)

stefie = @discretise(
    SL[k,j] == e[k],
    j ∈ X123, k ∈ X123)

Sst = assemble(SL, X123, X123)
bst = assemble(e, X123)
ust, chst = solve(BEAST.GMRESSolver(Sst),bst)

fcrst, geost = facecurrents(ust, X123)
Plotly.plot(patch(geost, norm.(fcrst)))


Φ, Θ = [0.0], range(0,stop=π,length=100)
pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for ϕ in Φ for θ in Θ]

ffd_mt    = potential(MWFarField3D(wavenumber=κ), pts, u1, X)
ffd_mt_pc = potential(MWFarField3D(wavenumber=κ), pts, u2, X)
ffd_err   = potential(MWFarField3D(wavenumber=κ), pts, u1-u2, X)
ffd_st    = potential(MWFarField3D(wavenumber=κ), pts, ust, X123)

using Plots
plot(norm.(ffd_mt))
scatter!(norm.(ffd_mt_pc))
scatter!(norm.(ffd_st))
scatter!(norm.(ffd_err))