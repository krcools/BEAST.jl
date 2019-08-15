using CompScienceMeshes
using BEAST
using LinearAlgebra

# G = meshsphere(1.0, 0.25)
G = meshsphere(1.0, 0.30)

c = 1.0
S = BEAST.HH3DSingleLayerTDBIO(c)
D = BEAST.HH3DDoubleLayerTDBIO(speed_of_light=c)
Id = BEAST.Identity()

# width, delay, scaling = 24.0, 36.0, 1.0
width, delay, scaling = 16.0, 24.0, 1.0
gaussian = creategaussian(width, delay, scaling)
fgaussian = fouriertransform(gaussian)
e = BEAST.planewave(point(0,0,1), c, gaussian)
de = BEAST.planewave(point(0,0,1), c, derive(gaussian))
h = BEAST.gradient(e)

# X = lagrangecxd0(G)
X = lagrangecxd0(G)
Y = duallagrangec0d1(G)
# Y = lagrangecxd0(G)

# X = duallagrangecxd0(G, boundary(G))
# Y = lagrangec0d1(G)

# Δt, Nt = 0.9032321, 301
# Δt, Nt = 1.2032321, 301
Δt, Nt = 0.16, 300
# T = timebasisc0d1(Δt, Nt)
P = timebasiscxd0(Δt, Nt)
H = timebasisc0d1(Δt, Nt)
δ = timebasisdelta(Δt, Nt)

# assemble the right hand side

bd  = assemble(n⋅h,     X ⊗ P)
Z1d = assemble(Id ⊗ Id, X ⊗ P, X ⊗ P)
Z0d = assemble(D,       X ⊗ P, X ⊗ P)
Zd = Z0d + (-0.5)*Z1d
u = marchonintime(inv(Zd[:,:,1]), Zd, bd, Nt)

bs = assemble(e, X ⊗ δ)
Zs = assemble(S, X ⊗ δ, X ⊗ P)
v = marchonintime(inv(Zs[:,:,1]), Zs, -bs, Nt)

U, Δω, ω0 = fouriertransform(u, Δt, 0.0, 2)
V, Δω, ω0 = fouriertransform(v, Δt, 0.0, 2)

ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmin(abs.(ω.-1.0)); ω1 = ω[i1]
U_td = U[:,i1] / fgaussian(ω1)
V_td = V[:,i1] / fgaussian(ω1)

fd_S = Helmholtz3D.singlelayer(wavenumber=1.0)
fd_e = Helmholtz3D.planewave(wavenumber=1.0, direction=ẑ)
fd_eh = assemble(strace(fd_e,G), X)
fd_Sh = assemble(fd_S, X, X)
U_fd = fd_Sh \ fd_eh

fcru_td, _ = facecurrents(U_td, X)
fcrv_td, _ = facecurrents(V_td, X)
fcru_fd, _ = facecurrents(U_fd, X)

using Plots
plot()
plot!(u[1,:])
plot!(v[1,:])

# import PlotlyJS
# ptch = patch(G, norm.(fcru_fd))
# PlotlyJS.plot(ptch)

Sstatic = Helmholtz3D.singlelayer(gamma=0.0)
Zstatic = assemble(Sstatic, X, X)
