using CompScienceMeshes, BEAST
o, x, y, z = euclidianbasis(3)

Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))
X = raviartthomas(Γ)

Δt, Nt = 0.6, 200

δ = timebasisdelta(Δt, Nt)
h = timebasisc0d1(Δt, Nt)
p = timebasiscxd0(Δt, Nt)

width, delay, scaling = 8.0, 12.0, 1.0
gaussian = creategaussian(width, delay, scaling)

direction, polarisation = ẑ, x̂
E = BEAST.planewave(polarisation, direction, derive(gaussian), 1.0)

@assert p(-Δt/2) ≈ 1
@assert h(0) ≈ 1

rhsd = assemble(E, X⊗δ)
rhsp = assemble(E, X⊗p)
rhsh = assemble(E, X⊗h)

imax = floor(Int,delay/Δt)+5

plot()
plot!(rhsd[95,:])
scatter!(rhsp[95,:]/Δt^1, markersize=2.0)
