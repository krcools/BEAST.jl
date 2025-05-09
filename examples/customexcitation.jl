macro testitem(name, tags, code) esc(code) end

using BEAST, Test

@testitem "example: custom excitation" tags=[:example] begin
    using CompScienceMeshes

    κ = 1.0

    # Custom excitations can be defined as functions of a point in Cartesian space,
    # complemented with a function specifying the return type.
    E(x) = point(ComplexF64, 1.0, 1.0im, 0.0) * exp(-im*κ*x[3])
    BEAST.scalartype(::typeof(E)) = ComplexF64

    # Such a custom field plays nicely with the tangential trace
    # operators defined as part of the BEAST framework
    e = (n × E) × n

    fn = joinpath(pathof(BEAST), "../../examples/assets/sphere45.in")
    Γ = readmesh(fn)
    # Γ = meshsphere(radius=1.0, h=0.1)
    RT = raviartthomas(Γ)

    t = Maxwell3D.singlelayer(wavenumber=κ)

    @hilbertspace j
    @hilbertspace k

    Txx = assemble(t[k,j], j∈RT, k∈RT)
    ex = assemble(e[k], k∈RT)

    iTXX = BEAST.GMRESSolver(Txx; maxiter=1000, reltol=1e-5)
    uX, ch = solve(iTXX, ex)
    @test ch.isconverged
end

using LinearAlgebra
import Plotly

fcr, geo = facecurrents(uX, RT)
Plotly.plot(patch(geo, norm.(fcr)))