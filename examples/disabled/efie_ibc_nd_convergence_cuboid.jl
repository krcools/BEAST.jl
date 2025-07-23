using CompScienceMeshes
using BEAST

using Makeitso
using BakerStreet
using DataFrames

using LinearAlgebra
BLAS.set_num_threads(32)

fn = splitext(splitdir(@__FILE__)[end])[1]

# const κ = 1.00
# const η = 1.00
# const ρ = 10.0

struct SingField{T} <: BEAST.Functional
    a::T
    b::T
end

function (f::SingField)(p)
    x = cartesian(p)
    x[3] ≈ 0.5 || return point(0,0,0)
    return sqrt(f.b - x[2]) * sqrt(x[2] - f.a) / sqrt(f.b-x[1]) / sqrt(x[1]-f.a) * point(1,0,0)
end

function BEAST.integrand(f::SingField, gx, ϕx)
    return dot(gx[1], ϕx)
end

function BEAST.scalartype(::SingField{T}) where {T}
    T
end

@target params () -> begin
    (;κ=1.0, η=1.0, ρ=1.0)
end

@target solutions () -> begin

    (;κ, η, ρ) = params
    
    function payload(;h)
        d = 1.0
        Γ = meshcuboid(d, d, d/2, h)
        X = raviartthomas(Γ)

        t = Maxwell3D.singlelayer(wavenumber=κ)
        δ = 0.1
        d = -2*Maxwell3D.singlelayer(gamma=1/δ)
        # E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
        # e = (n × E) × n
        e = SingField(0.0, 1.0)

        @hilbertspace j
        @hilbertspace k

        a = t - ρ*d
        A = assemble(a[k,j], X, X)
        b = assemble(e[k], X)

        Ai = BEAST.GMRESSolver(A; restart=1500, abstol=1e-8, maxiter=10_000)
        u = Ai * b
        return (;u, X)
    end

    α = 0.8
    h = collect(0.4 * α.^(0:17))
    runsims(payload, "$(fn)/solutions"; h)
end

@target refnormsquared (solutions) -> begin
    function payload(;h)
        (;u, X) = getrow(solutions; h=h)
        s = -Maxwell3D.singlelayer(gamma=1.0)
        S = assemble(s, X, X)
        (;uSu=real(dot(u, S*u)))
    end

    runsims(payload, "$fn/refnormsquared"; h=solutions.h)
end

@target errors (solutions, refnormsquared) -> begin

    hmin = minimum(solutions.h)
    uref, Xref = getrow(solutions; h=hmin)[[:u, :X]]
    uSuref = getrow(refnormsquared; h=hmin)[:uSu]
    
    s = -Maxwell3D.singlelayer(gamma=1.0)
    qstrat12 = BEAST.NonConformingIntegralOpQStrat(BEAST.DoubleNumSauterQstrat(3, 4, 6, 6, 6, 6))

    function payload(;h)
        @show h
        (;h, u, X) = getrow(solutions; h=h)
        (;uSu) = getrow(refnormsquared; h=h)

        S12 = assemble(s, Xref, X; quadstrat=[qstrat12])
        err = sqrt(uSuref - 2*real(dot(uref, S12*u)) + uSu)
        @show err
        return (;err)
    end

    runsims(payload, "$(fn)/errors", h=solutions.h[end:-1:2])
end


@make errors
# @make refnormsquared

using Plots
# rn = sqrt(refnormsquared)
plot(log10.(errors.h), log10.(errors.err), marker=:.);
plot!(log10.(errors.h), -0.27 .+ 0.75*log10.(errors.h), label="p = 0.75")


# @make solutions
# import Plotly
# using LinearAlgebra
# (;u, X) = solutions[1,:]
# fcr, geo = facecurrents(u, X)
# Plotly.plot([
#     CompScienceMeshes.patch(geo, log10.(norm.(fcr));
#         lighting_specular=0.0,
#         lighting_diffuse=0.0,
#         lighting_fresnel=0.0),
#     CompScienceMeshes.wireframe(skeleton(geo,1))    
# ])

