using CompScienceMeshes
using BEAST

using Makeitso
using BakerStreet
using DataFrames

fn = splitext(splitdir(@__FILE__)[end])[1]


@target solutions () -> begin
    function payload(;h)
        d = 1.0
        Γ = meshsphere(radius=d, h=h)
        X = raviartthomas(Γ)

        κ, η = 1.0, 1.0
        t = Maxwell3D.singlelayer(wavenumber=κ)
        s = -Maxwell3D.singlelayer(gamma=κ)
        E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
        e = (n × E) × n

        @hilbertspace j
        @hilbertspace k

        A = assemble(t[k,j], X, X)
        S = assemble(s[k,j], X, X)
        b = assemble(e[k], X)

        Ai = BEAST.GMRESSolver(A; restart=1500, abstol=1e-8, maxiter=10_000)
        u = Ai * b
        u_Snorm = real(sqrt(dot(u, S*u)))
        return (;u, X, u_Snorm)
    end

    α = 0.8
    h = collect(0.4 * α.^(0:15))
    return runsims(payload, "$fn/solutions"; h)
end


@target refnormsquared (solutions) -> begin
    (;u, X) = solutions[1,:]
    s = -Maxwell3D.singlelayer(gamma=1.0)
    S = assemble(s, X, X)
    return real(dot(u, S*u))
end


@target errors (solutions, refnormsquared) -> begin

    uref = solutions.u[1]
    Xref = solutions.X[1]

    s = -Maxwell3D.singlelayer(gamma=1.0)
    qstrat12 = BEAST.NonConformingIntegralOpQStrat(BEAST.DoubleNumSauterQstrat(3, 4, 6, 6, 6, 6))

    function payload(;h)
        r = getrow(solutions; h=h)
        (;h, u, X) = r
        Γ = geometry(X)

        S12 = assemble(s, Xref, X; quadstrat=[qstrat12])
        S22 = assemble(s, X, X)
        err = sqrt(refnormsquared - 2*real(dot(uref, S12*u)) + real(dot(u, S22*u)))
        return (;err)
    end

    return runsims(payload, "$(fn)/errors", h=solutions.h[end:-1:2])
end


@make errors

using Plots
plot(log10.(errors.h), log10.(errors.err), marker=:.)
plot!(log10.(errors.h), -0.05 .+ 0.75*log10.(errors.h), label="p = 0.75")
plot!(log10.(errors.h), 1.0 .+ 1.50*log10.(errors.h), label="p = 1.50")




