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
        # Γ = meshcuboid(d, d, d/2, h)
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
    runsims(payload, "$fn/solutions"; h)
    return loadsims("$fn/solutions")
end


@target errors (solutions) -> begin
    df = loadsims("$(fn)/solutions")
    numrows = size(df,1)

    uref = df.u[1]
    Xref = df.X[1]
    Γref = geometry(Xref)

    errs = zeros(numrows)
    s = -Maxwell3D.singlelayer(gamma=1.0)
    qstrat12 = BEAST.NonConformingIntegralOpQStrat(BEAST.DoubleNumSauterQstrat(3, 4, 6, 6, 6, 6))

    S11 = assemble(s, Xref, Xref)
    term11 = real(dot(uref, S11*uref))

    function payload(;h)
        r = getrow(df; h=h)
        (;h, u, X) = r
        Γ = geometry(X)

        S12 = assemble(s, Xref, X; quadstrat=[qstrat12])
        S22 = assemble(s, X, X)
        err = sqrt(term11 - 2*real(dot(uref, S12*u)) + real(dot(u, S22*u)))
        return (;err)
    end

    runsims(payload, "$(fn)/errors", h=df.h[end:-1:2])
    return loadsims("$(fn)/errors")
end

# @target compute_errors (solutions) -> begin
#     df = loadsims("sphere/solutions")
#     numrows = size(df,1)

#     uref = df.u[1]
#     Xref = df.X[1]
#     Γref = geometry(Xref)

#     errs = zeros(numrows)
#     s = -Maxwell3D.singlelayer(gamma=1.0)
#     qstrat12 = BEAST.NonConformingIntegralOpQStrat(BEAST.DoubleNumSauterQstrat(3, 4, 6, 6, 6, 6))

#     S11 = assemble(s, Xref, Xref)
#     term11 = real(dot(uref, S11*uref))

#     for i in reverse(2:numrows)
#         r = df[i,:]
#         (;h, u, X) = r
#         @show r, length(u)
#         Γ = geometry(X)

#         S12 = assemble(s, Xref, X; quadstrat=[qstrat12])
#         S22 = assemble(s, X, X)
#         errs[i] = sqrt(term11 - 2*real(dot(uref, S12*u)) + real(dot(u, S22*u)))
#         @show errs[i]
#     end

#     return errs
# end


@target weak_errors (solutions) -> begin
    df = loadsims("sphere/solutions")

    α = df.h[1] / df.h[2]
    W1 = reverse(df.u_Snorm)
    W2 = Iterators.drop(W1,1)
    W3 = Iterators.drop(W1,2)
    p = [log(abs((w2-w1)/(w3-w2))) / log(1/α) for (w1,w2,w3) in zip(W1,W2,W3)]
end


@make errors
# df = loadsims("solutions")
# error()

# using LinearAlgebra
# using Plots
# plot(df.h, real.(df.u_Snorm), marker=:.)
# plot(df.h[end:-1:2], p, marker=:.)

# refsol = df.u_Snorm[1]
# plot(log.(df.h), log.(abs.(df.u_Snorm .- refsol)), marker=:.)
# plot!(log.(df.h), 1.5 * log.(df.h))

# fcr, geo = facecurrents(df.u[1], df.X[1])
# import Plotly
# Plotly.plot(patch(geo, log10.(norm.(fcr))))



