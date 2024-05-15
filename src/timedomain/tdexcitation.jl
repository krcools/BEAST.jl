abstract type TDFunctional{T} end
Base.eltype(x::TDFunctional{T}) where {T} = T


function quaddata(exc::TDFunctional, testrefs, timerefs, testels, timeels)

    testqd = quadpoints(testrefs, testels, (2,))
    timeqd = quadpoints(timerefs, timeels, (10,))

    testqd, timeqd

end

function quaddata(excitation::TDFunctional,
        test_refspace, time_refspace::DiracBoundary,
        test_elements, time_elements)

    test_quad_data = quadpoints(test_refspace, test_elements, (2,))

    test_quad_data, nothing
end

function quadrule(exc::TDFunctional, testrefs, timerefs, p, τ, r, ρ, qd)

    MultiQuadStrategy(
        qd[1][1,p],
        SingleQuadStrategy2(
            qd[2][1,r]
        )
    )

end

function quadrule(exc::TDFunctional, testrefs, timerefs::DiracBoundary, p, τ, r, ρ, qd)

    MultiQuadStrategy(
        qd[1][1,p],
        nothing
    )

end

function assemble(exc::TDFunctional, testST; quaddata=quaddata, quadrule=quadrule)
    
    stagedtimestep = isa(temporalbasis(testST), BEAST.StagedTimeStep)
    if stagedtimestep
        return staged_assemble(exc, testST; quaddata=quaddata, quadrule=quadrule)
    end
    
    testfns = spatialbasis(testST)
    timefns = temporalbasis(testST)
    Z = zeros(eltype(exc), numfunctions(testfns), numfunctions(timefns))
    store(v,m,k) = (Z[m,k] += v)
    assemble!(exc, testST, store, quaddata=quaddata, quadrule=quadrule)
    return Z
end

function staged_assemble(exc::TDFunctional, testST::SpaceTimeBasis; 
    quaddata=quaddata, quadrule=quadrule)

    @warn "staged assemble of the right-hand side"
    testfns = spatialbasis(testST)
    timefns = temporalbasis(testST)
    stageCount = numstages(timefns)
    Nt = timefns.Nt
    Δt = timefns.Δt
    Z = zeros(eltype(exc), numfunctions(testfns) * stageCount, Nt)
    for i = 1:stageCount
        store(v,m,k) = (Z[(m-1)*stageCount+i,k] += v)
        tbsd = TimeBasisDeltaShifted(timebasisdelta(Δt, Nt), timefns.c[i])
        assemble!(exc, testfns ⊗ tbsd, store,
            quaddata=quaddata, quadrule=quadrule)
    end
    return Z
end

function assemble!(exc::TDFunctional, testST, store;
    quaddata=quaddata, quadrule=quadrule)
    testfns = spatialbasis(testST)
    timefns = temporalbasis(testST)

    testrefs = refspace(testfns)
    timerefs = refspace(timefns)

    testels, testad = assemblydata(testfns)
    timeels, timead = assemblydata(timefns)

    qd = quaddata(exc, testrefs, timerefs, testels, timeels)

    z = zeros(eltype(exc), numfunctions(testrefs), numfunctions(timerefs))
    for p in eachindex(testels)
        τ = testels[p]
        for r in eachindex(timeels)
            ρ = timeels[r]

            fill!(z, 0)
            qr = quadrule(exc, testrefs, timerefs, p, τ, r, ρ, qd)
            momintegrals!(z, exc, testrefs, timerefs, τ, ρ, qr)

            for i in 1 : numfunctions(testrefs)
                for d in 1 : numfunctions(timerefs)

                    v = z[i,d]
                    for (m,a) in testad[p,i]
                        for (k,c) in timead[r,d]
                            store(a*c*v, m, k)
                        end
                    end

                end
            end
        end
    end
end


abstract type NumQuadStrategy end

mutable struct MultiQuadStrategy{P,R} <: NumQuadStrategy
    quad_points::P
    inner_rule::R
end

# TODO: consolidate with the existing definition of SingleQuadStrategy
mutable struct SingleQuadStrategy2{P} <: NumQuadStrategy
    quad_points::P
end

timequadrule(qr::MultiQuadStrategy, p) = qr.inner_rule

function momintegrals!(z, exc::TDFunctional, testrefs, timerefs, τ, ρ, qr)

    for p in qr.quad_points
        x = p.point
        w = p.weight
        f = p.value
        dx = w

        # try
        #     @assert ρ.vertices[1][1] <= cartesian(x)[1] <= ρ.vertices[2][1]
        # catch
        #     @show ρ.vertices[1][1]
        #     @show cartesian(x)[1]
        #     @show ρ.vertices[2][1]
        #     error("")
        # end

        tqr = timequadrule(qr,p)
        timeintegrals!(z, exc, testrefs, timerefs, x, ρ, dx, tqr, f)

    end

end


function timeintegrals!(z, exc::TDFunctional, testrefs, timerefs, testpoint, timeelement, dx, qr, f)

    for p in qr.quad_points
        t = p.point
        w = p.weight
        U = p.value
        dt = w #* jacobian(t) # * volume(timeelement)

        for i in 1 : numfunctions(testrefs)
            for k in 1 : numfunctions(timerefs)
                z[i,k] += dot(f[i][1]*U[k], exc(testpoint,t)) * dt * dx
            end
        end
    end
end



function timeintegrals!(z, exc::TDFunctional,
        spacerefs, timerefs::DiracBoundary,
        testpoint, timeelement,
        dx, qr, testvals)

        # since timeelement uses barycentric coordinates,
        # the first/left vertex has coords u = 1.0!
        testtime = neighborhood(timeelement, point(0.0))
        @assert cartesian(testtime)[1] ≈ timeelement.vertices[2][1]

        for i in 1 : numfunctions(spacerefs)
            z[i,1] += dot(testvals[i][1], exc(testpoint, testtime)) * dx
        end
end
