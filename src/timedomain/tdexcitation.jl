abstract type TDFunctional{T} end
Base.eltype(x::TDFunctional{T}) where {T} = T

defaultquadstrat(exc::TDFunctional, testfns) = NumSpaceNumTimeQStrat(2, 10)

function assemble(exc::TDFunctional, testST; quadstrat=defaultquadstrat)
    
    stagedtimestep = isa(temporalbasis(testST), BEAST.StagedTimeStep)
    if stagedtimestep
        return staged_assemble(exc, testST; quadstrat)
    end
    
    testfns = spatialbasis(testST)
    timefns = temporalbasis(testST)
    Z = zeros(eltype(exc), numfunctions(testfns), numfunctions(timefns))
    store(v,m,k) = (Z[m,k] += v)
    assemble!(exc, testST, store; quadstrat)
    return Z
end

function staged_assemble(exc::TDFunctional, testST::SpaceTimeBasis; quadstrat=defaultquadstrat)

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
        assemble!(exc, testfns ⊗ tbsd, store; quadstrat)
    end
    return Z
end

function assemble!(exc::TDFunctional, testST, store; quadstrat=defaultquadstrat)

    testfns = spatialbasis(testST)
    timefns = temporalbasis(testST)

    testrefs = refspace(testfns)
    timerefs = refspace(timefns)

    testels, testad = assemblydata(testfns)
    timeels, timead = assemblydata(timefns)

    @show quadstrat
    @assert quadstrat != nothing
    qs = quadstrat(exc, testST)
    qd = quaddata(exc, testrefs, timerefs, testels, timeels, qs)
    
    num_testshapes = numfunctions(testrefs, domain(first(testels)))
    z = zeros(eltype(exc), num_testshapes, numfunctions(timerefs))
    for p in eachindex(testels)
        τ = testels[p]
        for r in eachindex(timeels)
            ρ = timeels[r]

            fill!(z, 0)
            qr = quadrule(exc, testrefs, timerefs, p, τ, r, ρ, qd, qs)
            momintegrals!(z, exc, testrefs, timerefs, τ, ρ, qr)

            for i in 1 : num_testshapes
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
