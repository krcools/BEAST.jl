#

#This function ignores Quaddata such that the program does not throw an error.
function quaddata(biop, tshapes, bshapes, 
    test_elements, bsis_elements, qs::GumerovKanekoDuraiswamiStrat)
    # println("Quaddata ignored.")
    return nothing
end


#This function reads the assemble information and do the assemble
function assemblechunk_body!(biop::Helmholtz3DOp{T,Val{0}},
    test_space::LagrangeBasis{0,-1}, test_elements, test_assembly_data, test_cell_ptrs,
    trial_space::LagrangeBasis{0,-1}, trial_elements, trial_assembly_data, trial_cell_ptrs,
    qd, zlocal, store, quadstrat::GumerovKanekoDuraiswamiStrat) where {T}

    myid = Threads.threadid()
    myid == 1 && print("dots out of 10: ")
    todo, done, pctg = length(test_elements), 0, 0
    for (p,(tcell,tptr)) in enumerate(zip(test_elements, test_cell_ptrs))
        for (q,(bcell,bptr)) in enumerate(zip(trial_elements, trial_cell_ptrs))

            A = momintegrals(biop, tcell, bcell, quadstrat)
            store(A/4/pi, tptr, bptr)
        end

        done += 1
        new_pctg = round(Int, done / todo * 100)
        if new_pctg > pctg + 9
            myid == 1 && print(".")
            pctg = new_pctg
        end 
    end
    myid == 1 && println("")
end


#This function looks at what operator is used and return the calculated result for two elements
function momintegrals(biop, tcell, bcell, qrule::GumerovKanekoDuraiswamiStrat)
    vertices2 = tcell.vertices
    vertices1 = bcell.vertices

    if biop isa HH3DSingleLayerFDBIO
        A, ~, ~, ~ = AnalytIntegralD0.GalerkinLaplaceTriGS(vertices1[1], vertices1[2], vertices1[3], vertices2[1], vertices2[2], vertices2[3])
    elseif biop isa HH3DDoubleLayerFDBIO
        ~, A, ~, ~ = AnalytIntegralD0.GalerkinLaplaceTriGS(vertices1[1], vertices1[2], vertices1[3], vertices2[1], vertices2[2], vertices2[3])
    elseif biop isa HH3DDoubleLayerTransposedFDBIO
        ~, A, ~, ~ = AnalytIntegralD0.GalerkinLaplaceTriGS(vertices2[1], vertices2[2], vertices2[3], vertices1[1], vertices1[2], vertices1[3])
    else
        error("Operator not supported. Try Helmholtz3D.doublelayer(;) or Helmholtz3D.singlelayer(;) or Helmholtz3D.doublelayer_transposed(;)")
    end

    return A
end