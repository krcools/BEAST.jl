

defaultquadstrat(op::IntegralOperator, tfs::RefSpace, bfs::RefSpace) =
    DoubleNumSauterQstrat(2,3,5,5,4,3)

"""
    blockassembler(operator, test_space, trial_space) -> assembler

Return a callable object for the creation of blocks within a BEM matrix.

This function performs all tasks common to the assembly of several blocks within
a single boundary element matrix. The return value can be used to generate blocks
by calling it as follows:

    assembler(I,J,storefn)

where `I` and `J` are arrays of indices in `test_space` and `trial_space`, respectively,
corresponding to the rows and columns of the desired block.

Note that the block will be constructed in compressed form, i.e. the rows and columns
of the store that are written into are the positions within `I` and `J` (as opposed
to the positions within `1:numfunctions(test_space)` and `1:numfunctions(trial_space)`).
In particular the size of the constructed block will be `(length(I), length(J))`.

This last property allows the assembly of permutations of the BEM matrix by supplying
for `I` and `J` permutations of `1:numfunctions(test_space)` and
`1:numfunctions(trial_space)`.
"""
function blockassembler end


# """
#     quadrule(operator,test_refspace,trial_refspace,p,test_element,q_trial_element, qd)

# Returns an object that contains all the dynamic (runtime) information that
# defines the integration strategy that will be used by `momintegrals!` to compute
# the interactions between the local test/trial functions defined on the specified
# geometric elements. The indices `p` and `q` refer to the position of the test
# and trial elements as encountered during iteration over the output of
# `geometry`.

# The last argument `qd` provides access to all precomputed data required for
# quadrature. For example it might be desirable to precompute all the quadrature
# points for all possible numerical quadrature schemes that can potentially be
# required during matrix assembly. This makes sense, since the number of point is
# order N (where N is the number of faces) but these points will appear in N^2
# computations. Precomputation requires some extra memory but can save a lot on
# computation time.
# """
# function quadrule end


"""
  elements(geo)

Create an iterable collection of the elements stored in `geo`. The order in which
this collection produces the elements determines the index used for lookup in the
data structures returned by `assemblydata` and `quaddata`.
"""
elements(geo) = [chart(geo,cl) for cl in geo]

elements(sp::Space) = elements(geometry(sp))

"""
    assemblechunk!(biop::IntegralOperator, tfs, bfs, store)

Computes the matrix of operator biop wrt the finite element spaces tfs and bfs
"""
function assemblechunk!(biop::IntegralOperator, tfs::Space, bfs::Space, store;
        quadstrat=defaultquadstrat(biop, tfs, bfs))

    tr = assemblydata(tfs); tr == nothing && return
    br = assemblydata(bfs); br == nothing && return

    test_elements, tad, tcells = tr
    bsis_elements, bad, bcells = br

    tgeo = geometry(tfs)
    bgeo = geometry(bfs)

    tdom = domain(chart(tgeo, first(tgeo)))
    bdom = domain(chart(bgeo, first(bgeo)))

    tshapes = refspace(tfs); num_tshapes = numfunctions(tshapes, tdom)
    bshapes = refspace(bfs); num_bshapes = numfunctions(bshapes, bdom)

    qs = if CompScienceMeshes.refines(tgeo, bgeo)
        TestRefinesTrialQStrat(quadstrat)
    elseif CompScienceMeshes.refines(bgeo, tgeo)
        TrialRefinesTestQStrat(quadstrat)
    else
        quadstrat
    end

    qd = quaddata(biop, tshapes, bshapes, test_elements, bsis_elements, qs)
    zlocal = zeros(scalartype(biop, tfs, bfs), 2num_tshapes, 2num_bshapes)
    assemblechunk_body!(biop,
        tfs, test_elements, tad, tcells,
        bfs, bsis_elements, bad, bcells,
        qd, zlocal, store; quadstrat=qs)

    # if CompScienceMeshes.refines(tgeo, bgeo)
    #     assemblechunk_body_test_refines_trial!(biop,
    #         tfs, test_elements, tad, tcells,
    #         bfs, bsis_elements, bad, bcells,
    #         qd, zlocal, store; quadstrat)
    #     qs = TestRefinesTrialQStrat(quadstrat)
    #     # assemblechunk_body!(biop,
    #     #     tfs, test_elements, tad, tcells,
    #     #     bfs, bsis_elements, bad, bcells,
    #     #     qd, zlocal, store; quadstrat=qs)
    # elseif CompScienceMeshes.refines(bgeo, tgeo)
    #     qs = TrialRefinesTestQStrat(quadstrat)
    #     assemblechunk_body!(biop,
    #         tfs, test_elements, tad, tcells,
    #         bfs, bsis_elements, bad, bcells,
    #         qd, zlocal, store; quadstrat=qs)
    #     # assemblechunk_body_trial_refines_test!(biop,
    #     #     tfs, test_elements, tad, tcells,
    #     #     bfs, bsis_elements, bad, bcells,
    #     #     qd, zlocal, store; quadstrat)
    # else
    #     assemblechunk_body!(biop,
    #         tfs, test_elements, tad, tcells,
    #         bfs, bsis_elements, bad, bcells,
    #         qd, zlocal, store; quadstrat)
    # end
end


function assemblechunk_body!(biop,
        test_space, test_elements, test_assembly_data, test_cell_ptrs,
        trial_space, trial_elements, trial_assembly_data, trial_cell_ptrs,
        qd, zlocal, store; quadstrat)

    test_shapes = refspace(test_space)
    trial_shapes = refspace(trial_space)

    myid = Threads.threadid()
    myid == 1 && print("dots out of 10: ")
    todo, done, pctg = length(test_elements), 0, 0
    for (p,(tcell,tptr)) in enumerate(zip(test_elements, test_cell_ptrs))
        for (q,(bcell,bptr)) in enumerate(zip(trial_elements, trial_cell_ptrs))

        fill!(zlocal, 0)
        qrule = quadrule(biop, test_shapes, trial_shapes, p, tcell, q, bcell, qd, quadstrat)
        momintegrals!(zlocal, biop,
            test_space,  tptr, tcell,
            trial_space, bptr, bcell, qrule)
        I = length(test_assembly_data[p])
        J = length(trial_assembly_data[q])
        for j in 1 : J, i in 1 : I
            zij = zlocal[i,j]
            for (n,b) in trial_assembly_data[q][j]
                zb = zij*b
                for (m,a) in test_assembly_data[p][i]
                    store(a*zb, m, n)
        end end end end

        done += 1
        new_pctg = round(Int, done / todo * 100)
        if new_pctg > pctg + 9
            myid == 1 && print(".")
            pctg = new_pctg
    end end
    myid == 1 && println("")
end


# function assemblechunk_body_test_refines_trial!(biop,
#     test_functions, test_charts, test_assembly_data, test_cells,
#     trial_functions, trial_charts, trial_assembly_data, trial_cells,
#     qd, zlocal, store; quadstrat)

#     test_shapes = refspace(test_functions)
#     trial_shapes = refspace(trial_functions)

#     myid = Threads.threadid()
#     myid == 1 && print("dots out of 10: ")
#     todo, done, pctg = length(test_charts), 0, 0
#     for (p,(tcell,tchart)) in enumerate(zip(test_cells, test_charts))
#         for (q,(bcell,bchart)) in enumerate(zip(trial_cells, trial_charts))

#             fill!(zlocal, 0)
#             qrule = quadrule(biop, test_shapes, trial_shapes, p, tchart, q, bchart, qd, quadstrat)
#             # @show ("1", qrule)
#             momintegrals_test_refines_trial!(zlocal, biop,
#                 test_functions, tcell, tchart,
#                 trial_functions, bcell, bchart,
#                 qrule, quadstrat)

#             I = length(test_assembly_data[p])
#             J = length(trial_assembly_data[q])
#             for j in 1 : J, i in 1 : I
#                 zij = zlocal[i,j]
#                 for (n,b) in trial_assembly_data[q][j]
#                     zb = zij*b
#                     for (m,a) in test_assembly_data[p][i]
#                         store(a*zb, m, n)
#         end end end end

#         done += 1
#         new_pctg = round(Int, done / todo * 100)
#         if new_pctg > pctg + 9
#             myid == 1 && print(".")
#             pctg = new_pctg
#         end end
#     myid == 1 && println("")
# end


# function assemblechunk_body_trial_refines_test!(biop,
#     test_functions, test_charts, test_assembly_data, test_cells,
#     trial_functions, trial_charts, trial_assembly_data, trial_cells,
#     qd, zlocal, store; quadstrat)

#     test_shapes = refspace(test_functions)
#     trial_shapes = refspace(trial_functions)

#     myid = Threads.threadid()
#     myid == 1 && print("dots out of 10: ")
#     todo, done, pctg = length(test_charts), 0, 0
#     for (p,(tcell,tchart)) in enumerate(zip(test_cells, test_charts))
#         for (q,(bcell,bchart)) in enumerate(zip(trial_cells, trial_charts))

#             fill!(zlocal, 0)
#             qrule = quadrule(biop, test_shapes, trial_shapes, p, tchart, q, bchart, qd, quadstrat)
#             momintegrals_trial_refines_test!(zlocal, biop,
#                 test_functions, tcell, tchart,
#                 trial_functions, bcell, bchart,
#                 qrule, quadstrat)
                
#             I = length(test_assembly_data[p])
#             J = length(trial_assembly_data[q])
#             for j in 1 : J, i in 1 : I
#                 zij = zlocal[i,j]
#                 for (n,b) in trial_assembly_data[q][j]
#                     zb = zij*b
#                     for (m,a) in test_assembly_data[p][i]
#                         store(a*zb, m, n)
#         end end end end

#         done += 1
#         new_pctg = round(Int, done / todo * 100)
#         if new_pctg > pctg + 9
#             myid == 1 && print(".")
#             pctg = new_pctg
#         end end
#     myid == 1 && println("")
# end



function blockassembler(biop::IntegralOperator, tfs::Space, bfs::Space;
        quadstrat=defaultquadstrat(biop, tfs, bfs))

    tgeo = geometry(tfs)
    bgeo = geometry(bfs)

    qs = if CompScienceMeshes.refines(tgeo, bgeo)
        TestRefinesTrialQStrat(quadstrat)
    elseif CompScienceMeshes.refines(bgeo, tgeo)
        TrialRefinesTestQStrat(quadstrat)
    else
        quadstrat
    end

    test_elements, test_assembly_data,
        trial_elements, trial_assembly_data,
        quadrature_data, zlocals = assembleblock_primer(biop, tfs, bfs; quadstrat=qs)

    return (test_ids, trial_ids, store) ->
        assembleblock_body!(biop,
            tfs, test_ids,   test_elements,  test_assembly_data,
            bfs, trial_ids, trial_elements, trial_assembly_data,
            quadrature_data, zlocals, store; quadstrat=qs)

    # if CompScienceMeshes.refines(tgeo, bgeo)
    #     return (test_ids, trial_ids, store) -> begin
    #         assembleblock_body_test_refines_trial!(biop,
    #             tfs, test_ids,   test_elements,  test_assembly_data,
    #             bfs, trial_ids, trial_elements, trial_assembly_data,
    #             quadrature_data, zlocals, store; quadstrat)
    #     end
    # elseif CompScienceMeshes.refines(bgeo, tgeo)
    #     return (test_ids, trial_ids, store) -> begin
    #         assembleblock_body_trial_refines_test!(biop,
    #             tfs, test_ids,   test_elements,  test_assembly_data,
    #             bfs, trial_ids, trial_elements, trial_assembly_data,
    #             quadrature_data, zlocals, store; quadstrat)
    #     end
    # else
    #     return (test_ids, trial_ids, store) -> begin
    #         assembleblock_body!(biop,
    #             tfs, test_ids,   test_elements,  test_assembly_data,
    #             bfs, trial_ids, trial_elements, trial_assembly_data,
    #             quadrature_data, zlocals, store; quadstrat)
    #     end
    # end
end


function assembleblock(operator::AbstractOperator, test_functions, trial_functions;
        quadstrat=defaultquadstrat(operator, test_functions, trial_functions))

    Z, store = allocatestorage(operator, test_functions, trial_functions)
    assembleblock!(operator, test_functions, trial_functions, store; quadstrat)

    sdata(Z)
end

function assembleblock!(biop::IntegralOperator, tfs::Space, bfs::Space, store;
        quadstrat=defaultquadstrat(biop, tfs, bfs))

    test_elements, tad, trial_elements, bad, quadrature_data, zlocals =
        assembleblock_primer(biop, tfs, bfs; quadstrat)

    active_test_dofs  = collect(1:numfunctions(tfs))
    active_trial_dofs = collect(1:numfunctions(bfs))

    assembleblock_body!(biop,
        tfs, active_test_dofs, test_elements, tad,
        bfs, active_trial_dofs, trial_elements, bad,
        quadrature_data, zlocals, store; quadstrat)
end


function assembleblock_primer(biop, tfs, bfs;
        quadstrat=defaultquadstrat(biop, tfs, bfs))

    test_elements, tad = assemblydata(tfs; onlyactives=false)
    bsis_elements, bad = assemblydata(bfs; onlyactives=false)

    tgeo = geometry(tfs)
    bgeo = geometry(bfs)

    tdom = domain(chart(tgeo, first(tgeo)))
    bdom = domain(chart(bgeo, first(bgeo)))

    tshapes = refspace(tfs); num_tshapes = numfunctions(tshapes, tdom)
    bshapes = refspace(bfs); num_bshapes = numfunctions(bshapes, bdom)

    qd = quaddata(biop, tshapes, bshapes, test_elements, bsis_elements, quadstrat)

    zlocals = Matrix{scalartype(biop, tfs, bfs)}[]

    for i in 1:Threads.nthreads()
        push!(zlocals, zeros(scalartype(biop, tfs, bfs), num_tshapes, num_bshapes))
    end

    return test_elements, tad, bsis_elements, bad, qd, zlocals
end

function assembleblock_body!(biop::IntegralOperator,
        tfs, test_ids, test_elements, test_assembly_data,
        bfs, trial_ids, bsis_elements, trial_assembly_data,
        quadrature_data, zlocals, store; quadstrat)

    test_shapes  = refspace(tfs)
    trial_shapes = refspace(bfs)

    # Enumerate all the active test elements
    active_test_el_ids  = Vector{Int}()
    active_trial_el_ids = Vector{Int}()

    test_id_in_blk  = Dict{Int,Int}()
    trial_id_in_blk = Dict{Int,Int}()

    for (i,m) in enumerate(test_ids);   test_id_in_blk[m] = i; end
    for (i,m) in enumerate(trial_ids); trial_id_in_blk[m] = i; end

    for m in test_ids,  sh in tfs.fns[m]; push!(active_test_el_ids,  sh.cellid); end
    for m in trial_ids, sh in bfs.fns[m]; push!(active_trial_el_ids, sh.cellid); end

    active_test_el_ids = unique!(sort!(active_test_el_ids))
    active_trial_el_ids = unique!(sort!(active_trial_el_ids))

    @assert length(active_test_el_ids) <= length(test_elements)
    @assert length(active_trial_el_ids) <= length(bsis_elements)

    @assert maximum(active_test_el_ids) <= length(test_elements) "$(maximum(active_test_el_ids)), $(length(test_elements))"
    @assert maximum(active_trial_el_ids) <= length(bsis_elements) "$(maximum(active_trial_el_ids)), $(length(bsis_elements))"

    for p in active_test_el_ids
        tcell = test_elements[p]
        for q in active_trial_el_ids
            bcell = bsis_elements[q]

            fill!(zlocals[Threads.threadid()], 0)
            qrule = quadrule(biop, test_shapes, trial_shapes, p, tcell, q, bcell, quadrature_data, quadstrat)
            momintegrals!(zlocals[Threads.threadid()], biop,
                tfs, p, tcell,
                bfs, q, bcell, qrule)

            for j in 1 : size(zlocals[Threads.threadid()],2)
                for i in 1 : size(zlocals[Threads.threadid()],1)
                    for (n,b) in trial_assembly_data[q,j]
                        n′ = get(trial_id_in_blk, n, 0)
                        n′ == 0 && continue
                        for (m,a) in test_assembly_data[p,i]
                            m′ = get(test_id_in_blk, m, 0)
                            m′ == 0 && continue
                            store(a*zlocals[Threads.threadid()][i,j]*b, m′, n′)
end end end end end end end

# function assembleblock_body_trial_refines_test!(biop::IntegralOperator,
#         tfs, test_ids, test_elements, test_assembly_data,
#         bfs, trial_ids, bsis_elements, trial_assembly_data,
#         quadrature_data, zlocals, store; quadstrat)

#     test_shapes  = refspace(tfs)
#     trial_shapes = refspace(bfs)

#     # Enumerate all the active test elements
#     active_test_el_ids  = Vector{Int}()
#     active_trial_el_ids = Vector{Int}()

#     test_id_in_blk  = Dict{Int,Int}()
#     trial_id_in_blk = Dict{Int,Int}()

#     for (i,m) in enumerate(test_ids);   test_id_in_blk[m] = i; end
#     for (i,m) in enumerate(trial_ids); trial_id_in_blk[m] = i; end

#     for m in test_ids,  sh in tfs.fns[m]; push!(active_test_el_ids,  sh.cellid); end
#     for m in trial_ids, sh in bfs.fns[m]; push!(active_trial_el_ids, sh.cellid); end

#     active_test_el_ids = unique!(sort!(active_test_el_ids))
#     active_trial_el_ids = unique!(sort!(active_trial_el_ids))

#     @assert length(active_test_el_ids) <= length(test_elements)
#     @assert length(active_trial_el_ids) <= length(bsis_elements)

#     @assert maximum(active_test_el_ids) <= length(test_elements) "$(maximum(active_test_el_ids)), $(length(test_elements))"
#     @assert maximum(active_trial_el_ids) <= length(bsis_elements) "$(maximum(active_trial_el_ids)), $(length(bsis_elements))"

#     for p in active_test_el_ids
#         tcell = test_elements[p]
#         for q in active_trial_el_ids
#             bcell = bsis_elements[q]

#             fill!(zlocals[Threads.threadid()], 0)
#             qrule = quadrule(biop, test_shapes, trial_shapes, p, tcell, q, bcell, quadrature_data, quadstrat)
#             momintegrals_trial_refines_test!(zlocals[Threads.threadid()], biop,
#                 tfs, p, tcell,
#                 bfs, q, bcell,
#                 qrule, quadstrat)
#             for j in 1 : size(zlocals[Threads.threadid()],2)
#                 for i in 1 : size(zlocals[Threads.threadid()],1)
#                     for (n,b) in trial_assembly_data[q,j]
#                         n′ = get(trial_id_in_blk, n, 0)
#                         n′ == 0 && continue
#                         for (m,a) in test_assembly_data[p,i]
#                             m′ = get(test_id_in_blk, m, 0)
#                             m′ == 0 && continue
#                             store(a*zlocals[Threads.threadid()][i,j]*b, m′, n′)
# end end end end end end end

# function assembleblock_body_test_refines_trial!(biop::IntegralOperator,
#     tfs, test_ids, test_elements, test_assembly_data,
#     bfs, trial_ids, bsis_elements, trial_assembly_data,
#     quadrature_data, zlocals, store; quadstrat)

#     test_shapes  = refspace(tfs)
#     trial_shapes = refspace(bfs)

#     # Enumerate all the active test elements
#     active_test_el_ids  = Vector{Int}()
#     active_trial_el_ids = Vector{Int}()

#     test_id_in_blk  = Dict{Int,Int}()
#     trial_id_in_blk = Dict{Int,Int}()

#     for (i,m) in enumerate(test_ids);   test_id_in_blk[m] = i; end
#     for (i,m) in enumerate(trial_ids); trial_id_in_blk[m] = i; end

#     for m in test_ids,  sh in tfs.fns[m]; push!(active_test_el_ids,  sh.cellid); end
#     for m in trial_ids, sh in bfs.fns[m]; push!(active_trial_el_ids, sh.cellid); end

#     active_test_el_ids = unique!(sort!(active_test_el_ids))
#     active_trial_el_ids = unique!(sort!(active_trial_el_ids))

#     @assert length(active_test_el_ids) <= length(test_elements)
#     @assert length(active_trial_el_ids) <= length(bsis_elements)

#     @assert maximum(active_test_el_ids) <= length(test_elements) "$(maximum(active_test_el_ids)), $(length(test_elements))"
#     @assert maximum(active_trial_el_ids) <= length(bsis_elements) "$(maximum(active_trial_el_ids)), $(length(bsis_elements))"

#     for p in active_test_el_ids
#         tcell = test_elements[p]
#         for q in active_trial_el_ids
#             bcell = bsis_elements[q]

#             fill!(zlocals[Threads.threadid()], 0)
#             qrule = quadrule(biop, test_shapes, trial_shapes, p, tcell, q, bcell, quadrature_data, quadstrat)
#             momintegrals_test_refines_trial!(zlocals[Threads.threadid()], biop,
#                 tfs, p, tcell,
#                 bfs, q, bcell,
#                 qrule, quadstrat)
#             for j in 1 : size(zlocals[Threads.threadid()],2)
#                 for i in 1 : size(zlocals[Threads.threadid()],1)
#                     for (n,b) in trial_assembly_data[q,j]
#                         n′ = get(trial_id_in_blk, n, 0)
#                         n′ == 0 && continue
#                         for (m,a) in test_assembly_data[p,i]
#                             m′ = get(test_id_in_blk, m, 0)
#                             m′ == 0 && continue
#                             store(a*zlocals[Threads.threadid()][i,j]*b, m′, n′)
# end end end end end end end


# function assembleblock_body_nested!(biop::IntegralOperator,
#     tfs, test_ids, test_elements, test_assembly_data,
#     bfs, trial_ids, bsis_elements, trial_assembly_data,
#     quadrature_data, zlocals, store; quadstrat)

#     test_shapes  = refspace(tfs)
#     trial_shapes = refspace(bfs)

#     # Enumerate all the active test elements
#     active_test_el_ids  = Vector{Int}()
#     active_trial_el_ids = Vector{Int}()

#     test_id_in_blk  = Dict{Int,Int}()
#     trial_id_in_blk = Dict{Int,Int}()

#     for (i,m) in enumerate(test_ids);   test_id_in_blk[m] = i; end
#     for (i,m) in enumerate(trial_ids); trial_id_in_blk[m] = i; end

#     for m in test_ids,  sh in tfs.fns[m]; push!(active_test_el_ids,  sh.cellid); end
#     for m in trial_ids, sh in bfs.fns[m]; push!(active_trial_el_ids, sh.cellid); end

#     active_test_el_ids = unique(sort(active_test_el_ids))
#     active_trial_el_ids = unique(sort(active_trial_el_ids))

#     for p in active_test_el_ids
#         tcell = test_elements[p]
#         for q in active_trial_el_ids
#             bcell = bsis_elements[q]

#             fill!(zlocals[Threads.threadid()], 0)
#             qrule = quadrule(biop, test_shapes, trial_shapes, p, tcell, q, bcell, quadrature_data, quadstrat)
#             momintegrals_test_refines_trial!(zlocals[Threads.threadid()], biop,
#                 tfs, p, tcell,
#                 bfs, q, bcell,
#                 qrule, quadstrat)
#             # momintegrals_test_refines_trial!(biop, test_shapes, trial_shapes, tcell, bcell, zlocals[Threads.threadid()], qrule, quadstrat)

#             for j in 1 : size(zlocals[Threads.threadid()],2)
#                 for i in 1 : size(zlocals[Threads.threadid()],1)
#                     for (n,b) in trial_assembly_data[q,j]
#                         n′ = get(trial_id_in_blk, n, 0)
#                         n′ == 0 && continue
#                         for (m,a) in test_assembly_data[p,i]
#                             m′ = get(test_id_in_blk, m, 0)
#                             m′ == 0 && continue
#                             store(a*zlocals[Threads.threadid()][i,j]*b, m′, n′)
# end end end end end end end


function assemblerow!(biop::IntegralOperator, test_functions::Space, trial_functions::Space, store;
        quadstrat=defaultquadstrat(biop, test_functions, trial_functions))

    tgeo = geometry(test_functions)
    bgeo = geometry(trial_functions)

    tdom = domain(chart(tgeo, first(tgeo)))
    bdom = domain(chart(bgeo, first(bgeo)))

    test_elements = elements(tgeo)
    trial_elements, trial_assembly_data = assemblydata(trial_functions)

    test_shapes  = refspace(test_functions)
    trial_shapes = refspace(trial_functions)

    num_test_shapes  = numfunctions(test_shapes, tdom)
    num_trial_shapes = numfunctions(trial_shapes, bdom)

    quadrature_data = quaddata(biop, test_shapes, trial_shapes, test_elements, trial_elements,
        quadstrat)
    zlocal = zeros(scalartype(biop, test_functions, trial_functions),
        num_test_shapes, num_trial_shapes)

    @assert length(trial_elements) == numcells(geometry(trial_functions))
    @assert numfunctions(test_functions) == 1

    assemblerow_body!(biop,
        test_functions, test_elements, test_shapes,
        trial_assembly_data, trial_functions, trial_elements, trial_shapes,
        zlocal, quadrature_data, store; quadstrat)
end


function assemblerow_body!(biop,
    test_functions, test_elements, test_shapes,
    trial_assembly_data, trial_functions, trial_elements, trial_shapes,
    zlocal, quadrature_data, store; quadstrat)

    test_function = test_functions.fns[1]
    for shape in test_function
        p = shape.cellid
        i = shape.refid
        a = shape.coeff
        tcell = test_elements[p]
        for (q,bcell) in enumerate(trial_elements)

            fill!(zlocal, 0)
            qrule = quadrule(biop, test_shapes, trial_shapes, p, tcell, q, bcell, quadrature_data, quadstrat)
            momintegrals!(zlocal, biop,
                test_functions, nothing, tcell,
                trial_functions, nothing, bcell,
                qrule)

            for j in 1:size(zlocal,2)
                for (n,b) in trial_assembly_data[q,j]
                    store(a*zlocal[i,j]*b, 1, n)
end end end end end


function assemblecol!(biop::IntegralOperator, test_functions::Space, trial_functions::Space, store;
        quadstrat=defaultquadstrat(biop, test_functions, trial_functions))

    test_elements, test_assembly_data = assemblydata(test_functions)
    trial_elements = elements(geometry(trial_functions))

    test_shapes  = refspace(test_functions)
    trial_shapes = refspace(trial_functions)

    num_test_shapes  = numfunctions(test_shapes)
    num_trial_shapes = numfunctions(trial_shapes)

    quadrature_data = quaddata(biop, test_shapes, trial_shapes, test_elements, trial_elements, quadstrat)
    zlocal = zeros(
        scalartype(biop, test_functions, trial_functions),
        num_test_shapes, num_trial_shapes)

    @assert length(test_elements) == numcells(geometry(test_functions))
    @assert numfunctions(trial_functions) == 1

    assemblecol_body!(biop,
        test_assembly_data, test_functions, test_elements,  test_shapes,
        trial_functions,   trial_elements, trial_shapes,
        zlocal, quadrature_data, store; quadstrat)
end


function assemblecol_body!(biop,
    test_assembly_data, test_functions, test_elements, test_shapes,
    trial_functions, trial_elements, trial_shapes,
    zlocal, quadrature_data, store; quadstrat)

    trial_function = trial_functions.fns[1]
    for shape in trial_function
        q = shape.cellid
        j = shape.refid
        b = shape.coeff

        bcell = trial_elements[q]
        for (p,tcell) in enumerate(test_elements)

            fill!(zlocal, 0)
            qrule = quadrule(biop, test_shapes, trial_shapes, p, tcell, q, bcell, quadrature_data, quadstrat)
            momintegrals!(zlocal, biop,
                test_functions, nothing, tcell,
                trial_functions, nothing, bcell, qrule)

            for i in 1:size(zlocal,1)
                for (m,a) in test_assembly_data[p,i]
                    store(a*zlocal[i,j]*b, m, 1)
end end end end end



