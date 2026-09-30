

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

"""
    AssemblyScratch{Z,T,Q}

Reusable local buffers for one block-assembly operation.

An `AssemblyScratch` must be owned by at most one assembly task at a time.
"""
struct AssemblyScratch{Z,T,Q}
    zlocal::Z
    tadjq::T
    qbuffer::Q
end

"""
    AbstractAssemblyScratchSource

Abstract selector for the source of block-assembly scratch buffers.
"""
abstract type AbstractAssemblyScratchSource end

"""
    NewAssemblyScratch()

Select fresh buffers for an assembly operation.

`assemblyscratch` allocates a new `AssemblyScratch` when this source is used.
"""
struct NewAssemblyScratch <: AbstractAssemblyScratchSource end

"""
    StoredAssemblyScratch(scratch)

Select an existing scratch bundle for an assembly operation.

The referenced bundle is returned without allocation. It must not be used by
another task until the current operation has finished.
"""
struct StoredAssemblyScratch{S} <: AbstractAssemblyScratchSource
    scratch::S
end

"""
    AssemblyScratchPool(available, available_condition)

Thread-safe pool of reusable scratch bundles for block assembly.

Each bundle is exclusively borrowed by one task. The condition protects the
available-bundle list and lets tasks wait when all bundles are in use.
"""
struct AssemblyScratchPool{S}
    available::Vector{S}
    available_condition::Threads.Condition
end

"""
    acquirescratch(pool::AssemblyScratchPool)

Acquire an available scratch bundle from `pool`.

The call waits if every bundle is currently in use. The caller owns the
returned bundle until it passes it to [`releasescratch!`](@ref).
"""
function acquirescratch(pool::AssemblyScratchPool)
    lock(pool.available_condition)
    try
        while isempty(pool.available)
            wait(pool.available_condition)
        end
        return pop!(pool.available)
    finally
        unlock(pool.available_condition)
    end
end

"""
    releasescratch!(pool::AssemblyScratchPool, scratch)

Return `scratch` to `pool` and notify one waiting task.

Every successful call to [`acquirescratch`](@ref) must be paired with one
call to this function, including when assembly exits through an exception.
"""
function releasescratch!(pool::AssemblyScratchPool, scratch)
    lock(pool.available_condition)
    try
        push!(pool.available, scratch)
        notify(pool.available_condition)
    finally
        unlock(pool.available_condition)
    end
    return nothing
end

"""
    assemblyscratch(
        ::NewAssemblyScratch,
        biop,
        test_space,
        trial_space,
        test_assembly_data,
        trial_assembly_data,
        quadstrat,
    )

Allocate an `AssemblyScratch` for the supplied assembly configuration.

The local matrix, trial-adjoint workspace, and quadrature buffer are sized
from the assembly data, spaces, and quadrature strategy.
"""
function assemblyscratch(
    ::NewAssemblyScratch,
    biop,
    test_space,
    trial_space,
    test_assembly_data,
    trial_assembly_data,
    quadstrat,
)
    return AssemblyScratch(
        zeros(
            scalartype(biop, test_space, trial_space),
            size(test_assembly_data.data, 2),
            size(trial_assembly_data.data, 2),
        ),
        Vector{eltype(trial_assembly_data.data)}(
            undef, size(trial_assembly_data.data, 1)
        ),
        quadraturebuffer(quadstrat, test_space, trial_space),
    )
end

"""
    assemblyscratch(
        source::StoredAssemblyScratch,
        biop,
        test_space,
        trial_space,
        test_assembly_data,
        trial_assembly_data,
        quadstrat,
    )

Return the stored scratch bundle without allocating new buffers.

The caller is responsible for ensuring that the bundle is not used
concurrently by another assembly task.
"""
function assemblyscratch(
    source::StoredAssemblyScratch,
    biop,
    test_space,
    trial_space,
    test_assembly_data,
    trial_assembly_data,
    quadstrat,
)
    return source.scratch
end


# """
#     quadrule(operator,test_refspace,trial_refspace,p,test_element,q_trial_element, qd)

# Returns an object that contains all the dynamic (runtime) information that
# defines the integration strategy that will be used by `integrate!` to compute
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
    integrate!(operator, test_refspace, trial_refspace, test_index, test_chart,
        trial_index, trial_chart, quad_data, quadstrat,
        out, test_space, test_ptr, trial_space, trial_ptr; action::QuadRuleAction)
    integrate!(out, operator, test_space, test_ptr, test_chart,
        trial_space, trial_ptr, trial_chart, qrule)

For `IntegralOperator` assembly, `integrate!` is overloaded twice over (folded
together in BEAST 2.10 from what used to be the separate `quadrule` and
`momintegrals!` functions):

- One family of methods, dispatching on `quadstrat`, builds the quadrature rule
  appropriate for the given pair of elements (the role `quadrule` used to play
  on its own), then either evaluates it into `out` or returns it unevaluated, depending on
  `action` (see [`ApplyIntegrate`](@ref), [`ReturnQRule`](@ref),
  [`ApplyIntegrateNonConforming`](@ref)). Doing this in one step, rather than
  returning the rule to a separate caller for it to dispatch on, avoids a
  dynamic dispatch on the wide union of rule types a given `quadstrat` can
  produce.
- The other family of methods, dispatching on the concrete type of an
  already-built `qrule` (the role `momintegrals!` used to play), performs the
  actual numerical integration into `out`.
"""
function integrate! end

# Called by the `integrate!` methods above that dispatch on `quadstrat`, once one
# of them has built a concrete `qrule`, from within the same branch/method that
# constructed it, so the type of `qrule` is still known to the compiler and this
# dispatch resolves statically instead of at runtime.
function integrate!(action::ApplyIntegrate, out, op, test_space, tptr, tcell, trial_space, bptr, bcell, qrule)
    return integrate!(
        out, op, test_space, tptr, tcell, trial_space, bptr, bcell, qrule, action.qbuffer
    )
end

integrate!(::ReturnQRule, out, op, test_space, tptr, tcell, trial_space, bptr, bcell, qrule) = qrule

# Same idea as the ApplyIntegrate case above, but for the non-conforming-mesh rules,
# whose `qrule` (test_space/trial_space here are local refspaces, not a full Space)
# is reached only after the generic dispatcher in momintegrals.jl has already
# stripped Space down to refspace. Skips straight to the refspace-level `integrate!`
# instead of trying to re-derive a Space that isn't there to begin with.
function integrate!(action::ApplyIntegrateNonConforming, out, op, test_space, tptr, tcell, trial_space, bptr, bcell, qrule)
    return integrate!(
        op, test_space, trial_space, tcell, bcell, out, qrule, action.qbuffer
    )
end

function integrate!(op, test_local_space, trial_local_space, test_chart, trial_chart, out, qrule, qbuffer)
    return integrate!(op, test_local_space, trial_local_space, test_chart, trial_chart, out, qrule)
end


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

    numfunctions(tfs) == 0 && return
    numfunctions(bfs) == 0 && return

    tr = assemblydata(tfs); tr == nothing && return
    br = assemblydata(bfs); br == nothing && return

    test_elements, tad, tcells = tr
    bsis_elements, bad, bcells = br

    tgeo = geometry(tfs)
    bgeo = geometry(bfs)

    # tdom = domain(chart(tgeo, first(tgeo)))
    # bdom = domain(chart(bgeo, first(bgeo)))

    tshapes = refspace(tfs); #num_tshapes = numfunctions(tshapes, tdom)
    bshapes = refspace(bfs); #num_bshapes = numfunctions(bshapes, bdom)

    qs = if CompScienceMeshes.refines(tgeo, bgeo)
        TestRefinesTrialQStrat(quadstrat)
    elseif CompScienceMeshes.refines(bgeo, tgeo)
        TrialRefinesTestQStrat(quadstrat)
    else
        quadstrat
    end

    qd = quaddata(biop, tshapes, bshapes, test_elements, bsis_elements, qs)
    # @show "after" qs
    # assemblechunk_body!(biop,
    #     tfs, test_elements, tad, tcells,
    #     bfs, bsis_elements, bad, bcells,
    #     qd, zlocal, store; quadstrat=qs)

    # @show length(qd.tpo/ints), length(qd.bpoints), maximum(tcells), maximum(bcells)

    assemblechunk_body!(biop, tfs, bfs,
        test_elements, tcells, tad, eachindex(tcells),
        bsis_elements, bcells, bad, eachindex(bcells),
        qd, store; quadstrat=qs)

end

@testitem "assemble!: zero sized block" begin
    using CompScienceMeshes
    import BEAST.BlockArrays

    fn = joinpath(dirname(pathof(BEAST)), "../examples/assets/sphere45.in")
    m1 = readmesh(fn)
    m2 = m1[Int[]]

    X = BEAST.DirectProductSpace([raviartthomas(m) for m in [m1, m2]])
    T = Maxwell3D.singlelayer(gamma=1.0)

    @hilbertspace j[1:2]
    @hilbertspace k[1:2]
    a = T[k[1],j[1]] + T[k[1],j[2]] + T[k[2],j[2]] + T[k[2],j[1]]

    A = assemble(a, X, X)
    M = AbstractMatrix(A)

    n1 = numfunctions(X[1])
    n2 = numfunctions(X[2])

    @test n2 == 0

    @test BlockArrays.blocksize(M) == (2,2)
    @test BlockArrays.blocksizes(M) == [(n1,n1) (n1,n2); (n2,n1) (n2,n2)]
end


function assemblechunk_body!(biop, test_space, trial_space,
    test_elements, test_element_ptrs, test_assembly_data, active_test_els,
    trial_elements, trial_element_ptrs, trial_assembly_data, active_trial_els,
    qd, store;
    quadstrat,
    scheduler=:serial,
    scratchsource=NewAssemblyScratch())

    num_tshapes = size(test_assembly_data.data,2)
    num_bshapes = size(trial_assembly_data.data,2)

    @tasks for p in eachindex(active_test_els)
        @set scheduler = scheduler
        @local scratch = assemblyscratch(
            scratchsource,
            biop,
            test_space,
            trial_space,
            test_assembly_data,
            trial_assembly_data,
            quadstrat,
        )
        zlocal = scratch.zlocal
        tadjq = scratch.tadjq
        qaction = ApplyIntegrate(scratch.qbuffer)
        P = active_test_els[p]
        tcell = test_elements[P]
        tptr = test_element_ptrs[P]

        for (q,Q) in enumerate(active_trial_els)
            bcell = trial_elements[Q]
            bptr = trial_element_ptrs[Q]

            fill!(zlocal, 0)

            integrate!(biop, refspace(test_space), refspace(trial_space),
                P, tcell, Q, bcell, qd, quadstrat,
                zlocal, test_space, tptr, trial_space, bptr; action=qaction)
            for j in 1 : num_bshapes
                tadjq .= @view trial_assembly_data.data[:,j,q]
                for i in 1 : num_tshapes
                    zij = zlocal[i,j]
                    badip = @view test_assembly_data.data[:,i,p]
                    for (n,b) in tadjq
                        (n < 1 || iszero(b)) && continue
                        zb = zij*b
                        for (m,a) in badip
                            (m < 1 || iszero(a)) && continue
                            store(a*zb, m, n)
end end end end end end end


# function assemblechunk_body_colored!(biop,
#         test_space, testad, testelementids,
#         trial_space, trialad, trialelementids,
#         qd, store, scheduler; quadstrat)

#     test_shapes = refspace(test_space)
#     trial_shapes = refspace(trial_space)

#     test_elements, test_assembly_data, test_cell_ptrs = testad
#     trial_elements, trial_assembly_data, trial_cell_ptrs = trialad

#     num_tshapes = numfunctions(test_shapes, domain(chart(geometry(test_space), first(geometry(test_space)))))
#     num_bshapes = numfunctions(trial_shapes, domain(chart(geometry(trial_space), first(geometry(trial_space)))))

#     @tasks for p in testelementids
#         @set scheduler = scheduler
#         @local zlocal = zeros(scalartype(biop, test_space, trial_space), num_tshapes, num_bshapes)
#         tcell, tptr = test_elements[p], test_cell_ptrs[p]

#         for q in trialelementids
#             bcell, bptr = trial_elements[q], trial_cell_ptrs[q]
#             fill!(zlocal, 0)
#             @inline qrule = quadrule(biop, test_shapes, trial_shapes, p, tcell, q, bcell, qd, quadstrat)
#             integrate!(zlocal, biop,
#                 test_space,  tptr, tcell,
#                 trial_space, bptr, bcell, qrule)
#             for j in 1:length(trial_assembly_data[q])
#                 QQ = @view trial_assembly_data.data[:,j,q]
#                 for i in 1:length(test_assembly_data[p])
#                 zij = zlocal[i,j]
#                 PP = @view test_assembly_data.data[:,i,p]
#                 for (n,b) in QQ
#                     n < 1 && break
#                     iszero(b) && continue
#                     zb = zij*b
#                     for (m,a) in PP
#                         m < 1 && break
#                         iszero(a) && continue
#                         store(a*zb, m, n)
# end end end end end end end

struct AssembleblockbodyFunctor{B,T1,T2,T3,T4,T5,T6,T7,T8,T9}
    biop::B
    tfs::T1
    testelements::T2
    testassemblydata::T3
    bfs::T4
    trialelements::T5
    trialassemblydata::T6
    quadraturedata::T7
    scratchpool::T8
    quadstrat::T9
end

function active_element_ids(space, ids)
    element_ids = collect(Int,
        shape.cellid for id in ids for shape in basisfunction(space, id))
    return unique!(sort!(element_ids))
end

function (f::AssembleblockbodyFunctor)(testids, trialids, store)

    tad = f.testassemblydata
    bad = f.trialassemblydata

    # The following code assumed that the testelements cache was build for
    # all the elements in geoemtry(tfs). Similar for the trial side.
    active_test_els = active_element_ids(f.tfs, testids)
    active_trial_els = active_element_ids(f.bfs, trialids)

    tad1 = reduce_assembly_data(tad, testids, active_test_els)
    bad1 = reduce_assembly_data(bad, trialids, active_trial_els)

    test_element_ptrs = eachindex(f.testelements)
    trial_element_ptrs = eachindex(f.trialelements)

    scratch = acquirescratch(f.scratchpool)
    try
        assemblechunk_body!(f.biop, f.tfs, f.bfs,
            f.testelements, test_element_ptrs, tad1, active_test_els,
            f.trialelements, trial_element_ptrs, bad1, active_trial_els,
            f.quadraturedata, store;
            quadstrat=f.quadstrat,
            scheduler=:serial,
            scratchsource=StoredAssemblyScratch(scratch))
    finally
        releasescratch!(f.scratchpool, scratch)
    end
    # assembleblock_body!(f.biop, f.tfs, testids, f.testelements, f.testassemblydata,
    #     f.bfs, trialids, f.trialelements, f.trialassemblydata,
    #     f.quadraturedata, f.zlocals, store; quadstrat=f.quadstrat)
end

function blockassembler(biop::IntegralOperator, tfs::Space, bfs::Space;
        quadstrat=defaultquadstrat(biop, tfs, bfs), primer=nothing)

    tgeo = geometry(tfs)
    bgeo = geometry(bfs)

    qs = if CompScienceMeshes.refines(tgeo, bgeo)
        TestRefinesTrialQStrat(quadstrat)
    elseif CompScienceMeshes.refines(bgeo, tgeo)
        TrialRefinesTestQStrat(quadstrat)
    else
        quadstrat
    end

    if primer === nothing
        primer = assembleblock_primer(biop, tfs, bfs; quadstrat=qs)
    end

    test_elements, test_assembly_data,
        trial_elements, trial_assembly_data,
        quadrature_data, scratchpool = primer

    return AssembleblockbodyFunctor(
        biop,
        tfs,
        test_elements,
        test_assembly_data,
        bfs,
        trial_elements,
        trial_assembly_data,
        quadrature_data,
        scratchpool,
        qs,
    )

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


# function assembleblock(operator::AbstractOperator, test_functions, trial_functions;
#         quadstrat=defaultquadstrat(operator, test_functions, trial_functions))

#     Z, store = allocatestorage(operator, test_functions, trial_functions)
#     assembleblock!(operator, test_functions, trial_functions, store; quadstrat)

#     sdata(Z)
# end

# function assembleblock!(biop::IntegralOperator, tfs::Space, bfs::Space, store;
#         quadstrat=defaultquadstrat(biop, tfs, bfs))

#     test_elements, tad, trial_elements, bad, quadrature_data, zlocals =
#         assembleblock_primer(biop, tfs, bfs; quadstrat)

#     active_test_dofs  = collect(1:numfunctions(tfs))
#     active_trial_dofs = collect(1:numfunctions(bfs))

#     assembleblock_body!(biop,
#         tfs, active_test_dofs, test_elements, tad,
#         bfs, active_trial_dofs, trial_elements, bad,
#         quadrature_data, zlocals, store; quadstrat)
# end


function assembleblock_primer(biop, tfs, bfs;
        quadstrat=defaultquadstrat(biop, tfs, bfs))

    test_elements, tad = assemblydata(tfs; onlyactives=false)
    bsis_elements, bad = assemblydata(bfs; onlyactives=false)

    tshapes = refspace(tfs)
    bshapes = refspace(bfs)

    qd = quaddata(biop, tshapes, bshapes, test_elements, bsis_elements, quadstrat)

    scratches = [
        assemblyscratch(NewAssemblyScratch(), biop, tfs, bfs, tad, bad, quadstrat)
        for _ in 1:(2 * Threads.nthreads())
    ]
    scratchpool = AssemblyScratchPool(scratches, Threads.Condition())

    return test_elements, tad, bsis_elements, bad, qd, scratchpool
end

# function assembleblock_body!(biop::IntegralOperator,
#         tfs, test_ids, test_elements, test_assembly_data,
#         bfs, trial_ids, bsis_elements, trial_assembly_data,
#         quadrature_data, zlocals, store; quadstrat)

#     tgeo = geometry(tfs)
#     bgeo = geometry(bfs)

#     tdom = domain(chart(tgeo, first(tgeo)))
#     bdom = domain(chart(bgeo, first(bgeo)))

#     tshapes = refspace(tfs); num_tshapes = numfunctions(tshapes, tdom)
#     bshapes = refspace(bfs); num_bshapes = numfunctions(bshapes, bdom)

#         # test_shapes  = refspace(tfs)
#         # trial_shapes = refspace(bfs)

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

#     # zlocal = take!(zlocals)
#     zlocal = zeros(scalartype(biop, tfs, bfs), num_tshapes, num_bshapes)
#     for p in active_test_el_ids
#         tcell = test_elements[p]
#         for q in active_trial_el_ids
#             bcell = bsis_elements[q]

#             fill!(zlocal, 0)
#             qrule = quadrule(biop, tshapes, bshapes, p, tcell, q, bcell, quadrature_data, quadstrat)
#             integrate!(zlocal, biop,
#                 tfs, p, tcell,
#                 bfs, q, bcell, qrule)

#             for j in 1 : size(zlocal,2)
#                 for i in 1 : size(zlocal,1)
#                     for (n,b) in trial_assembly_data[q,j]
#                         n′ = get(trial_id_in_blk, n, 0)
#                         n′ == 0 && continue
#                         for (m,a) in test_assembly_data[p,i]
#                             m′ = get(test_id_in_blk, m, 0)
#                             m′ == 0 && continue
#                             store(a*zlocal[i,j]*b, m′, n′)
#     end end end end end end
#     # put!(zlocals, zlocal)
# end

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
    qbuffer = quadraturebuffer(quadstrat, test_functions, trial_functions)
    qaction = ApplyIntegrate(qbuffer)
    for shape in test_function
        p = shape.cellid
        i = shape.refid
        a = shape.coeff
        tcell = test_elements[p]
        for (q,bcell) in enumerate(trial_elements)

            fill!(zlocal, 0)
            integrate!(biop, test_shapes, trial_shapes, p, tcell, q, bcell, quadrature_data, quadstrat,
                zlocal, test_functions, nothing, trial_functions, nothing; action=qaction)

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
    qbuffer = quadraturebuffer(quadstrat, test_functions, trial_functions)
    qaction = ApplyIntegrate(qbuffer)
    for shape in trial_function
        q = shape.cellid
        j = shape.refid
        b = shape.coeff

        bcell = trial_elements[q]
        for (p,tcell) in enumerate(test_elements)

            fill!(zlocal, 0)
            integrate!(biop, test_shapes, trial_shapes, p, tcell, q, bcell, quadrature_data, quadstrat,
                zlocal, test_functions, nothing, trial_functions, nothing; action=qaction)

            for i in 1:size(zlocal,1)
                for (m,a) in test_assembly_data[p,i]
                    store(a*zlocal[i,j]*b, m, 1)
end end end end end
