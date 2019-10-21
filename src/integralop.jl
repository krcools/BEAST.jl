abstract type IntegralOperator <: Operator end



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
    quadrule(operator,test_refspace,trial_refspace,p,test_element,q_trial_element, qd)

Returns an object that contains all the dynamic (runtime) information that
defines the integration strategy that will be used by `momintegrals!` to compute
the interactions between the local test/trial functions defined on the specified
geometric elements. The indices `p` and `q` refer to the position of the test
and trial elements as encountered during iteration over the output of
`geometry`.

The last argument `qd` provides access to all precomputed data required for
quadrature. For example it might be desirable to precompute all the quadrature
points for all possible numerical quadrature schemes that can potentially be
required during matrix assembly. This makes sense, since the number of point is
order N (where N is the number of faces) but these points will appear in N^2
computations. Precomputation requires some extra memory but can save a lot on
computation time.
"""
function quadrule end


"""
  elements(geo)

Create an iterable collection of the elements stored in `geo`. The order in which
this collection produces the elements determines the index used for lookup in the
data structures returned by `assemblydata` and `quaddata`.
"""
elements(geo) = [chart(geo,cl) for cl in cells(geo)]

elements(sp::Space) = elements(geometry(sp))

"""
    assemblechunk!(biop::IntegralOperator, tfs, bfs, store)

Computes the matrix of operator biop wrt the finite element spaces tfs and bfs
"""
function assemblechunk!(biop::IntegralOperator, tfs::Space, bfs::Space, store)

    test_elements, tad = assemblydata(tfs)
    bsis_elements, bad = assemblydata(bfs)

    tshapes = refspace(tfs); num_tshapes = numfunctions(tshapes)
    bshapes = refspace(bfs); num_bshapes = numfunctions(bshapes)

    qd = quaddata(biop, tshapes, bshapes, test_elements, bsis_elements)
    zlocal = zeros(scalartype(biop, tfs, bfs), 2num_tshapes, 2num_bshapes)

    if !CompScienceMeshes.refines(tfs.geo, bfs.geo)
        assemblechunk_body!(biop,
            tshapes, test_elements, tad,
            bshapes, bsis_elements, bad,
            qd, zlocal, store)
    else
        assemblechunk_body_nested_meshes!(biop,
            tshapes, test_elements, tad,
            bshapes, bsis_elements, bad,
            qd, zlocal, store)
    end
end


function assemblechunk_body!(biop,
        test_shapes, test_elements, test_assembly_data,
        trial_shapes, trial_elements, trial_assembly_data,
        qd, zlocal, store)

    print("dots out of 10: ")
    todo, done, pctg = length(test_elements), 0, 0
    for (p,tcell) in enumerate(test_elements)
        for (q,bcell) in enumerate(trial_elements)

        fill!(zlocal, 0)
        strat = quadrule(biop, test_shapes, trial_shapes, p, tcell, q, bcell, qd)
        momintegrals!(biop, test_shapes, trial_shapes, tcell, bcell, zlocal, strat)
        I = length(test_assembly_data[p])
        J = length(trial_assembly_data[q])
        for j in 1 : J, i in 1 : I
            zij = zlocal[i,j]
            for (n,b) in trial_assembly_data[q][j]
                zb = zij*b
                for (m,a) in test_assembly_data[p][i]
                    azb = a*zb
                    store(azb, m, n)
        end end end end

        done += 1
        new_pctg = round(Int, done / todo * 100)
        if new_pctg > pctg + 9
            print(".")
            pctg = new_pctg
    end end
    println("")
end


function assemblechunk_body_nested_meshes!(biop,
        test_shapes, test_elements, test_assembly_data,
        trial_shapes, trial_elements, trial_assembly_data,
        qd, zlocal, store)

    print("dots out of 10: ")
    todo, done, pctg = length(test_elements), 0, 0
    for (p,tcell) in enumerate(test_elements)
        for (q,bcell) in enumerate(trial_elements)

        fill!(zlocal, 0)
        strat = quadrule(biop, test_shapes, trial_shapes, p, tcell, q, bcell, qd)
        momintegrals_nested!(biop, test_shapes, trial_shapes, tcell, bcell, zlocal, strat)
        I = length(test_assembly_data[p])
        J = length(trial_assembly_data[q])
        for j in 1 : J, i in 1 : I
            for (n,b) in trial_assembly_data[q][j], (m,a) in test_assembly_data[p][i]
                store(a*zlocal[i,j]*b, m, n)
        end end end

        done += 1
        new_pctg = round(Int, done / todo * 100)
        if new_pctg > pctg + 9
            print(".")
            pctg = new_pctg
    end end
    println("")
end



function blockassembler(biop::IntegralOperator, tfs::Space, bfs::Space)

    test_elements, test_assembly_data,
        trial_elements, trial_assembly_data,
        quadrature_data, zlocal = assembleblock_primer(biop, tfs, bfs)

    return function f(test_ids, trial_ids, store)
        assembleblock_body!(biop,
            tfs, test_ids,   test_elements,  test_assembly_data,
            bfs, trial_ids, trial_elements, trial_assembly_data,
            quadrature_data, zlocal, store)
    end
end


function assembleblock(operator::AbstractOperator, test_functions, trial_functions)
    Z, store = allocatestorage(operator, test_functions, trial_functions)
    assembleblock!(operator, test_functions, trial_functions, store)
    sdata(Z)
end

function assembleblock!(biop::IntegralOperator, tfs::Space, bfs::Space, store)

    test_elements, tad, trial_elements, bad, quadrature_data, zlocal =
        assembleblock_primer(biop, tfs, bfs)

    active_test_dofs  = collect(1:numfunctions(tfs))
    active_trial_dofs = collect(1:numfunctions(bfs))

    assembleblock_body!(biop,
        tfs, active_test_dofs, test_elements, tad,
        bfs, active_trial_dofs, trial_elements, bad,
        quadrature_data, zlocal, store)
end


function assembleblock_primer(biop, tfs, bfs)

    test_elements, tad = assemblydata(tfs)
    bsis_elements, bad = assemblydata(bfs)

    tshapes = refspace(tfs); num_tshapes = numfunctions(tshapes)
    bshapes = refspace(bfs); num_bshapes = numfunctions(bshapes)

    qd = quaddata(biop, tshapes, bshapes, test_elements, bsis_elements)
    zlocal = zeros(scalartype(biop, tfs, bfs), num_tshapes, num_bshapes)

    return test_elements, tad, bsis_elements, bad, qd, zlocal

end

function assembleblock_body!(biop::IntegralOperator,
        tfs, test_ids, test_elements, test_assembly_data,
        bfs, trial_ids, bsis_elements, trial_assembly_data,
        quadrature_data, zlocal, store)

    test_shapes  = refspace(tfs)
    trial_shapes = refspace(bfs)

    # Enumerate all the active test elements
    active_test_el_ids  = Vector{Int}()
    active_trial_el_ids = Vector{Int}()

    #test_id_in_blk  = zeros(Int, numfunctions(tfs))
    #trial_id_in_blk = zeros(Int, numfunctions(bfs))
    test_id_in_blk  = Dict{Int,Int}()
    trial_id_in_blk = Dict{Int,Int}()

    for (i,m) in enumerate(test_ids);   test_id_in_blk[m] = i; end
    for (i,m) in enumerate(trial_ids); trial_id_in_blk[m] = i; end

    for m in test_ids,  sh in tfs.fns[m]; push!(active_test_el_ids,  sh.cellid); end
    for m in trial_ids, sh in bfs.fns[m]; push!(active_trial_el_ids, sh.cellid); end

    active_test_el_ids = unique(sort(active_test_el_ids))
    active_trial_el_ids = unique(sort(active_trial_el_ids))

    for p in active_test_el_ids
        tcell = test_elements[p]
        for q in active_trial_el_ids
            bcell = bsis_elements[q]

            fill!(zlocal, 0)
            strat = quadrule(biop, test_shapes, trial_shapes, p, tcell, q, bcell, quadrature_data)
            momintegrals!(biop, test_shapes, trial_shapes, tcell, bcell, zlocal, strat)

            for j in 1 : size(zlocal,2)
                for i in 1 : size(zlocal,1)
                    for (n,b) in trial_assembly_data[q,j]
                        n′ = get(trial_id_in_blk, n, 0)
                        n′ == 0 && continue
                        for (m,a) in test_assembly_data[p,i]
                            m′ = get(test_id_in_blk, m, 0)
                            m′ == 0 && continue
                            store(a*zlocal[i,j]*b, m′, n′)
end end end end end end end


function assemblerow!(biop::IntegralOperator, test_functions::Space, trial_functions::Space, store)

    test_elements = elements(geometry(test_functions))
    trial_elements, trial_assembly_data = assemblydata(trial_functions)

    test_shapes  = refspace(test_functions)
    trial_shapes = refspace(trial_functions)

    num_test_shapes  = numfunctions(test_shapes)
    num_trial_shapes = numfunctions(trial_shapes)

    quadrature_data = quaddata(biop, test_shapes, trial_shapes, test_elements, trial_elements)
    zlocal = zeros(scalartype(biop, test_functions, trial_functions),
        num_test_shapes, num_trial_shapes)

    @assert length(trial_elements) == numcells(geometry(trial_functions))
    @assert numfunctions(test_functions) == 1

    assemblerow_body!(biop,
        test_functions, test_elements, test_shapes,
        trial_assembly_data, trial_elements, trial_shapes,
        zlocal, quadrature_data, store)
end


function assemblerow_body!(biop,
    test_functions, test_elements, test_shapes,
    trial_assembly_data, trial_elements, trial_shapes,
    zlocal, quadrature_data, store)

    test_function = test_functions.fns[1]
    for shape in test_function
        p = shape.cellid
        i = shape.refid
        a = shape.coeff
        tcell = test_elements[p]
        for (q,bcell) in enumerate(trial_elements)

            fill!(zlocal, 0)
            strat = quadrule(biop, test_shapes, trial_shapes, p, tcell, q, bcell, quadrature_data)
            momintegrals!(biop, test_shapes, trial_shapes, tcell, bcell, zlocal, strat)

            for j in 1:size(zlocal,2)
                for (n,b) in trial_assembly_data[q,j]
                    store(a*zlocal[i,j]*b, 1, n)
end end end end end


function assemblecol!(biop::IntegralOperator, test_functions::Space, trial_functions::Space, store)

    test_elements, test_assembly_data = assemblydata(test_functions)
    trial_elements = elements(geometry(trial_functions))

    test_shapes  = refspace(test_functions)
    trial_shapes = refspace(trial_functions)

    num_test_shapes  = numfunctions(test_shapes)
    num_trial_shapes = numfunctions(trial_shapes)

    quadrature_data = quaddata(biop, test_shapes, trial_shapes, test_elements, trial_elements)
    zlocal = zeros(
        scalartype(biop, test_functions, trial_functions),
        num_test_shapes, num_trial_shapes)

    @assert length(test_elements) == numcells(geometry(test_functions))
    @assert numfunctions(trial_functions) == 1

    assemblecol_body!(biop,
        test_assembly_data, test_elements,  test_shapes,
        trial_functions,   trial_elements, trial_shapes,
        zlocal, quadrature_data, store)
end


function assemblecol_body!(biop,
    test_assembly_data, test_elements, test_shapes,
    trial_functions, trial_elements, trial_shapes,
    zlocal, quadrature_data, store)

    trial_function = trial_functions.fns[1]
    for shape in trial_function
        q = shape.cellid
        j = shape.refid
        b = shape.coeff

        bcell = trial_elements[q]
        for (p,tcell) in enumerate(test_elements)

            fill!(zlocal, 0)
            strat = quadrule(biop, test_shapes, trial_shapes, p, tcell, q, bcell, quadrature_data)
            momintegrals!(biop, test_shapes, trial_shapes, tcell, bcell, zlocal, strat)

            for i in 1:size(zlocal,1)
                for (m,a) in test_assembly_data[p,i]
                    store(a*zlocal[i,j]*b, m, 1)
end end end end end




# mutable struct SauterSchwabStrategy
#     hits::Int64
# end
# function momintegrals!(biop, tshs::subReferenceSpace, bshs::subReferenceSpace, tcell, bcell, z, strat::SauterSchwabStrategy)
#
#     A = biop.alpha
#     k = biop.gamma
#
#     M, N = size(z)
#
#     hits = strat.hits
#     print(hits)
#     acc = 3
#
#     if hits == 0
#         ssm = PositiveDistance(acc)
#     elseif hits == 1
#         ssm = CommonVertex(acc)
#     elseif hits == 2
#         ssm = CommonEdge(acc)
#     elseif hits == 3
#         ssm = CommonFace(acc)
#     else
#         error("hits can not exceed 3")
#     end
#     M= tcell.N
#     N= bcell.N
#     print("M = $(M) and N = $(N) \n")
#     z = Array{Complex{Float64},2}(M,N)
#     for i = 1:M
#         for j = 1:N
#             # print("i = $(i) and j = $(j) \n")
#             function integrand(u,v)
#                 upt = neighborhood(tcell,u)
#                 vpt = neighborhood(bcell,v)
#                 tshape = shapefuns(upt)
#                 bshape = shapefuns(vpt)
#                 y = cartesian(upt)
#                 x = cartesian(vpt)
#                 kernel = A * exp(-complex(0,1)*k*norm(x-y))/(4.0*π*norm(x-y))
#                 # kernel = kernelvals(biop, upt, vpt)
#                 ujac = jacobian(upt)
#                 vjac = jacobian(vpt)
#                 return tshape[i]*bshape[j]*kernel* ujac * vjac
#                 # return kernel* ujac * vjac
#             end
#             z[i,j] = (sauterschwab_parameterized(tcell,bcell,integrand,ssm))
#             end
#     end
#
#     return z
# end



# mutable struct QuadData{WPV1,WPV2}
#   tpoints::Matrix{Vector{WPV1}}
#   bpoints::Matrix{Vector{WPV2}}
# end
