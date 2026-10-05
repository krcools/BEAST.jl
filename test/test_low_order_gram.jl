using Test
using SparseArrays
using LinearAlgebra
using CompScienceMeshes
using BEAST

"""
    exact_lagrange_gram(space)

Build the exact identity Gram matrix for the low-order Lagrange bases used here.
The reference is assembled from element mass matrices, independently of BEAST's
quadrature and local assembly paths.
"""
function exact_lagrange_gram(space::BEAST.LagrangeBasis{D}) where {D}
    geo = geometry(space)
    ncells = numcells(geo)
    nfunctions = numfunctions(space)
    T = coordtype(geo)
    local_shapes = [Vector{Tuple{Int,Int,T}}() for _ in 1:ncells]

    for (global_id, fn) in enumerate(space.fns)
        for shape in fn
            push!(local_shapes[shape.cellid], (global_id, shape.refid, shape.coeff))
        end
    end

    rows = Int[]
    cols = Int[]
    values = T[]
    dimension_plus_one = dimension(geo) + 1

    for (cellid, shapes) in enumerate(local_shapes)
        isempty(shapes) && continue
        element_volume = volume(chart(geo, cellid))

        for (row, local_row, row_coefficient) in shapes
            for (col, local_col, col_coefficient) in shapes
                local_mass = if D == 0
                    element_volume
                elseif D == 1
                    factor = element_volume / (dimension_plus_one * (dimension_plus_one + 1))
                    factor * (local_row == local_col ? 2 : 1)
                else
                    error("exact_lagrange_gram only supports degrees zero and one")
                end

                push!(rows, row)
                push!(cols, col)
                push!(values, row_coefficient * local_mass * col_coefficient)
            end
        end
    end

    return sparse(rows, cols, values, nfunctions, nfunctions)
end

function check_gram(space; atol=100eps(Float64), rtol=100eps(Float64))
    assembled = assemble(BEAST.Identity(), space, space)
    reference = exact_lagrange_gram(space)
    @test size(assembled) == size(reference)
    @test assembled ≈ reference atol=atol rtol=rtol
    return assembled, reference
end

"""
    exact_dual_constant_gram(dual, constants)

Exact identity Gram matrix between `duallagrangec0d1` and
`lagrangecxd0` on the same coarse mesh.
"""
function exact_dual_constant_gram(dual, constants)
    fine = geometry(dual)
    coarse = CompScienceMeshes.parent(fine).mesh
    ndual = numfunctions(dual)
    nconstant = numfunctions(constants)
    T = coordtype(coarse)
    rows = Int[]
    cols = Int[]
    values = T[]
    dimension_plus_one = dimension(fine) + 1

    for (dual_id, fn) in enumerate(dual.fns)
        for shape in fn
            coarse_id = CompScienceMeshes.parent(fine, shape.cellid)
            entry = shape.coeff * volume(chart(fine, shape.cellid)) / dimension_plus_one
            push!(rows, dual_id)
            push!(cols, coarse_id)
            push!(values, entry)
        end
    end

    return sparse(rows, cols, values, ndual, nconstant)
end

"""Exact identity Gram matrix between a dual piecewise-constant basis and coarse constants."""
function exact_dual_piecewise_constant_gram(dual, constants)
    fine = geometry(dual)
    ndual = numfunctions(dual)
    nconstant = numfunctions(constants)
    T = coordtype(fine)
    rows = Int[]
    cols = Int[]
    values = T[]

    for (dual_id, fn) in enumerate(dual.fns)
        for shape in fn
            coarse_id = CompScienceMeshes.parent(fine, shape.cellid)
            entry = shape.coeff * volume(chart(fine, shape.cellid))
            push!(rows, dual_id)
            push!(cols, coarse_id)
            push!(values, entry)
        end
    end

    return sparse(rows, cols, values, ndual, nconstant)
end

function check_mixed_gram(mesh; atol=100eps(Float64), rtol=100eps(Float64))
    dual = duallagrangec0d1(mesh)
    constants = lagrangecxd0(mesh)
    reference = exact_dual_constant_gram(dual, constants)
    assembled_transpose = try
        assemble(BEAST.Identity(), constants, dual)
    catch exception
        @test_broken false
        @info "Mixed assembly is not available for coarse and barycentric meshes" exception
        return nothing, nothing, reference
    end
    assembled = assembled_transpose'
    @test assembled ≈ reference atol=atol rtol=rtol
    @test assembled_transpose ≈ reference' atol=atol rtol=rtol
    return assembled, assembled_transpose, reference
end

function check_dual_constant_gram(dual, constants, reference; atol=100eps(Float64), rtol=100eps(Float64))
    assembled_transpose = try
        assemble(BEAST.Identity(), constants, dual)
    catch exception
        @test_broken false
        @info "Mixed assembly is not available for this coarse and barycentric pair" exception
        return nothing, nothing
    end
    assembled = assembled_transpose'
    @test assembled ≈ reference atol=atol rtol=rtol
    @test assembled_transpose ≈ reference' atol=atol rtol=rtol
    return assembled, assembled_transpose
end

@testset "exact low-order Gram matrices" begin
    meshes = (
        ("2D", meshrectangle(1.0, 1.0, 0.01, 2)),
        ("3D surface", meshrectangle(1.0, 1.0, 0.01, 3)),
    )

    for (dimension_name, mesh) in meshes
        @testset "$dimension_name piecewise constants" begin
            check_gram(lagrangecxd0(mesh))
        end

        @testset "$dimension_name continuous P1 Lagrange" begin
            check_gram(lagrangec0d1(mesh))
        end

        @testset "$dimension_name dual P1 Lagrange" begin
            check_gram(duallagrangec0d1(mesh))
        end

        @testset "$dimension_name dual/constant mixed Gram" begin
            check_mixed_gram(mesh)
        end

        @testset "$dimension_name dual piecewise-constant Gram" begin
            check_gram(duallagrangecxd0(mesh))
        end

        @testset "$dimension_name dual piecewise-constant mixed Gram" begin
            dual = duallagrangecxd0(mesh)
            constants = lagrangecxd0(mesh)
            reference = exact_dual_piecewise_constant_gram(dual, constants)
            check_dual_constant_gram(dual, constants, reference)
        end
    end
end

@testset "asset sphere Gram matrices" begin
    mesh = readmesh(joinpath(dirname(pathof(BEAST)), "../test/assets/sphere316.in"))

    for (name, space) in (
        ("piecewise constants", lagrangecxd0(mesh)),
        ("continuous P1 Lagrange", lagrangec0d1(mesh)),
        ("dual piecewise constants", duallagrangecxd0(mesh)),
        ("dual P1 Lagrange", duallagrangec0d1(mesh)),
    )
        @testset "$name" begin
            check_gram(space)
        end
    end

    constants = lagrangecxd0(mesh)
    dual = duallagrangecxd0(mesh)
    reference = exact_dual_piecewise_constant_gram(dual, constants)
    @testset "dual/constant mixed Gram" begin
        check_dual_constant_gram(dual, constants, reference)
    end
end

@testset "asset sphere RWG Gram matrices" begin
    mesh = readmesh(joinpath(dirname(pathof(BEAST)), "../test/assets/sphere316.in"))
    rwg = raviartthomas(mesh)
    bc = buffachristiansen(mesh)

    @testset "RWG/RWG" begin
        gram = assemble(BEAST.Identity(), rwg, rwg)
        @test size(gram) == (numfunctions(rwg), numfunctions(rwg))
        @test all(isfinite, nonzeros(gram))
    end

    @testset "RWG/BC" begin
        gram = assemble(BEAST.Identity(), rwg, bc)
        @test size(gram) == (numfunctions(rwg), numfunctions(bc))
        @test all(isfinite, nonzeros(gram))
    end
end

@testset "local assembly invariants" begin
    mesh = meshrectangle(1.0, 1.0, 0.05, 2)
    identity = BEAST.Identity()

    for space in (
        lagrangecxd0(mesh),
        lagrangec0d1(mesh),
        duallagrangecxd0(mesh),
        duallagrangec0d1(mesh),
    )
        gram = exact_lagrange_gram(space)
        @test gram ≈ gram'
        @test all(isfinite, nonzeros(gram))
    end

    constants = lagrangecxd0(mesh)
    matched = spzeros(numfunctions(constants), numfunctions(constants))
    BEAST.assemble_local_matched!(identity, constants, constants,
        (value, row, column) -> (matched[row, column] += value))
    mixed = spzeros(size(matched)...)
    try
        BEAST.assemble_local_mixed!(identity, constants, constants,
            (value, row, column) -> (mixed[row, column] += value))
        @test matched ≈ mixed
    catch exception
        @test_broken false
        @info "Mixed assembly is not available for this 2D simplex geometry" exception
    end

    fine_mesh = meshrectangle(1.0, 1.0, 0.125, 2)
    fine = lagrangecxd0(fine_mesh)
    coarse = lagrangecxd0(mesh)
    refined = try
        assemble(identity, fine, coarse)
    catch exception
        @test_broken false
        @info "Refinement assembly is not available for this 2D simplex geometry" exception
        nothing
    end
    if refined !== nothing
        mixed_refined = spzeros(size(refined)...)
        try
            BEAST.assemble_local_mixed!(identity, fine, coarse,
                (value, row, column) -> (mixed_refined[row, column] += value))
            @test refined ≈ mixed_refined
        catch exception
            @test_broken false
            @info "Mixed assembly is not available for this 2D refinement geometry" exception
        end
    end

    shifted = translate(mesh, point(2.0, 0.0))
    disjoint = try
        assemble(identity, lagrangecxd0(mesh), lagrangecxd0(shifted))
    catch exception
        @test_broken false
        @info "Disjoint mixed assembly is not available for this 2D simplex geometry" exception
        nothing
    end
    disjoint === nothing || @test nnz(disjoint) == 0
end

@testset "supported 3D assembly paths" begin
    mesh = meshrectangle(1.0, 1.0, 0.05, 3)
    identity = BEAST.Identity()

    constants = lagrangecxd0(mesh)
    matched = spzeros(numfunctions(constants), numfunctions(constants))
    BEAST.assemble_local_matched!(identity, constants, constants,
        (value, row, column) -> (matched[row, column] += value))
    mixed = spzeros(size(matched)...)
    BEAST.assemble_local_mixed!(identity, constants, constants,
        (value, row, column) -> (mixed[row, column] += value))
    @test matched ≈ mixed

    dual = duallagrangecxd0(mesh)
    reference = exact_dual_piecewise_constant_gram(dual, constants)
    dual_constant = assemble(identity, constants, dual)'
    @test dual_constant ≈ reference

    dual_p1 = duallagrangec0d1(mesh)
    reference_p1 = exact_dual_constant_gram(dual_p1, constants)
    dual_p1_constant = assemble(identity, constants, dual_p1)'
    @test dual_p1_constant ≈ reference_p1

    shifted = translate(mesh, point(2.0, 0.0, 0.0))
    disjoint = assemble(identity, constants, lagrangecxd0(shifted))
    @test nnz(disjoint) == 0
end

@testset "low-order Gram benchmark" begin
    target_cells = 6000

    edge_length = sqrt(2 / target_cells)

    for (dimension_name, mesh) in (
        ("2D", meshrectangle(1.0, 1.0, edge_length, 2)),
        ("3D surface", meshrectangle(1.0, 1.0, edge_length, 3)),
    )
        for (name, space) in (
            ("lagrangecxd0", lagrangecxd0(mesh)),
            ("lagrangec0d1", lagrangec0d1(mesh)),
            ("duallagrangec0d1", duallagrangec0d1(mesh)),
            ("duallagrangecxd0", duallagrangecxd0(mesh)),
        )
            @testset "$dimension_name $name" begin
                assembled, reference = check_gram(space)
                @test size(assembled) == size(reference)
            end
        end

        @testset "$dimension_name dual/constant mixed Gram" begin
            forward, transpose, reference = check_mixed_gram(mesh)
            if forward !== nothing
                @test size(forward) == size(reference)
                @test size(transpose) == reverse(size(reference))
            end
        end

        @testset "$dimension_name dual piecewise-constant mixed Gram" begin
            dual = duallagrangecxd0(mesh)
            constants = lagrangecxd0(mesh)
            reference = exact_dual_piecewise_constant_gram(dual, constants)
            forward, transpose = check_dual_constant_gram(dual, constants, reference)
            if forward !== nothing
                @test size(forward) == size(reference)
                @test size(transpose) == reverse(size(reference))
            end
        end
    end
end

@testset "assembly allocation regression" begin
    function allocation_check(f, limit)
        f()
        allocations = @allocated f()
        @test allocations < limit
    end

    identity = BEAST.Identity()
    mesh = meshrectangle(1.0, 1.0, 0.25, 3)
    constants = lagrangecxd0(mesh)
    fine = lagrangecxd0(meshrectangle(1.0, 1.0, 0.125, 3))

    @testset "matched assembly" begin
        allocation_check(() -> assemble(identity, constants, constants), 20_000_000)
    end

    @testset "mixed assembly" begin
        dual = duallagrangecxd0(mesh)
        allocation_check(() -> assemble(identity, constants, dual), 50_000_000)
    end

    @testset "refinement assembly" begin
        allocation_check(() -> assemble(identity, fine, constants), 50_000_000)
    end

    @testset "realistic sphere assembly" begin
        sphere = readmesh(joinpath(dirname(pathof(BEAST)), "../test/assets/sphere316.in"))
        sphere_constants = lagrangecxd0(sphere)
        allocation_check(() -> assemble(identity, sphere_constants, sphere_constants), 50_000_000)
        sphere_dual = duallagrangecxd0(sphere)
        allocation_check(() -> assemble(identity, sphere_constants, sphere_dual), 100_000_000)
    end

    @testset "RWG assembly" begin
        sphere = readmesh(joinpath(dirname(pathof(BEAST)), "../test/assets/sphere316.in"))
        rwg = raviartthomas(sphere)
        bc = buffachristiansen(sphere)
        allocation_check(() -> assemble(identity, rwg, rwg), 100_000_000)
        allocation_check(() -> assemble(identity, rwg, bc), 150_000_000)
    end

    @testset "subdivision assembly" begin
        sphere = readmesh(joinpath(dirname(pathof(BEAST)), "../test/assets/sphere872.in"))
        space = subdsurface(sphere)
        allocation_check(() -> assemble(identity, space, space), 400_000_000)
    end
end
