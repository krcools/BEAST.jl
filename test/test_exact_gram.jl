using Test
using SparseArrays
using LinearAlgebra
using CompScienceMeshes
using CollisionDetection
using BEAST

function exact_lagrange_gram(space::BEAST.LagrangeBasis{D}) where {D}
    geo = geometry(space)
    T = coordtype(geo)
    local_shapes = [Vector{Tuple{Int,Int,T}}() for _ in 1:numcells(geo)]
    for (global_id, fn) in enumerate(space.fns), shape in fn
        push!(local_shapes[shape.cellid], (global_id, shape.refid, shape.coeff))
    end

    rows, cols, values = Int[], Int[], T[]
    nlocal = dimension(geo) + 1
    for (cellid, shapes) in enumerate(local_shapes)
        element_volume = volume(chart(geo, cellid))
        for (row, local_row, row_coefficient) in shapes,
            (col, local_col, col_coefficient) in shapes

            local_mass = if D == 0
                element_volume
            elseif D == 1
                factor = element_volume / (nlocal * (nlocal + 1))
                factor * (local_row == local_col ? 2 : 1)
            else
                error("Only degree-zero and degree-one references are supported")
            end
            push!(rows, row)
            push!(cols, col)
            push!(values, row_coefficient * local_mass * col_coefficient)
        end
    end
    sparse(rows, cols, values, numfunctions(space), numfunctions(space))
end

function exact_dual_constant_gram(dual, constants; degree_one)
    fine = geometry(dual)
    T = coordtype(fine)
    rows, cols, values = Int[], Int[], T[]
    divisor = degree_one ? dimension(fine) + 1 : 1
    for (dual_id, fn) in enumerate(dual.fns), shape in fn
        coarse_id = CompScienceMeshes.parent(fine, shape.cellid)
        value = shape.coeff * volume(chart(fine, shape.cellid)) / divisor
        push!(rows, dual_id)
        push!(cols, coarse_id)
        push!(values, value)
    end
    sparse(rows, cols, values, numfunctions(dual), numfunctions(constants))
end

@testset "exact low-order Gram matrices" begin
    for mesh in (meshrectangle(1.0, 1.0, 0.05, 2), meshrectangle(1.0, 1.0, 0.05, 3))
        for space in (
            lagrangecxd0(mesh),
            lagrangec0d1(mesh),
            duallagrangecxd0(mesh),
            duallagrangec0d1(mesh),
        )
            @test assemble(Identity(), space, space) ≈ exact_lagrange_gram(space)
        end

        if universedimension(mesh) == 3
            constants = lagrangecxd0(mesh)
            for (dual, degree_one) in (
                (duallagrangecxd0(mesh), false),
                (duallagrangec0d1(mesh), true),
            )
                reference = exact_dual_constant_gram(dual, constants; degree_one)
                @test assemble(Identity(), constants, dual)' ≈ reference
                @test assemble(Identity(), constants, dual) ≈ reference'
            end
        end
    end
end

@testset "local assembly edge cases" begin
    mesh = meshrectangle(1.0, 1.0, 0.01, 3)
    identity = Identity()
    constants = lagrangecxd0(mesh)

    matched = spzeros(numfunctions(constants), numfunctions(constants))
    BEAST.assemble_local_matched!(identity, constants, constants,
        (value, row, column) -> (matched[row, column] += value))
    mixed = spzeros(size(matched)...)
    BEAST.assemble_local_mixed!(identity, constants, constants,
        (value, row, column) -> (mixed[row, column] += value))
    @test matched ≈ mixed

    shifted = translate(mesh, point(2.0, 0.0, 0.0))
    @test nnz(assemble(identity, constants, lagrangecxd0(shifted))) == 0

    coarse = meshrectangle(1.0, 1.0, 0.05, 3)
    fine = meshrectangle(1.0, 1.0, 0.025, 3)
    fine_space = lagrangecxd0(fine)
    coarse_space = lagrangecxd0(coarse)
    refined = assemble(identity, fine_space, coarse_space)
    mixed_refined = spzeros(size(refined)...)
    BEAST.assemble_local_mixed!(identity, fine_space, coarse_space,
        (value, row, column) -> (mixed_refined[row, column] += value))
    @test refined ≈ mixed_refined
end

@testset "subdivision basis assembly" begin
    mesh = readmesh(joinpath(dirname(@__FILE__), "assets", "sphere872.in"))
    space = subdsurface(mesh)
    matrix = assemble(Identity(), space, space)
    @test matrix ≈ matrix'
    @test all(isfinite, nonzeros(matrix))
end

@testset "nonuniform mesh assembly" begin
    # The final boundary cells have a different size when the edge length
    # does not divide the rectangle dimensions.
    mesh = meshrectangle(1.0, 0.7, 0.01, 3)
    @test assemble(Identity(), lagrangecxd0(mesh), lagrangecxd0(mesh)) ≈
        exact_lagrange_gram(lagrangecxd0(mesh))
    @test assemble(Identity(), lagrangec0d1(mesh), lagrangec0d1(mesh)) ≈
        exact_lagrange_gram(lagrangec0d1(mesh))
end

@testset "local assembly inference" begin
    mesh = meshrectangle(1.0, 1.0, 0.01, 3)
    identity = Identity()
    constants = lagrangecxd0(mesh)
    store = (value, row, column) -> nothing
    elements, _ = assemblydata(constants)

    @test @inferred(BEAST.assemble_local_matched!(identity, constants, constants, store)) === nothing
    @test @inferred(BEAST.assemble_local_mixed!(identity, constants, constants, store)) === nothing
    @test @inferred(BEAST.elementstree(elements)) isa CollisionDetection.Octree
end

@testset "local assembly allocation smoke" begin
    mesh = meshrectangle(1.0, 1.0, 0.02, 3)
    space = lagrangecxd0(mesh)
    assemble(Identity(), space, space)
    allocations = @allocated assemble(Identity(), space, space)
    @test allocations < 10_000_000
end
