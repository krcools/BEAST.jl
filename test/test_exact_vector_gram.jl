using Test
using LinearAlgebra
using CompScienceMeshes
using BEAST

function rt_shape_values(ch, refid, coeff, points)
    ref = BEAST.RTRefSpace{coordtype(ch)}()
    reference_points = (
        point(coordtype(ch), 1, 0),
        point(coordtype(ch), 0, 1),
        point(coordtype(ch), 0, 0),
    )
    values = ntuple(length(points)) do i
        value = ref(neighborhood(ch, reference_points[i]))[refid].value
        coeff * value
    end
    return values
end

function rt_shape_value(ch, refid, coeff, x)
    values = rt_shape_values(ch, refid, coeff, ch.vertices)
    barycentric = carttobary(ch, x)
    weights = (barycentric[1], barycentric[2], 1 - sum(barycentric))
    return coeff * zero(values[1]) + sum(weights[i] * values[i] for i in 1:3)
end

function exact_rt_gram(test_space, trial_space)
    test_geometry = geometry(test_space)
    trial_geometry = geometry(trial_space)
    result = zeros(Float64, numfunctions(test_space), numfunctions(trial_space))

    for (trial_cell_id, trial_cell) in enumerate(trial_geometry)
        trial_chart = chart(trial_geometry, trial_cell)
        test_cell_id = trial_cell_id
        if test_geometry !== trial_geometry
            test_cell_id = CompScienceMeshes.parent(trial_geometry, trial_cell_id)
        end
        test_chart = chart(test_geometry, test_cell_id)

        test_shapes = [(i, shape) for (i, fn) in enumerate(test_space.fns)
                       for shape in fn if shape.cellid == test_cell_id]
        trial_shapes = [(j, shape) for (j, fn) in enumerate(trial_space.fns)
                        for shape in fn if shape.cellid == trial_cell_id]

        isempty(test_shapes) && continue
        isempty(trial_shapes) && continue

        test_values = [
            rt_shape_value(test_chart, shape.refid, shape.coeff, vertex)
            for (_, shape) in test_shapes, vertex in trial_chart.vertices
        ]
        trial_values = [
            rt_shape_values(trial_chart, shape.refid, shape.coeff, trial_chart.vertices)
            for (_, shape) in trial_shapes
        ]
        area = volume(trial_chart)
        for (test_index, (i, _)) in enumerate(test_shapes)
            for (trial_index, (j, _)) in enumerate(trial_shapes)
                value = zero(area)
                for a in 1:3, b in 1:3
                    value += dot(test_values[test_index, a], trial_values[trial_index][b])
                end
                for a in 1:3
                    value += dot(test_values[test_index, a], trial_values[trial_index][a])
                end
                result[i, j] += area * value / 12
            end
        end
    end
    return result
end

@testset "exact vector Grams on larger surface meshes" begin
    meshes = (
        ("sphere", meshsphere(1.0, 0.1)),
        ("rectangle", meshrectangle(1.0, 0.8, 0.1, 3)),
        ("cube", meshcuboid(1.0, 0.9, 0.8, 0.1; generator=:compsciencemeshes)),
    )

    for (name, mesh) in meshes
        @testset "$name RT/RT" begin
            rt = raviartthomas(mesh)
            assembled = Matrix(assemble(Identity(), rt, rt))
            reference = exact_rt_gram(rt, rt)
            @test assembled ≈ reference atol=sqrt(eps()) rtol=sqrt(eps())
        end

        @testset "$name RT/BC" begin
            rt = raviartthomas(mesh)
            bc = buffachristiansen(mesh)
            # BC is defined on the barycentric refinement, so BEAST's
            # refinement path assembles BC x RT. Transpose it for RT x BC.
            assembled = Matrix(assemble(Identity(), rt, bc))
            reference = exact_rt_gram(rt, bc)
            @test assembled ≈ reference atol=sqrt(eps()) rtol=sqrt(eps())
        end
    end
end
