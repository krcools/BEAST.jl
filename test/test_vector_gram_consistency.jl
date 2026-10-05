using Test
using SparseArrays
using CompScienceMeshes
using BEAST

function assemble_mixed_reference(operator, test_space, trial_space)
    factory, store = BEAST.allocatestorage(
        operator,
        test_space,
        trial_space,
        Val{:bandedstorage},
        BEAST.LongDelays{:ignore},
    )
    BEAST.assemble_local_mixed!(operator, test_space, trial_space, store)
    return factory()
end

@testset "RT/BC assembly consistency" begin
    meshes = (
        ("sphere", readmesh(joinpath(dirname(pathof(BEAST)), "../test/assets/sphere316.in"))),
        ("rectangle", meshrectangle(1.0, 0.8, 0.2, 3)),
        ("cube", meshcuboid(1.0, 0.9, 0.8, 0.3; generator=:compsciencemeshes)),
    )
    identity = BEAST.Identity()

    for (name, mesh) in meshes
        @testset "$name" begin
            rt = raviartthomas(mesh)
            bc = buffachristiansen(mesh)

            assembled_bc_rt = assemble(identity, bc, rt)
            assembled_rt_bc = assemble(identity, rt, bc)
            mixed_bc_rt = assemble_mixed_reference(identity, bc, rt)
            mixed_rt_bc = assemble_mixed_reference(identity, rt, bc)

            @test size(assembled_bc_rt) == (numfunctions(bc), numfunctions(rt))
            @test size(assembled_rt_bc) == (numfunctions(rt), numfunctions(bc))
            @test size(assembled_bc_rt) == size(mixed_bc_rt)
            @test size(assembled_rt_bc) == size(mixed_rt_bc)
            @test nnz(assembled_bc_rt) > 0
            @test nnz(assembled_rt_bc) > 0
            @test all(isfinite, nonzeros(assembled_bc_rt))
            @test all(isfinite, nonzeros(assembled_rt_bc))
            @test assembled_bc_rt ≈ mixed_bc_rt
            @test assembled_rt_bc ≈ mixed_rt_bc
            @test assembled_rt_bc ≈ assembled_bc_rt'
        end
    end
end
