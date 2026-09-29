using BEAST
using CUDA
using CompScienceMeshes
using LinearAlgebra

ext = Base.get_extension(BEAST, :BEASTCUDAExt)
Γ = meshsphere(1.0, 0.1)
X = raviartthomas(Γ)
operator = Maxwell3D.singlelayer(wavenumber=2π)
quadrature = BEAST.DoubleNumSauterQstrat(2, 2, 2, 2, 2, 2)

assembler = ext.gpu_blockassembler(operator, X, X;
    quadstrat=quadrature,
    # device=CUDA.device(),  # Select the GPU that owns this assembler.
    # gpu_blocksize=256,     # Set the number of threads used by its kernels.
)

test_ids = collect(1:32)
trial_ids = collect((numfunctions(X)-31):numfunctions(X))
block = CUDA.zeros(scalartype(operator, X, X), length(test_ids), length(trial_ids))
assembler(block, test_ids, trial_ids)

test_blocks = [test_ids, trial_ids]
trial_blocks = [trial_ids, test_ids]
blocks = ext.gpu_batched_blockassemble(
    assembler,
    test_blocks,
    trial_blocks;
    # budget=256 * 2^20,     # Lower this limit to process the blocks in more chunks.
)

# devices = collect(CUDA.devices())
# assemblers = ext.gpu_blockassemblers(operator, X, X; quadstrat=quadrature, devices)
# blocks = ext.gpu_batched_blockassemble(assemblers, test_blocks, trial_blocks)

# reference = assemble(operator, X, X; threading=:single, quadstrat=quadrature)
# norm(Array(block) - reference[test_ids, trial_ids]) / norm(reference[test_ids, trial_ids])
