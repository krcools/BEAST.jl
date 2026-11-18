using BEAST
using CUDA
using CompScienceMeshes
using LinearAlgebra

ext = Base.get_extension(BEAST, :BEASTCUDAExt)
Γ = meshsphere(1.0, 0.1)
X = raviartthomas(Γ)
operator = Maxwell3D.singlelayer(wavenumber=2π)
quadrature = BEAST.DoubleNumSauterQstrat(2, 2, 2, 2, 2, 2)

devices = collect(CUDA.devices())
tiling = ext.TilingStrategy(ext.EqualTiling(4), ext.EqualTiling(4))

matrix = assemble(operator, X, X;
    threading=:gpu,
    quadstrat=quadrature,
    # devices=devices,       # Use all visible GPUs instead of the current one.
    # nstreams=2,            # Assemble two tile pairs concurrently per GPU.
    # tilingstrat=tiling,    # Split the matrix explicitly to limit peak memory.
)

# reference = assemble(operator, X, X; threading=:single, quadstrat=quadrature)
# norm(matrix - reference) / norm(reference)
