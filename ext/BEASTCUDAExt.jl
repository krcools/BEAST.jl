module BEASTCUDAExt

using CUDA
using CUDA.CUSPARSE

using BEAST
import BEAST: assemble!, Threading, Space, IntegralOperator
import BEAST: _integrands
import BEAST: LagrangeRefSpace, RTRefSpace
using BEAST.CompScienceMeshes
using BEAST.SparseArrays
using BEAST.StaticArrays

include("BEASTCUDAExt/gpu_utils.jl")
include("BEASTCUDAExt/gpu_basis.jl")
include("BEASTCUDAExt/gpu_integrals.jl")
include("BEASTCUDAExt/gpu_assemble_integralop.jl")

end
