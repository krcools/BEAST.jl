module BEASTCUDAExt

using CUDA
using CUDA.Adapt
using CUDA.CUSPARSE

using BEAST
import BEAST: assemble!, Threading, Space, IntegralOperator
import BEAST: _integrands
import BEAST: LagrangeRefSpace, RTRefSpace
using BEAST.CompScienceMeshes
using BEAST.LinearAlgebra
using BEAST.SauterSchwabQuadrature
using BEAST.SparseArrays
using BEAST.StaticArrays

Adapt.@adapt_structure CommonVertex
Adapt.@adapt_structure CommonEdge
Adapt.@adapt_structure CommonFace

include("BEASTCUDAExt/tiling.jl")
include("BEASTCUDAExt/gpu_utils.jl")
include("BEASTCUDAExt/gpu_basis.jl")
include("BEASTCUDAExt/gpu_integrals.jl")
include("BEASTCUDAExt/gpu_assemble_integralop.jl")

end
