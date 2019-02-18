# write your own tests here

using StaticArrays

module PkgTests

using Distributed
using LinearAlgebra
using SparseArrays
using Test
using Pkg

import BEAST

include("test_fourier.jl")
include("test_specials.jl")

include("test_basis.jl")
include("test_directproduct.jl")
include("test_raviartthomas.jl")
include("test_rt.jl")
include("test_rtx.jl")

include("test_dvg.jl")
include("test_bcspace.jl")
include("test_trace.jl")
include("test_timebasis.jl")
include("test_rtports.jl")

include("test_gram.jl")
include("test_vector_gram.jl")

include("test_assemblerow.jl")

include("test_wiltonints.jl")
include("test_sauterschwabints.jl")
include("test_ss_nested_meshes.jl")
include("test_nitsche.jl")
include("test_nitschehh3d.jl")

include("test_tdop_scaling.jl")
include("test_tdrhs_scaling.jl")

include("test_farfield.jl")

#include("test_assemblerow.jl")

try
    Pkg.installed("BogaertInts10")
    @info "`BogaertInts10` detected. Including relevant tests."
    include("test_bogaertints.jl")
catch
    @info "`Could not load BogaertInts10`. Related tests skipped."
end


end
