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
include("test_subd_basis.jl")
include("test_rt2.jl")
include("test_nd2.jl")

include("test_dvg.jl")
include("test_bcspace.jl")
include("test_trace.jl")
include("test_ttrace.jl")
include("test_timebasis.jl")
include("test_rtports.jl")
include("test_ndjunction.jl")
include("test_ndspace.jl")
include("test_restrict.jl")
include("test_ndlcd_restrict.jl")
include("test_interpolate_and_restrict.jl")
include("test_rt3d.jl")
include("test_gradient.jl")
include("test_mult.jl")

include("test_gram.jl")
include("test_vector_gram.jl")
include("test_local_storage.jl")
include("test_embedding.jl")

include("test_assemblerow.jl")
include("test_mixed_blkassm.jl")
include("test_local_assembly.jl")
include("test_assemble_refinements.jl")

include("test_dipole.jl")

include("test_wiltonints.jl")
include("test_sauterschwabints.jl")
include("test_hh3dints.jl")
include("test_ss_nested_meshes.jl")
include("test_nitsche.jl")
include("test_nitschehh3d.jl")

include("test_curlcurlgreen.jl")
include("test_hh3dtd_exc.jl")
include("test_hh3dexc.jl")
include("test_hh3d_nearfield.jl")
include("test_tdassembly.jl")
include("test_tdhhdbl.jl")
include("test_tdmwdbl.jl")
include("test_compressed_storage.jl")
include("test_tdefie_irk.jl")
include("test_dyadicop.jl")
# include("test_matrixconv.jl")

include("test_tdop_scaling.jl")
include("test_tdrhs_scaling.jl")
include("test_td_tensoroperator.jl")

include("test_variational.jl")

include("test_handlers.jl")
include("test_gridfunction.jl")

include("test_hh_lsvie.jl")

using TestItemRunner
@run_package_tests


try
    Pkg.installed("BogaertInts10")
    @info "`BogaertInts10` detected. Including relevant tests."
    include("test_bogaertints.jl")
catch
    @info "`Could not load BogaertInts10`. Related tests skipped."
end


end
