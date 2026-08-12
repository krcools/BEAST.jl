# write your own tests here

using StaticArrays

module PkgTests

using TestItemRunner

using Distributed
using LinearAlgebra
using SparseArrays
using Test
using Pkg

import BEAST

@testitem "fourier" begin
    include("test_fourier.jl")
end

@testitem "specials" begin
    include("test_specials.jl")
end

@testitem "basis" begin
    include("test_basis.jl")
end

@testitem "lagrange" begin
    include("test_lagrange.jl")
end

@testitem "directproduct" begin
    include("test_directproduct.jl")
end

@testitem "raviartthomas" begin
    include("test_raviartthomas.jl")
end

@testitem "rt" begin
    include("test_rt.jl")
end

@testitem "rtx" begin
    include("test_rtx.jl")
end

@testitem "subd_basis" begin
    include("test_subd_basis.jl")
end

@testitem "rt2" begin
    include("test_rt2.jl")
end

@testitem "nd2" begin
    include("test_nd2.jl")
end

@testitem "dvg" begin
    include("test_dvg.jl")
end

@testitem "bcspace" begin
    include("test_bcspace.jl")
end

@testitem "trace" begin
    include("test_trace.jl")
end

@testitem "ttrace" begin
    include("test_ttrace.jl")
end

@testitem "timebasis" begin
    include("test_timebasis.jl")
end

@testitem "rtports" begin
    include("test_rtports.jl")
end

@testitem "ndjunction" begin
    include("test_ndjunction.jl")
end

@testitem "ndspace" begin
    include("test_ndspace.jl")
end

@testitem "restrict" begin
    include("test_restrict.jl")
end

@testitem "ndlcd_restrict" begin
    include("test_ndlcd_restrict.jl")
end

@testitem "interpolate_and_restrict" begin
    include("test_interpolate_and_restrict.jl")
end

@testitem "rt3d" begin
    include("test_rt3d.jl")
end

@testitem "gradient" begin
    include("test_gradient.jl")
end

@testitem "mult" begin
    include("test_mult.jl")
end

@testitem "gram" begin
    include("test_gram.jl")
end

@testitem "vector_gram" begin
    include("test_vector_gram.jl")
end

@testitem "local_storage" begin
    include("test_local_storage.jl")
end

@testitem "embedding" begin
    include("test_embedding.jl")
end

@testitem "assemblerow" begin
    include("test_assemblerow.jl")
end

# include("test_mixed_blkassm.jl")
@testitem "local_assembly" begin
    include("test_local_assembly.jl")
end

@testitem "assemble_refinements" begin
    include("test_assemble_refinements.jl")
end

@testitem "dipole" begin
    include("test_dipole.jl")
end

@testitem "sauterschwabints1D" begin
    include("test_sauterschwabints1D.jl")
end

@testitem "telles" begin
    include("test_telles.jl")
end

@testitem "hh2d_exec" begin
    include("test_hh2d_exec.jl")
end

@testitem "hh2d_ops" begin
    include("test_hh2d_ops.jl")
end

@testitem "hh2d_mie" begin
    include("test_hh2d_mie.jl")
end

@testitem "hh2d_mie_higher_order" begin
    include("test_hh2d_mie_higher_order.jl")
end

@testitem "hh2d_nearfield" begin
    include("test_hh2d_nearfield.jl")
end

@testitem "wiltonints" begin
    include("test_wiltonints.jl")
end

@testitem "sauterschwabints" begin
    include("test_sauterschwabints.jl")
end

@testitem "hh3dints" begin
    include("test_hh3dints.jl")
end

@testitem "ss_nested_meshes" begin
    include("test_ss_nested_meshes.jl")
end

@testitem "nitsche" begin
    include("test_nitsche.jl")
end

@testitem "nitschehh3d" begin
    include("test_nitschehh3d.jl")
end

@testitem "curlcurlgreen" begin
    include("test_curlcurlgreen.jl")
end

@testitem "hh3dtd_exc" begin
    include("test_hh3dtd_exc.jl")
end

# include("test_hh3dexc.jl")
@testitem "hh3d_nearfield" begin
    include("test_hh3d_nearfield.jl")
end

@testitem "tdassembly" begin
    include("test_tdassembly.jl")
end

@testitem "tdhhdbl" begin
    include("test_tdhhdbl.jl")
end

@testitem "tdmwdbl" begin
    include("test_tdmwdbl.jl")
end

@testitem "compressed_storage" begin
    include("test_compressed_storage.jl")
end

@testitem "tdefie_irk" begin
    include("test_tdefie_irk.jl")
end

@testitem "dyadicop" begin
    include("test_dyadicop.jl")
end

@testitem "tdop_scaling" begin
    include("test_tdop_scaling.jl")
end

@testitem "tdrhs_scaling" begin
    include("test_tdrhs_scaling.jl")
end

@testitem "td_tensoroperator" begin
    include("test_td_tensoroperator.jl")
end

@testitem "variational" begin
    include("test_variational.jl")
end

@testitem "handlers" begin
    include("test_handlers.jl")
end

@testitem "ncrossbdm" begin
    include("test_ncrossbdm.jl")
end

@testitem "gridfunction" begin
    include("test_gridfunction.jl")
end

@testitem "itsolver" begin
    include("test_itsolver.jl")
end

@testitem "hh_lsvie" begin
    include("test_hh_lsvie.jl")
end

@testitem "composed_basis" begin
    include("test_composed_basis.jl")
end

@testitem "composed_operator" begin
    include("test_composed_operator.jl")
end

@testitem "analytic_excitation" begin
    include("test_analytic_excitation.jl")
end

@testitem "vie" begin
    include("test_vie.jl")
end

@testitem "evie_dvie" begin
    include("test_evie_dvie.jl")
end

@testitem "coloring" begin
    include("test_coloring.jl")
end

@run_package_tests filter = ti -> begin
    # @show ti.tags
    # @show isempty(intersect([:example, :diagnostics], ti.tags))
    isempty(intersect([:example, :diagnostics], ti.tags))
end verbose = true

try
    Pkg.installed("BogaertInts10")
    @info "`BogaertInts10` detected. Including relevant tests."
    include("test_bogaertints.jl")
catch
    @info "`Could not load BogaertInts10`. Related tests skipped."
end


end
