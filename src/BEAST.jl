module BEAST

using Distributed
using LinearAlgebra
using Pkg
using SharedArrays
using SparseArrays

using SauterSchwabQuadrature
using FastGaussQuadrature

import LinearAlgebra: cross, dot
import LinearAlgebra: ×, ⋅

import SharedArrays: sdata

export dot

export planewave
export RefSpace, numfunctions, coordtype, scalartype, assemblydata, geometry, refspace, valuetype
export lagrangecxd0, lagrangec0d1, duallagrangec0d1
export duallagrangecxd0
export lagdimension
export restrict
export raviartthomas, raowiltonglisson, positions
export brezzidouglasmarini
export portcells, rt_ports, getindex_rtg, subset
export StagedTimeStep
export subdsurface,subdBasis,assemblydata,refspace
export spatialbasis, temporalbasis
export ⊗
export timebasisc0d1
export timebasiscxd0
export timebasisdelta
export timebasisshiftedlagrange
export TimeBasisDeltaShifted
export ntrace
export strace
export SingleLayer
export DoubleLayer
export DoubleLayerTransposed
export HyperSingular
export HH3DSingleLayerTDBIO
export HH3DDoubleLayerTDBIO
export ∂n
export HH3DHyperSingularFDBIO
export NitscheHH3
export MWSingleLayerTDIO
export MWDoubleLayerTDIO
export MWFarField3DTD
export MWFarField3D
export MWSingleLayer3D, MWHyperSingular, MWWeaklySingular
export MWDoubleLayer3D
export PlaneWaveMW
export TangTraceMW, CrossTraceMW
export curl
export MWSingleLayerField3D
export SingleLayerTrace
export DoubleLayerRotatedMW3D
export MWSingleLayerPotential3D
export gmres
export @hilbertspace, @varform, @discretise
export solve
export convolve
export marchonintime
export RungeKuttaConvolutionQuadrature
export laplace_to_z, inverse_z_transform, real_inverse_z_transform
export butcher_tableau_radau_2stages
export butcher_tableau_radau_3stages
export creategaussian
export derive, integrate
export fouriertransform
export assemble
export Identity
export shapevals
export NCross
export quadrule, elements
export blockassembler
export localoperator
export localoperator2
export assemble
export facecurrents
export potential
export get_scatter_parameters
export quaddata



export kernelvals
export integrand
export shapevals
export ScalarTrace
export PlaneWaveDirichlet
export PlaneWaveNeumann

struct NormalVector end

using CompScienceMeshes
using Combinatorics
using FFTW
using SparseArrays

include("utils/polynomial.jl")
include("utils/sparsend.jl")
include("utils/specialfns.jl")
include("utils/combinatorics.jl")
include("utils/linearspace.jl")
include("utils/matrixconv.jl")
include("utils/polyeig.jl")

include("bases/basis.jl")
include("bases/lincomb.jl")
include("bases/trace.jl")
include("bases/restrict.jl")
include("bases/divergence.jl")

include("bases/local/laglocal.jl")
include("bases/local/rtlocal.jl")
include("bases/local/ndlocal.jl")
include("bases/local/bdmlocal.jl")

include("bases/lagrange.jl")
include("bases/rtspace.jl")
include("bases/rtxspace.jl")
include("bases/bcspace.jl")
include("bases/ndspace.jl")
include("bases/bdmdiv.jl")


include("bases/subdbasis.jl")
include("bases/stagedtimestep.jl")

include("bases/timebasis.jl")
include("bases/tensorbasis.jl")

include("excitation.jl")
include("operator.jl")
include("localop.jl")
include("identityop.jl")
include("integralop.jl")
include("quaddata.jl")
include("postproc.jl")

include("quadrature/double_quadrature.jl")
include("quadrature/singularity_extraction.jl")

include("timedomain/tdintegralop.jl")
include("timedomain/tdexcitation.jl")
include("timedomain/motlu.jl")
include("timedomain/tdtimeops.jl")
include("timedomain/rkcq.jl")
include("timedomain/zdomain.jl")

# Support for Maxwell equations
include("maxwell/mwexc.jl")
include("maxwell/mwops.jl")
include("maxwell/wiltonints.jl")
include("maxwell/sauterschwabints_rt.jl")
include("maxwell/sauterschwabints_bdm.jl")
include("maxwell/nxdbllayer.jl")
include("maxwell/nitsche.jl")
include("maxwell/farfield.jl")
include("maxwell/nearfield.jl")
include("maxwell/spotential.jl")
include("maxwell/maxwell.jl")

# Support for the Helmholtz equation
include("helmholtz2d/helmholtzop.jl")

include("helmholtz3d/hh3dexc.jl")
include("helmholtz3d/hh3dops.jl")
include("helmholtz3d/nitsche.jl")

include("decoupled/dpops.jl")

include("helmholtz3d/timedomain/tdhh3dops.jl")
include("helmholtz3d/timedomain/tdhh3dexc.jl")

include("maxwell/timedomain/mwtdops.jl")
include("maxwell/timedomain/mwtdexc.jl")
include("maxwell/timedomain/tdfarfield.jl")

include("utils/butchertableau.jl")
include("utils/variational.jl")

include("solvers/solver.jl")
include("solvers/lusolver.jl")
include("solvers/itsolver.jl")




const x̂ = point(1,0,0)
const ŷ = point(0,1,0)
const ẑ = point(0,0,1)
export x̂, ŷ, ẑ

const n = NormalVector()
export n

end # module
