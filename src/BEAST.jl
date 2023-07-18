module BEAST

using Distributed
using LinearAlgebra
using SharedArrays
using SparseArrays
using FillArrays
using BlockArrays
using SparseMatrixDicts

using ConvolutionOperators
using SauterSchwabQuadrature
using SauterSchwab3D
using FastGaussQuadrature
using LinearMaps
using LiftedMaps

using AbstractTrees
using NestedUnitRanges

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
export brezzidouglasmarini3d
export nedelecd3d
export nedelecc3d
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
export ttrace 
export SingleLayer
export DoubleLayer
export DoubleLayerTransposed
export HyperSingular
export HH3DSingleLayerTDBIO
export HH3DDoubleLayerTDBIO
export ∂n, grad
export HH3DHyperSingularFDBIO
export DirichletTrace
export HH3DMonopole, gradHH3DMonopole
export HH3DLinearPotential

export HH3DSingleLayerNear
export HH3DDoubleLayerNear
export HH3DDoubleLayerTransposedNear
export HH3DHyperSingularNear

export NitscheHH3
export MWSingleLayerTDIO
export MWDoubleLayerTDIO
export MWFarField3DTD
export MWFarField3D
export MWSingleLayer3D, MWHyperSingular, MWWeaklySingular
export MWDoubleLayer3D
export PlaneWaveMW
export dipolemw3d, DipoleMW
export TangTraceMW, CrossTraceMW
export curl
export gradient
export MWSingleLayerField3D
export SingleLayerTrace
export DoubleLayerRotatedMW3D
export MWSingleLayerPotential3D

export VIEOperator

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

export DofInterpolate


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

function convolve end

include("utils/polynomial.jl")
include("utils/specialfns.jl")
include("utils/combinatorics.jl")
include("utils/linearspace.jl")
include("utils/zeromap.jl")
include("utils/rank1map.jl")

include("bases/basis.jl")
include("bases/lincomb.jl")
include("bases/trace.jl")
include("bases/restrict.jl")
include("bases/divergence.jl")

include("bases/local/laglocal.jl")
include("bases/local/rtlocal.jl")
include("bases/local/ndlocal.jl")
include("bases/local/bdmlocal.jl")
include("bases/local/ncrossbdmlocal.jl")
include("bases/local/ndlcclocal.jl")
include("bases/local/ndlcdlocal.jl")
include("bases/local/bdm3dlocal.jl")

include("bases/lagrange.jl")
include("bases/rtspace.jl")
include("bases/rtxspace.jl")
include("bases/bcspace.jl")
include("bases/ndspace.jl")
include("bases/bdmdiv.jl")
include("bases/ncrossbdmspace.jl")
include("bases/ndlccspace.jl")
include("bases/ndlcdspace.jl")
include("bases/dual3d.jl")
include("bases/bdm3dspace.jl")


include("bases/subdbasis.jl")
include("bases/stagedtimestep.jl")

include("bases/timebasis.jl")
include("bases/tensorbasis.jl")

include("operator.jl")

include("quadrature/quadstrats.jl")

include("excitation.jl")
include("localop.jl")
include("multiplicativeop.jl")
include("identityop.jl")
include("integralop.jl")
include("interpolation.jl")
include("quaddata.jl")

include("quadrature/quaddata.jl")
include("quadrature/quadrule.jl")

include("quadrature/double_quadrature.jl")
include("quadrature/singularity_extraction.jl")
include("quadrature/sauterschwabints.jl")

include("postproc.jl")
include("postproc/segcurrents.jl")
include("postproc/farfield.jl")

include("timedomain/tdintegralop.jl")
include("timedomain/tdexcitation.jl")
include("timedomain/motlu.jl")
include("timedomain/tdtimeops.jl")
include("timedomain/rkcq.jl")
include("timedomain/zdomain.jl")


# Support for Maxwell equations
include("maxwell/mwexc.jl")
include("maxwell/mwops.jl")
include("maxwell/nxdbllayer.jl")
include("maxwell/wiltonints.jl")
include("maxwell/nitsche.jl")
include("maxwell/farfield.jl")
include("maxwell/nearfield.jl")
include("maxwell/spotential.jl")
include("maxwell/maxwell.jl")
include("maxwell/sourcefield.jl")

# Support for the Helmholtz equation
include("helmholtz2d/helmholtzop.jl")

include("helmholtz3d/hh3dexc.jl")
include("helmholtz3d/hh3dops.jl")
include("helmholtz3d/nitsche.jl")
include("helmholtz3d/hh3dnear.jl")
include("helmholtz3d/hh3dfar.jl")
include("helmholtz3d/hh3d_sauterschwabqr.jl")
include("helmholtz3d/wiltonints.jl")

#suport for Volume Integral equation
include("volumeintegral/vie.jl")
include("volumeintegral/vieexc.jl")
include("volumeintegral/vieops.jl")
include("volumeintegral/farfield.jl")
include("volumeintegral/sauterschwab_ints.jl")

include("decoupled/dpops.jl")
include("decoupled/potentials.jl")

include("helmholtz3d/timedomain/tdhh3dops.jl")
include("helmholtz3d/timedomain/tdhh3dexc.jl")
include("helmholtz3d/timedomain/tdhh3dpp.jl")

include("maxwell/timedomain/mwtdops.jl")
include("maxwell/timedomain/mwtdexc.jl")
include("maxwell/timedomain/tdfarfield.jl")

include("utils/butchertableau.jl")
include("utils/variational.jl")

include("solvers/solver.jl")
include("solvers/lusolver.jl")
include("solvers/itsolver.jl")

include("utils/plotlyglue.jl")




const x̂ = point(1,0,0)
const ŷ = point(0,1,0)
const ẑ = point(0,0,1)
export x̂, ŷ, ẑ

const n = NormalVector()
export n

end # module
