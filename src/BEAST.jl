module BEAST

export planewave
export planewavemw3d
export n, ι

using Combinatorics

import Base.*


immutable NormalVector end
const n = NormalVector()

include("utils/polynomial.jl")
include("utils/sparsend.jl")
include("utils/specialfns.jl")
include("utils/combinatorics.jl")

include("bases/basis.jl")
include("bases/trace.jl")
include("bases/restrict.jl")
include("bases/divergence.jl")

include("bases/laglocal.jl")
include("bases/rtlocal.jl")
include("bases/ndlocal.jl")

include("bases/lagrange.jl")
include("bases/rtspace.jl")
include("bases/bcspace.jl")
include("bases/ndspace.jl")

include("bases/timebasis.jl")
include("bases/tensorbasis.jl")

include("quaddata.jl")
include("excitation.jl")
include("operator.jl")
include("localop.jl")
include("identityop.jl")
include("integralop.jl")
include("postproc.jl")

include("timedomain/tdintegralop.jl")
include("timedomain/tdexcitation.jl")
include("timedomain/motlu.jl")
include("timedomain/tdtimeops.jl")

# Support for Maxwell equations
include("maxwell/wiltonints.jl")
include("maxwell/mwops.jl")
include("maxwell/nxdbllayer.jl")
include("maxwell/mwexc.jl")
#include("maxwell/bogaertints.jl")
include("maxwell/nitsche.jl")
include("maxwell/farfield.jl")
include("maxwell/nearfield.jl")
include("maxwell/spotential.jl")

# Support for the Helmholtz equation
include("helmholtz2d/helmholtzop.jl")

include("helmholtz3d/hh3dops.jl")
include("helmholtz3d/hh3dexc.jl")
include("helmholtz3d/nitsche.jl")

include("helmholtz3d/timedomain/tdhh3dops.jl")
include("helmholtz3d/timedomain/tdhh3dexc.jl")

include("maxwell/timedomain/mwtdops.jl")
include("maxwell/timedomain/mwtdexc.jl")

include("utils/variational.jl")
include("lusolver.jl")

# try
#     Pkg.installed("LinearForms")
#     info("`LinearForms` detected: form compiler support enabled.")
#     @eval using LinearForms
#     @eval include("lusolver.jl")
# catch
#     #error("Please install prerequisite 'LinearForms':\n\n    Pkg.clone(\"https://github.com/krcools/LinearForms.jl\")")
#     warn("`LinearForms` not installed: form compiler support disabled.")
# end

using CompScienceMeshes

end # module
