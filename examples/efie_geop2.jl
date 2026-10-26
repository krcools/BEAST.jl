# # EFIE on Curvilinear Geometries
#
# This tutorial demonstrates how to solve the Electric Field Integral Equation
# (EFIE) for scattering by a perfectly conducting sphere using higher order
# geometry descriptions and higher order basis functions. Two discretisations are
# compared:
#
# 1. a flat-faceted mesh carrying lowest order divergence conforming
#    Graglia-Wilton-Peterson (GWP) functions, and
# 2. a curvilinear mesh of second order triangles carrying second order GWP
#    functions.
#

using CompScienceMeshes
using BEAST
using LinearAlgebra
import PlotlyDocumenter #hide
const CSM = CompScienceMeshes;

# ## Problem definition
#
# The EFIE is formulated in terms of the Maxwell single layer operator. The
# excitation is the tangential trace of a plane wave travelling in the
# ``z``-direction and polarised along ``x``.

κ, η = 3.0, 1.0
t = Maxwell3D.singlelayer(wavenumber=κ)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n;

# Placeholders for the trial and test functions are introduced. These are used
# to write down the bilinear and linear forms making up the variational problem.

@hilbertspace j
@hilbertspace k;

# ## Flat facets and lowest order basis functions
#
# GMsh is used to create a mesh for the unit sphere. With `order=1`, the mesh
# consists of flat triangles. On this mesh, the lowest order GWP space is
# constructed. This space is equivalent to the space of Raviart-Thomas functions.

radius = 1.0
h = 0.3

Γ1 = CSM.gmshsphere(radius, h; order=1)
X1 = BEAST.gwpdiv(Γ1; order=0)
@show numfunctions(X1)

# For accurate results with higher order basis functions on curved elements, the
# default quadrature strategy is replaced by one with more quadrature points. For
# the operator, the first two arguments of `DoubleNumSauterQstrat` set the
# orders of the outer and inner rules used for well separated triangle pairs.
# The remaining arguments set the orders of the Sauter-Schwab rules used for
# pairs of triangles that coincide, or share an edge or a vertex. The excitation
# is integrated with a rule of order 12.
#
# The macro `@defaultquadstrat` defines a method that dispatches on the *types*
# of the operator and the spaces, so a space must be available when it is
# invoked. The new default applies to all spaces of the same type as `X1`.

BEAST.@defaultquadstrat (t,X1,X1) BEAST.DoubleNumSauterQstrat(4,5,6,6,6,6)
BEAST.@defaultquadstrat (e,X1) BEAST.SingleNumQStrat(12);

# The system matrix and right hand side are assembled. The system matrix is
# assembled in parallel, using a colouring of the mesh cells to avoid race
# conditions. The resulting dense system is solved by LU decomposition.

A1 = assemble(t[k,j], ∏(X1), ∏(X1); threading=:cellcoloring)
b1 = assemble(e[k], ∏(X1))
u1 = BEAST.lu(AbstractMatrix(A1)) * b1;

# ## Curvilinear geometry and second order basis functions
#
# Next, the sphere is approximated by a mesh of second order triangles. These
# elements are parametrised by quadratic polynomials and conform much more
# closely to the true geometry. Second order GWP functions are defined on this
# mesh.

Γ2 = CSM.gmshsphere(radius, h; order=2)
X2 = BEAST.gwpdiv(Γ2; order=2)
@show numfunctions(X2)

# The curvilinear space is of a different type than `X1`, so the quadrature
# strategy has to be set for it as well. Assembly and solution then proceed
# exactly as before.

BEAST.@defaultquadstrat (t,X2,X2) BEAST.DoubleNumSauterQstrat(4,5,6,6,6,6)
BEAST.@defaultquadstrat (e,X2) BEAST.SingleNumQStrat(12);

A2 = assemble(t[k,j], ∏(X2), ∏(X2); threading=:cellcoloring)
b2 = assemble(e[k], ∏(X2))
u2 = BEAST.lu(AbstractMatrix(A2)) * b2;

# ## Visualisation
#
# The function `facecurrents` evaluates the solution on a refinement of the
# mesh. The keyword argument `order` controls the number of subdivisions per
# element, which allows the variation of the higher order current within each
# element to be visualised. For curvilinear meshes, the refinement follows the
# curved geometry.

fcr1, geo1 = facecurrents(u1, X1; order=3)
fcr2, geo2 = facecurrents(u2, X2; order=3);

# The magnitude of the surface current is plotted for both discretisations,
# using the same colour scale.

using PlotlyBase
cmin = 0.0
cmax = maximum(norm.(fcr2))
plt1 = Plot([patch(geo1, norm.(fcr1); caxis=(cmin, cmax))],
    Layout(title="Flat facets, GWP order 0"))
PlotlyDocumenter.to_documenter(plt1) #hide

#-

plt2 = Plot([patch(geo2, norm.(fcr2); caxis=(cmin, cmax))],
    Layout(title="Curvilinear mesh, GWP order 2"))
PlotlyDocumenter.to_documenter(plt2) #hide

# Comparing both plots shows how much more detail the higher order solution on
# the curvilinear mesh resolves within each element.
