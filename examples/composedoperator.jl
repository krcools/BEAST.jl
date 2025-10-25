
# # Green function
# The green function or more generaly the kernel function in the double integral, depending on both the coordinates x and y, is supplied by the user as a function.
# There are standard implementations of the Helmholtz 3D Green function and the gradient of the Helmholtz 3D Green function:
# ```math
# \begin{aligned}
# G(x,y) &= \frac{e^{-γ|x-y|}}{4\pi|x-y|} \\
# \nabla_x G(x,y) &= -\frac{e^{-γ|x-y|}}{4\pi|x-y|^2} \left( γ + \frac{1}{|x-y|} \right) (x-y)
# \end{aligned}
# ```
# Those can be called as 
using BEAST
using CompScienceMeshes
γ = 2*pi*1im*10^8
G = BEAST.HH3DGreen(γ)
∇ₓG = BEAST.HH3DGradGreen(γ);
# writing your own kernelfunction is possible by defining the following structure
struct MyGreenFunction{T} <: BEAST.Kernel{T}

end
# and the function 
function (g::MyGreenFunction)(x,y)
    return "# the value of the kernel function #"
end
# **_NOTE:_** Trace operators will not know how to handle your custom kernel function, the principle value will be used in that case. To implement the correct treatment of your kernel function by the trace operator, please read the manual on the trace operator.

# # Potential operator
# A Potential operator maps a function defined on the mesh to a function defined in the domain. Such a map is typically defined as a single integral over your mesh as 
# ```math
# u(x) = \int_{\Gamma} G(x,y) f(y) dy
# ```
# where ``G`` is the kernel function and ``f`` is the function defined on the mesh.
# Potentials are used in both the definition of boundary or volume integral operators and in the postprocessing routines to evaluate the fields in the domain.
# A potential operator can be constructed by calling
KernelFunction = G
Operation = *
BasisOperation = B -> B
D = 2
BEAST.PotentialIntegralOperator{D}(KernelFunction, Operation, BasisOperation);
# where `D` is the dimension of the mesh, `KernelFunction` is the kernel function, for example the green function defined before, `operation` is the operation to be used between the kernelfunction and the basis function in the integral (e.g. `*` for scalar multiplication, `×` for the cross product), and `BasisOperation` is a function that is applied on the basis only, for example gradient, divergence, curl,...
# The double-layer potential ``u(x) = \int_{\Gamma} \nabla G(x,y) \times f(y) dy`` can for example be written as
K = BEAST.PotentialIntegralOperator{2}(BEAST.HH3DGradGreen(γ), ×, B -> B);
# The hyper singular part of the single-layer potential ``u(x) = \int_{\Gamma} \nabla G(x,y) \nabla \cdot f(y) dy`` can be written as
L = BEAST.PotentialIntegralOperator{2}(BEAST.HH3DGradGreen(γ), *, B -> divergence(B));
# **_NOTE:_** Be carefull when using the dot product as operation, this will take the complex conjugate of the first argument which is not always desired. Use f(x,y) -> transpose(x)*y instead.
# ## Potential operators in the post-processing
# The potential operator can be used in the post-processing routines to evaluate the fields in the domain.
# For example, the double-layer potential is evaluated as
#-
# ````
# potential(K, points, coefficients, basis)
# ````
#-
# # Trace operator
# The trace operator can be applied on a potential operator. This switches to a principal value interpretation of the integral by adding extra nessesary (+/-)1/2 contributions. The oprators yielded by this trace operator include the second duality paring integral and can be cast directly into the `assemble` functionality.
# The Trace operator can be defined as
TestFunctionMap = T->T
Interior = true
PrincipalValue=false
U = 2
T = BEAST.TraceOperator{U}(TestFunctionMap, Interior, PrincipalValue);
# with `U` the dimension of the testmesh, `TestFunctionMap` a function applied on the testspace, `Interior` a boolean if the trace is taken from the interior or exterior of the testdomain, `PrincipalValue` a boolean, true if only the principal value is computed, if false, also the (local) trace term is computed.
# The following standard implementations of the trace operator are available:\
# BEAST.extttrace: ``\lim_{\boldsymbol{r}\notin\Omega \to \partial \Omega} - n\times \left(n\times \boldsymbol{V}\left(\boldsymbol{r}\right)\right)  ``\
# BEAST.intttrace: ``\lim_{\boldsymbol{r}\in\Omega \to \partial \Omega} - n\times \left( n\times \boldsymbol{V}\left(\boldsymbol{r}\right)\right)  ``\
# BEAST.extntrace: ``\lim_{\boldsymbol{r}\notin\Omega \to \partial \Omega}  n\cdot \boldsymbol{V}\left(\boldsymbol{r}\right)  ``\
# BEAST.intntrace: ``\lim_{\boldsymbol{r}\notin\Omega \to \partial \Omega}  n\cdot \boldsymbol{V}\left(\boldsymbol{r}\right)  ``\
# BEAST.extrtrace: ``\lim_{\boldsymbol{r}\notin\Omega \to \partial \Omega} n\times \boldsymbol{V}\left(\boldsymbol{r}\right)  ``\
# BEAST.intrtrace: ``\lim_{\boldsymbol{r}\in\Omega \to \partial \Omega}  n\times \boldsymbol{V}\left(\boldsymbol{r}\right)  ``\
# BEAST.exttrace: ``\lim_{\boldsymbol{r}\notin\Omega \to \partial \Omega}   V\left(\boldsymbol{r}\right)  ``\
# BEAST.inttrace: ``\lim_{\boldsymbol{r}\notin\Omega \to \partial \Omega}   V\left(\boldsymbol{r}\right)  ``

# ## Conventions and default behavior
# Meshes suplied to the trace operator (those meshes are supplied trough the supplied basis) should have outward pointing normal vectors. The volumes enclosed by the different meshes should be disjoint. 
# ##  Displaced Meshes
# For an intuitive idea of this feature we refer to the "gap" idea in the global multi-trace method. The displacement indicates for two overlapping triangles from within which domain the trace was taken. Intiutively, this feature contains the information about the limit taken and the interior or exterior information. 
# Default behavior: Domains are expected to be disjoint with outward pointing normalvector. The direction of the limit is the direction in which all gaps are closed. This is achieved by displacing over -1.0*normalvector.
# If different behavior is expected form the mesh, one should wrap it in a displaced mesh attribute as:
eps = 1.0
Γ = meshcuboid(1.0,1.0,1.0,0.3)
Γ2 = BEAST.GlobalDisplacementMesh(Γ,eps);
# where `eps` is the distance over which a chart is displaced in the direction of the normal on this chart. Be carefull if multiple traingles from differen meshes overlap that the size of eps is chosen in such a way that the desired limits are taken correctly in all operators.
# More advanced displacement meshes can be constructed as a subtype of BEAST.DisplacementMesh{T},  with T the scalartype. The new structure should implement the function displacementchart, which defines for a givenen chart index the simplex together with the local displacement.

# ## Principal value
# The principal value of the integral can be called by setting `PrincipalValue` = true in 
TestFunctionMap = T->T
Interior = true
PrincipalValue=true
T = BEAST.TraceOperator{U}(TestFunctionMap, Interior, PrincipalValue);
# The following standard implementations are available:

# BEAST.pvttrace: ``- n\times n\times \boldsymbol{V}\left(\boldsymbol{r}\right)  `` \
# BEAST.pvntrace: ``  n\cdot \boldsymbol{V}\left(\boldsymbol{r}\right)  `` \
# BEAST.pvrtrace: `` n\times \boldsymbol{V}\left(\boldsymbol{r}\right)  `` \
# BEAST.pvtrace: ``   V\left(\boldsymbol{r}\right)  `` \


# ## Volume Trace 
# A potential operator can be mapped to a volumetric operator with:
BEAST.volumetrace(K) ;

# # Hilbert Space Framework
# The hilbertspace frameweork is implemented for potentials to implement potential operators on directproduct spaces. For example, with the above defined potentials. The potential function can be called on this general potential operator in the post-processing.
@hilbertspace j m 
p = K[j] + L[m] ;

# The Trace operator mapping a potential to a direct product of trace spaces can be called as 
@hilbertspace j m 
t = T[j] + T[m];

# This generarl trace operator can be called on a general potential, yielding all the nessesary boundary element operators. This system can be assembled using the general assemble routines.
lhs = t(p);

# # Global Multi-Trace Example 
# Below, two examples are sketched of the global multi-trace pmchwt integral equation. Example one contains two cubes next to eachother where cube one is cut out of cube two (see plot), for which the default routines can be used, example two shows cube two inside of cube one, here we have to define the displacement ourself. Both examples should yield the same result.
# ## The common part
using BEAST
using CompScienceMeshes
using PlotlyJS
using PlotlyDocumenter #hide

import BEAST.BlockArrays.BlockedVector
import BEAST.NestedUnitRanges.nestedrange
import Plots.heatmap as heatmap
# physical constants

ω = 2*pi*10.0^8
ϵ0, ϵ1, ϵ2 = [1.0, 2.0, 3.0]*8.854*10^(-12)
μ0, μ1, μ2 = [2.0, 2.0, 1.0]*4*pi*10^(-7)
κ0, κ1, κ2 = ω*sqrt(ϵ0*μ0), ω*sqrt(ϵ1*μ1), ω*sqrt(ϵ2*μ2);


# Define meshes
Γ1 = meshcuboid(2.0,2.0,2.0,0.125)
Γ2 = translate!(meshcuboid(1.0,1.0,1.0,0.125),(point(1.0,0.0,0.0)))

function remove_both_if_double(a)
    c = sort.(a)
    b = [i for i in a if length(findall(==(sort(i)),c))==1]
    return b
end
Γ3 = weld(Γ1,Γ2;glueop=remove_both_if_double)
CompScienceMeshes.orient(Γ3)

# Plot of structures
plt = Plot([wireframe(Γ3),patch(Γ2;color=:blue)])
PlotlyDocumenter.to_documenter(plt) #hide
#-
plt = Plot([wireframe(Γ1),patch(Γ2;color=:blue)])
PlotlyDocumenter.to_documenter(plt) #hide

# Function Spaces. The mesh Γ1 is artificailly expanded to acheave correct traces in the second example.
X1 = raviartthomas(BEAST.GlobalDisplacementMesh(Γ1,1))
X2 = raviartthomas(Γ2)
X3 = raviartthomas(Γ3)

Xtot1 = BEAST.DirectProductSpace([X3,X3,X2,X2])
Xtot2 = BEAST.DirectProductSpace([X1,X1,X2,X2]);

# Define Green functions in each domain
G0 = BEAST.HH3DGreen(1im*κ0)
G1 = BEAST.HH3DGreen(1im*κ1)
G2 = BEAST.HH3DGreen(1im*κ2)
dG0 = BEAST.HH3DGradGreen(1im*κ0)
dG1 = BEAST.HH3DGradGreen(1im*κ1)
dG2 = BEAST.HH3DGradGreen(1im*κ2);

# Define Single-Layer Potentials in each domain
S0 = BEAST.PotentialIntegralOperator{2}(G0,*,B->B) +1/κ0^2* BEAST.PotentialIntegralOperator{2}(dG0,*,B->divergence(B))
S1 = BEAST.PotentialIntegralOperator{2}(G1,*,B->B) +1/κ1^2* BEAST.PotentialIntegralOperator{2}(dG1,*,B->divergence(B))
S2 = BEAST.PotentialIntegralOperator{2}(G2,*,B->B) +1/κ2^2* BEAST.PotentialIntegralOperator{2}(dG2,*,B->divergence(B));


# Define the Traces in each domain
Tin = BEAST.intttrace
Tex = BEAST.extttrace;


# Define Double-Layer Potentials in each domain
K0 = BEAST.PotentialIntegralOperator{2}(dG0,×,B->B)
K1 = BEAST.PotentialIntegralOperator{2}(dG1,×,B->B)
K2 = BEAST.PotentialIntegralOperator{2}(dG2,×,B->B);

# Define the incident field
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ0)
H = -1/(im*ω*μ0)*curl(E)

e = (n × E) × n
h = (n × H) × n;

# Construct the hilbert spaces
@hilbertspace j1 m1 j2 m2 ;

# Post-Processing definitions
xrange = range(-2.0, stop=4.0, length=100)
zrange = range(-2.0, stop=4.0, length=100)
pts = [point(x,1.0,z) for x ∈ xrange, z ∈ zrange];


# ## 2 Cubes next to eachother (with cutting cube one out of cube two)

# Scattered Electric field inside the domains
E0 = 1im*ω*μ0* S0[j1] + K0[m1] + 1im*ω*μ0* S0[j2] + K0[m2]
E1 = -1im*ω*μ1* S1[j1] - K1[m1]
E2 = -1im*ω*μ2* S2[j2] - K2[m2];

# Scattered Magnetic field inside the domains
H0 = 1im*ω*ϵ0* S0[m1] - K0[j1] + 1im*ω*ϵ0* S0[m2] - K0[j2]
H1 = -1im*ω*ϵ1* S1[m1] + K1[j1]
H2 = -1im*ω*ϵ2* S2[m2] + K2[j2];

# Construct Trace Operators
TH0 = Tex[j1] + Tex[j2] 
TH1 = -Tin[j1] 
TH2 = -Tin[j2] 

TE0 = Tex[m1] + Tex[m2]
TE1 = -Tin[m1]
TE2 = -Tin[m2];

# Design System matrix

Z = TH0(H0) + TH1(H1) +  TH2(H2) + TE0(E0) + TE1(E1) + TE2(E2);

# rhs

rhs = e[m1]+e[m2]+h[j1]+h[j2];


# Assemble and solve
M1 = assemble(Z,Xtot1,Xtot1);
R1 = assemble(rhs,Xtot1);

u1 = Matrix(M1)\Vector(R1);
ax = nestedrange(Xtot1, 1, numfunctions);
u1 = BlockedVector(u1, (ax,));
@show "system solved" 
# Plot Facecurrents

using LinearAlgebra
fcr1, geo = facecurrents(u1[j2], X2)
fcr2, geo = facecurrents(u1[j1], X3)
plt = Plot([patch(Γ2,norm.(fcr1); cmin = 0.0,cmax=0.005), patch(Γ3,norm.(fcr2);opacity=0.5, cmin = 0.0,cmax=0.005, showscale=false)])
PlotlyDocumenter.to_documenter(plt) #hide
# Plot the Electric Field 

Es1 = potential(E0+E1+E2, pts, u1, Xtot1)
Ein1 = E.(pts)
Etot1 = Es1+Ein1

heatmap(xrange,zrange,real.(getindex.(Etot1,1)))

# ## 2 Cubes inside eachother

# Scattered Electric field inside the domains
E0 = 1im*ω*μ0* S0[j1] + K0[m1] 
E1 = -1im*ω*μ1* S1[j1] - K1[m1]+1im*ω*μ1* S1[j2] + K1[m2]
E2 = -1im*ω*μ2* S2[j2] - K2[m2];

# Scattered Magnetic field inside the domains
H0 = 1im*ω*ϵ0* S0[m1] - K0[j1] 
H1 = -1im*ω*ϵ1* S1[m1] + K1[j1]+ 1im*ω*ϵ1* S1[m2] - K1[j2]
H2 = -1im*ω*ϵ2* S2[m2] + K2[j2];

# Construct Trace Operators
TH0 = Tex[j1] 
TH1 = -Tin[j1] + Tex[j2]
TH2 = -Tin[j2] 

TE0 = Tex[m1] 
TE1 = -Tin[m1] + Tex[m2]
TE2 = -Tin[m2];

# Design System matrix

Z = TH0(H0) + TH1(H1) +  TH2(H2) + TE0(E0) + TE1(E1) + TE2(E2) ;

#rhs

rhs = e[m1]+h[j1];

# Assemble and solve

M2 = assemble(Z,Xtot2,Xtot2);
R2 = assemble(rhs,Xtot2);

u2 = Matrix(M2)\Vector(R2);
ax = nestedrange(Xtot2, 1, numfunctions);
u2 = BlockedVector(u2, (ax,));
@show "system solved" 
# Plot Facecurrents

using LinearAlgebra
fcr1, geo = facecurrents(u2[j2], X2)
fcr2, geo = facecurrents(u2[j1], X1)
plt = Plot([patch(Γ2,norm.(fcr1);cmin = 0.0,cmax=0.005), patch(Γ1,norm.(fcr2);opacity=0.5, cmin = 0.0,cmax=0.005, showscale=false)])
PlotlyDocumenter.to_documenter(plt) #hide
# Plot the Electric Field 

Es2 = potential(E0+E1+E2, pts, u2, Xtot2)
Ein2 = E.(pts)
Etot2 = Es2+Ein2

heatmap(xrange,zrange,real.(getindex.(Etot2,1)))