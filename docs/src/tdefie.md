# Solving the Time Domain EFIE using Marching-on-in-Time

If broadband information is required or if the system under study will be coupled to non-linear components, the scattering problem should be solved directly in the time domain, i.e. as a hyperbolic evolution problem.

## Building the geometry

Building the geometry and defining the spatial finite elements happens in completely the same manner as for frequency domain simulations:

```julia
using CompScienceMeshes
using BEAST
D, dx = 1.0, 0.3
Γ = meshsphere(1.0, dx)
X = raviartthomas(Γ)
nothing # hide
```

Time domain currents are approximated in the tensor product of a spatial and temporal finite element space:

$j(x) \approx \sum_{i=1}^{N_T} \sum_{m=1}^{N_S} u_{i,m} T_i(t) f_m(x)$

This package only supports translation invariant temopral basis functions, i.e.

$T_i(t) = T(t - i \Delta t),$

where $\Delta t$ is the time step used to solve the problem. This time step depends on the bandwidth of the incident field and the desired accuracy. Space-Time Galerkin solvers of the type used here are not subject to stability conditions linking spatial and temporal discretisation resolutions.

We need a temporal trial space and test space. Common examples of temporal trial spaces are the shifted quadratic spline and shifted lagrange basis functions. In this example, a shifted quadratic spline `S` is used for the trial space, while a delta function `U` is used as the temporal test space to obtain a time-stepping solution.  

```julia
Δt, Nt = 0.25, 300
S = BEAST.timebasisspline2(Δt, Nt)
U = BEAST.timebasisdelta(Δt, Nt)
nothing # hide
```

We want to solve the EFIE, i.e. we want to find the current $j$ such that

$Tj = -e^i,$

where the incident electric field can be any Maxwell solution in the background medium. To describe this problem in Julia we create a retarded potential operator objects and a functional representing the incident field:

```julia
x = point(1.0,0.0,0.0)
y = point(0.0,1.0,0.0)
z = point(0.0,0.0,1.0)
gaussian = BEAST.creategaussian(30Δt, 60Δt)
E = BEAST.planewave(x, z, BEAST.derive(gaussian), 1.0)
T = BEAST.MWSingleLayerTDIO(1.0,-1.0,-1.0,2,0)
nothing; # hide
```

Using the finite element spaces defined above this retarded potential equation can be discretized.

```julia
V = X ⊗ S #Space and time trial basis
W = X ⊗ U #Space and time test basis

B = assemble(E, W)
Z = assemble(T, W, V)
nothing # hide
```

The variable `Z` can efficiently store the matrices corresponding to different delays (`assemble` knows about the specific sparsity pattern of such matrices and returns a sparse array of rank three fit for purpose).

The algorithm below solves the discrete convolution problem by marching on in time:

```julia
Z0 = Z[:,:,1]
W0 = inv(Z0)
x = BEAST.marchonintime(W0,Z,-B,Nt)
nothing # hide
```
Computing the values of the induced current is now possible in the same manner as for frequency domain simulations by first converting our MOT solution back to the frequency domain using the fourier transform, along with some adjustments.

```julia
Xefie, Δω, ω0 = fouriertransform(x, Δt, 0.0, 2)
ω = collect(ω0 + (0:Nt-1)*Δω)
_, i1 = findmin(abs(ω-1.0))
ω1 = ω[i1]

ue = Xefie[:,i1]
ue /= fouriertransform(gaussian)(ω1)
fcr, geo = facecurrents(ue, X)
nothing # hide
```

For now the package still relies upon Matlab for some of its visualisation. This dependency will be removed in the near future:

```julia
include(Pkg.dir("CompScienceMeshes","examples","matlab_patches.jl"))
mat"clf"
patch(geo, real.(norm.(fcr)))
mat"cd($(pwd()))"
mat"print('current.png', '-dpng')"
```
