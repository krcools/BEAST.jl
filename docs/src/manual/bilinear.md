# Systems of boundary integral equations and bilinear forms

For the most simple variational formulations such as the EFIE and MFIE, the boundary element matrices and right hand sided can be manually assembled.

For more complex formulations, such as those encountered in solving the transmission problem and in problems involving composite systems, this approach quickly becomes unwieldy. To facilitate the formulation and solution of these problems, BEAST.jl provides the ability to define quite general linear and bilinear forms.

## The single body transmission problem

THe transmission problem is defined by the material properties of the exterior and interior domain and the incident field:

```@example transmission
using CompScienceMeshes
using BEAST

κ1 = 1.0
κ2 = 2.0

E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ1)
H = -1/(im*κ1)*curl(E)
nothing # hide
```

The PMCHWT can be used to model this transmission problem. To write down the variational formulation, pairs of placeholders for the trial functions and test functions are declared by the `@hilbertspace` macro.

```@example transmission
T1  = Maxwell3D.singlelayer(wavenumber=κ1)
K1  = Maxwell3D.doublelayer(wavenumber=κ1)
T2  = Maxwell3D.singlelayer(wavenumber=κ2)
K2  = Maxwell3D.doublelayer(wavenumber=κ2)

e = (n × E) × n
h = (n × H) × n

@hilbertspace j m
@hilbertspace k l

a =
    (T1+T2)[k,j] + (-K1-K2)[k,m] +
    (K1+K2)[l,j] + (T1+T2)[l,m]
b =
    -e[k] -h[l]
```

Note that we are able to define the formulation without specifying or constructing the mesh, the boundary element spaces, or the system matrices. The above constitutes a lazy representation of the formulation that takes neglible time to construct! 

Now this formulation is applied to a concrete geometry.

```@example transmission
Γ = meshsphere(;radius=1.0, h=0.35)
RT = raviartthomas(Γ)

X = BEAST.DirectProductSpace([RT,RT])

bx = assemble(b, X)
Axx = assemble(a, X, X)
```

The type of `Axx` reveals that the structure of the underlying Hilbert space (with two generators in this case) is preserved. Also noteworthy is that only non-zero blocks in the system matrix are stored. This is important both for to limit memory consumption and to avoid unnecessary computations in performing the matrix-vector product.

BEAST.jl provides wrappers for iterative solvers that conform to the `LinearMaps.jl` interface. For example, `SXX = GMRESSolver(Axx)` acts like a `LinearMap` representing the inverse of `Axx`. Computing the action `uX = SXX * bx` runs the Krylov iterative process with right hand side `bx` and returns the solution to the linear system:

```@example transmission
SXX = BEAST.GMRESSolver(Axx)
uX = SXX * bx
typeof(uX)
```

The solution of this iterative process *remembers* its block structure. This makes it easy to select the electric and magnetic components:

```@example transmission
import PlotlyBase
import PlotlyDocumenter # hide
using LinearAlgebra

fcrm, geom = facecurrents(uX[m], RT)
fcrj, geoj = facecurrents(uX[j], RT)

ptm = CompScienceMeshes.patch(geom, norm.(fcrm); caxis=(0,1.2) , showscale=false)
ptj= CompScienceMeshes.patch(geoj, norm.(fcrj); caxis=(0,1.2))

pl = [ PlotlyBase.Plot(ptm) PlotlyBase.Plot(ptj) ]
PlotlyDocumenter.to_documenter(pl) # hide
```

## Calderon preconditioning for the PMCHWT

Let's try to speed up convergence of the iterative solver by constructing a Calderon preconditioner. We opt for a block diagonal preconditioner containing only single layer contributions. The wavenumber used for construction of the preconditioner is purely imaginary to avoid the introduction of resonances.

```@example transmission
Ty  = Maxwell3D.singlelayer(gamma=κ1)

c = Ty[k,j] + Ty[l,m]
```

We also require discrete versions of the duality pairing to map primal test coefficients to dual expansion coefficients

```@example transmission
Nx = BEAST.NCross()
d = Nx[k,j] + Nx[l,m]
```

The corresponding block matrices can be built as before by appealing to the `assemble` function

```@example transmission
BC = buffachristiansen(Γ)
Y = BEAST.DirectProductSpace([BC,BC])

Cyy = assemble(c, Y, Y)
Dxy = assemble(d, X, Y)
```

A *lazy* inverse for the discrete duality can be constructed by wrapping the block matrix in a Krylov solver object:

```@example transmission
DYX = BEAST.GMRESSolver(Dxy, verbose=false)
DXY = BEAST.GMRESSolver(Dxy', verbose=false)
```

We have now all the ingredients to write down the system we want to solve:

```@example transmission
PAXx = DXY * Cyy * DYX * Axx
PbX = DXY * Cyy * DYX * bx

PSXx = BEAST.GMRESSolver(PAXx)

u1, ch1 = solve(SXX, bx)
u2, ch2 = solve(PSXx, PbX)

using LinearAlgebra
norm(u1-u2), ch1.iters, ch2.iters
```

Even for this small example, the solution was reconstructed in a much smaller number of iterations. Because objects of type `GMRESSolver` effectively behave as the inverse of the wrapped `LinearMap`, an internal iterative solver to compute the action of the inverse of the duality pairing is triggered during each outer iteration without introducing notational overhead!