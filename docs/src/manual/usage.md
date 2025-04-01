
# [General Usage](@id usageRef)

!!! info
    The fundamental approach, which applies in most cases is:
    1. Define trial and test functions.
    2. Define an operator and an excitation.
    3. Assemble the system matrix and the right-hand side

The available basis functions and corresponding geometry representations, available [operators](@ref operator), and excitations are defined in the corresponding sections of this documentation.


---
## Introductory Example: EFIE

The fundamental procedure is exemplified for the electric field integral equation (EFIE); further common steps are discussed afterwards:

```@setup introductory
import PlotlyBase
import PlotlyDocumenter
```

```@example introductory
using CompScienceMeshes
using BEAST

# --- 1. basis functions
Γ  = meshsphere(radius=1.0, h=0.4)   # triangulate sphere of radius one
RT = raviartthomas(Γ)                # define basis functions

# --- 2. operators & excitation
𝑇 = Maxwell3D.singlelayer(wavenumber=2.0)                             # integral operator
𝐸 = Maxwell3D.planewave(direction=x̂, polarization=ẑ, wavenumber=2.0)  # excitation
𝑒 = (n × 𝐸) × n # tangential part

# --- 3. compute the RHS and system matrix
e = assemble(𝑒, RT)         # assemble RHS
T = assemble(𝑇, RT, RT)     # assemble system matrix
nothing #hide
```

---
### Explanation

The example follows the 3 steps.
Specifically, in the example trial and test functions are the same ([Raviart-Thomas](@ref raviartthomasDef)), the excitation is a [plane wave](@ref planewaveEx), and the operator is the [Maxwell single layer operator](@ref MWsinglelayerDef).

!!! tip
    The [`assemble`](@ref assemble) function is the key function of this package, it accepts either *excitation + test function* or *operator + test + trial function*.

    The operator can also be a linear combination of several operators.

```@docs
assemble
```

---
### Further Common Steps

The linear system of equations can now, for example, be solved via the iterative GMRES solver of the [Krylov.jl](https://github.com/JuliaSmoothOptimizers/Krylov.jl) package.
However, other solver could be used.
Subsequently, different post-processing steps can be conducted, such as computing the scattered field from the determined expansion coefficients.

!!! tip
    Key functions for the post-processing are the [`potential`](@ref potential) and [`facecurrents`](@ref facecurrents) functions.

This is shown in the following:


```@example introductory
using Krylov, LinearAlgebra

# --- solve linear system iteratively
u, ch = Krylov.gmres(T, -e, rtol=1e-5)

fcr, geo = facecurrents(u, RT)
pt = CompScienceMeshes.patch(Γ, norm.(fcr))
pl = PlotlyBase.Plot(pt)
PlotlyDocumenter.to_documenter(pl) # hide
```

```@example introductory
# --- post processing: compute scattered electric field at two Cartesian points
points = [[3.0, 4.0, 2.0], [3.0, 4.0, 3.0]]
EF = potential(MWSingleLayerField3D(gamma=im*2.0), points, u, RT)
```


!!! warning
    more details

```@raw html
<figure>
  <img
    src="../../assets/currentREADME.png"
    alt="Setup"
    width="250" />
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  <img
    src="../../assets/currentRealREADME.png"
    alt="Setup"
    width="250" />
</figure>
<br/>
```

---
## Plotting & Exporting

Plotting details.

Export VTK.

...