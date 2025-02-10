
# [Maxwell Single Layer Operator](@id MWsinglelayerDef)

The Maxwell single layer operator, also known als electric field integral operator (EFIO) is encountered in many time-harmonic BEM scattering formulations in electromagnetics.
So far, only the 3D variant is implemented.

---
## Definition

The operator is defined as (see, e.g., [raoElectromagneticScatteringSurfaces1982](@cite))
```math
\bm{\mathcal{T}} \bm b = α \bm{\mathcal{T}}_{\!\!s} \bm b + β \bm{\mathcal{T}}_{\!\!h} \bm b
```
for a vector field ``\bm{b}`` and parameters ``α`` and ``β``, as well as, the weakly singular operator (also known as vector potential operator)
```math
\bm{\mathcal{T}}_{\!\!s} \bm b = \int_\Gamma g_γ(\bm x,\bm y) \, \bm b(\bm y) \,\mathrm{d}\bm y
```
and the hyper singular operator (also known as scalar potential operator)
```math
\bm{\mathcal{T}}_{\!\!h} \bm b = ∇\int_\Gamma g_γ(\bm x,\bm y) \, ∇_Γ⋅\bm b(\bm y) \,\mathrm{d}\bm y
```
with the free-space Green's function
```math
g_{γ}(\bm x,\bm y) = \dfrac{\mathrm{e}^{-γ|x-y|}}{4π|x-y|} \,.
```
The parameters are typically ``α=-\mathrm{j}k``, ``β=-1/(\mathrm{j}k)``, and ``γ = \mathrm{j}k`` with ``k`` denoting the wavenumber and ``\mathrm{j}`` the imaginary unit.


---
## As Bilinear Form

When handed to the [`assemble`](@ref assemble) function, the operators are interpreted as the corresponding bilinear forms
```math
a(\bm t, \bm b) = α ∬_{\Gamma \times \Gamma} \bm t(\bm x) ⋅ \bm b(\bm y) \, g_γ(\bm x,\bm y) \,\mathrm{d}\bm y \mathrm{d}\bm x + β ∬_{Γ×Γ} ∇_Γ⋅\bm t(\bm x) \,  ∇_Γ⋅\bm b(\bm y) \, g_γ(\bm x,\bm y) \,\mathrm{d}\bm y \mathrm{d}\bm x
```
for the complete Maxwell single layer operator,
```math
a_s(\bm t, \bm b) = α ∬_{\Gamma \times \Gamma} \bm t(\bm x) ⋅ \bm b(\bm y) \, g_γ(\bm x,\bm y) \,\mathrm{d}\bm y \mathrm{d}\bm x
```
for the weakly singular part, and
```math
a_h(\bm t, \bm b) = β ∬_{Γ×Γ} ∇_Γ⋅\bm t(\bm x) \,  ∇_Γ⋅\bm b(\bm y) \, g_γ(\bm x,\bm y) \,\mathrm{d}\bm y \mathrm{d}\bm x
```
for the hyper singular part.
Note that the gradient in the hypersingular operator has been moved to the test function using, e.g., a Stokes identity [nedelecWaveEquations2001; p. 73](@cite). 
The corresponding matrices contain the entries
```math
[\bm A]_{mn} =  a_x(\bm t_m, \bm b_n) \,.
```

### API

The weakly singular and the hyper singular operator can be used as such or combined in the single layer operator.

!!! tip
    The qualifier `Maxwell3D` has to be used. 


```@docs
Maxwell3D.singlelayer
Maxwell3D.weaklysingular
Maxwell3D.hypersingular
```

---
## As Linear Map

When handed to the [`potential`](@ref potential) function, the operator can be evaluated at provided points in space.
Commonly, this is used in post-processing.


### Far-Field

In the limit that the observation point ``\bm x \rightarrow \infty``, the operator simplifies to the far-field (FF) version
```math
(\bm{\mathcal{T}}_\mathrm{FF} \bm b)(\bm x) = α  \bm{u}_r \times \int_\Gamma \bm{b}(\bm y) \mathrm{e}^{\mathrm{j}\bm{u}_r \cdot\, \bm{y}} \,\mathrm{d}\bm{y} \times \bm{u}_r 
```
where ``\bm{u}_r`` is the unit vector in the direction of the evaluation point.


### API

```@docs
MWSingleLayerField3D
MWFarField3D
```

!!! tip
    The provided points for the [`potential`](@ref potential) should be in Cartesian coordinates. The returned fields are also in Cartesian from.

!!! warning
    Singularities are not addressed: the case when the evaluation point is close to the surface is not treated properly, so far.