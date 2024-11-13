
# [Maxwell Double Layer Operator](@id MWdoublelayerDef)

The Maxwell double layer operator is encountered in many time-harmonic BEM scattering formulations in electromagnetics.
So far, only the 3D variant is implemented.

---
## Definition

The operator is defined as (see, e.g., ...)
```math
\bm{\mathcal{K}} \bm b = α \int_\Gamma ∇_{\!x} g_γ(\bm x,\bm y) \times \bm b(\bm y) \,\mathrm{d}\bm y
```
for a vector field ``\bm{b}`` and a parameter ``α`` with the free-space Green's function
```math
g_{γ}(\bm x,\bm y) = \dfrac{\mathrm{e}^{-γ|x-y|}}{4π|x-y|} \,.
```
The parameters are typically ``α=1`` and ``γ = \mathrm{j}k`` with ``k`` denoting the wavenumber and ``\mathrm{j}`` the imaginary unit.
As variation, the rotaded double layer operator
```math
\bar{\bm{\mathcal{K}}} = \bm{n} \times \bm{\mathcal{K}} 
```
can be used which simply involves the cross product with the normal vector ``\bm{n}`` of the surface ``\Gamma``.


---
## As Bilinear Form

When handed to the [`assemble`](@ref assemble) function, the operators are interpreted as the corresponding bilinear forms
```math
a(\bm t, \bm b) = α ∬_{\Gamma \times \Gamma} \bm t(\bm x) ⋅  \,( ∇_{\!x} g_γ(\bm x,\bm y) \times \bm b(\bm y) ) \,\mathrm{d}\bm y \mathrm{d}\bm x
```
and
```math
\bar{a}(\bm t, \bm b) = α ∬_{\Gamma \times \Gamma} \bm t(\bm x) ⋅  \, (\bm{n} \times ( ∇_{\!x} g_γ(\bm x,\bm y) \times \bm b(\bm y) )) \,\mathrm{d}\bm y \mathrm{d}\bm x
```
resulting in the matrix
```math
[\bm A]_{mn} =  a(\bm t_m, \bm b_n) \,.
```


### API

```@docs
Maxwell3D.doublelayer
```


---
## As Linear Map

When handed to the [`potential`](@ref potential) function, the operator can be evaluated at provided points in space.
Commonly, this is used in post-processing.


### Far-Field

In the limit that the observation point ``\bm x \rightarrow \infty``, the operator simplifies to the far-field (FF) version
```math
(\bm{\mathcal{K}}_\mathrm{FF} \bm b)(\bm x) = α \bm{u}_r \times \int_\Gamma \bm{b}(\bm y) \mathrm{e}^{\mathrm{j}\bm{u}_r \cdot\, \bm{y}} \,\mathrm{d}\bm{y}
```

### API

```@docs
BEAST.MWDoubleLayerField3D
BEAST.MWDoubleLayerFarField3D
MWDoubleLayerRotatedFarField3D
```

!!! tip
    The provided points for the [`potential`](@ref potential) should be in Cartesian coordinates. The returned fields are also in Cartesian from.

!!! warning
    Singularities are not addressed: the case when the evaluation point is close to the surface is not treated properly, so far.