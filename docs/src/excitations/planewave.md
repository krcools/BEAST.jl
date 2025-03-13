
# [Plane Wave Excitation](@id planewaveEx)

A plane wave can be used as excitation in time-harmonic and time-domain scenarios.

---
## Time-Harmonic

A time-harmonic plane wave with amplitude ``a``, wave vector ``\bm k = k \hat{\bm k}``, and polarization ``\hat{\bm p}`` (vectors with a hat denote unit vectors) is defined by the field
```math
\bm e_\mathrm{PW}(\bm x) = a \hat{\bm p}  \, \mathrm{e}^{-\mathrm{j} \bm k \cdot \bm x}  \,,
```
where the polarization and wave vector are orthogonal, that is,
```math
\bm k \cdot \hat{\bm p} = 0
```
holds.

### API

!!! warning
    There are currently two APIs for the BEM plane wave. Fix in future.

```@docs
BEAST.planewavemw3d
Maxwell3D.planewave
BEAST.planewavevie
```

---
## Time-Domain

A plane wave in the time-domain is defined as ...


### API

Some more details would be helpful here.

```@docs
BEAST.planewave
```