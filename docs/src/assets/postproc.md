# Post-processing and Visualisation

## Field computation

The main API for the computation of fields in a vector of points is the potential function. This function is capable of computing

```math
F(x) = \int_{\Gamma} K(x,y) u(y) dy
```

with ``u(y) = \Sigma_{i=1}^N u_i f_i(y)`` and ``K(x,y)`` the integation kernel defining the type of potential. For example, the far field for a vector valued surface density is

```math
F(x) = \int_{\Gamma} e^{ik \frac{x \cdot y}{|x|}} u(y) dy
```

For a far field potential such as this, the value only depends on the direction, not the magnitude, of ``x``, as can be read off from the normalisation in the exponent.

The following script computes the far field along a semi-circle in the xz-plane.

```julia
Θ, Φ = range(0.0,stop=2π,length=100), 0.0
dirs = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ for ϕ in Φ]
farfield = potential(MWFarField3D(wavenumber=κ), dirs, u, X)
```
