
# [Identity Operator](@id identityDef)

The identity operator is implemented in two flavors `BEAST.Identity()` and `BEAST.NCross()`.

---
## Definition

The identity operator
```math
\bm{\mathcal{I}} \bm b = \bm b
```
returns the function it is provided unchanged.


### As bilinear form

When handed to the [`assemble`](@ref assemble) function, the operator is interpreted as the bilinear form
```math
a(\bm t, \bm b) = \int_\Gamma \bm t(\bm x)  ⋅ \bm{\mathcal{I}} \bm b(\bm x) \,\mathrm{d}\bm x = \int_\Gamma \bm t(\bm x)  ⋅ \bm b(\bm x) \,\mathrm{d}\bm x \,.
```
Hence, the resulting matrix ``\bm A`` contains the entries
```math
[\bm A]_{mn} =  \int_\Gamma \bm t_m(\bm x)  ⋅ \bm b_n(\bm x) \,\mathrm{d}\bm x \,.
```

### API

```@docs; canonical=false
BEAST.Identity
```


---
## Variant with ``\bm n \times``

As a variation, the identity operator is provided with a cross product with the normal vector ``\bm n`` of the surface ``\Gamma``.

### The bilinear form

The [`assemble`](@ref assemble) function, interprets the operator as the bilinear form
```math
a(\bm t, \bm b) = \int_\Gamma \bm t(\bm x)  ⋅ \bm{\mathcal{I}} (\bm n(\bm x) \times \bm b(\bm x)) \,\mathrm{d}\bm x = \int_\Gamma \bm t(\bm x)  ⋅ (\bm n(\bm x) \times \bm b(\bm x)) \,\mathrm{d}\bm x \,.
```
Hence, the resulting matrix ``\bm A`` contains the entries
```math
[\bm A]_{mn} =  \int_\Gamma \bm t_m(\bm x)  ⋅ (\bm n(\bm x) \times \bm b_n(\bm x)) \,\mathrm{d}\bm x \,.
```

### API

```@docs; canonical=false
BEAST.NCross
```