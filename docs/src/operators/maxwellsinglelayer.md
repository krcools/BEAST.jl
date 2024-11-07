
# [Maxwell Single Layer Operator](@id MWsinglelayerDef)

## Definition

```math
\bm{\mathcal{T}}
```


When applied to 
```math
a(\bm t, \bm b) = α ∬_{\Gamma \times \Gamma} \bm t(\bm x) ⋅ \bm b(\bm y) \, g_γ(\bm x,\bm y) \,\mathrm{d}\bm y \mathrm{d}\bm x + β ∬_{Γ×Γ} ∇_Γ⋅\bm t(\bm x) \,  ∇_Γ⋅\bm b(\bm y) \, g_γ(\bm x,\bm y) \,\mathrm{d}\bm y \mathrm{d}\bm x
```

with the free-space Green's function

```math
g_{γ}(\bm x,\bm y) = \dfrac{\mathrm{e}^{-γ|x-y|}}{4π|x-y|}
```

## API

```@docs
Maxwell3D.singlelayer
```