# Geometry

```@setup geo
import PlotlyBase
import PlotlyDocumenter
```

```@example geo
using CompScienceMeshes
Γ = meshsphere(radius=1.0, h=0.35)
pt = CompScienceMeshes.patch(Γ)
pl = PlotlyBase.Plot(pt)
PlotlyDocumenter.to_documenter(pl) # hide
```
