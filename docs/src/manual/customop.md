# Custom Operators

BEAST.jl provides the single and double layer operators for the 3D Laplace operator, the 2D and 3D Helmholtz equations, and the 3D Maxwell equations. Nevertheless there are situations where the user may want to define their own integral operator. BEAST.jl makes this process as easy as possible.

Defining a new integral operator can be as simple as defining a type representing the operator, a method for the evaluation of the integral operator's integrand, and a method for `BEAST.scalartype` so the framework can determine the most appropriate storage type.

```@example customop
using LinearAlgebra
using CompScienceMeshes
using BEAST

struct BellCurveOp{T} <: BEAST.IntegralOperator
    width::T
end
```

By default, assembly of new operators is dealt using the `DoubleNumSauterQstrat` quadrature strategy. This strategy selects either a simplex tensorial quadrature rule when triangles are well-separated, or a bespoke Sauter-Schwab rule when triangles have vertices in common. These rules are to a large extend independent of the integrand. When defining operators exhibiting highly irregular integrands, or when using panels for which the mentioned quadrature rules are not applicable, additional implementation may be required.

At the heart of the operator defintion lies the routine that allows evaluation of the integrand. This should be provided in the following format:

```@example customop
function (igd::BEAST.Integrand{<:BellCurveOp})(p,q,f,g)
    α = igd.operator.width

    x = CompScienceMeshes.cartesian(p)
    y = CompScienceMeshes.cartesian(q)
    R = LinearAlgebra.norm(x-y)

    BEAST._integrands(f,g) do fi,gj
        dot(fi.value, gj.value) * exp(-(R/α)^2)
    end
end
```

Note that `(p,q)` are neighborhoods (see CompScienceMeshes.jl). The convenience function `BEAST._integrands` populates a `StaticArrays.SMatrix` with all posible combinations of the operator integrans with trial and test functions without allocating memory.

Finally, a method declaring the scalar type used by the operator needs to be defined so that BEAST.jl can provide the most efficient storage for any combination of real or complex values operators and finite element spaces.

```@example customop
BEAST.scalartype(op::BellCurveOp{T}) where {T} = T
```

```@example customop
Γ = CompScienceMeshes.meshsphere(radius=1.0, h=0.3)
X = BEAST.raviartthomas(Γ)

op = BellCurveOp(2.0)
Z = assemble(op, X, X)

sv = svdvals(Z)
import PlotlyBase
t1 = PlotlyBase.scatter(y=sv)
plt = PlotlyBase.Plot([t1])
import PlotlyDocumenter # hide
PlotlyDocumenter.to_documenter(plt) #hide
```

Being a Hilbert-Schmidt operator, the quickly decaying spectrum is to be expected.