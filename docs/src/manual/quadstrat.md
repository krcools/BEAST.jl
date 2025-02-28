# Customising the choice of quadrature rules

At the heart of any boundary element code are routines that can compute matrix entries that take on the form of double integrals over pairs of panels of an integrand that is the product of a test function, a trial function, and the fundamental solution of the equation under study. The fundamental solution is singular at the origin, implying that simple quadrature methods are not sufficient to compute the entries to satisfactory accuracy. A mechanism to activate the most appropriate rule is required. BEAST.jl attempts to offer a reasonable default, but the advanced user may wish to take control beyond this.

This page describes how the user can intervene in this system.

## List of implemented quadrature strategies

The algorithm that determines which quadrature rule is used for any given geometric constellation of interacting panels is called the quadrature strategy. The list of strategies currently implemented in BEAST is:

```@example introductory
using TypeTree
using BEAST

print(join(tt(BEAST.AbstractQuadStrat), ""))
nothing # hide
```

The type of the quadrature strategy object determines the quadrature rule selection algorithm. The value of the fields contained by the quadrature strategy object typically selects the order or number of quadrature points used by the various possible quadrature rules.

The methods responsible for caching quadrature related data and quadrature rule selection are named `quaddata` and `quadrule`, respectively. For a given combination of a discrete boundary integral operator and quadrature strategy, the methods that will be dispatched to can be queried as follows:

```@example introductory
using CompScienceMeshes
using BEAST

Γ  = meshsphere(radius=1.0, h=0.45) # triangulate sphere of radius one
RT = raviartthomas(Γ)
𝑇 = Maxwell3D.singlelayer(wavenumber=2.0)
qs = BEAST.DoubleNumWiltonSauterQStrat(6, 7, 7, 8, 6, 6, 6, 6)
```

```@example introductory
BEAST.quadinfo(𝑇, RT, RT; quadstrat=qs)
nothing # hide
```

## Explicitly providing the quadrature strategy

The assembly function takes a keyword argument that allows to specify a specific quadrature strategy to use.

```@example introductory
using CompScienceMeshes
using BEAST

Γ  = meshsphere(radius=1.0, h=0.45)   # triangulate sphere of radius one
RT = raviartthomas(Γ)
𝑇 = Maxwell3D.singlelayer(wavenumber=2.0)

Z1 = assemble(𝑇, RT, RT; quadstrat=BEAST.DoubleNumWiltonSauterQStrat(2, 3, 6, 7, 5, 5, 5, 5))
Z2 = assemble(𝑇, RT, RT; quadstrat=BEAST.DoubleNumWiltonSauterQStrat(6, 7, 7, 8, 6, 6, 6, 6))

using LinearAlgebra
2 * norm(Z1-Z2) / norm(Z1+Z2)
```

This method works well when only a single integral operator appears in the boundary integral equation. For systems containing multiple equations or when linear combinations of different operators appear, specifying a single quadrature strategy at the callsite of `assemble` will likely not be appropriate.


## Setting the default quadrature rule

To query the set default quadrature strategy for a triple `(op, testfns, trialfns)`,

```@example introductory
using CompScienceMeshes
using BEAST

Γ  = meshsphere(radius=1.0, h=0.45)   # triangulate sphere of radius one
RT = raviartthomas(Γ)
𝑇 = Maxwell3D.singlelayer(wavenumber=2.0)

BEAST.defaultquadstrat(𝑇, RT, RT)
```

A new default can be set using the `@defaultquadstrat` macro. This creates a new method for the function `defaultquadstrat` that will be dispatched to for the arguments of the provided types.

```@example introductory
BEAST.@defaultquadstrat (𝑇, RT, RT) BEAST.DoubleNumWiltonSauterQStrat(6, 7, 7, 8, 6, 6, 6, 6) 
Z2 = assemble(𝑇, RT, RT)

BEAST.@defaultquadstrat (𝑇, RT, RT) BEAST.DoubleNumWiltonSauterQStrat(2, 3, 6, 7, 5, 5, 5, 5) 
Z1 = assemble(𝑇, RT, RT)

using LinearAlgebra
2 * norm(Z1-Z2) / norm(Z1+Z2)
```

The above number provides some insight into the accuracy of the selected quadrature strategy. The advantage of setting the default quadrature strategy like this is that it is global and will automatically affect all subsequent calls to assembly. This means the set strategy will be used even if assembly is called as part of the assembly of a more complicated larger system, potentially containing linear combinations of integral operators.

The downside is that the default can only be set for concrete types of `(op, testfns, trialfns)` and that the call to the macro needs to be at Module scope. This makes it less appealing for complicated simulations and sweeps over the quadrature accuracy.

# Specifying a quadrature strategy selection method

This method is somehwat more involved, but is the most general and should allow for full control of the quadrature rules used in assembly.

```@example introductory
function myquadstrat1(op, testfns, trialfns)
    if op isa BEAST.MWSingleLayer3D
        return BEAST.DoubleNumWiltonSauterQStrat(2, 3, 6, 7, 5, 5, 5, 5)
    end
    return BEAST.defaultquadstrat(op, testfns, trialfns)
end

function myquadstrat2(op, testfns, trialfns)
    if op isa BEAST.MWSingleLayer3D
        return BEAST.DoubleNumWiltonSauterQStrat(6, 7, 7, 8, 6, 6, 6, 6)
    end
    return BEAST.defaultquadstrat(op, testfns, trialfns)
end


Z1 = assemble(𝑇, RT, RT; quadstrat=myquadstrat1)
Z3 = assemble(𝑇, RT, RT; quadstrat=myquadstrat2)
2 * norm(Z1-Z2) / norm(Z1+Z2)
```

Best practice is too return `BEAST.defaultquadstrat(op, testnfs, trialfns)` by default to ensure that all operators are supported, also those for which no explicit overwrite is specified.


