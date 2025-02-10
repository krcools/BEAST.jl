# Quadrature strategies

There are many ways to approximately compute the singular integrals that appear in boundary element discretisations of surface and volume integral euqations.

BEAST.jl is configured to select reasonable defaults, but advanced users may want to select their own quadrature rules. This section provides information on how to do this.

!!! warning
    TODO: discuss directly available/implemented singularity treatment

## quaddata and quadrule

Numerical quadrature is governed by a pair of functions that need to be designed to work together:

- `quaddata`: this function is executed before the assembly loop is entered. It's job is to compute all data needed for quadrature that the developer wants to be cached. Typically this is all geometric information such as the parametric and cartesian coordinates of all quadratures rules for all elemenents. It makes sense to cache this data as it will be used many times over in the double for loop that governs assembly. Typically, near singular interactions require more careful treatment than far interactions. This means that multiple quadrature rules per elements can be required. In such cases, the developer may want to opt to sture quadrature points and weights for all these rules. The function returns a quaddata object that holds all the cached data.
- `quadrule`: quadrule is executed inside the assembly hotloop. It receives a pair of elements and the quaddata object as its arguments. Based on this, the relevant cached data is extracted and stored in a quadrule object. The type of this object will determine the actual quadrature routined that will be called upon to do the numerical quadrature.

## quadstrat

The pair of quaddata and quadrule methods that is used is determined by the type of the operator and finite elements, and a `quadstrat` object. This object is passed to the assembly routine and passed on to `quaddata` and `quadstrat`, so it can be considered during dispatch.

Parameters, such as those that determine the accuracy of the numerical quadrature, are part of the runtime payload of the quadstrat object. This is usefull when the user is interested on the impact of these parameters on the performance and the accuracy of the solver without having to supply a new pair of `quadstrat`/`quaddata` methods for each possible value of these parameters.

Roughly this leads to the following (simplified) assembly routine:

```julia
function assemble(op, tfs, bfs, store; quadstrat=QS)

    tad, tels = assemblydata(op, tfs)
    bad, bels = assemblydata(op, bfs)
    
    qd = quadata(op,tels,bels,quadstrat)
    for tel in tels
        for bel in bels
            qr = quadrule(op,tel,bel,qd,quadstrat)
            zlocal = momintegrals(op,tel,bel,qr)

            for i in axes(zlocal,1)
                for j in axes(zlocal,2)
                    m, a = tad[tel,i]
                    n, b = bad[bel,j]
                    store(a*zlocal[i,j]*b,m,n)
end end end end end
```

It is conceivable that the types and functions described above look like this:

```julia
struct DoubleNumQS
    test_precision
    trial_precision
end

function quaddata(op, tels, bels, quadstrat::DoubleNumQS)
    tqps = [quadpoints(tel,precision=quadstrat.test_precision) for tel in tels]
    bqps = [quadpoints(bel,precision=quadstrat.basis_precision) for bel in bels]
    return (test_quadpoints=tqps, basis_quadpoints=bqps)
end


struct DoubleNumQR
    test_quadpoints
    trial_quadpoints
end

struct HighPrecisionQR end

function quadrule(op, tel, bel, qd, quadstrat::DoubleNumQs)
    if wellseparated(tel, bel)
        return DoubleNumQR(qd.test_quadpoints[tel], qd.basis_quadpoints[bel])
    else
        return HighPrecisionQR(tel, bel)
    end
end

function momintegrals(op, tel, bel, qr::DoubleNumQR)
    ...
end

function momintegrals(op, tel, bel, qr::HighPrecisionQR)
    ...
end
```