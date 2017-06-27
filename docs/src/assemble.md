# The Matrix Assemble Routine

A lot of the design of this package derives from the need to express boundary element and finite element matrix assembly in a concise but general manner that is compatible with a wide range of linear and bilinear forms as encountered in the solution of variational problems.

In this section the matrix assembly routine at the center of this package will be discussed. As a case study we will go over the steps required to extend support for new kernels and new finite element spaces.

The matrix assemble routine is surprisingly short:

```julia
function assemblechunk!(biop::IntegralOperator, tfs::Space, bfs::Space, store)

    test_elements, tad = assemblydata(tfs)
    bsis_elements, bad = assemblydata(bfs)

    tshapes = refspace(tfs); num_tshapes = numfunctions(tshapes)
    bshapes = refspace(bfs); num_bshapes = numfunctions(bshapes)

    T = promote_type(scalartype(biop), scalartype(tfs), scalartype(bfs))
    zlocal = zeros(T, num_tshapes, num_bshapes)

    qd = quaddata(biop, tshapes, bshapes, test_elements, bsis_elements)
    for (p,tcell) in enumerate(test_elements), (q,bcell) in enumerate(bsis_elements)

        fill!(zlocal, 0)
        strat = quadrule(biop, tshapes, bshapes, p, tcell, q, bcell, qd)
        momintegrals!(biop, tshapes, bshapes, tcell, bcell, zlocal, strat)

        for j in 1 : num_bshapes, i in 1 : num_tshapes
            z = zlocal[i,j]
            for (m,a) in tad[p,i], (n,b) in bad[q,j]
                store(a*z*b, m, n)
end end end end
```

Support for direct product spaces, linear combinations of kernels, non-standard storage of matrix elements, and parallel execution is provided by layers on top of this assembly routine. In this section we will focus on discussing the design and implementation of this inner building block that lies at the basis of more general functionality.

Finite element spaces are usually stored as a collection of functions that in turn each comprise contributions from a limited number of geometric elements that make up the support of the function. In FEM and BEM matrix assembly, however, we need the *tranposed* information: given a geometric cell, and a local shape function, we need the ability to efficiently retrieve the list of basis functions whose definition contains the given local shape function on the given cell and the weight (aka coefficient) by which it contributes. The data structure that contains this information is referred to as the assembly data `ad`. In particular, `ad[e,s]`, where `e` is the index of a geometric cell and `s` is the index of a local shape function, returns an iterable collection of pairs `(m,w)` where `m` is an index into the iterable collection of basis functions making up the finite element space and `w` is a weight, such that shape function `s` on geometric element `e` contributes with weight `w` to basis function `n`.

```julia
test_elements, tad = assemblydata(tfs)
bsis_elements, bad = assemblydata(bfs)
```

One final note on the `assemblydata` function: as you can see from the above snippet, the function returns, in addition to the actual assembly data, an iterable collection of geometric elements. This collection is a subset of the collection of elements making up the geometry on which the finite element space is defined. The elements returned are those that actually appear in the domain of one or more of the functions that span the finite element space. The double for loop that iterates over pairs of trial and testing functions will only visit those used elements. Elements that are part of the geometry but do not appear as part of the support of a function are skipped. This behaviour is required to guarantee scalability when using multiple threads in assembling the matrix: each thread is assigned a subset of the basis functions; visiting unused elements in all threads is harmful for the overall efficiency.

With this assembly data in hand, matrix assembly can be done by iterating over geometric cells, rather than over basis functions. Doing this avoids visiting a given geometric cell more than once. When computing matrices resulting from discretisation with e.g. Raviart-Thomas elements, this can speed up assembly time with a factor 9.

The problem of matrix assembly is now reduced to the computation of interactions between local shape functions defined on all possible pairs of geometric cells. The space of local shape functions can be retrieved by calling

```julia
tshapes = refspace(tfs); num_tshapes = numfunctions(tshapes)
bshapes = refspace(bfs); num_bshapes = numfunctions(bshapes)
```

Here, `num_tshapes` and `num_bshapes` are the number of local shape functions. For example, when using Raviart-Thomas elements, the number of local shape functions equals three (one for every edge of the reference triangle).

Based on this dimension, and based on the types used to represent numbers in the fields over which the spaces and the kernel are defined, the storage for local shape function interaction is pre-allocated:

```julia
T = promote_type(scalartype(biop), scalartype(tfs), scalartype(bfs))
zlocal = zeros(T, num_tshapes, num_bshapes)
```

Note that the computation of the storage type ensures that high precision or complex data types are only used when required. At all times the minimal storage type is selected. Not only does this keep memory use down, it also results in faster linear algebra computations such as matrix-vector multiplication.

Before entering the double for loop that is responsible for the enumeration of all pairs of geometric cells (a trial cell pairs with a test cell), the implementer is given the opportunity to precompute data for use in the integration kernels. For example when using numerical quadrature rules to compute the double integral in the expression of the matrix entries, it is likely that a set of quadrature points for any given trial cells will be reused in interactions with a large number of test cells. To avoid computing these points and weights over and over, the client developer is given the opportunity to compute and store them by providing an appropriate method for `quaddata` . If memory use is more important the runtime, the client programmer is perfectly allowed to compute points and weights on the fly without storing them.


```julia
fill!(zlocal, 0)
strat = quadrule(biop, tshapes, bshapes, p, tcell, q, bcell, qd)
momintegrals!(biop, tshapes, bshapes, tcell, bcell, zlocal, strat)
```

For a given pair `(tcell,bcell)` of test cell and trial cell (with respective indices `p` and `q` in collections `test_elements` and `bsis_elements`), all possible interactions between local shape functions are computed. After resetting the buffer used to store these interactions, the quadrature strategy is determined. The quadrature strategy in general could depend on:

- the kernel `biop` defining the integral operator,
- the local test and trial shape functions `tshapes` and `bshapes` (functions of high polynomial degree and functions that are highly oscillatory typically require bespoke integration methods),
- and the geometric test and trial cells `tcell` and `bcell` (cells that touch or are near to each other lead to quickly varying or even singular integrands requiring dedicated integration rules).

The method returns an object `strat` that: (i) describes (by its type and its data fields) the integration strategy that is appropriate to compute the current set of local interactions, (ii) contains all data precomputed and stored in `qd` that is relevant to this particular integration (for example a set of quadrature points and weights). This explains why the indices `p` and `q` where passed too `quadrule`: they allow for the quick retrieval of relevant pre-stored data from `qd`.

The routing that is responsible for the actual computation of the interactions between the local shape functions takes the quadrule object `strat` as one of its arguments. The idea is that `momintegrals!` has many methods, not only for different types of kernel and shape functions, but also for different types of `strat`. For example, there are implementations of `momintegrals!` for the computation of the Maxwellian single layer operator w.r.t. spaces of Raviart-Thomas elements that employ double numerical quadrature, singularity extraction, and even more advanced integration routines.

*Note*: the type of `strat` depends on the orientation of the two interacting geometric cells. This information is only available at runtime. In other words, there will be a slight type instability at this point in the code. This is by design however, and not different from the use of virtual functions in an c++ implementation. Numerical experiments show that this form of runtime polymorphism results in negligible runtime overhead.

When all possible interactions between local shape functions have been computed, they need to be stored in the global system matrix. This is done in the matrix assembly loop:

```julia
for j in 1 : num_bshapes
    for i in 1 : num_tshapes
        z = zlocal[i,j]
        for (m,a) in tad[p,i]
            for (n,b) in bad[q,j]
                store(a*z*b, m, n)
            end
        end
    end
end
```
For both the test and trial local shape functions, the global indices at which they appear in the finite element space (and the corresponding weights) are retrieved from the assembly data objects. The contributing value `v = a*z*b` is constructed and its storage is delegated to the `store` method, which we received as one of the arguments passed to `assemble_chunk!`. In the simplest case, `assemble_chunk!` can be used like this:

```julia
Z = zeros(Complex128, numfunctions(tfs), numfunctions(bfs))
store(v, m, n) = (Z[m,n] += v)
assemble_chunk!(kernel, tfs, bfs, store)
```

In other words `store` will simply add the computed value to the specified entry in the global system matrix. Allowing the caller to specify `store` as an argument allows for more flexibility than hardcoding this behaviour in the assembly routine. Indeed, when computing blocks of a larger system, or when e.g. the transposed or a multiple of a given operator is desired, a fairly simple redefinition of `store` can provide this functionality. This is also the reason why `assemble_chunk!` ends in an exclamation mark: even though strictly speaking none of the arguments are modified, the function clearly has an effect on variables defined outside of its scope!

## Case Study: Implementation of the Nitsche Operator Assembly

In the Nitsche method for the Maxwell system, penalty terms are added to the classic discretisation of the EFIE. When discretized using a non-conforming finite elements space (typically because the underlying geometric mesh is not conforming), the penalty term will force the solution to be divergence conforming in some weak sense. The penalty term derives from the following bilinear form:

```math
p(v,u) = \int_{\gamma} v(x) \int_{Γ} \frac{e^{-ik|x-y|}}{4π|x-y|} u(y) dy dx
```

Note that ``u(x)`` is supported by a 2D surface ``Γ`` whereas ``v(y)`` is supported by a 1D curve ``γ``. The complete implementation of this operator could look like

```julia
type SingleLayerTrace{T} <: MaxwellOperator3D
    gamma::T
end

function quaddata(operator::SingleLayerTrace,
    localtestbasis::LagrangeRefSpace,  localtrialbasis::LagrangeRefSpace,
    testelements,  trialelements)

  tqd = quadpoints(localtestbasis,  testelements,  (10,))
  bqd = quadpoints(localtrialbasis, trialelements, (8,))

  return QuadData(tqd, bqd)
end

function quadrule(op::SingleLayerTrace, g::LagrangeRefSpace, f::LagrangeRefSpace, i, τ, j, σ, qd)
    DoubleQuadStrategy(
        qd.tpoints[1,i],
        qd.bpoints[1,j]
    )
end

integrand(op::SingleLayerTrace, kernel, g, τ, f, σ) = f[1]*g[1]*kernel.green
```

Every kernel corresponds with a type. Kernels can potentially depend on a set of parameters; these appear as fields in the type. Here our Nitsche kernel depends on the wavenumber. In quaddata we precompute quadrature points for all geometric cells in the supports of test and trial elements. This is fairly sloppy: only one rule for test and trial integration is considered. A high accuracy implementation would typically compute points for both low quality and high quality quadrature rules.

Also `quadrule` is sloppy: we always select a `DoubleQuadStrategy` to perform the computation of interactions between local shape functions. No singularity extraction or other advanced technique is considered for nearby interactions. Clearly amateurs at work here!

`BEAST` provides a default implementation of an integration routine using double numerical quadrature. All that is required to tap into that implementation is a method overloading `integrand`. From the above formula it is clear what this method should look like.

That's it!
