```@meta
CurrentModule = BEAST
```

# BEAST.jl documentation

BEAST provides a number of types modelling concepts and a number of algorithms for the efficient and simple implementation of boundary and finite element solvers. It provides full implementations of these concepts for the LU based solution of boundary integral equations for the Maxwell and Helmholtz systems.

Because Julia only compiles code at execution time, users of this library can hook into the code provided in this package at any level. In the extreme case it suffices to provide overwrites of the `assemble` functions. In that case, only the LU solution will be performed by the code here.

At the other end it suffices that users only supply integration kernels that act on the element-element interaction level. This package will manage all required steps for matrix assembly.

For the Helmholtz 2D and Maxwell 3D systems, complete implementations are supplied. These models will be discussed in detail to give a more concrete idea of the APIs provides and how to extend them.

Central to the solution of boundary integral equations is the assembly of the system matrix. The system matrix is fully determined by specifying a kernel G, a set of trial functions, and a set of test functions.

## Basis

Sets of both trial and testing functions are implemented by models following the basis concept. The term basis is somewhat misleading as it is nowhere required nor enforced that these functions are linearly independent. Models implementing the Basis concept need to comply to the following semantics.


- [`numfunctions(basis)`](@ref): number of functions in the Basis.
- [`coordtype(basis)`](@ref): type of (the components of) the values taken on by the functions in the Basis.
- [`scalartype(d)`](@ref): the scalar field underlying the vector space the basis functions take value in.
- [`refspace(basis)`](@ref): returns the ReferenceSpace of local shape functions on which the Basis is built.
- [`assemblydata(basis)`](@ref): `assemblydata` returns an iterable collection `elements` of geometric elements and a look table `ad` for use in assembly of interaction matrices. In particular, for an index `element_idx` into `elements` and an index `local_shape_idx` in basis of local shape functions `refspace(basis)`, `ad[element_idx, local_shape_idx]` returns the iterable collection of `(global_idx, weight)` tuples such that the local shape function at `local_shape_idx` defined on the element at `element_idx` contributes to the basis function at `global_idx` with a weight of `weight`.
- [`geometry(basis)`](@ref): returns an iterable collection of Elements. The order in which these Elements are encountered corresponds to the indices used in the assembly data structure.


## Reference Space

The *reference space* concept defines an API for working with spaces of local shape functions. The main role of objects implementing this concept is to allow specialization of the functions that depend on the precise reference space used.

The functions that depend on the type and value of arguments modeling *reference space* are:

- [`numfunctions(refspace)`](@ref): returns the number of shape functions on each element.

## Kernel

A kernel is a fairly simple concept that mainly exists as part of the definition of a Discrete Operator. A kernel should obey the following semantics:

In many function definitions the kernel object is referenced by `operator` or something similar. This is a misleading name as an operator definition should always be accompanied by the domain and range space.

## Discrete Operator

Informally speaking, a Discrete Operator is a concept that allows for the computation of an interaction matrix. It is a kernel together with a test and trial basis. A Discrete Operator can be passed to `assemble` and friends to compute its matrix representation.

A discrete operator is a triple `(kernel, test_basis, trial_basis)`, where `kernel` is a Kernel, and `test_basis` and `trial_basis` are Bases. In addition, the following expressions should be implemented and behave according to the correct semantics:

- [`quaddata(operator,test_refspace,trial_refspace,test_elements,trial_elements)`](@ref): create the data required for the computation of element-element interactions during assembly of discrete operator matrices.
- [`quadrule(operator,test_refspace,trial_refspace,p,test_element,q_trial_element,qd)`](@ref): returns an integration strategy object that will be passed to `momintegrals!` to select an integration strategy. This rule can depend on the test/trial reference spaces and interacting elements. The indices `p` and `q` refer to the position of the interacting elements in the enumeration defined by `geometry(basis)` and allow for fast retrieval of any element specific data stored in the quadrature data object `qd`.
- [`momintegrals!(operator,test_refspace,trial_refspace,test_element,trial_element,zlocal,qr)`](@ref): this function computes the local interaction matrix between the set of local test and trial shape functions and a specific pair of elements. The target matrix `zlocal` is provided as an argument to minimise memory allocations over subsequent calls. `qr` is an object returned by `quadrule` and contains all static and dynamic data defining the integration strategy used.

In the context of fast methods such as the Fast Multipole Method other algorithms on Discrete Operators will typically be defined to compute matrix vector products. These algorithms do not explicitly compute and store the interaction matrix (this would lead to unacceptable computational and memory complexity).

```@docs
elements
```

```@docs
numfunctions
coordtype
scalartype
assemblydata
geometry
refspace
```

```@docs
quaddata
quadrule
momintegrals!
```
