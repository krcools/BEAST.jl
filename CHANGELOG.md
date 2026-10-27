- Support for higher order geometries. For an example, see `examples/efie_geop2.jl` and the corresponding page in the documentation.

# New in 2.10.0

- `momintegrals!` has been renamed to `integrate!`. The family of `quadrule` methods used by frequency-domain `IntegralOperator` assembly (the ones dispatching on an `AbstractQuadStrat` for a pair of test/trial elements) has also been folded into `integrate!`: each such method now builds its quadrature rule and, from within the same method/branch, immediately calls `integrate!` on it, instead of returning it to a separate call site. This avoids a dynamic dispatch on the possible quadrature-rule types in the assembly hot loop. Whether the rule is applied or just handed back is controlled by a new `action` keyword: `action=ApplyIntegrate()` (the default) evaluates the rule into the output buffer, `action=ReturnQRule()` recovers the old behaviour of returning the rule object unevaluated, and `action=ApplyIntegrateNonConforming()` is used internally by the non-conforming-mesh rules (`NonConformingOverlapQRule`, `NonConformingTouchQRule`), whose own `integrate!` methods only have local refspaces in scope rather than a full test/trial `Space`. `quadrule` methods outside this family (e.g. `LocalOperator`, excitations, farfield/nearfield postprocessing) are unaffected and keep their name. `quadinfo` now returns the built quadrature rule instead of `nothing`.
- `solve` and `solve!`  on `GMRESSolver` allow for the defaults preconditioners to be overwritten
- `materialize` keyword to assemble for bilinear forms allows for fine-grained control over the assembly of operators (e.g. dense representations vs H-matrix representations)
- Restart for GMRES defaults to the maximum number of iterations
- `assemble` of bilinear forms checks for recurring blocks and reuses them.
- Support for visualisation through `PlotlyJS` through an extension

# New in 2.9.0

- Cell coloring based lock-free multi-threaded assembly for frequency domain integral operators
- Higher order Lagrange elements (cx and c0) on segments
- All LinearMaps can be cast into bilinear forms. Assembly of the latter trivially returns the underlying LinearMap. This enables for example the construction inverses that explicitly exploit block diagonal structure.
- Use CompScienceMeshes 0.11 providing typed cells.

# New in 2.8.0
- Improved support of for the Helmholtz equation.
- Higher-order quasi-Helmholtz projectors.
- Bourhis and Nedelec Theta operators.
- Curl of Lagrange basis now based on GWP basis to support curl for any order of basis.
- Direct solvers (Cholesky,LU) wrapped as LinearMap.
- Better integration of GMRES from Krylov.jl

# New in 2.7.3
- Generic way of introducing potentials added.
- Generic Trace operators of those potentials introduced.
- Composed user defined operators can be constructed. More information can be found in the documentation.

## Fixed bugs:
- Minus sign in the normalvector cross RT and ND space.
- Construction of the duallagrangec0d1 space.
- Quadstrat in the LinearCombinationOfOperators.

# New in 2.7.2
- GMRES from Krylov.jl wrapped as a LinearMap
- Function `union` for a vector of spaces of the same type and defined on the same mesh returns a space containing the union of the functions in the input spaces.
- Function `extend` computes the extension by zero of functions defined on a submesh. Currently works only when the corresponding cellls in the submesh and supermesh are ordered identically (oriented identically is not sufficient).

## New in 2.7.0

- PR #157: Improved support for 2D Helmholtz problems
