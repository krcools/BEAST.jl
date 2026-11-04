# Unreleased

- Interpolators onto finite element spaces. `LagrangeInterpolator` interpolates scalar fields onto continuous Lagrange spaces of any order. `GWPInterpolator` interpolates tangential fields onto GWP(curl) spaces and, when applied to a space instead of a field, returns the sparse matrix that interpolates each function of that space onto the GWP space.
- `InterpolatorOp` assembles the interpolation matrix between two spaces that carry an interpolator. Such spaces are wrapped in `SpaceWithDoFs`; `gwpcurl(mesh; order, interpolator=true)` constructs one directly. For the lowest order this reproduces `CompScienceMeshes.embedding`, but it applies to any order.
- `examples/pmchwt_single_trace.jl` uses `InterpolatorOp` and GWP spaces, and now works for any order.
- **Breaking:** the higher-order quasi-Helmholtz projectors are now actual projectors: the inverse Gram matrix is included, so it no longer needs to be applied by the caller. `PΛ` and `ℙΣ` have been renamed to `PΛH` and `ℙΣH`.
- `verbose` keyword for `assemble` to suppress progress output, e.g. when several matrices are assembled concurrently.
- Example `examples/efie_bunny.jl`: EFIE on a cubic Stanford bunny with face currents.

## Fixed bugs:
- Quasi-local assembly threw a `BoundsError` when called from the interactive thread pool (the default since Julia 1.12, e.g. `-t 8`).

# New in 2.11.0

- Support for higher order geometries. For an example, see `examples/efie_geop2.jl` and the corresponding page in the documentation. `facecurrents` can generate values on curve aware refinements.
- `momintegrals!` has been renamed to `integrate!`. The family of `quadrule` methods used by frequency-domain `IntegralOperator` assembly (the ones dispatching on an `AbstractQuadStrat` for a pair of test/trial elements) has also been folded into `integrate!`: each such method now builds its quadrature rule and, from within the same method/branch, immediately calls `integrate!` on it, instead of returning it to a separate call site. This avoids a dynamic dispatch on the possible quadrature-rule types in the assembly hot loop. Whether the rule is applied or just handed back is controlled by a new `action` keyword: `action=ApplyIntegrate()` (the default) evaluates the rule into the output buffer, `action=ReturnQRule()` recovers the old behaviour of returning the rule object unevaluated, and `action=ApplyIntegrateNonConforming()` is used internally by the non-conforming-mesh rules (`NonConformingOverlapQRule`, `NonConformingTouchQRule`), whose own `integrate!` methods only have local refspaces in scope rather than a full test/trial `Space`. `quadrule` methods outside this family (e.g. `LocalOperator`, excitations, farfield/nearfield postprocessing) are unaffected and keep their name. `quadinfo` now returns the built quadrature rule instead of `nothing`.
- Telles near-singularity quadrature for 2D Helmholtz operators (`DoubleNumSauterTellesQstrat`) and for near-field post-processing.
- `MWStaticExtractedHyperSingular`: the Maxwell hypersingular operator with its static kernel extracted (Müller), which can be tested with functions that are not divergence conforming.
- Singularity extraction (Wilton) for the Helmholtz 3D Nitsche operator with `lagrangec0d1` trial functions.
- Assembly of systems exploits symmetry in addition to repetition: blocks of operators for which `issymmetric` returns `true` (currently `MWSingleLayer3D` and `MWDoubleLayer3D`) are reused in transposed form.
- Quadrature and scratch buffers are reused across integration calls, reducing allocations during assembly.
- Assembly of Lagrange and Nédélec elements on planar triangle meshes (`udim = 2`).
- `gradient` for `DirectProductSpace`s.
- Local operators and excitations can be assembled for bases with no functions.
- Interpolation (and hence `restrict`) of Lagrange elements on segments and tetrahedra.
- Local assembly uses a visitor-based octree search.

## Fixed bugs:
- Dual linear Lagrange elements on surfaces (`duallagrangec0d1`) now sum to one on open surfaces too. Their counterparts on curves have been rewritten and are now continuous across coarse vertices.
- `duallagrangecxd0` on curves.
- `unitfunctioncxd0` for meshes with two-dimensional vertices.
- The broken three-argument `restrict` for quadratic Lagrange elements was removed; the generic fallback through interpolation is used instead.
- Integral operators with a test mesh that refines the trial mesh (`TestRefinesTrialQRule`) now work on curves, not only on triangles.
- Local operators between a refined test space and a trial space supported on part of the coarse mesh (e.g. dual Lagrange test functions) no longer throw a `BoundsError`.

## Dependencies:
- CompScienceMeshes 1, CollisionDetection 0.2, PlotlyBase 0.10, SauterSchwabQuadrature 2.5.

# New in 2.10.0

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
