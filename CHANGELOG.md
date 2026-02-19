# Changelog

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
