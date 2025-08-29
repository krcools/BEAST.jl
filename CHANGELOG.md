# Changelog

# New in 2.7.2

- GMRES from Krylov.jl wrapped as a LinearMap
- Function `union` for a vector of spaces of the same type and defined on the same mesh returns a space containing the union of the functions in the input spaces.
- Function `extend` computes the extension by zero of functions defined on a submesh. Currently works only when the corresponding cellls in the submesh and supermesh are ordered identically (oriented identically is not sufficient).

## New in 2.7.0

- PR #157: Improved support for 2D Helmholtz problems
