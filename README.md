# BEAST

Boundary Element Analysis and Simulation Toolkit

[![Build Status](https://travis-ci.org/krcools/BEAST.jl.svg?branch=master)](https://travis-ci.org/krcools/BEAST.jl)

[![Coverage Status](https://coveralls.io/repos/krcools/BEAST.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/krcools/BEAST.jl?branch=master)

[![codecov.io](http://codecov.io/github/krcools/BEAST.jl/coverage.svg?branch=master)](http://codecov.io/github/krcools/BEAST.jl?branch=master)

## Introduction

This package contains common basis functions and assembly routines for the implementation of
boundary element methods. Examples are included for the 2D and 3D Helmholtz equations and for
the 3D Maxwell equations.

Support for the space-time Galerkin based solution of time domain integral equations is in
place for the 3D Helmholtz and Maxwell equations.

## Prerequisites

In addition to the Julia packages in REQUIRE.jl, the following packages are required:

* `krcools/CompScienceMeshes`
* `krcools/LinearForms`

In addition, some functionality requires `gmsh` to be installed and on the system path. `LinearForms`
is only required for syntactic sugar; all assembly routines can be directly called.

## Example script

To solve scattering by a time harmonic electromagnetic plane wave by a perfectly conducting
sphere:

```julia
using BEAST
using CompScienceMeshes
using LinearForms

o, x, y, z = euclidianbasis(Float64, 3)
n = BEAST.n

Γ = meshsphere(1.0, 0.2)
RT = raviartthomas(Γ)

κ = 1.0
t = MWSingleLayer3D(im*κ)
E = PlaneWaveMW(z, x, κ)
e = (n × E) × n

j, = hilbertspace(:j)
k, = hilbertspace(:k)

EFIE = @varform t[k,j] == e[k]
efie = @discretise EFIE j∈RT k∈RT

u = solve(efie)
fcr, geo = facecurrents(u, RT)
```
