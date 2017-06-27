# BEAST

Boundary Element Analysis and Simulation Toolkit

[![Build Status](https://travis-ci.org/krcools/BEAST.jl.svg?branch=master)](https://travis-ci.org/krcools/BEAST.jl)
[![Coverage Status](https://coveralls.io/repos/krcools/BEAST.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/krcools/BEAST.jl?branch=master)
[![codecov.io](http://codecov.io/github/krcools/BEAST.jl/coverage.svg?branch=master)](http://codecov.io/github/krcools/BEAST.jl?branch=master)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://krcools.github.io/BEAST.jl/latest/)

## Introduction

This package contains common basis functions and assembly routines for the implementation of
boundary element methods. Examples are included for the 2D and 3D Helmholtz equations and for
the 3D Maxwell equations.

Support for the space-time Galerkin based solution of time domain integral equations is in
place for the 3D Helmholtz and Maxwell equations.

## Installation

Installation is done simply by cloning this repo from within Julia (v0.6):

* `Pkg.clone("https://github.com/krcools/BEAST.jl")`

Prequisite packages will be pulled in automatically. In addition, some functionality requires `gmsh` to
be installed and on the system path.

## Hello World

To solve scattering by a time harmonic electromagnetic plane wave by a perfectly conducting
sphere:

```julia
using BEAST, CompScienceMeshes
o, x, y, z = euclidianbasis(3)

Γ = readmesh(Pkg.dir("BEAST","examples","sphere2.in"))
RT = raviartthomas(Γ)

κ = 1.0
t = MWSingleLayer3D(im*κ)
E = planewavemw3d(direction=z, polarization=x, wavenumber=κ)
e = (n × E) × n

@hilbertspace j; @hilbertspace k
efie = @discretise t[k,j]==e[k]  j∈RT k∈RT

u = solve(efie)
```
