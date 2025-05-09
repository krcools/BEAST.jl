
# BEAST.jl

This Julia package, the *boundary element analysis and simulation toolkit (BEAST)*, provides routines to convert integral and differential equations to linear systems of equations
via the boundary element method (BEM) and the finite element method (FEM). 
To this end, the (Petrov-) **Galerkin method** is employed.

Currently, the focus is on equations encountered in **classical electromagnetism**, where frequency and time domain equations are covered.
Several [operators](@ref operator), [basis functions](@ref basesRef), and geometry representations are implemented.

!!! note
    SI units and for time-harmonic simulations a time convention of ``\mathrm{e}^{\,\mathrm{j}\omega t}`` are used everywhere.

!!! tip
    To use the code have a look at the [general usage](@ref usageRef).

    However, the code is designed such that users can easily hook into the code at any level and implement new features.
    To do so, have a look at the [internals documentation](@ref InternalsRef) and the [contribution guidelines](@ref contribRef).
    Design goals are extendability and a performant execution. 


---
## Installation

Installing BEAST is done by entering the package manager (enter `]` at the julia REPL) and issuing:

```
pkg> add BEAST 
```


---
## Overview

The following [operators](@ref operator), [basis functions](@ref basesRef), and geometry representations are implemented.
To see details and all variations, have a look at the corresponding sections of this documentation.

### Operators

- **Boundary Integral operators**
    + Maxwell (3D)
        - Single Layer (time-harmonic & time-domain)
        - Double Layer (time-harmonic & time-domain)
    + Helmholtz (2D & 3D)
        - Single Layer
        - Double Layer

- **Volume Integral operators**
    + Maxwell
        - Single Layer (time-harmonic)
        - Double Layer (time-harmonic)

- **Local Operators**
    + Identity (+ variations thereof)

```@raw html
<br/>
```

### Basis functions

- **Spatial**
    + Low Order 
        - Raviart-Thomas / Rao-Wilton-Glisson
        - Buffa-Christiansen
        - Brezzi-Douglas-Marini
    + High Order
        - Graglia-Wilton-Peterson (GWP)
        - B-Spline based


- **Temporal**
    + Lagrange?
    + ...

```@raw html
<br/>
```

### Geometry Representations

- **Low Order** ("flat")
    + Triangular
    + Quadrilaterals
    + Tetrahedra


- **High Order** (curvilinear)
    + NURBS-surfaces

