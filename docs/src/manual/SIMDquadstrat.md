# SIMDDoubleNum Quadstrat
This quadrature strategy uses the SIMD properties of the hardware to accelerate the double-num interactions. It should be used in combination with a quadrature strategy for the near interacting triangels, for example with the sacuterschwab quadrature.
The SIMD quadtarture is implemented for the Maxwell3D.singlelayer operator assembled in a raviart-thomas basis on triangles and for the composed operator structures at any basis embedded in a 3D space.
The Maxwell3D.doublelayer operator is assembled by casting the structure to the correpsonding compsoed operator.
Example
```julia
    op = Maxwell3D.singlelayer(wavenumber=1.0)
    assemble(op,X,X,quadstrat = BEAST.SIMDDoubleNumSauterQstrat(6,6,6,6,6,6))
    op = Maxwell3D.doublelayer(wavenumber=1.0)
    assemble(op,X,X,quadstrat = BEAST.SIMDDoubleNumSauterQstrat(6,6,6,6,6,6))
    op = Maxwell3D.doublelayer(wavenumber=1.0)
    assemble(op,X,X,quadstrat = BEAST.SIMDDoubleNumSauterQstrat(6,6,6,6,6,6))
    op = BEAST.CompDoubleInt(x->x, BEAST.Dot(),BEAST.HH3DGreen(1.0*im),BEAST.STimesV(), x->x)
    assemble(op,X,X,quadstrat = BEAST.SIMDDoubleNumSauterQstrat(6,6,6,6,6,6))
```

Note: The operators used in the operator fields of the CompDoubleInt should be of the specific types shown below
vector scalar multiplication: BEAST.VTimesS
scalar vector multiplication: BEAST.STimesV
scalar scalar multiplication: BEAST.STimesS
dot product: BEAST.Dot
cross product: BEAST.Cross
