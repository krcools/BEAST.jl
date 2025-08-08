# Defining Custom Excitation

BEAST.jl aims to use your own fields as excitations as easy as possible. Specifying your own field can be done by defining a function of the Cartesian coordinates and supplementing this definition with a method of `BEAST.scalartype` to announce the scalar type of the value returned by the field evaluation. Scalar type is, as the name suggest, always a scalar type such as `ComplexF64`, even if the field is vector valued.

Upon supplying this data, the field function happily interacts with the rest of the BEAST.jl infrastructure. For example, taking the trace of the field is in no way different from doing so for the supplied exciations:

```julia
    using CompScienceMeshes

    κ = 1.0

    # Custom excitations can be defined as functions of a point in Cartesian space,
    # complemented with a function specifying the return type.
    E(x) = point(ComplexF64, 1.0, 1.0im, 0.0) * exp(-im*κ*x[3])
    BEAST.scalartype(::typeof(E)) = ComplexF64

    # Such a custom field plays nicely with the tangential trace
    # operators defined as part of the BEAST framework
    e = (n × E) × n

    fn = joinpath(pathof(BEAST), "../../examples/assets/sphere45.in")
    Γ = readmesh(fn)
    RT = raviartthomas(Γ)

    t = Maxwell3D.singlelayer(wavenumber=κ)

    @hilbertspace j
    @hilbertspace k

    Txx = assemble(t[k,j], j∈RT, k∈RT)
    ex = assemble(e[k], k∈RT)

    iTXX = BEAST.GMRESSolver(Txx; maxiter=1000, reltol=1e-5)
    uX = iTXX * ex
```

