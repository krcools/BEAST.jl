

"""
    quaddata(operator, test_refspace, trial_refspace, test_elements, trial_elements)

Returns an object cashing data required for the computation of boundary element
interactions. It is up to the client programmer to decide what (if any) data is
cached. For double numberical quadrature, storing the integration points for
example can significantly speed up matrix assembly.

- `operator` is an integration kernel.
- `test_refspace` and `trial_refspace` are reference space objects. `quadata`
is typically overloaded on the type of these local spaces of shape functions.
(See the implementation in `maxwell.jl` for an example).
- `test_elements` and `trial_elements` are iterable collections of the geometric
elements on which the finite element space are defined. These are provided to
allow computation of the actual integrations points - as opposed to only their
coordinates.
"""
function quaddata end
