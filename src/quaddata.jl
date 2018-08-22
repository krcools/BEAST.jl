

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


"""
    quadrule(operator, test_refspace, trial_refspace, test_index, test_chart, trial_index, trial_chart, quad_data)

Based on the operator kernel and the test and trial elements, this function builds
an object whose type and data fields specify the quadrature rule that needs to be
used to accurately compute the interaction integrals. The `quad_data` object created
by `quaddata` is passed to allow reuse of any precomputed data such as quadrature
points and weights, geometric quantities, etc.

The type of the returned quadrature rule will help in deciding which method of
`momintegrals` to dispatch to.
"""
function quadrule(op::IntegralOperator, g::RefSpace, f::RefSpace, i, τ, j, σ, qd)
    # defines coincidence of points
    dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))

    # decides on whether to use singularity extraction
    xtol = 0.2
    k = norm(op.gamma)

    hits = 0
    xmin = xtol
    for t in τ.vertices
      for s in σ.vertices
        d = norm(t-s)
        xmin = min(xmin, k*d)
        if d < dtol
          hits +=1
          break
        end
      end
    end

  hits == 3   && return SauterSchwabQuadrature.CommonFace(qd.gausslegendre[3])
  hits == 2   && return SauterSchwabQuadrature.CommonEdge(qd.gausslegendre[2])
  hits == 1   && return SauterSchwabQuadrature.CommonVertex(qd.gausslegendre[1])
  xmin < xtol && return DoubleQuadStrategy(
      qd.tpoints[2,i],
      qd.bpoints[2,j],)
  return DoubleQuadStrategy(
    qd.tpoints[1,i],
    qd.bpoints[1,j],)
end
