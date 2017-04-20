export HH3DHyperSingularFDBIO
#export HH3DDoubleLayerFDBIO

"""
```
∫_Γ dx ∫_Γ dy \left(α G g(x) n_x ⋅ n_y f(y) + β G \mbox{curl} g(x) ⋅ \mbox{curl} f(y) \right)
```

with ``G(x,y) = \frac{e^{-γ |x-y|}}{4 π |x-y|}``
"""
type HH3DHyperSingularFDBIO{T,K} <: MaxwellOperator3D
    "coefficient of the weakly singular term"
    alpha::T
    "coefficient of the hyper singular term"
    beta::T
    "`im*κ` with `κ` the wave number"
    gamma::K
end

scalartype(op::HH3DHyperSingularFDBIO) = promote_type(typeof(op.alpha), typeof(op.beta), typeof(op.gamma))

HH3DHyperSingularFDBIO(gamma) = HH3DHyperSingularFDBIO(gamma^2, one(gamma), gamma)
HH3DSingleLayerFDBIO(gamma) = HH3DHyperSingularFDBIO(one(gamma), zero(gamma), gamma)


function quaddata(op::HH3DHyperSingularFDBIO,
        test_refspace::LagrangeRefSpace, trial_refspace::LagrangeRefSpace,
        test_elements, trial_elements)

        test_eval(x)  = test_refspace(x,  Val{:withcurl})
        trial_eval(x) = trial_refspace(x, Val{:withcurl})

        qp = quadpoints(test_eval,  test_elements,  (3,)), quadpoints(trial_eval, trial_elements, (4,))
end


function quadrule(op::HH3DHyperSingularFDBIO,
        test_refspace, trial_refspace,
        i, test_element, j, trial_element,
        quadrature_data)

        DoubleQuadStrategy(
            quadrature_data[1][1,i],
            quadrature_data[2][1,j]
        )
end


function integrand(op::HH3DHyperSingularFDBIO,
        kernel, test_values, test_element, trial_values, trial_element)

    α = op.alpha
    β = op.beta

    G = kernel.green

    if length(test_values) != 2
        @show test_values
    end
    g, curlg = test_values
    f, curlf = trial_values

    nx = normal(test_element)
    ny = normal(trial_element)

    α*dot(nx,ny)*g*f*G + β*dot(curlg,curlf)*G
end
