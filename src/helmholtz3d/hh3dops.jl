export HH3DHyperSingularFDBIO

abstract type Helmholtz3DOp <: MaxwellOperator3D end


"""
```
∫_Γ dx ∫_Γ dy \left(α G g(x) n_x ⋅ n_y f(y) + β G \mbox{curl} g(x) ⋅ \mbox{curl} f(y) \right)
```

with ``G(x,y) = \frac{e^{-γ |x-y|}}{4 π |x-y|}``
"""
struct HH3DHyperSingularFDBIO{T,K} <: Helmholtz3DOp
    "coefficient of the weakly singular term"
    alpha::T
    "coefficient of the hyper singular term"
    beta::T
    "`im*κ` with `κ` the wave number"
    gamma::K
end

HH3DHyperSingularFDBIO(gamma) = HH3DHyperSingularFDBIO(gamma^2, one(gamma), gamma)
#HH3DSingleLayerFDBIO(gamma) = HH3DHyperSingularFDBIO(one(gamma), zero(gamma), gamma)

scalartype(op::HH3DHyperSingularFDBIO) = promote_type(typeof(op.alpha), typeof(op.beta), typeof(op.gamma))


"""
```math
a(u,v) = α ∬_{Γ×Γ} u(x) G_{γ}(|x-y|) v(y)
```

with ``G_{γ}(r) = \frac{e^{-γr}}{4πr}``.
"""
struct HH3DSingleLayer{T,K} <: Helmholtz3DOp
    alpha::T
    gamma::K
end

HH3DSingleLayerFDBIO(gamma) = HH3DSingleLayer(one(gamma), gamma)


function quaddata(op::Helmholtz3DOp, test_refspace::LagrangeRefSpace,
        trial_refspace::LagrangeRefSpace, test_elements, trial_elements)

    test_eval(x)  = test_refspace(x,  Val{:withcurl})
    trial_eval(x) = trial_refspace(x, Val{:withcurl})

    qp = quadpoints(test_eval,  test_elements,  (3,)), quadpoints(trial_eval, trial_elements, (4,))
end


function quadrule(op::Helmholtz3DOp, test_refspace, trial_refspace, i,
        test_element, j, trial_element, quadrature_data)

    DoubleQuadStrategy(
        quadrature_data[1][1,i],
        quadrature_data[2][1,j])
end


function integrand(op::HH3DHyperSingularFDBIO,
        kernel, test_values, test_element, trial_values, trial_element)

    α = op.alpha
    β = op.beta

    G = kernel.green

    g, curlg = test_values
    f, curlf = trial_values

    nx = normal(test_element)
    ny = normal(trial_element)

    α*dot(nx,ny)*g*f*G + β*dot(curlg,curlf)*G
end


function integrand(op::HH3DSingleLayer, kernel, test_values,
        test_element, trial_values, trial_element)

    α = op.alpha
    G = kernel.green
    g, curlg = test_values
    f, curlf = trial_values

    α*g*f*G
end


module Helmholtz3D
    using ..BEAST

    singlelayer(;
        gamma=error("propagation constant is a required argument"),
        alpha=one(gamma)) = HH3DSingleLayer(alpha,gamma)

    hypersingular(;
        gamma=error("propagation constant is a required argument"),
        alpha=gamma^2,
        beta=ones(gamma)) = HH3DHyperSingularFDBIO(alpha, beta, gamma)

end

export Helmholtz3D
