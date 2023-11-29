abstract type Helmholtz3DOp{T,K} <: MaxwellOperator3D{T,K} end
abstract type Helmholtz3DOpReg{T,K} <: MaxwellOperator3DReg{T,K} end
"""
```
∫_Γ dx ∫_Γ dy \\left(α G g(x) n_x ⋅ n_y f(y) + β G \\mbox{curl} g(x) ⋅ \\mbox{curl} f(y) \\right)
```
with ``G(x,y) = \\frac{e^{-γ |x-y|}}{4 π |x-y|}``
"""
struct HH3DHyperSingularFDBIO{T,K} <: Helmholtz3DOp{T,K}
    "coefficient of the weakly singular term"
    alpha::T
    "coefficient of the hyper singular term"
    beta::T
    "`im*κ` with `κ` the wave number"
    gamma::K
end

function sign_upon_permutation(op::HH3DHyperSingularFDBIO, I, J)
    return Combinatorics.levicivita(I) * Combinatorics.levicivita(J)
end

struct HH3DHyperSingularReg{T,K} <: Helmholtz3DOpReg{T,K}
    "coefficient of the weakly singular term"
    alpha::T
    "coefficient of the hyper singular term"
    beta::T
    "`im*κ` with `κ` the wave number"
    gamma::K
end

struct HH3DHyperSingularSng{T,K} <: Helmholtz3DOp{T,K}
    "coefficient of the weakly singular term"
    alpha::T
    "coefficient of the hyper singular term"
    beta::T
    "`im*κ` with `κ` the wave number"
    gamma::K
end
HH3DHyperSingularFDBIO(gamma) = HH3DHyperSingularFDBIO(gamma^2, one(gamma), gamma)

"""
```math
a(u,v) = α ∬_{Γ×Γ} u(x) G_{γ}(|x-y|) v(y)
```
with ``G_{γ}(r) = \\frac{e^{-γr}}{4πr}``.
"""
struct HH3DSingleLayerFDBIO{T,K} <: Helmholtz3DOp{T,K}
    alpha::T
    gamma::K
end

struct HH3DSingleLayerReg{T,K} <: Helmholtz3DOpReg{T,K}
    alpha::T
    gamma::K
end

struct HH3DSingleLayerSng{T,K} <: Helmholtz3DOp{T,K}
    alpha::T
    gamma::K
end

struct HH3DDoubleLayerFDBIO{T,K} <: Helmholtz3DOp{T,K}
    alpha::T
    gamma::K
end

struct HH3DDoubleLayerReg{T,K} <: Helmholtz3DOpReg{T,K}
    alpha::T
    gamma::K
end

struct HH3DDoubleLayerSng{T,K} <: Helmholtz3DOp{T,K}
    alpha::T
    gamma::K
end

struct HH3DDoubleLayerTransposedFDBIO{T,K} <: Helmholtz3DOp{T,K}
    alpha::T
    gamma::K
end

struct HH3DDoubleLayerTransposedReg{T,K} <: Helmholtz3DOpReg{T,K}
    alpha::T
    gamma::K
end

struct HH3DDoubleLayerTransposedSng{T,K} <: Helmholtz3DOp{T,K}
    alpha::T
    gamma::K
end

defaultquadstrat(::Helmholtz3DOp, ::LagrangeRefSpace, ::LagrangeRefSpace) =
    DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)

defaultquadstrat(::Helmholtz3DOp, ::subReferenceSpace, ::subReferenceSpace) = 
    DoubleNumWiltonSauterQStrat(4,7,4,7,4,4,4,4)

regularpart(op::HH3DHyperSingularFDBIO) = HH3DHyperSingularReg(op.alpha, op.beta, op.gamma)
singularpart(op::HH3DHyperSingularFDBIO) = HH3DHyperSingularSng(op.alpha, op.beta, op.gamma)

function (igd::Integrand{<:HH3DHyperSingularFDBIO})(x,y,f,g)
    α = igd.operator.alpha
    β = igd.operator.beta
    γ = gamma(igd.operator)

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1 / R
    green = exp(-γ*R)*(i4pi*iR)
    nx = normal(x)
    ny = normal(y)

    _integrands(f,g) do fi, gi
        α*dot(nx,ny)*gi.value*fi.value*green + β*dot(gi.curl,fi.curl)*green
    end
end


HH3DSingleLayerFDBIO(gamma) = HH3DSingleLayerFDBIO(one(gamma), gamma)

regularpart(op::HH3DSingleLayerFDBIO) = HH3DSingleLayerReg(op.alpha, op.gamma)
singularpart(op::HH3DSingleLayerFDBIO) = HH3DSingleLayerSng(op.alpha, op.gamma)

function (igd::Integrand{<:HH3DSingleLayerFDBIO})(x,y,f,g)
    α = igd.operator.alpha
    γ = gamma(igd.operator)

   r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1 / R
    green = exp(-γ*R)*(i4pi*iR)

    αG = α * green

    _integrands(f,g) do fi, gi
        dot(gi.value, αG*fi.value)
    end
end


HH3DDoubleLayerFDBIO(gamma) = HH3DDoubleLayerFDBIO(one(gamma), gamma)
regularpart(op::HH3DDoubleLayerFDBIO) = HH3DDoubleLayerReg(op.alpha, op.gamma)
singularpart(op::HH3DDoubleLayerFDBIO) = HH3DDoubleLayerSng(op.alpha, op.gamma)

function (igd::Integrand{<:HH3DDoubleLayerFDBIO})(x,y,f,g)
    γ = gamma(igd.operator)
    α = igd.operator.alpha

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-γ*R)*(iR*i4pi)
    gradgreen = -(γ + iR) * green * (iR * r)
    αgradgreen = α * gradgreen
    n = normal(y)
    fvalue = getvalue(f)
    gvalue = getvalue(g)

    return _krondot(fvalue,gvalue) * dot(n, -αgradgreen)
end


HH3DDoubleLayerTransposedFDBIO(gamma) = HH3DDoubleLayerTransposedFDBIO(one(gamma), gamma)
regularpart(op::HH3DDoubleLayerTransposedFDBIO) = HH3DDoubleLayerTransposedReg(op.alpha, op.gamma)
singularpart(op::HH3DDoubleLayerTransposedFDBIO) = HH3DDoubleLayerTransposedSng(op.alpha, op.gamma)

function (igd::Integrand{<:HH3DDoubleLayerTransposedFDBIO})(x,y,f,g)
    γ = gamma(igd.operator)
    α = igd.operator.alpha

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-γ*R)*(iR*i4pi)
    gradgreen = -(γ + iR) * green * (iR * r)
    αgradgreen = α * gradgreen
    n = normal(x)
    fvalue = getvalue(f)
    gvalue = getvalue(g)

    return _krondot(fvalue,gvalue) * dot(n, αgradgreen)
end
