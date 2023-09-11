
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

function sign_upon_permutation(op::HH3DSingleLayerFDBIO, I, J)
    return 1
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

function sign_upon_permutation(op::HH3DDoubleLayerFDBIO, I, J)
    return Combinatorics.levicivita(J)
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

function sign_upon_permutation(op::HH3DDoubleLayerTransposedFDBIO, I, J)
    return Combinatorics.levicivita(I)
end

defaultquadstrat(::Helmholtz3DOp, ::LagrangeRefSpace, ::LagrangeRefSpace) = DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)

function quaddata(op::Helmholtz3DOp, test_refspace::LagrangeRefSpace,
        trial_refspace::LagrangeRefSpace, test_elements, trial_elements,
        qs::DoubleNumWiltonSauterQStrat)

    test_eval(x)  = test_refspace(x,  Val{:withcurl})
    trial_eval(x) = trial_refspace(x, Val{:withcurl})
    
    # The combinations of rules (6,7) and (5,7 are) BAAAADDDD
    # they result in many near singularity evaluations with any
    # resemblence of accuracy going down the drain! Simply don't!
    # (same for (5,7) btw...).
    # test_qp = quadpoints(test_eval,  test_elements,  (6,))
    # bssi_qp = quadpoints(trial_eval, trial_elements, (7,))

    test_qp = quadpoints(test_eval,  test_elements,  (qs.outer_rule_far,))
    bsis_qp = quadpoints(trial_eval, trial_elements, (qs.inner_rule_far,))

    gausslegendre = (
      _legendre(qs.sauter_schwab_common_vert,0,1),
      _legendre(qs.sauter_schwab_common_edge,0,1),
      _legendre(qs.sauter_schwab_common_face,0,1),)

    return (;test_qp, bsis_qp, gausslegendre)
end

function quaddata(op::Helmholtz3DOp, test_refspace::LagrangeRefSpace,
    trial_refspace::LagrangeRefSpace, test_elements, trial_elements,
    qs::DoubleNumQStrat)

    test_eval(x)  = test_refspace(x,  Val{:withcurl})
    trial_eval(x) = trial_refspace(x, Val{:withcurl})

    # The combinations of rules (6,7) and (5,7 are) BAAAADDDD
    # they result in many near singularity evaluations with any
    # resemblence of accuracy going down the drain! Simply don't!
    # (same for (5,7) btw...).
    # test_qp = quadpoints(test_eval,  test_elements,  (6,))
    # bssi_qp = quadpoints(trial_eval, trial_elements, (7,))

    test_qp = quadpoints(test_eval,  test_elements,  (qs.outer_rule,))
    bsis_qp = quadpoints(trial_eval, trial_elements, (qs.inner_rule,))

    # gausslegendre = (
    # _legendre(qs.sauter_schwab_common_vert,0,1),
    # _legendre(qs.sauter_schwab_common_edge,0,1),
    # _legendre(qs.sauter_schwab_common_face,0,1),)

    return (;test_qp, bsis_qp)
end

defaultquadstrat(::Helmholtz3DOp, ::subReferenceSpace, ::subReferenceSpace) = DoubleNumWiltonSauterQStrat(4,7,4,7,4,4,4,4)

function quaddata(op::Helmholtz3DOp, test_refspace::subReferenceSpace,
        trial_refspace::subReferenceSpace, test_elements, trial_elements,
        qs::DoubleNumWiltonSauterQStrat)

    test_qp = quadpoints(test_refspace,  test_elements,  (qs.outer_rule_far,))
    bsis_qp = quadpoints(trial_refspace, trial_elements, (qs.inner_rule_far,))

    return test_qp, bsis_qp
end



function quadrule(op::Helmholtz3DOp,
        test_refspace::LagrangeRefSpace,
        trial_refspace::LagrangeRefSpace,
        i, test_element, j, trial_element, qd,
        qs::DoubleNumWiltonSauterQStrat)


    tol, hits = sqrt(eps(eltype(eltype(test_element.vertices)))), 0
    dmin2 = floatmax(eltype(eltype(test_element.vertices)))

    for t in test_element.vertices
        for s in trial_element.vertices
            d2 = LinearAlgebra.norm_sqr(t-s)
            dmin2 = min(dmin2, d2)
            hits += (d2 < tol)
    end end

    hits == 3 && return SauterSchwabQuadrature.CommonFace(qd.gausslegendre[3])
    hits == 2 && return SauterSchwabQuadrature.CommonEdge(qd.gausslegendre[2])
    hits == 1 && return SauterSchwabQuadrature.CommonVertex(qd.gausslegendre[1])

    test_quadpoints  = qd.test_qp
    trial_quadpoints = qd.bsis_qp
    h2 = volume(trial_element)
    xtol2 = 0.2 * 0.2
    k2 = abs2(gamma(op))
    max(dmin2*k2, dmin2/16h2) < xtol2 && return WiltonSERule(
        test_quadpoints[1,i],
        DoubleQuadRule(
            test_quadpoints[1,i],
            trial_quadpoints[1,j]))

    return DoubleQuadRule(
        qd.test_qp[1,i],
        qd.bsis_qp[1,j])
end


function quadrule(op::HH3DSingleLayerFDBIO,
        test_refspace::LagrangeRefSpace{T,0} where T,
        trial_refspace::LagrangeRefSpace{T,0} where T,
        i, test_element, j, trial_element, qd,
        qs::DoubleNumQStrat)

    return DoubleQuadRule(
        qd.test_qp[1,i],
        qd.bsis_qp[1,j])
end

regularpart(op::HH3DHyperSingularFDBIO) = HH3DHyperSingularReg(op.alpha, op.beta, op.gamma)
singularpart(op::HH3DHyperSingularFDBIO) = HH3DHyperSingularSng(op.alpha, op.beta, op.gamma)

#= function quadrule(op::HH3DHyperSingularFDBIO,
    test_refspace::LagrangeRefSpace,
    trial_refspace::LagrangeRefSpace,
    i, test_element, j, trial_element, qd,
    qs::DoubleNumWiltonSauterQStrat)

    tol2, hits = eps(eltype(eltype(test_element.vertices))), 0
    for t in test_element.vertices
        for s in trial_element.vertices
            hits += (LinearAlgebra.norm_sqr(t-s) < tol2)
    end end

    hits == 3 && return SauterSchwabQuadrature.CommonFace(qd.gausslegendre[3])
    hits == 2 && return SauterSchwabQuadrature.CommonEdge(qd.gausslegendre[2])
    hits == 1 && return SauterSchwabQuadrature.CommonVertex(qd.gausslegendre[1])

    test_quadpoints  = qd.test_qp
    trial_quadpoints = qd.bsis_qp
    h2 = volume(trial_element)
    xtol2 = 0.2 * 0.2
    k2 = abs2(op.gamma)
    max(dmin2*k2, dmin2/16h2) < xtol2 && return WiltonSERule(
        test_quadpoints[1,i],
        DoubleQuadRule(
            test_quadpoints[1,i],
            trial_quadpoints[1,j]))

    return DoubleQuadRule(
        qd.test_qp[1,i],
        qd.bsis_qp[1,j])
end =#


function quadrule(op::HH3DHyperSingularFDBIO,
    test_refspace::LagrangeRefSpace{T,1} where T,
    trial_refspace::LagrangeRefSpace{T,1} where T,
    i, test_element, j, trial_element, qd,
    qs::DoubleNumQStrat)

    tol, hits = sqrt(eps(eltype(eltype(test_element.vertices)))), 0
    for t in test_element.vertices
        for s in trial_element.vertices
            norm(t-s) < tol && (hits +=1)
    end end

    # hits == 3 && return SauterSchwabQuadrature.CommonFace(qd.gausslegendre[3])
    # hits == 2 && return SauterSchwabQuadrature.CommonEdge(qd.gausslegendre[2])
    # hits == 1 && return SauterSchwabQuadrature.CommonVertex(qd.gausslegendre[1])

    test_quadpoints  = qd.test_qp
    trial_quadpoints = qd.bsis_qp

    # hits != 0 && return WiltonSERule(
    #     test_quadpoints[1,i],
    #     DoubleQuadRule(
    #         test_quadpoints[1,i],
    #         trial_quadpoints[1,j]))

    return DoubleQuadRule(
        qd.test_qp[1,i],
        qd.bsis_qp[1,j])
end


function quadrule(op::HH3DSingleLayerFDBIO, test_refspace::subReferenceSpace,
    trial_refspace::subReferenceSpace, i, test_element, j, trial_element, quadrature_data,
    qs::DoubleNumWiltonSauterQStrat)

    # tol, hits = 1e-10, 0
    # for t in getelementVertices(test_element)
    #     for s in getelementVertices(trial_element)
    #         norm(t-s) < tol && (hits +=1; break)
    # end end
    #
    # test_quadpoints  = quadrature_data[1]
    # trial_quadpoints = quadrature_data[2]

    # hits != 0 && return WiltonSERule(
    #     test_quadpoints[1,i],
    #     DoubleQuadRule(
    #         test_quadpoints[1,i],
    #         trial_quadpoints[1,j]))

    return DoubleQuadRule(
        quadrature_data[1][1,i],
        quadrature_data[2][1,j])
    # return SauterSchwabStrategy(hits)
     # test_qd = qd[1]
     # trial_qd = qd[2]
     #
     # DoubleQuadRule(
     #    test_qd[1,i], # rule 1 on test element j
     #    trial_qd[1,j] # rule 1 on trial element i
     # )
end


function quadrule(op::Helmholtz3DOp,
    test_refspace::RefSpace, trial_refspace::RefSpace,
    i, test_element, j, trial_element, quadrature_data,
    qs::DoubleNumWiltonSauterQStrat)

    return DoubleQuadRule(
        quadrature_data[1][1,i],
        quadrature_data[2][1,j])
end

#= function quadrule(op::HH3DDoubleLayerTransposedFDBIO,
    test_refspace::LagrangeRefSpace{T,1} where T,
    trial_refspace::LagrangeRefSpace{T,0} where T,
    i, test_element, j, trial_element, quadrature_data,
    qs::DoubleNumWiltonSauterQStrat)

    tol, hits = sqrt(eps(eltype(eltype(test_element.vertices)))), 0
    dmin2 = floatmax(eltype(eltype(test_element.vertices)))

    for t in test_element.vertices
        for s in trial_element.vertices
            d2 = LinearAlgebra.norm_sqr(t-s)
            dmin2 = min(dmin2, d2)
            hits += (d2 < tol)
    end end

    hits == 3 && return SauterSchwabQuadrature.CommonFace(quadrature_data.gausslegendre[3])
    hits == 2 && return SauterSchwabQuadrature.CommonEdge(quadrature_data.gausslegendre[2])
    hits == 1 && return SauterSchwabQuadrature.CommonVertex(quadrature_data.gausslegendre[1])

    test_quadpoints  = quadrature_data[1]
    trial_quadpoints = quadrature_data[2]
    test_quadpoints  = quadrature_data.test_qp
    trial_quadpoints = quadrature_data.bsis_qp
    h2 = volume(trial_element)
    xtol2 = 0.2 * 0.2
    k2 = abs2(gamma(op))

    max(dmin2*k2, dmin2/16h2) < xtol2 && return WiltonSERule(
        test_quadpoints[1,i],
        DoubleQuadRule(
            test_quadpoints[1,i],
            trial_quadpoints[1,j]))

    return DoubleQuadRule(
        quadrature_data[1][1,i],
        quadrature_data[2][1,j])
end =#

#= function quadrule(op::HH3DDoubleLayerFDBIO,
    test_refspace::LagrangeRefSpace{T,0} where T,
    trial_refspace::LagrangeRefSpace{T,1} where T,
    i, test_element, j, trial_element, quadrature_data,
    qs::DoubleNumWiltonSauterQStrat)

    tol, hits = sqrt(eps(eltype(eltype(test_element.vertices)))), 0
    dmin2 = floatmax(eltype(eltype(test_element.vertices)))

    for t in test_element.vertices
        for s in trial_element.vertices
            d2 = LinearAlgebra.norm_sqr(t-s)
            dmin2 = min(dmin2, d2)
            hits += (d2 < tol)
    end end

    hits == 3 && return SauterSchwabQuadrature.CommonFace(quadrature_data.gausslegendre[3])
    hits == 2 && return SauterSchwabQuadrature.CommonEdge(quadrature_data.gausslegendre[2])
    hits == 1 && return SauterSchwabQuadrature.CommonVertex(quadrature_data.gausslegendre[1])
    test_quadpoints  = quadrature_data[1]
    trial_quadpoints = quadrature_data[2]

    test_quadpoints  = quadrature_data.test_qp
    trial_quadpoints = quadrature_data.bsis_qp
    h2 = volume(trial_element)
    xtol2 = 0.2 * 0.2
    k2 = abs2(gamma(op))

#=     max(dmin2*k2, dmin2/16h2) < xtol2 && return WiltonSERule(
        test_quadpoints[1,i],
        DoubleQuadRule(
            test_quadpoints[1,i],
            trial_quadpoints[1,j])) =#
    return DoubleQuadRule(
        quadrature_data[1][1,i],
        quadrature_data[2][1,j])
end =#


function quadrule(op::Helmholtz3DOp,
    test_refspace::RefSpace, trial_refspace::RefSpace,
    i, test_element, j, trial_element, quadrature_data,
    qs::DoubleNumQStrat)

    return DoubleQuadRule(
        quadrature_data[1][1,i],
        quadrature_data[2][1,j])
end

function quadrule(op::Helmholtz3DOp,
    test_refspace::RTRefSpace, trial_refspace::RTRefSpace,
    i, test_element, j, trial_element, quadrature_data,
    qs::DoubleNumWiltonSauterQStrat)

    test_quadpoints  = quadrature_data[1]
    trial_quadpoints = quadrature_data[2]

    return DoubleQuadRule(
        quadrature_data[1][1,i],
        quadrature_data[2][1,j])
end



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





HH3DSingleLayerFDBIO(gamma) = HH3DSingleLayerFDBIO(one(gamma), gamma)

regularpart(op::HH3DSingleLayerFDBIO) = HH3DSingleLayerReg(op.alpha, gamma(op))
singularpart(op::HH3DSingleLayerFDBIO) = HH3DSingleLayerSng(op.alpha, gamma(op))

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


function integrand(op::Union{HH3DSingleLayerFDBIO,HH3DSingleLayerReg},
    kernel, test_values, test_element, trial_values, trial_element)

α = op.alpha
G = kernel.green

g = test_values.value
f = trial_values.value

α*dot(g, G*f)
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


function integrand(biop::HH3DDoubleLayerFDBIO,
        kernel, fp, mp, fq, mq)

    nq = normal(mq)
    fp[1] * dot(nq, -kernel.gradgreen) * fq[1]
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

function integrand(biop::HH3DDoubleLayerTransposedFDBIO,
        kernel, fp, mp, fq, mq)

    np = normal(mp)
    fp[1] * dot(np, kernel.gradgreen) * fq[1]
end
