
abstract type Helmholtz3DOp <: MaxwellOperator3D end
abstract type Helmholtz3DOpReg <: MaxwellOperator3DReg end
"""
```
∫_Γ dx ∫_Γ dy \\left(α G g(x) n_x ⋅ n_y f(y) + β G \\mbox{curl} g(x) ⋅ \\mbox{curl} f(y) \\right)
```
with ``G(x,y) = \\frac{e^{-γ |x-y|}}{4 π |x-y|}``
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
scalartype(op::HH3DHyperSingularFDBIO) = promote_type(typeof(op.alpha), typeof(op.beta), typeof(op.gamma))

"""
```math
a(u,v) = α ∬_{Γ×Γ} u(x) G_{γ}(|x-y|) v(y)
```
with ``G_{γ}(r) = \\frac{e^{-γr}}{4πr}``.
"""
struct HH3DSingleLayerFDBIO{T,K} <: Helmholtz3DOp
    alpha::T
    gamma::K
end

struct HH3DSingleLayerReg{T,K} <: Helmholtz3DOpReg
    alpha::T
    gamma::K
end

struct HH3DSingleLayerSng{T,K} <: Helmholtz3DOp
    alpha::T
    gamma::K
end

struct HH3DDoubleLayer{T,K} <: Helmholtz3DOp
    alpha::T
    gamma::K
end

struct HH3DDoubleLayerTransposed{T,K} <: Helmholtz3DOp
    alpha::T
    gamma::K
end

function quaddata(op::Helmholtz3DOp, test_refspace::LagrangeRefSpace,
        trial_refspace::LagrangeRefSpace, test_elements, trial_elements)

    test_eval(x)  = test_refspace(x,  Val{:withcurl})
    trial_eval(x) = trial_refspace(x, Val{:withcurl})

    # The combinations of rules (6,7) and (5,7 are) BAAAADDDD
    # they result in many near singularity evaluations with any
    # resemblence of accuracy going down the drain! Simply don't!
    # (same for (5,7) btw...).
    # test_qp = quadpoints(test_eval,  test_elements,  (6,))
    # bssi_qp = quadpoints(trial_eval, trial_elements, (7,))

    test_qp = quadpoints(test_eval,  test_elements,  (4,))
    bsis_qp = quadpoints(trial_eval, trial_elements, (7,))

    return test_qp, bsis_qp
end

function quaddata(op::Helmholtz3DOp, test_refspace::subReferenceSpace,
        trial_refspace::subReferenceSpace, test_elements, trial_elements)

    test_qp = quadpoints(test_refspace,  test_elements,  (4,))
    bsis_qp = quadpoints(trial_refspace, trial_elements, (7,))

    return test_qp, bsis_qp
end



function quadrule(op::HH3DSingleLayerFDBIO,
        test_refspace::LagrangeRefSpace{T,0} where T,
        trial_refspace::LagrangeRefSpace{T,0} where T,
        i, test_element, j, trial_element, quadrature_data)

    tol, hits = sqrt(eps(eltype(eltype(test_element.vertices)))), 0
    for t in test_element.vertices
        for s in trial_element.vertices
            norm(t-s) < tol && (hits +=1; break)
    end end

    test_quadpoints  = quadrature_data[1]
    trial_quadpoints = quadrature_data[2]

    hits != 0 && return WiltonSEStrategy(
        test_quadpoints[1,i],
        DoubleQuadStrategy(
            test_quadpoints[1,i],
            trial_quadpoints[1,j]))

    return DoubleQuadStrategy(
        quadrature_data[1][1,i],
        quadrature_data[2][1,j])
end


function quadrule(op::HH3DSingleLayerFDBIO, test_refspace::subReferenceSpace,
    trial_refspace::subReferenceSpace, i, test_element, j, trial_element, quadrature_data)

    # tol, hits = 1e-10, 0
    # for t in getelementVertices(test_element)
    #     for s in getelementVertices(trial_element)
    #         norm(t-s) < tol && (hits +=1; break)
    # end end
    #
    # test_quadpoints  = quadrature_data[1]
    # trial_quadpoints = quadrature_data[2]

    # hits != 0 && return WiltonSEStrategy(
    #     test_quadpoints[1,i],
    #     DoubleQuadStrategy(
    #         test_quadpoints[1,i],
    #         trial_quadpoints[1,j]))

    return DoubleQuadStrategy(
        quadrature_data[1][1,i],
        quadrature_data[2][1,j])
    # return SauterSchwabStrategy(hits)
     # test_qd = qd[1]
     # trial_qd = qd[2]
     #
     # DoubleQuadStrategy(
     #    test_qd[1,i], # rule 1 on test element j
     #    trial_qd[1,j] # rule 1 on trial element i
     # )
end


function quadrule(op::Helmholtz3DOp,
    test_refspace::RefSpace, trial_refspace::RefSpace,
    i, test_element, j, trial_element, quadrature_data)

    test_quadpoints  = quadrature_data[1]
    trial_quadpoints = quadrature_data[2]

    return DoubleQuadStrategy(
        quadrature_data[1][1,i],
        quadrature_data[2][1,j])
end

function quadrule(op::Helmholtz3DOp,
    test_refspace::RTRefSpace, trial_refspace::RTRefSpace,
    i, test_element, j, trial_element, quadrature_data)

    test_quadpoints  = quadrature_data[1]
    trial_quadpoints = quadrature_data[2]

    return DoubleQuadStrategy(
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





HH3DSingleLayerFDBIO(gamma) = HH3DSingleLayerFDBIO(one(gamma), gamma)

regularpart(op::HH3DSingleLayerFDBIO) = HH3DSingleLayerReg(op.alpha, op.gamma)
singularpart(op::HH3DSingleLayerFDBIO) = HH3DSingleLayerSng(op.alpha, op.gamma)

function integrand(op::Union{HH3DSingleLayerFDBIO,HH3DSingleLayerReg},
        kernel, test_values, test_element, trial_values, trial_element)

    α = op.alpha
    G = kernel.green

    g = test_values.value
    f = trial_values.value

    α*dot(g, G*f)
end


function innerintegrals!(op::HH3DSingleLayerSng, test_neighborhood,
        test_refspace::LagrangeRefSpace{T,0} where {T},
        trial_refspace::LagrangeRefSpace{T,0} where {T},
        test_elements, trial_element, zlocal, quadrature_rule::WiltonSEStrategy, dx)

    γ = op.gamma
    α = op.alpha

    s1, s2, s3 = trial_element.vertices

    x = cartesian(test_neighborhood)
    n = normalize((s1-s3)×(s2-s3))
    ρ = x - dot(x - s1, n) * n

    scal, vec = WiltonInts84.wiltonints(s1, s2, s3, x, Val{1})
    ∫G = (scal[2] - γ*scal[3] + 0.5*γ^2*scal[4]) / (4π)

    zlocal[1,1] += α * ∫G * dx
    return nothing
end



function integrand(biop::HH3DDoubleLayer,
        kernel, fp, mp, fq, mq)

    nq = normal(mq)
    fp[1] * dot(nq, -kernel.gradgreen) * fq[1]
end



function integrand(biop::HH3DDoubleLayerTransposed,
        kernel, fp, mp, fq, mq)

    np = normal(mp)
    fp[1] * dot(np, kernel.gradgreen) * fq[1]
end


module Helmholtz3D

    using ..BEAST
    const Mod = BEAST

    function singlelayer(;
        gamma=nothing,
        wavenumber=nothing,
        alpha=nothing)

        if (gamma == nothing) && (wavenumber == nothing)
            error("Supply one of (not both) gamma or wavenumber")
        end

        if (gamma != nothing) && (wavenumber != nothing)
            error("Supply one of (not both) gamma or wavenumber")
        end

        if gamma == nothing
            if iszero(real(wavenumber))
                gamma = -imag(wavenumber)
            else
                gamma = im*wavenumber
            end
        end

        @assert gamma != nothing
        alpha == nothing && (alpha = one(real(typeof(gamma))))

        Mod.HH3DSingleLayerFDBIO(alpha,gamma)
    end

    hypersingular(;
            gamma=error("propagation constant is a required argument"),
            alpha=gamma^2,
            beta=one(gamma)) =
        Mod.HH3DHyperSingularFDBIO(alpha, beta, gamma)

    planewave(;
            direction=error("direction is a required argument"),
            wavenumber=error("wavenumber is a required arguement"),
            amplitude=one(eltype(direction))) =
        Mod.HH3DPlaneWave(direction, wavenumber, amplitude)

    doublelayer(;gamma=error("gamma missing"), alpha=one(gamma)) =
        Mod.HH3DDoubleLayer(alpha, gamma)

    doublelayer_transposed(;gamma=error("gamma missing"), alpha=one(gamma)) =
        Mod.HH3DDoubleLayerTransposed(alpha, gamma)
end

export Helmholtz3D
