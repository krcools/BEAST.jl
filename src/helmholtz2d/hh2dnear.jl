
struct HH2DSingleLayerNear{T, K}
    alpha::T
    gamma::K

    function HH2DSingleLayerNear(alpha, gamma)
        gamma = hh2d_makegammacomplexifneeded(gamma)
        return new{typeof(alpha),typeof(gamma)}(alpha, gamma)
    end

    function HH2DSingleLayerNear(op::HH2DSingleLayerFDBIO)
        return HH2DSingleLayerNear(op.alpha, op.gamma)
    end 
end

struct HH2DDoubleLayerNear{T, K}
    alpha::T
    gamma::K

    function HH2DDoubleLayerNear(alpha, gamma)
        gamma = hh2d_makegammacomplexifneeded(gamma)
        return new{typeof(alpha),typeof(gamma)}(alpha, gamma)
    end

    function HH2DDoubleLayerNear(op::HH2DDoubleLayerFDBIO)
        return HH2DDoubleLayerNear(op.alpha, op.gamma)
    end 
end

struct HH2DDoubleLayerTransposedNear{T, K}
    alpha::T
    gamma::K

    function HH2DDoubleLayerTransposedNear(alpha, gamma)
        gamma = hh2d_makegammacomplexifneeded(gamma)
        return new{typeof(alpha),typeof(gamma)}(alpha, gamma)
    end

    function HH2DDoubleLayerTransposedNear(op::HH2DDoubleLayerTransposedFDBIO)
        return HH2DDoubleLayerTransposedNear(op.alpha, op.gamma)
    end 
end

struct HH2DHyperSingularNear{T, K}
    alpha::T
    beta::T
    gamma::K

    function HH2DHyperSingularNear(alpha, beta, gamma)
        gamma = hh2d_makegammacomplexifneeded(gamma)
        return new{typeof(alpha),typeof(gamma)}(alpha, beta, gamma)
    end

    function HH2DHyperSingularNear(op::HH2DHyperSingularFDBIO)
        return HH2DHyperSingularNear(op.alpha, op.beta, op.gamma)
    end
end


HH2DNear = Union{HH2DSingleLayerNear, HH2DDoubleLayerNear, HH2DDoubleLayerTransposedNear, HH2DHyperSingularNear}
defaultquadstrat(op::HH2DNear, basis) = SingleNumQStrat(20)
quaddata(op::HH2DNear,rs,els,qs::SingleNumQStrat) = quadpoints(rs,els,(qs.quad_rule,))
quadrule(op::HH2DNear,refspace,p,y,q,el,qdata,qs::SingleNumQStrat) = qdata[1,q]

# TODO: This function is useful for changing the gradient to a curl
# However, it is a case of type piracy.
# (Evidently, the case is not supported by cross function due to mismatch of length)

@inline function LinearAlgebra.cross(z::SVector{3,T1}, v::SVector{2,T2}) where {T1,T2}
    z == ẑ || throw(ArgumentError("Only ẑ = (0,0,1) supported"))
    return SVector{2,T2}(-v[2],  v[1])
end

@inline function LinearAlgebra.cross(v::SVector{2,T2}, z::SVector{3,T1}) where {T1,T2}
    z == ẑ || throw(ArgumentError("Only ẑ = (0,0,1) supported"))
    return SVector{2,T2}( v[2], -v[1])
end

# WARNING: Matching the HH3D we use x for the sources
# TODO: Discuss if we change this?
function kernelvals(op::HH2DNear,y,p)

    γ = op.gamma
    if iszero(real(γ))
        k = imag(γ)
    else
        k = -im*γ
    end
    x = cartesian(p)
    r = y - x
    R = norm(r)

    kR = k * R
    hankels = hankelh2.([0 1], kR)
    #green = - op.alpha * im / 4 * hankels[1]
    #gradgreen = op.alpha * k * im / 4 * hankels[2] * r / R
    green = - im / 4 * hankels[1]
    gradgreen =  k * im / 4 * hankels[2] * r / R

    ddhankel = hankels[1] - 1/(kR) * hankels[2]

    #txty = dot(normal(tgeo), normal(bgeo))

    nx = normal(p)
    tx = tangents(p, 1) / norm(tangents(p, 1))

    # Needed for hypersingular operator
    gradnxgradgreen = ddhankel * (k^2/R^2 * dot(r, nx) * r)
    gradnxgradgreen += k * hankels[2] * (nx / R - r/R^3 * dot(r, nx)) 
    gradnxgradgreen *= - im / 4

    (;γ, r, R, green, gradgreen, nx, tx, gradnxgradgreen)
end

function integrand(op::HH2DSingleLayerNear, krn, y, f, p)

    α = op.alpha
    G = krn.green
    fx = f.value

    return α * G * fx
end

function integrand(op::HH2DDoubleLayerNear, krn, y, f, p)

    α = op.alpha
    ∇G = -krn.gradgreen
    nx = krn.nx

    fx = f.value

    ∂G∂n = nx ⋅ ∇G

    return α * ∂G∂n * fx
end

function integrand(op::HH2DDoubleLayerTransposedNear, krn, y, f, p)

    α = op.alpha
    ∇G = krn.gradgreen
    fx = f.value

    return α * ∇G * fx
end

function integrand(op::HH2DHyperSingularNear, krn, y, f, p)
    α = op.alpha
    β = op.beta

    G = krn.green
    ∇G = -krn.gradgreen

    #∇nₓ∇ₓG = krn.gradnxgradgreen

    fx = f.value
    dfx = f.derivative
    tx = krn.tx

    # Based on (2.86) in Kumagai et al, “Integral equation methods for electromagnetics”
    # The formula returns the electric field, but we like to match the behavior of the
    # dyadic hypersingular operator and this requires an additional rotation
    return cross(ẑ, -α * G * tx * fx + β * ∇G * dfx)
end
