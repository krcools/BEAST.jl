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

HH2DNear{T, K} = Union{
    HH2DSingleLayerNear{T, K},
    HH2DDoubleLayerNear{T, K},
    HH2DDoubleLayerTransposedNear{T, K},
    HH2DHyperSingularNear{T, K}
} where {T, K}

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

mutable struct KernelValsHH2DNear{K1,K2,F}
    gamma::K1
    vect::SVector{2,F}
    dist::F
    green::K2
    gradgreen::SVector{2,K2}
    n::SVector{2,F}
    t::SVector{2,F}
end
# WARNING: Matching the HH3D we use x for the sources
# TODO: Discuss if we change this?

# Kernel values for 2D Laplace BIE operators
function kernelvals(op::HH2DNear{T, K}, y, p) where {T, K <: Val{0}}

    xc = cartesian(p)
    yc = cartesian(y)
    r = yc - xc
    R = norm(r)

    #Logarithmic singularity of Hankel function: Gk = -1/(2π) * log(kR)
    green = -1/(2π) * log(R)
    gradgreen = -1/(2π) * (1/R) * (r/R)

    n = normal(p)
    t = tangents(p, 1) / norm(tangents(p, 1))

    return KernelValsHH2DNear(nothing, r, R, green, gradgreen, n, t)
end

# Kernel values for 2D Helmholtz equation BIE operators
function kernelvals(op::HH2DNear{T, K}, y, p) where {T, K <: Real}
    γ = op.gamma

    xc = cartesian(p)
    yc = cartesian(y)
    r = yc - xc
    R = norm(r)

    green = 1/(2π) * besselk(0, γ * R)
    gradgreen = -γ/(2π) * besselk(1, γ * R) * (r/R)

    n = normal(p)
    t = tangents(p, 1) / norm(tangents(p, 1))

    return KernelValsHH2DNear(γ, r, R, green, gradgreen, n, t)
end

# Kernel values for 2D Helmholtz equation BIE operators
function kernelvals(op::HH2DNear{T, K}, y, p) where {T, K <: Complex}

    γ = op.gamma
    # Even though the evaluation of the Hankel function delivers
    # in general a complex output (sole exception if wavenumber is purely imaginary)
    # the Hankel function is evaluated much faster if the wavenumber, and thus
    # the argument of the Hankelfunction, is purely real
    if iszero(real(γ))
        k = imag(γ)
    else
        k = -im*γ
    end
    xc = cartesian(p)
    yc = cartesian(y)
    r = yc - xc
    R = norm(r)

    kR = k * R
    hankels = hankelh2.([0 1], kR)
    green = - im / 4 * hankels[1]
    gradgreen = k * im / 4 * hankels[2] * r / R

    n = normal(p)
    t = tangents(p, 1) / norm(tangents(p, 1))

    return KernelValsHH2DNear(γ, r, R, green, gradgreen, n, t)
end


function integrand(op::HH2DSingleLayerNear, krn, x, g, y)

    α = op.alpha
    G = krn.green
    ϕ = g.value

    return α * G * ϕ
end

function integrand(op::HH2DDoubleLayerNear, krn, x, g, y)

    α = op.alpha
    ∇G = -krn.gradgreen
    ny = krn.n
    ϕ = g.value
    ∂G∂n = dot(ny, ∇G)

    return α * ∂G∂n * ϕ
end

function integrand(op::HH2DDoubleLayerTransposedNear, krn, x, g, y)

    α = op.alpha
    ∇G = krn.gradgreen
    ϕ = g.value

    return α * ∇G * ϕ
end

function integrand(op::HH2DHyperSingularNear, krn, x, g, y)

    α = op.alpha
    β = op.beta
    G = krn.green
    ∇G = -krn.gradgreen
    ϕ = g.value
    dϕ = g.derivative
    tx = krn.t

    # Based on (2.86) in Kumagai et al, “Integral equation methods for electromagnetics”
    # The formula returns the electric field, but we like to match the behavior of the
    # dyadic hypersingular operator and this requires an additional rotation
    return cross(ẑ, (-α * G * tx * ϕ + β * ∇G * dϕ))
end

function integrand(op::HH2DHyperSingularNear{T, K}, krn, x, g, y) where {T, K <: Val{0}}

    β = op.beta
    ∇G = -krn.gradgreen
    dϕ = g.derivative

    return cross(ẑ, β * ∇G * dϕ)
end

function defaultquadstrat(op::BEAST.HH2DNear, X::BEAST.Space)
    return defaultquadstrat(op, X, X)
end
defaultquadstrat(op::HH2DNear, tfs, bfs) = SingleNumQStrat(20)
quaddata(op::HH2DNear,rs,els,qs::SingleNumQStrat) = quadpoints(rs,els,(qs.quad_rule,))
quadrule(op::HH2DNear,refspace,p,y,q,el,qdata,qs::SingleNumQStrat) = qdata[1,q]