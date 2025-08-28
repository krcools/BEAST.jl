
struct HH2DPlaneWave{P,K,T} <: Functional{T}
    direction::P
    gamma::K
    amplitude::T
end

function (f::HH2DPlaneWave)(r)
    d = f.direction
    γ = f.gamma
    a = f.amplitude
    return a * exp(-γ*dot(d,r))
end

scalartype(f::HH2DPlaneWave{P,K,T}) where {P,K,T} = promote_type(eltype(P), K, T)

struct gradHH2DPlaneWave{P,K,T} <: Functional{T}
    direction::P
    gamma::K
    amplitude::T
end

function (f::gradHH2DPlaneWave)(r)
    d = f.direction
    γ = f.gamma
    a = f.amplitude

    return - a * γ * d * exp(-γ * dot(d, r))
end

scalartype(f::gradHH2DPlaneWave{P,K,T}) where {P,K,T} = promote_type(eltype(P), K, T)

struct curlHH2DPlaneWave{P,K,T} <: Functional{T}
    direction::P
    polarization::P
    gamma::K
    amplitude::T
end

function (f::curlHH2DPlaneWave)(r)
    d = f.direction
    p = f.polarization
    γ = f.gamma
    a = f.amplitude
    return a * p * exp(-γ*dot(d,r))
end

scalartype(f::curlHH2DPlaneWave{P,K,T}) where {P,K,T} = promote_type(eltype(P), K, T)

function curl(m::HH2DPlaneWave)
    d = m.direction
    polarization = -SVector(+d[2], -d[1]) # By using a minus sign here, the amplitude stays positive
    return curlHH2DPlaneWave(d, polarization, m.gamma, m.amplitude * (m.gamma))
end

*(a::Number, m::HH2DPlaneWave) = HH2DPlaneWave(m.direction, m.gamma, a * m.amplitude)
*(a::Number, m::gradHH2DPlaneWave) = gradHH2DPlaneWave(m.direction, m.gamma, a * m.amplitude)
*(a::Number, m::curlHH2DPlaneWave) = curlHH2DPlaneWave(m.direction, m.polarization, m.gamma, a * m.amplitude)


"""
    HH2DMonopole

Potential of a monopole-type point source (e.g., of an electric charge)
"""
struct HH2DMonopole{P,K,T}
    position::P
    gamma::K
    amplitude::T
end

scalartype(x::HH2DMonopole{P,K,T}) where {P,K,T} = promote_type(eltype(P), K, T)

function (f::HH2DMonopole)(r::CompScienceMeshes.MeshPointNM)
    return f(cartesian(r))
end

function (f::HH2DMonopole)(r)
    γ = f.gamma
    p = f.position
    a = f.amplitude

    return a / (4 * im) * hankelh2(0, -im*γ * norm(r - p))
end

struct curlHH2DMonopole{P,K,T}
    position::P
    gamma::K
    amplitude::T
end

scalartype(x::curlHH2DMonopole{P,K,T}) where {P,K,T} = promote_type(eltype(P), K, T)

function (f::curlHH2DMonopole)(r)
    a = f.amplitude
    γ = f.gamma
    p = f.position
    vecR = r-p
    R = norm(vecR)

    return a / (4 * im) * 1/R * (-im*γ) * (-hankelh2(1,(-im*γ*R))) * (SVector(-vecR[2],vecR[1]))
end

function (f::curlHH2DMonopole)(mp::CompScienceMeshes.MeshPointNM)
    fieldval = f(cartesian(mp))
    t = tangents(mp, 1)
    return dot(t, fieldval)
end

function curl(m::HH2DMonopole)
    return curlHH2DMonopole(m.position, m.gamma, m.amplitude)
end

struct gradHH2DMonopole{P,K,T}
    position::P
    gamma::K
    amplitude::T
end

scalartype(x::gradHH2DMonopole{P,K,T}) where {P,K,T} = promote_type(eltype(P), K, T)

function (f::gradHH2DMonopole)(r)
    a = f.amplitude
    γ = f.gamma
    p = f.position
    vecR = r - p
    R = norm(vecR)

    return a / (4 * im) * -hankelh2(1, -im*γ * R) * (-im*γ) * vecR / R
end

function grad(m::HH2DMonopole)
    return gradHH2DMonopole(m.position, m.gamma, m.amplitude)
end

*(a::Number, m::HH2DMonopole) = HH2DMonopole(m.position, m.gamma, a * m.amplitude)
*(a::Number, m::curlHH2DMonopole) = curlHH2DMonopole(m.position, m.gamma, a * m.amplitude)
*(a::Number, m::gradHH2DMonopole) = gradHH2DMonopole(m.position, m.gamma, a * m.amplitude)

dot(::NormalVector, m::gradHH2DMonopole) = NormalDerivative(HH2DMonopole(m.position, m.gamma, m.amplitude))

function (f::NormalDerivative{T,F})(manipoint) where {T,F<:HH2DMonopole}
    m = f.field
    grad_m = grad(m)
    n = normal(manipoint)
    r = cartesian(manipoint)
    return dot(n, grad_m(r))
end

"""
    HH2DDirectedMonopole

Potential of a monopole-type point source (e.g., of an electric charge)
"""
struct HH2DDirectedMonopole{P,K,T}
    position::P
    direction::P
    gamma::K
    amplitude::T
end

scalartype(x::HH2DDirectedMonopole{P,K,T}) where {P,K,T} = promote_type(eltype(P), K, T)

function (f::HH2DDirectedMonopole)(r)
    γ = f.gamma
    p = f.position
    d = f.direction
    a = f.amplitude

    x = r[1] - p[1]
    y = r[2] - p[2]
    
    Ix = d[1]
    Iy = d[2]

    R = norm(r - p)

    k = -im*γ

    dhankelh2(x) = -hankelh2(1, x) # H₀^(2)'

    return a / (4 * im) * k * dhankelh2(k*R) * (Iy * x - Ix * y) / R
end

struct curlHH2DDirectedMonopole{P,K,T} <: Functional{T}
    position::P
    direction::P
    gamma::K
    amplitude::T
end

function (f::curlHH2DDirectedMonopole)(r)
    γ = f.gamma
    p = f.position
    d = f.direction
    a = f.amplitude

    x = r[1] - p[1]
    y = r[2] - p[2]
    
    Ix = d[1]
    Iy = d[2]

    R = norm(r - p)

    k = -im*γ

    dhankelh2(x) = -hankelh2(1, x)  # H₀^(2)'
    ddhankelh2(x) = hankelh2(1, x)/x - hankelh2(0, x)  # H₀^(2)''

    X = ddhankelh2(k * R) * k * (Iy * x - Ix * y) / R^2 * y + dhankelh2(k * R) * (-Ix / R  - (Iy * x - Ix * y) * y / R^3)
    Y = ddhankelh2(k * R) * k * (Iy * x - Ix * y) / R^2 * x + dhankelh2(k * R) * (+Iy / R  - (Iy * x - Ix * y) * x / R^3)
 
    return a / (4 * im) * k * SVector(X, -Y)
end

scalartype(f::curlHH2DDirectedMonopole{P,K,T}) where {P,K,T} = promote_type(eltype(P), K, T)

function curl(m::HH2DDirectedMonopole)
    return curlHH2DDirectedMonopole(m.position, m.direction, m.gamma, m.amplitude)
end

*(a::Number, m::HH2DDirectedMonopole) = HH2DDirectedMonopole(
    m.position,
    m.direction,
    m.gamma,
    a * m.amplitude
)

#=
*(a::Number, m::gradHH2DDirectedMonopole) = gradHH2DDirectedMonopole(
    m.position,
    m.direction,
    m.gamma,
    a * m.amplitude
)
=#

*(a::Number, m::curlHH2DDirectedMonopole) = curlHH2DDirectedMonopole(
    m.position,
    m.direction,
    m.gamma,
    a * m.amplitude
)

struct ScalarTrace{T,F} <: Functional{T}
    field::F
end

ScalarTrace(f::F) where {F} = ScalarTrace{scalartype(f), F}(f)
ScalarTrace{T}(f::F) where {T,F} = ScalarTrace{T,F}(f)

strace(f::F) where {F} = ScalarTrace{scalartype(f), F}(f)
strace(f, mesh::Mesh) = ScalarTrace(f)

@testitem "scalar trace" begin
    using CompScienceMeshes
    fn = joinpath(dirname(pathof(BEAST)), "../test/assets/sphere45.in")
    Γ = readmesh(fn)
    X = lagrangecxd0(Γ)

    U = Helmholtz3D.planewave(wavenumber=1.0, direction=ẑ)
    u = strace(U)

    ux = assemble(u, X)
    @test u isa BEAST.ScalarTrace
    @test BEAST.scalartype(u) == ComplexF64
end

(s::ScalarTrace)(x) = s.field(cartesian(x))
integrand(s::ScalarTrace, tx, fx) = dot(tx.value, fx)
scalartype(s::ScalarTrace{T}) where {T} = T

# We assume an orthogonal system (t, n, z)
# This deviates from Moritas book, where they assume
# (n, t, z)
mutable struct TangentTrace{T,F} <: Functional{T}
    field::F
end

TangentTrace(f::F) where {F} = TangentTrace{scalartype(f), F}(f)
TangentTrace{T}(f::F) where {T,F} = TangentTrace{T,F}(f)
scalartype(s::TangentTrace{T}) where {T} = T

function (ϕ::TangentTrace)(p)
    F = ϕ.field
    x = cartesian(p)
    t = tangents(p,1)
    t = t / norm(t)

    return dot(t, F(x))
end

shapevals(f::Functional, ϕ, ts) = shapevals(ValOnly(), ϕ, ts)

function (f::NormalDerivative{T,F})(manipoint) where {T,F<:HH2DPlaneWave}
    d = f.field.direction
    γ = f.field.gamma
    a = f.field.amplitude
    n = normal(manipoint)
    r = cartesian(manipoint)
    return -γ*a * dot(d,n) * exp(-γ*dot(d,r))
end

dot(::NormalVector, m::gradHH2DPlaneWave) = NormalDerivative(HH2DPlaneWave(m.direction, m.gamma, m.amplitude))

struct PlaneWaveDirichlet{T,P} <: Functional{T}
    wavenumber::T
    direction::P
end

scalartype(x::PlaneWaveDirichlet{T}) where {T} = complex(T)

struct PlaneWaveNeumann{T,P} <: Functional{T}
    wavenumber::T
    direction::P
end

scalartype(x::PlaneWaveNeumann{T}) where {T} = complex(T)


function (field::PlaneWaveDirichlet)(mp)

    wavenumber = field.wavenumber
    direction  = field.direction

    cart = cartesian(mp)
    exp(-im * wavenumber * dot(direction, cart))
end


function(field::PlaneWaveNeumann)(mp)

    wavenumber = field.wavenumber
    direction  = field.direction

    cart = cartesian(mp)
    norm = normal(mp)

    wave = exp(-im * wavenumber * dot(direction, cart))
    grad = -im * wavenumber * direction * wave

    d = norm[1] * grad[1]
    for i in 2:length(norm)  d += norm[i] * grad[i]  end
    return d
end

function integrand(pw::PlaneWaveDirichlet, sv, fx)
    tx = sv[1]
    return dot(tx, fx)
end

function integrand(pw::PlaneWaveNeumann, sv, fx)
    tx = sv[1]
    d = tx[1] * fx[1]
    for i in 2:length(tx) d += tx[i]*fx[i] end
    return d
end
