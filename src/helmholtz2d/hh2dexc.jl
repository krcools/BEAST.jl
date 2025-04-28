
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

function (f::HH2DPlaneWave)(r::CompScienceMeshes.MeshPointNM)
    return f(cartesian(r))
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

    return -γ * d * exp(-γ * dot(d, r))
end

function (f::gradHH2DPlaneWave)(mp::CompScienceMeshes.MeshPointNM)
    r = cartesian(mp)
    return dot(normal(mp), f(r))
end

scalartype(f::gradHH2DPlaneWave{P,K,T}) where {P,K,T} = promote_type(eltype(P), K, T)

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

shapevals(f::Functional, ϕ, ts) = shapevals(ValOnly(), ϕ, ts)

function (f::NormalDerivative{T,F})(manipoint) where {T,F<:HH2DPlaneWave}
    d = f.field.direction
    γ = f.field.gamma
    a = f.field.amplitude
    n = normal(manipoint)
    r = cartesian(manipoint)
    return -γ*a * dot(d,n) * exp(-γ*dot(d,r))
end

*(a::Number, m::HH2DPlaneWave) = HH2DPlaneWave(m.direction, m.gamma, a * m.amplitude)
*(a::Number, m::gradHH2DPlaneWave) = gradHH2DPlaneWave(m.direction, m.gamma, a * m.amplitude)

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
