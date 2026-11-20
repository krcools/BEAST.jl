abstract type OsrcEpsilonParameter end

"""
    OsrcStandardEpsilon(curvature)

Standard complex wavenumber damping parameter.
Designed for optimal performance of the OSRC operator on the sphere at high-frequency, defined by:
    ϵ = 0.39κ^(1/3)H^(2/3),
where `H` is the mean curvature of the surface.
"""

struct OsrcStandardEpsilon{T <: Real} <: OsrcEpsilonParameter
    curvature::T
end

function getepsilon(wavenumber::Real, complexdamping::OsrcStandardEpsilon)
    return 0.39*wavenumber^(1/3)*complexdamping.curvature^(2/3)
end

struct RationalFunction
    R0::ComplexF64
    Aj::Vector{ComplexF64}
    Bj::Vector{ComplexF64}
    Np::Int
end

function scalarapproximation(rational::RationalFunction)
    return z -> rational.R0 - sum(rational.Aj[i] / (rational.Bj[i]*(1 + rational.Bj[i] * z)) for i in 1:rational.Np)
end

abstract type RationalApproximant end
"""
    construct(approximant)

Construct a rational function representation of an approximant.
A `RationalApproximant` specifies a procedure for approximating a function by a rational function. 
"""
function construct end
scalarapproximation(approximant::RationalApproximant, args...;kwargs...) = scalarapproximation(construct(approximant, args...;kwargs...))

"""
    PadeApproximant(Np, θ_p)

Rotated branch cut Padé rational approximation of √(1+z), 
with `Np` number of terms and a rotation angle `θ_p` in the complex plane.
"""
struct PadeApproximant{T <: Real} <: RationalApproximant
    Np::Int
    θ_p::T
end

function construct(approximant::PadeApproximant)
    (Np, θ_p) = (approximant.Np, approximant.θ_p)
    imag_conv = -im 
    # get the real and rotated pade coefficients
    aj, bj = ComplexF64[], ComplexF64[]
    Aj, Bj = ComplexF64[], ComplexF64[]
    for j in 1:Np
        a =  2/(2*Np+1) * (sin(j*pi/(2*Np+1)))^2
        b = cos(j*pi/(2*Np+1))^2
        A = exp(-imag_conv*θ_p/2) * a / (1 + b*(exp(-imag_conv * θ_p) - 1))^2
        B = exp(-imag_conv*θ_p) * b / (1 + b*(exp(-imag_conv * θ_p) - 1))
        push!(aj, a)
        push!(bj, b)
        push!(Aj, A)
        push!(Bj, B)
    end
    # compute the constant term
    C0 = exp(imag_conv*θ_p/2) * (1 + sum(aj[j]*(exp(-imag_conv*θ_p)-1)/(1+bj[j]*(exp(-imag_conv*θ_p)-1)) for j in 1:Np))
    R0 = C0 + sum(Aj[j]/Bj[j] for j in 1:Np)
    return RationalFunction(R0, Aj, Bj, Np)
end

struct InverseApproximant{A <: RationalApproximant} <: RationalApproximant
    approximant::A
end

function construct(inverseapproximant::InverseApproximant)
    rational = construct(inverseapproximant.approximant)
    return inverserational(rational)
end

"""
    inverserational(rational)

Obtain the rational representation of the inverse of a rational function.
"""
function inverserational(rational::RationalFunction)
    (R0, Aj, Bj, Np) = (rational.R0, rational.Aj, rational.Bj, rational.Np)
    # compute the roots/zeros
    poles = -1.0 ./ Bj
    residues = Aj ./ (Bj .^ 2)
    M = LinearAlgebra.Diagonal(poles) + ones(ComplexF64, Np) * transpose(residues ./ R0)
    roots = LinearAlgebra.eigvals(Matrix(M))
    # compute new residues
    dR = [sum(residues[j] / (roots[k] - poles[j])^2 for j in 1:Np) for k in 1:Np]

    # convert to rational function representation
    Binv = -1.0 ./ roots
    Ainv = -(Binv .^ 2) ./ dR
    R0inv = 1.0 / R0
    return RationalFunction(R0inv, Ainv, Binv, Np)
end

struct SliceLinearMap
    P::LinearMap
    Q::LinearMap
    rows::Vector{Int}
    cols::Vector{Int}
end

struct SliceRow
    rows::Vector{Int}
    n::Int
end

function (f::SliceRow)(vector::AbstractVector)
    return vector[f.rows]
end

struct SliceColumn
    cols::Vector{Int}
    n::Int
end

function (f::SliceColumn)(vector::AbstractVector)
    y = zeros(eltype(vector), f.n)
    y[f.cols] = vector
    return y
end

"""
    SliceLinearMap(n, m, rows, cols)

Obtain a left and right multiplier LinearMap which efficiently slice out the rows: `rows` and columns: `cols` of a `n × m` LinearMap.
"""
function SliceLinearMap(n::Int, m::Int, rows::AbstractVector{<:Int}, cols::AbstractVector{<:Int})
    n_rows, n_cols = length(rows), length(cols)
    P = LinearMap(SliceRow(rows, n), SliceColumn(rows, n), n_rows, n; ismutating=false)
    Q =  LinearMap(SliceColumn(cols, m), SliceRow(cols, m), m, n_cols; ismutating=false)
    return SliceLinearMap(P, Q, rows, cols)
end

# OSRC operators
abstract type OsrcOperator <: Operator end

struct MteOsrc{T <: Complex, A <: RationalApproximant} <: OsrcOperator
    κ_ϵ::T
    approximant::A
    order::Int
end

struct EtmOsrc{T <: Complex, A <: RationalApproximant} <: OsrcOperator
    κ_ϵ::T
    approximant::A
    order::Int
end

function constructosrc(::Type{T}, wavenumber::Real, ϵ::Real, approximant::RationalApproximant; order=0) where {T <: OsrcOperator}
    κ_ϵ = wavenumber - im*ϵ
    return T(κ_ϵ, approximant, order)
end

function constructosrc(::Type{T}, wavenumber::Real, approximant::RationalApproximant, complexdamping::OsrcEpsilonParameter; kwargs...) where {T <: OsrcOperator}
    ϵ = getepsilon(wavenumber, complexdamping)
    return constructosrc(T, wavenumber, ϵ, approximant; kwargs...)
end

MteOsrc(wavenumber::Real, approximant::RationalApproximant, args...; kwargs...) = constructosrc(MteOsrc, wavenumber, approximant, args...; kwargs...)
EtmOsrc(wavenumber::Real, approximant::RationalApproximant, args...; kwargs...) = constructosrc(EtmOsrc, wavenumber, InverseApproximant(approximant), args...; kwargs...)
# default constructors
MteOsrc(wavenumber::Real, Np::Int, θ_p::Real, curvature::Real; kwargs...) = MteOsrc(wavenumber, PadeApproximant(Np, θ_p), OsrcStandardEpsilon(curvature); kwargs...)
EtmOsrc(wavenumber::Real, Np::Int, θ_p::Real, curvature::Real; kwargs...) = EtmOsrc(wavenumber, PadeApproximant(Np, θ_p), OsrcStandardEpsilon(curvature); kwargs...)

function scalartype(op::OsrcOperator)
    T = scalartype(op.κ_ϵ)
    return Complex{T}
end

function assemble(op::MteOsrc, X::Space, Y::Space; quadstrat=defaultquadstrat)
    rot_X = n × X
    rot_Y = n × Y
    # create auxiliary basis functions
    Lagrangeq = BEAST.lagrangec0(X.geo; order=op.order+1)
    grad_Lagrangeq = n×BEAST.curl(Lagrangeq)
    N_X, N_Y = numfunctions(X), numfunctions(Y)

    # assemble all submatrices
    Id = BEAST.Identity()
    G = assemble(Id, rot_X, rot_Y; quadstrat)
    A = assemble(Id, BEAST.divergence(X), BEAST.divergence(Y); quadstrat)
    N_ϵ = (1/op.κ_ϵ)^2 * A
    M = assemble(Id, Lagrangeq, Lagrangeq; quadstrat)
    K_ϵ = (op.κ_ϵ)^2 * M
    L = -assemble(Id, rot_X, grad_Lagrangeq; quadstrat)
    L_transpose = -assemble(Id, grad_Lagrangeq, rot_Y; quadstrat)

    # construct rational function and obtain its parameters
    rational = construct(op.approximant)
    (R0, Aj, Bj, Np) = (rational.R0, rational.Aj, rational.Bj, rational.Np)

    # construct the auxiliary block system and extract the inverse Schur complement
    function constructΠinverse(j)
        AXY = [G-Bj[j]*N_ϵ       Bj[j]*L
                L_transpose     K_ϵ]
        SliceAXY = SliceLinearMap(size(AXY)..., collect(1:N_X), collect(1:N_Y))
        return SliceAXY.P * BEAST.lu(AXY) * SliceAXY.Q
    end
    MtE_map = - BEAST.lu(G - N_ϵ) * (R0 * LinearAlgebra.I - G * sum(Aj[j]/Bj[j] * constructΠinverse(j) for j in 1:Np))     
    return MtE_map
end

function assemble(op::EtmOsrc, X::Space, Y::Space; quadstrat=defaultquadstrat)
    rot_X = n × X
    rot_Y = n × Y
    # create auxiliary basis functions
    Lagrangeq = BEAST.lagrangec0(X.geo; order=op.order+1)
    grad_Lagrangeq = n×BEAST.curl(Lagrangeq)
    N_X, N_Y = numfunctions(X), numfunctions(Y)

    # assemble all submatrices
    Id = BEAST.Identity()
    G = assemble(Id, rot_X, rot_Y; quadstrat)
    A = assemble(Id, BEAST.divergence(X), BEAST.divergence(Y); quadstrat)
    N_ϵ = (1/op.κ_ϵ)^2 * A
    M = assemble(Id, Lagrangeq, Lagrangeq; quadstrat)
    K_ϵ = (op.κ_ϵ)^2 * M
    L = -assemble(Id, rot_X, grad_Lagrangeq; quadstrat)
    L_transpose = -assemble(Id, grad_Lagrangeq, rot_Y; quadstrat)

    # construct rational function and obtain its parameters
    rational = construct(op.approximant)
    (R0, Aj, Bj, Np) = (rational.R0, rational.Aj, rational.Bj, rational.Np)

    # construct the auxiliary block system and extract the inverse Schur complement
    function constructΠinverse(j)
        AXY = [G-Bj[j]*N_ϵ       Bj[j]*L
                L_transpose     K_ϵ]
        SliceAXY = SliceLinearMap(size(AXY)..., collect(1:N_X), collect(1:N_Y))
        return SliceAXY.P * BEAST.lu(AXY) * SliceAXY.Q
    end
    iG = BEAST.lu(G)
    EtM_map = iG * (R0 * LinearAlgebra.I - G * sum(Aj[j]/Bj[j] * constructΠinverse(j) for j in 1:Np)) * (G - N_ϵ)
    return EtM_map
end