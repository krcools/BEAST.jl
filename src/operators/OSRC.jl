using SparseArrays
using LinearAlgebra
using BEAST
import Polynomials

struct MtE_OSRC_op <: Operator
    wavenumber::Float64
    Np::Int
    θ_p::Float64
    curvature::Float64
    aj::Vector{ComplexF64}
    bj::Vector{ComplexF64}
    Aj::Vector{ComplexF64}
    Bj::Vector{ComplexF64}
end

# First a rotating branch-cut rational Padé approximation of the square root function ``\sqrt{1+z^2}`` is implemented.
imag_conv = -im     # The paper uses mathematical time-harmonic convention ``e^{-i \omega t}`` -> swap to match BEAST time-harmonic convention ``e^{i \omega t}``.
function MtE_OSRC_op(wavenumber::Float64, Np::Int, θ_p::Float64, curvature::Float64)
    # get the real and rotated pade coefficients
    aj = ComplexF64[]
    bj = ComplexF64[]
    Aj = ComplexF64[]
    Bj = ComplexF64[]
    for j in 1:Np
        a =  2/(2*Np+1)*(sin(j*pi/(2*Np+1)))^2
        b = (cos(j*pi/(2*Np+1)))^2
        A = exp(-imag_conv*θ_p/2) * a / (1 + b*(exp(-imag_conv * θ_p) - 1))^2
        B = exp(-imag_conv*θ_p) * b / (1 + b*(exp(-imag_conv * θ_p) - 1))

        push!(aj, a)
        push!(bj, b)
        push!(Aj, A)
        push!(Bj, B)
    end
    return MtE_OSRC_op(wavenumber, Np, θ_p, curvature, aj, bj, Aj, Bj)
end

function scalartype(op::MtE_OSRC_op)
    T = scalartype(op.wavenumber)
    return Complex{T}
end

function get_RNp(z, OSRC)
    R_Np = 1 + sum(OSRC.aj[j]*z/(1+OSRC.bj[j]*z) for j in 1:OSRC.Np)
    return R_Np
end

function get_C0(OSRC)
    C0 = exp(imag_conv*OSRC.θ_p/2) * get_RNp((exp(-imag_conv*OSRC.θ_p)-1), OSRC)
    return C0
end

function get_R0(OSRC)
    C0 = get_C0(OSRC)
    R0 = C0 + sum(OSRC.Aj[j]/OSRC.Bj[j] for j in 1:OSRC.Np)
    return R0
end

function rotated_pade(z, OSRC)
    return get_R0(OSRC)*I - sum(OSRC.Aj[j]*I/(OSRC.Bj[j]*(I + OSRC.Bj[j]*z)) for j in 1:OSRC.Np)
end

function rotated_implicit_pade(get_Π, OSRC)
    return get_R0(OSRC)*I - sum(OSRC.Aj[j]*I/(OSRC.Bj[j]*(get_Π(j, OSRC))) for j in 1:OSRC.Np)
end


# # Projection and embedding of LinearMaps
# Construct LinearMaps which efficiently slice out a relevant submatrix of a LinearMap (without constructing the actual LinearMap).

struct SlicedLinearMap
    A::LinearMap
    P::LinearMap
    Q::LinearMap
    rows::UnitRange{Int}
    cols::UnitRange{Int}
end

function SlicedLinearMap(A::LinearMap, rows::UnitRange{Int}, cols::UnitRange{Int})
    # selected rows and columns inside the matrix
    n_rows = length(rows)
    n_cols = length(cols)
    n, m = size(A)

    # functions used for the construction of the matrices
    function slice_row(vector::AbstractVector)
        return vector[rows]
    end

    function slice_column(vector::AbstractVector)
        y = zeros(ComplexF64, n)
        y[cols] = vector
        return y
    end

    P = LinearMap(slice_row, n_rows, n; ismutating=false)
    Q =  LinearMap(slice_column, n, n_cols; ismutating=false)
    return SlicedLinearMap(A, P, Q, rows, cols)
end

# The square root operator is regularized by adding a small imaginary component ``\epsilon`` to the wavenumber: ``k_{\epsilon} = k + i \epsilon``.
function MtE_damping(op)
    return 0.39*op.wavenumber^(1/3)*op.curvature^(2/3)
end

function assemble(op::MtE_OSRC_op,X::Space,Y::Space; quadstrat=defaultquadstrat)
    R_0 = get_R0(op)

    ϵ = MtE_damping(op)
    κ = op.wavenumber
    κ_ϵ = κ + imag_conv*ϵ

    #create auxilary basis functions
    L0_int = BEAST.lagrangec0d1(X.geo)
    grad_L0_int = BEAST.gradient(L0_int)
    # Define the relevant function spaces
    curl_X = BEAST.curl(X)
    curl_Y = BEAST.curl(Y)

    N_L0 = numfunctions(L0_int)
    N_X = numfunctions(X)
    N_Y = numfunctions(Y)

    # Assemble the submatrices of the blockmatrix of the system
    Id = BEAST.Identity();
    G = assemble(Id, X, Y; quadstrat)
    N_ϵ = (1/κ_ϵ)^2 * assemble(Id, curl_X, curl_Y; quadstrat)
    K_ϵ = κ_ϵ^2 * assemble(Id, L0_int, L0_int; quadstrat)
    L = assemble(Id, X, grad_L0_int; quadstrat)
    L_transpose = assemble(Id, grad_L0_int, Y; quadstrat)

    # construct the sparse system matrix and invert
    function create_j_phi_matrix17(j)
        B_j = op.Bj[j]
        # blockmatrix of sparse matrices
        AXY = [G-B_j*N_ϵ       B_j*L
                L_transpose     K_ϵ]
        SXY = BEAST.lu(AXY)
        Sliced_SXY = SlicedLinearMap(SXY, 1:N_X, 1:N_Y)
        P = Sliced_SXY.P
        Q = Sliced_SXY.Q
        SXY_sliced = P*SXY*Q
        return SXY_sliced
    end

    sum_Π_inv_matrix = sum(op.Aj[j]/op.Bj[j] * create_j_phi_matrix17(j) for j in 1:op.Np)
    G_N_ϵ_inv = BEAST.lu(G - N_ϵ)

    MtE_map = - (G_N_ϵ_inv * R_0 - G_N_ϵ_inv * G * sum_Π_inv_matrix)
    return MtE_map
end

## The Electric-to-Magnetic (EtM) OSRC operator is also implemented.
# More details for this OSRC operator are given in the paper <a href="https://www.researchgate.net/publication/261636307_Approximate_local_magnetic-to-electric_surface_operators_for_time-harmonic_Maxwell's_equations" target='new'> Approximate local magnetic-to-electric surface operators for time-harmonic Maxwell’s equations (C. Geuzaine et al. (2014))</a>.

# get Polynomial representation of numerator and denumerator Padé series, to obtain an inversion of the fraction.

function prod_poles(op::MtE_OSRC_op, range)
    return prod(Polynomials.Polynomial([1, op.Bj[j]]) for j in range)
end

function A_j_polynomial(op::MtE_OSRC_op, j)
    range = filter(x -> x != j, 1:op.Np)
    return  Polynomials.Polynomial([0, op.Aj[j]]) * prod_poles(op, range)
end

function get_numerator(op::MtE_OSRC_op)
    return BEAST.get_C0(op) *  prod_poles(op, 1:op.Np) + sum(A_j_polynomial(op, j) for j in 1:op.Np)
end

function get_denumerator(op::MtE_OSRC_op)
    return prod_poles(op, 1:op.Np)
end

struct EtM_OSRC_op
    wavenumber::Float64
    Np::Int
    θ_p::Float64
    curvature::Float64
    r0::ComplexF64
    rj::Vector{ComplexF64}
    qj::Vector{ComplexF64}
end

# construct inverse pade rational
function EtM_OSRC_op(wavenumber::Float64, Np::Int, θ_p::Float64, curvature::Float64)
    # construct MtE operator with regular rotational branch cut Padé approximation
    MtE_OSRC_op = BEAST.MtE_OSRC_op(wavenumber, Np, θ_p, curvature)
    pade_inv_num = get_denumerator(MtE_OSRC_op)
    pade_inv_denum = get_numerator(MtE_OSRC_op)
    qj = Polynomials.roots(pade_inv_denum)

    # polynmial division to get constant term r0
    r0, remainder = Polynomials.divrem(pade_inv_num, pade_inv_denum)  
    P_deriv = Polynomials.derivative(pade_inv_denum)
    residues = [remainder(p)/P_deriv(p) for p in qj] # compute residues=poles at every point
    rj = residues
    return EtM_OSRC_op(wavenumber, Np, θ_p, curvature, r0, rj, qj)
end

function rational_pade_inv(z, op::EtM_OSRC_op)
    s = op.r0
    for (r, q) in zip(op.rj, op.qj)
        s += r / (z - q)
    end
    return s
end


function  assemble(op::EtM_OSRC_op,X::Space,Y::Space; quadstrat=defaultquadstrat)
    ϵ = BEAST.MtE_damping(op)
    κ = op.wavenumber
    κ_ϵ = κ + imag_conv*ϵ

    # assemble M2 matrix

    curl_X = BEAST.curl(X);
    curl_Y = BEAST.curl(Y);
    P1 = BEAST.lagrangec0d1(X.geo);
    P1_grad = BEAST.gradient(P1);

    N_X = numfunctions(X)
    N_Y = numfunctions(Y)

    Id = BEAST.Identity();
    A = assemble(Id, X, Y; quadstrat);
    N = (1/κ_ϵ)^2 * assemble(Id, curl_X, curl_Y; quadstrat);
    K = κ_ϵ^2 * assemble(Id, P1, P1; quadstrat);
    L = assemble(Id, X, P1_grad; quadstrat);
    L_transpose = assemble(Id, P1_grad, Y; quadstrat)

    function get_M2(j)
        Bq = - (op.qj[j]*A + N);
        M2 = [Bq L
              L_transpose K];
        M2_inv = BEAST.lu(M2)
        SliceMap = BEAST.SlicedLinearMap(M2_inv, 1:N_X, 1:N_Y)
        P = SliceMap.P
        Q = SliceMap.Q
        M2_inv_sliced = P*M2_inv*Q
        return M2_inv_sliced
    end

    A_inv = BEAST.lu(A)
    EtM_map = A_inv * (op.r0*I + A * sum(op.rj[j]*get_M2(j) for j in 1:op.Np))*(A - N)
    return EtM_map
end

function scalartype(op::EtM_OSRC_op)
    T = scalartype(op.wavenumber)
    return Complex{T}
end