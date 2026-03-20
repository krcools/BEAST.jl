using SparseArrays, LinearAlgebra, BEAST

struct MtE_OSRC_op <: Operator
    wavenumber::ComplexF64
    Np::Int
    θ_p::Float64
    curvature::Float64
    aj::Vector{ComplexF64}
    bj::Vector{ComplexF64}
    Aj::Vector{ComplexF64}
    Bj::Vector{ComplexF64}
end

imag_conv = -im    
function MtE_OSRC_op(wavenumber::ComplexF64, Np::Int, θ_p::Float64, curvature::Float64)
    # get the real and rotated pade coefficients
    aj = ComplexF64[]
    bj = ComplexF64[]
    Aj = ComplexF64[]
    Bj = ComplexF64[]
    for j in 1:Np
        a =  2/(2*Np+1)*(sin(j*pi/(2*Np+1)))^2
        b = (cos(j*pi/(2*Np+1)))^2
        A = exp(-imag_conv*θ_p/2) * a / (1 + b*(exp(-imag_conv * θ_p) - 1))^2                               # check convention sign of imag_conv
        B = exp(-imag_conv*θ_p) * b / (1 + b*(exp(-imag_conv * θ_p) - 1))                                   # similar here

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

# # Projection and embedding of LinearMaps
# Construct LinearMaps which slice out a relevant submatrix of a LinearMap (without constructing the actual LinearMap).
# This is done very efficiently (see Benchmark)

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

# Next, the MtE operator is assembled as a linear map. This is implemented according to the discretization in the paper.

function assemble(op::MtE_OSRC_op,X::Space,Y::Space; quadstrat=defaultquadstrat)
    R_0 = get_R0(op)

    κ = op.wavenumber
    κ_ϵ = κ

    #create auxilary basis functions
    L0_int = BEAST.lagrangec0d1(X.geo)
    grad_L0_int = BEAST.gradient(L0_int)
    # Define the relevant function spaces
    curl_X = BEAST.curl(X)
    curl_Y = BEAST.curl(Y)

    N_X = numfunctions(X)
    N_Y = numfunctions(Y)

    # Assemble the submatrices of the blockmatrix of the system
    Id = BEAST.Identity()
    G = assemble(Id, X, Y)
    N_ϵ = (1/κ_ϵ)^2 * assemble(Id, curl_X, curl_Y)
    K_ϵ = κ_ϵ^2 * assemble(Id, L0_int, L0_int)
    L = assemble(Id, X, grad_L0_int)
    L_transpose = L'

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

    MtE_map = - G * (G_N_ϵ_inv * R_0 * (I - 1/R_0 * G * sum_Π_inv_matrix))
    return MtE_map
end

### Simplified version

struct MtE_OSRC_simp_op <: Operator
    wavenumber::ComplexF64
end

function MtE_OSRC_simp_op(wavenumber::ComplexF64)
    return MtE_OSRC_simp_op(wavenumber)
end

function scalartype(op::MtE_OSRC_simp_op)
    T = scalartype(op.wavenumber)
    return Complex{T}
end

function assemble(op::MtE_OSRC_simp_op, X::Space, Y::Space; quadstrat=defaultquadstrat)
    # Define the relevant function spaces
    curl_X = BEAST.curl(X)
    curl_Y = BEAST.curl(Y)

    # Assemble the submatrices of the blockmatrix of the system
    Id = BEAST.Identity()
    G = assemble(Id, X, Y)
    N_ϵ = (1/op.wavenumber)^2 * assemble(Id, curl_X, curl_Y)

    G_N_ϵ_inv = BEAST.lu(G - N_ϵ)
    return - G * G_N_ϵ_inv
end
