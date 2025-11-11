mutable struct LUFactorization{T,M<:Matrix{T},X} <: LinearMap{T}
    A::M
    F::Any
    factorized::Bool
    axes::X
end

function LUFactorization(A)
    M = typeof(A) 
    T = eltype(A)
    axs = reverse(axes(A))
    X = typeof(axs)
    B = Matrix(A)
    LUFactorization{T,Matrix{T},X}(B, nothing, false, axs)
end

# function LUFactorization(A::SparseMatrixCSC) 
#     T= eltype(A)
#     fac = LinearAlgebra.lu(A)
#     LUFactorization{T,typeof(fac)}(fac, false)
# end
    

Base.axes(A::LUFactorization) = A.axes
Base.size(A::LUFactorization) = reverse(size(A.A))

LinearAlgebra.adjoint(A::LUFactorization) = LUFactorization(adjoint(A.A))
LinearAlgebra.transpose(A::LUFactorization) = LUFactorization(transpose(A.A))

function lu(A)
    LUFactorization(A)
end

function Base.:*(A::LUFactorization, b::AbstractVector)
    T = promote_type(eltype(A), eltype(b))
    y = similar(Array{T}, axes(A,2))
    mul!(y, A, b)
    return y
end

function LinearAlgebra.mul!(y::AbstractVector, L::LUFactorization, b::AbstractVector)
    fill!(y,false)
    if L.factorized == false
        L.F = LinearAlgebra.lu!(L.A)
        L.factorized = true
    end
    y[:] = L.F \ Vector(b)
end


@testitem "LUFactorization" begin
    using LinearAlgebra

    using CompScienceMeshes, BEAST

    h = 0.5
    M = meshrectangle(1.0,1.0,h)
    Y = BEAST.gwpdiv(M;order=1)

    b = rand(numfunctions(Y))

    G = assemble(BEAST.Identity(),Y,Y)

    iG = BEAST.lu(G)

    x=iG*b
    x2 = inv(Matrix(G))*b

    @test norm(x - x2) < 1e-12

end

@testitem "LUFactorization blocked" begin
    using LinearAlgebra
    using BEAST.BlockArrays

    n = 3
    J = Matrix{Float64}(I,n,n)
    Z = zeros(Float64,n,n)
    A = BlockedArray([Z J; -J Z], [n,n], [n,n])

    Ai = BEAST.lu(A)
    b = collect(1:2n)

    x = Ai * b
    @test blocksize(x) == (2,)
    @test blocksizes(x) == [(n,), (n,)]

    @test x[1:n] ≈ -b[n+1:end]
    @test x[n+1:end] ≈ b[1:n]
end