mutable struct LUFactorization{T,M<:AbstractMatrix{T},X} <: LinearMap{T}
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
    LUFactorization{T,M,X}(A, nothing, false, axs)
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
        L.F = LinearAlgebra.lu(L.A)
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


@testitem "LUFactorization Sparse/Dense" begin
    using LinearAlgebra, SparseArrays
    using CompScienceMeshes, BEAST

    h = 0.5
    M = meshsphere(1.0,h)
    X = BEAST.gwpdiv(M;order=1)
    
    G = assemble(BEAST.Identity(),X,X)
    dG = Matrix(G)

    b = rand(numfunctions(X))

    Gi = BEAST.lu(G)
    dGi = BEAST.lu(dG)

    #First time computing the factorization
    @time x1 = Gi*b
    @time x2 = dGi*b
    #Second time only substitutions
    @time x1 = Gi*b
    @time x2 = dGi*b

    @test typeof(Gi.F) == SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}
    @test typeof(dGi.F) ==  LinearAlgebra.LU{Float64, Matrix{Float64}, Vector{Int64}}

end