struct LUFactorization{T,M} <: LinearMap{T}
    fac::M
end

Base.axes(A::LUFactorization) = reverse(axes(A.fac))
Base.size(A::LUFactorization) = reverse(size(A.fac))

LinearAlgebra.adjoint(A::LUFactorization) = LUFactorization(adjoint(A.fac))
LinearAlgebra.transpose(A::LUFactorization) = LUFactorization(transpose(A.fac))

function LUFactorization(A::SparseMatrixCSC) 
    T= eltype(A)
    fac = LinearAlgebra.lu(A)
    LUFactorization{T,typeof(fac)}(fac)
end

lu(A::SparseMatrixCSC) = LinearMap(LUFactorization(A))

function LUFactorization(fac::M) where M 
    T= eltype(fac)
    LUFactorization{T,M}(fac)
end

function Base.:*(A::LUFactorization, b::AbstractVector)
    T = promote_type(eltype(A), eltype(b))
    y = BlockedVector{T}(undef, (axes(A,2),))
    mul!(y, A, b)
    return y
end

function LinearAlgebra.mul!(y::AbstractVector, lu::LUFactorization, b::AbstractVector)
    fill!(y,zero(eltype(y)))
    y[:] = lu.fac \ Vector(b)
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