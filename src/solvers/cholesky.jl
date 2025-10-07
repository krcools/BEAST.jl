struct CholeskyFactorization{T,M} <: LinearMap{T}
    fac::M
end

Base.axes(A::CholeskyFactorization) = reverse(axes(A.fac))
Base.size(A::CholeskyFactorization) = reverse(size(A.fac))

LinearAlgebra.adjoint(A::CholeskyFactorization) = CholeskyFactorization(adjoint(A.fac))
LinearAlgebra.transpose(A::CholeskyFactorization) = CholeskyFactorization(transpose(A.fac))

function CholeskyFactorization(A::SparseMatrixCSC) 
    T= eltype(A)
    fac = LinearAlgebra.cholesky(Hermitian(A))
    CholeskyFactorization{T,typeof(fac)}(fac)
end

cholesky(A::SparseMatrixCSC) = LinearMap(CholeskyFactorization(A))

function CholeskyFactorization(fac::M) where M 
    T= eltype(fac)
    CholeskyFactorization{T,M}(fac)
end

function Base.:*(A::CholeskyFactorization, b::AbstractVector)
    T = promote_type(eltype(A), eltype(b))
    y = PseudoBlockVector{T}(undef, (axes(A,2),))
    mul!(y, A, b)
    return y
end

function LinearAlgebra.mul!(y::AbstractVector, chol::CholeskyFactorization, b::AbstractVector)
    fill!(y,zero(eltype(y)))
    y[:] = chol.fac \ Vector(b)
end