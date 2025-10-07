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
    y = PseudoBlockVector{T}(undef, (axes(A,2),))
    mul!(y, A, b)
    return y
end

function LinearAlgebra.mul!(y::AbstractVector, lu::LUFactorization, b::AbstractVector)
    fill!(y,zero(eltype(y)))
    y[:] = lu.fac \ Vector(b)
end
