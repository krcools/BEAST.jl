module SparseND

type Banded3D{T} <: AbstractArray{T,3}
    k0::Array{Int,2}
    k1::Array{Int,2}
    data::Array{T,3}
    function Banded3D{T}(k0,k1,data) where T
        @assert size(k0) == size(k1)
        @assert size(k0) == size(data,1,2)
        return new(k0,k1,data)
    end
end

Banded3D{T}(k0,k1,data::Array{T}) = Banded3D{T}(k0,k1,data)

import Base: size, getindex, setindex!, linearindexing

size(A::Banded3D) = size(A.data)
#linearindexing{T<:Banded3D}(::Type{T}) = Base.LinearSlow()
function getindex(A::Banded3D, m::Int, n::Int, k::Int)
    (A.k0[m,n] <= k <= A.k1[m,n]) || return zero(eltype(A.data))
    A.data[m,n,k-A.k0[m,n]+1]
end
function setindex!(A::Banded3D, v, m, n, k)
    @assert A.k0[m,n] <= k <= A.k1[m,n]
    A.data[m,n,k-A.k0[m,n]+1] = v
end

"""
    convolve{T}(A::Banded3D, x::Array{T,2})

Compute the *space-time* convolution of A and x.

Returns an array `y` of size `(size(A,1),size(x,2))` such that
`y[m,i] = sum([A[]])``
"""
function convolve{T}(A::Banded3D, x::Array{T,2})
    V = promote_type(eltype(A),eltype(x))
    y = zeros(V,size(A,1),size(x,2))
    for i in 1:size(A,1)
        for j in 1:size(A,2)
            for k in A.k0[i,j] : A.k1[i,j]
                for l in 1 : size(y,2)
                    l-k+1 < 1 && continue
                    y[i,l] += A[i,j,k] * x[j,l-k+1]
                end
            end
        end
    end
    return y
end

end # module
