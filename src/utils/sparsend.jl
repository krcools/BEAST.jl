module SparseND

struct Banded3D{T} <: AbstractArray{T,3}
    k0::Array{Int,2}
    data::Array{T,3}
    maxk0::Int
    maxk1::Int
    function Banded3D{T}(k0,data,maxk1) where T
        @assert size(k0) == size(data)[2:3]
        return new(k0,data,maximum(k0),maxk1)
    end
end

Banded3D(k0,data::Array{T},maxk1) where {T} = Banded3D{T}(k0,data,maxk1)

bandwidth(A::Banded3D) = size(A.data,1)

import Base: size, getindex, setindex!

# size(A::Banded3D) = size(A.data)
size(A::Banded3D) = tuple(size(A.k0)..., A.maxk1)

function getindex(A::Banded3D, m::Int, n::Int, k::Int)
    k0 = A.k0[m,n]
    k0 == 0 && return zero(eltype(A))
    l = k - k0 + 1
    l < 1 && return zero(eltype(A))
    l > bandwidth(A) && return zero(eltype(A))
    A.data[l,m,n]
end

function setindex!(A::Banded3D, v, m, n, k)
    k0 = A.k0[m,n]
    @assert k0 != 0
    @assert A.k0[m,n] <= k <= A.k0[m,n] + bandwidth(A) - 1
    A.data[k-A.k0[m,n]+1,m,n] = v
end

# """
#     convolve{T}(A::Banded3D, x::Array{T,2})
#
# Compute the *space-time* convolution of A and x.
#
# Returns an array `y` of size `(size(A,1),size(x,2))` such that
# `y[m,i] = sum([A[]])``
# """
# function convolve(A::Banded3D, x::Array{T,2}) where T
#     V = promote_type(eltype(A),eltype(x))
#     y = zeros(V,size(A,1),size(x,2))
#     bw = bandwidth(A)
#     for i in 1:size(A,1)
#         for l in 1 : size(y,2)
#             for j in 1:size(A,2)
#                 k0 = A.k0[i,j]
#                 for k in k0 : k0 + bw - 1
#                     l-k+1 < 1 && continue
#                     y[i,l] += A.data[k,i,j] * x[j,l-k+1]
#                 end
#             end
#         end
#     end
#     return y
# end

function Base.:+(A::Banded3D{T}, B::Banded3D{T}) where {T}

    M = size(A,1)
    N = size(A,2)

    @assert M == size(B,1)
    @assert N == size(B,2)

    Abw = size(A.data,1)
    Bbw = size(B.data,1)

    # keep track of empty columns
    Az = findall(A.k0 .== 0)
    Bz = findall(B.k0 .== 0)

    K0 = min.(A.k0, B.k0)
    K0[Az] = B.k0[Az]
    K0[Bz] = A.k0[Bz]

    AK1 = A.k0 .+ (Abw - 1); AK1[Az] .= -1
    BK1 = B.k0 .+ (Bbw - 1); BK1[Bz] .= -1
    K1 = max.(AK1, BK1)

    bw = maximum(K1 - K0 .+ 1)
    @assert bw > 0
    @assert bw >= Abw
    @assert bw >= Bbw
    data = zeros(T, bw, M, N)
    for m in axes(A.data,2)
        for n in axes(A.data,3)
            k0 = K0[m,n]
            @assert k0 != 0

            Ak0 = A.k0[m,n]
            if Ak0 != 0
                for Al in 1:Abw
                    k = Ak0 + Al - 1
                    l = k - k0 + 1
                    data[l,m,n] += A.data[Al,m,n]
                end
            end

            Bk0 = B.k0[m,n]
            if Bk0 != 0
                for Bl in 1:Bbw
                    k = Bk0 + Bl - 1
                    l = k - k0 + 1
                    l < 1 && @show m, n, k0, k, Bk0, Bl
                    d = data[l,m,n]
                    data[l,m,n] = d + B.data[Bl,m,n]
                end
            end
        end
    end

    Banded3D(K0, data, max(A.maxk1, B.maxk1))
end

end # module
