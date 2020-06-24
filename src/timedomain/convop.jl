struct ConvOp{T} <: AbstractArray{T,3}
    data::Array{T,3}
    k0::Array{Int,2}
    k1::Array{Int,2}
    tail::Array{T,2}
    length::Int
end

function Base.size(obj)
    return (size(obj.data)[2:3]...,obj.length)
end

function Base.getindex(obj::ConvOp, m::Int, n::Int, k::Int)
    if k < obj.k0[m,n]
        return zero(eltype(obj.data))
    end

    if k > obj.k1[m,n]
        return obj.tail[m,n]
    end

    return obj.data[k-obj.k0[m,n]+1,m,n]
end

#
# function convolve!(y, Z::ConvOp, x, j, k_start, k_stop=size(Z,3))
#     for m in axes(y,1)
#         for n in axes(x,1)
#             k0 = Z.k0[m,n]
#             k1 = Z.k1[m,n]
#             for k in max(k0,k_start):min(k1,k_stop)
#                 p = k - k0 + 1
#                 j-k < 0 && continue
#                 y[m] += Z.data[p,m,n] * x[n,j-k+1]
#             end
#
#             X = sum(x[n,1:j-k1])
#             y[m] += Z.tail[m,n] * X
#         end
#     end
# end


function convolve!(y, Z::ConvOp, x, X, j, k_start, k_stop=size(Z,3))
    # @info "The corrrect convolve!"
    for n in axes(x,1)
        for m in axes(y,1)
            k0 = Z.k0[m,n]
            k1 = Z.k1[m,n]
            for k in max(k0,k_start):min(k1,k_stop)
                p = k - k0 + 1
                j-k < 0 && continue
                y[m] += Z.data[p,m,n] * x[n,j-k+1]
            end

            j-k1 > 0 || continue
            y[m] += Z.tail[m,n] * X[n,j-k1]
        end
    end
end


# function convolve(Z::ConvOp, x, j, k_start)
#     y = zeros(eltype(Z), size(Z)[1])
#     convolve!(y, Z, x, j, k_start, size(Z)[3])
#     return y
# end
function polyeig(Z::ConvOp)
    kmax = maximum(Z.k1)
    M = size(Z,1)
    Q = zeros(eltype(Z), 2M, 2M, kmax+1)
    for k in 1:kmax
        Q[1:M,1:M,k] .= Z[:,:,k]
    end
    Id = Matrix{eltype(Z)}(LinearAlgebra.I,M,M)
    Q[M+1:2M,1:M,1] .= Id
    Q[M+1:2M,M+1:2M,1] .= -Id
    Q[M+1:2M,M+1:2M,1] .= Id
    Q[1:M,M+1:2M,kmax+1] .= Z[:,:,kmax+1]
    return eigvals(companion(Q))
    # return Q
end