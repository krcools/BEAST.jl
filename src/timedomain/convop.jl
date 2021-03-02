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


function convolve!(y, Z::ConvOp, x, X, j, k_start, k_stop=size(Z,3))
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


function polyeig(Z::ConvOp)
    kmax = maximum(Z.k1)
    M = size(Z,1)
    Q = zeros(eltype(Z), 2M, 2M, kmax+1)
    for k in 1:kmax
        Q[1:M,1:M,k] .= Z[:,:,k]
    end
    Id = Matrix{eltype(Z)}(LinearAlgebra.I,M,M)
    Q[M+1:2M,1:M,1] .= -Id
    Q[M+1:2M,M+1:2M,1] .= Id
    Q[M+1:2M,M+1:2M,2] .= -Id
    Q[1:M,M+1:2M,kmax+1] .= Z[:,:,kmax+1]
    return eigvals(companion(Q)), Q
    # return Q
end


function polyeig(Z)
    return eigvals(companion(Z))
end


function Base.:+(a::ConvOp, b::ConvOp)

    @assert size(a.data) == size(b.data)
    M,N = size(a.data)
    T = promote_type(eltype(a.data), eltype(b.data))

    k0 = min.(a.k0, b.k0)
    k1 = max.(a.k1, b.k1)

    bandwidth = maximum(k1 - k0) + 1
    data = zeros(T, bandwidth, M, N)

    k1max = maximum(k1)
    tail = zeros(T,M,N)
    for m in 1:M, n in 1:N
        tail[m,n] = a[m,n,k1max+1] + b[m,n,k1max+1]
    end

    for m in 1:M
        for n in 1:N
            for k in k0[m,n]:k1[m,n]
                data[k,m,n] = a[m,n,k] + b[m,n,k]
            end
        end
    end

    lgt = max(a.length, b.length)
    return ConvOp(data, tail, k0, k1, lgt)
end


