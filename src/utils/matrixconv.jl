struct MatrixConvolution{T} <: AbstractArray{T,3}
    arr::Array{T,3}
end

function Base.:+(x::MatrixConvolution, y::MatrixConvolution)

    A = x.arr
    B = y.arr

    T = promote_type(eltype(A), eltype(B))
    @assert axes(A)[1:2] == axes(B)[1:2]

    kmax = max(size(A,3), size(B,3))
    C = zeros(T,size(A)[1:2]...,kmax)

    for k in axes(A,3)
        C[:,:,k] += A[:,:,k]
    end

    for k in axes(B,3)
        C[:,:,k] += B[:,:,k]
    end

    return MatrixConvolution(C)
end


Base.size(x::MatrixConvolution) = size(x.arr)
Base.getindex(x::MatrixConvolution, i::Int) = x.arr[i]
Base.IndexStyle(x::MatrixConvolution) = IndexLinear()


function Base.:*(x::MatrixConvolution, y::Matrix)

    A = x.arr

    T = promote_type(eltype(A), eltype(y))
    @assert axes(A,2) == axes(y,1)

    C = zeros(T, size(A)...)
    for k in axes(C,3)
        C[:,:,k] = A[:,:,k]*y
    end

    return MatrixConvolution(C)
end

function Base.:*(y::Matrix, x::MatrixConvolution)

    A = x.arr

    T = promote_type(eltype(A), eltype(y))
    @assert axes(A,2) == axes(y,1)

    C = zeros(T, size(A)...)
    for k in axes(C,3)
        C[:,:,k] = y*A[:,:,k]
    end

    return MatrixConvolution(C)
end


Base.:*(a::Number, x::MatrixConvolution) = MatrixConvolution(a * x.arr)
Base.:*(x::MatrixConvolution, a::Number) = MatrixConvolution(x.arr * a)
Base.:/(x::MatrixConvolution, a::Number) = MatrixConvolution(x.arr / a)
Base.:-(x::MatrixConvolution) = (-1) * x

function convolve(Z::Array,x,j,k0)
    M,N,K = size(Z)
    y = similar(Z,M)
    fill!(y,0)
    for k ∈ k0 : min(j,K)
        i = j - k + 1
        y += Z[:,:,k] * x[:,i]
    end
    return y
end

convolve(x::MatrixConvolution, y::Matrix, i, j) = convolve(x.arr, y, i, j)


function Base.hvcat((M,N)::Tuple{Int,Int}, as::MatrixConvolution...)
    kmax = maximum(size(a,3) for a in as)

    @assert length(as) == M*N

    li = LinearIndices((1:N,1:M))
    for m in 1:M
        a = as[li[1,m]]
        M1 = size(a,1)
        for n in 2:N
            a = as[li[n,m]]
            @assert size(a,1) == M1
        end
    end

    for n in 1:N
        a = as[li[n,1]]
        N1 = size(a,2)
        for m in 2:M
            a = as[li[n,m]]
            @assert size(a,2) == N1
        end
    end

    Ms = [size(as[li[1,i]],1) for i in 1:M]
    Ns = [size(as[li[j,1]],2) for j in 1:N]

    cMs = pushfirst!(cumsum(Ms),0)
    cNs = pushfirst!(cumsum(Ns),0)
    T = promote_type(eltype.(as)...)
    data = zeros(T, last(cMs), last(cNs), kmax)

    @show size(data)
    @show eltype(data)

    for m in 1:M
        I = cMs[m]+1 : cMs[m+1]
        for n in 1:N
            J = cNs[n]+1 : cNs[n+1]
            a = as[li[n,m]]
            K = 1:size(a,3)
            data[I,J,K] .= a
        end
    end

    return MatrixConvolution(data)
end