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


convolve(x::MatrixConvolution, y::Matrix, i, j) = convolve(x.arr, y, i, j)
