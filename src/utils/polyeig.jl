export companion

function companion(Z)

    T = eltype(Z)
    K = size(Z,3)
    @assert K > 1

    M, N = size(Z)[1:2]
    C = similar(Z, M*(K-1), N*(K-1))
    fill!(C,0)

    @assert M == N
    Id = Matrix{T}(I, M, N)
    for m in 2:K-1
        n = m-1
        C[(m-1)*M+1:m*M, (n-1)*N+1:n*N] = Id
    end

    W = -inv(Z[:,:,1])
    for n in 1:K-1
        m = 1
        C[(m-1)*M+1:m*M, (n-1)*N+1:n*N] = W*Z[:,:,n+1]
    end

    return C
end
