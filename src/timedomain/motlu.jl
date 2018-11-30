

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

function convolve(Z::SparseND.Banded3D,x,j,k0)
    T = promote_type(eltype(Z), eltype(x))
    M,N,K = size(Z)
    @assert M == size(x,1)
y = zeros(T,M)
    for m in 1:M
        for n in 1:N
            for k in k0:min(j,K)
                i = j - k + 1
                y[m] += Z[m,n,k] * x[n,i]
            end
        end
    end
    return y
end

"""
    marchonintime(W0,Z,B,I)

Solve by marching-on-in-time the causal convolution problem defined by `(W0,Z,B)`
up to timestep `I`. Here, `Z` is an array of order 3 that contains a discretisation
of a time translation invariant retarded potential operator. `W0` is the inverse of
the slice `Z[:,:,1]`.
"""
function marchonintime(W0,Z,B,I)
    T = eltype(Z)
    M,N,K = size(Z)
    @assert M == size(B,1)
    x = zeros(T,N,I)
    for i in 1:I
        b = B[:,i] - convolve(Z,x,i,2)
        x[:,i] = W0 * b
        (i % 10 == 0) && print(i, "[", I, "] - ")
    end
    return x
end
