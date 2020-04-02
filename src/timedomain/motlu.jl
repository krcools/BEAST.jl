





"""
    marchonintime(W0,Z,B,I)

Solve by marching-on-in-time the causal convolution problem defined by `(W0,Z,B)`
up to timestep `I`. Here, `Z` is an array of order 3 that contains a discretisation
of a time translation invariant retarded potential operator. `W0` is the inverse of
the slice `Z[:,:,1]`.
"""
function marchonintime(W0,Z,B,I)
    T = eltype(W0)
    M,N = size(Z)
    @assert M == size(B,1)
    x = zeros(T,N,I)
    for i in 1:I
        b = B[:,i] - convolve(Z,x,i,2)
        x[:,i] += W0 * b
        (i % 10 == 0) && print(i, "[", I, "] - ")
    end
    return x
end

using BlockArrays

function convolve(Z::BlockArray, x, i, j_start)
    # ax1 = axes(Z,1)
    ax2 = axes(Z,2)
    T = eltype(eltype(Z))
    y = PseudoBlockVector{T}(undef,blocklengths(axes(Z,1)))
    fill!(y,0)
    for I in blockaxes(Z,1)
        for J in blockaxes(Z,2)
            xJ = view(x, ax2[J], :)
            try
                ZIJ = Z[I,J].banded
                y[I] .+= convolve(ZIJ, xJ, i, j_start)
            catch
                @info "Skipping unassigned block."
                continue
            end
        end
    end
    return y
end

function marchonintime(W0,Z::BlockArray,B,I)
    T = eltype(W0)
    M,N = size(W0)
    @assert M == size(B,1)
    x = zeros(T,N,I)
    for i in 1:I
        R = [ B[j][i] for j in 1:N ]
        S = convolve(Z,x,i,2)
        # @show size(R)
        # @show size(S)
        # b = R - convolve(Z,x,i,2)
        b = R - S
        x[:,i] += W0 * b
        (i % 10 == 0) && print(i, "[", I, "] - ")
    end
    return x
end
