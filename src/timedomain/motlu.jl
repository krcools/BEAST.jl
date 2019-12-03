





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
    cs = BlockArrays.cumulsizes(Z)
    bs = [blocksize(Z, (i,1))[1] for i in 1:nblocks(Z,1)]
    # @show bs
    T = eltype(eltype(Z))
    y = PseudoBlockVector{T}(undef,bs)
    fill!(y,0)
    for I in 1:nblocks(Z,1)
        # xI = view(x, cs[1][I] : cs[1][I+1]-1, :)
        for J in 1:nblocks(Z,2)
            xJ = view(x, cs[2][J] : cs[2][J+1]-1, :)
            isassigned(Z.blocks, I, J) || continue
            ZIJ = Z[Block(I,J)].banded
            # @show size(xJ) size(ZIJ)
            # @show size(y[Block(I)])
            y[Block(I)] .+= convolve(ZIJ, xJ, i, j_start)
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
