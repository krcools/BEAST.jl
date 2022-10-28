# function timeslice(Z::ConvolutionOperators.ConvOp, k)
#     T = eltype(Z.data)
#     Zk = zeros(T, size(Z)[1:2])
#     for n in axes(Z,2)
#         for m in axes(Z,1)
#             Zk[m,n] = Z[m,n,k]
#     end end
#     return Zk
# end






# function marchonintime(W0,Z,B,I)
#     T = eltype(W0)
#     M,N = size(Z)
#     @assert M == size(B,1)
#     x = zeros(T,N,I)
#     for i in 1:I
#         b = B[:,i] - convolve(Z,x,i,2)
#         x[:,i] += W0 * b
#         (i % 10 == 0) && print(i, "[", I, "] - ")
#     end
#     return x
# end

# function marchonintime(iZ0, Z::ConvolutionOperators.AbstractConvOp, B, Nt)

#     T = eltype(iZ0)
#     Ns = size(Z,1)
#     x = zeros(T,Ns,Nt)
#     csx = zeros(T,Ns,Nt)
#     y = zeros(T,Ns)

#     todo, done, pct = Nt, 0, 0
#     for i in 1:Nt
#         fill!(y,0)
#         ConvolutionOperators.convolve!(y, Z, x, csx, i, 2, Nt)
#         y .*= -1
#         y .+= B[:,i]
#         # @show norm(B[:,i])

#         x[:,i] .+= iZ0 * y
#         if i > 1
#             csx[:,i] .= csx[:,i-1] .+ x[:,i]
#         else
#             csx[:,i] .= x[:,i]
#         end

#         done += 1
#         new_pct = round(Int, done / todo * 100)
#         new_pct > pct+9 && (println("[$new_pct]"); pct=new_pct)
#     end
#     x
# end

# function convolve(Z::BlockArray, x, i, j_start)
#     # ax1 = axes(Z,1)
#     ax2 = axes(Z,2)
#     T = eltype(eltype(Z))
#     y = PseudoBlockVector{T}(undef,blocklengths(axes(Z,1)))
#     fill!(y,0)
#     for I in blockaxes(Z,1)
#         for J in blockaxes(Z,2)
#             xJ = view(x, ax2[J], :)
#             try
#                 ZIJ = Z[I,J].banded
#                 y[I] .+= convolve(ZIJ, xJ, i, j_start)
#             catch
#                 @info "Skipping unassigned block."
#                 continue
#             end
#         end
#     end
#     return y
# end

# function convolve!(y,Z::BlockArray, x, csx, i, j_start, j_stop)
#     ax1 = axes(Z,1)
#     ax2 = axes(Z,2)
#     T = eltype(eltype(Z))
#     fill!(y,0)
#     for I in blockaxes(Z,1)
#         for J in blockaxes(Z,2)
#             xJ = view(x, ax2[J], :)
#             csxJ = view(csx, ax2[J], :)
#             yI = view(y, ax1[I])
#             ConvolutionOperators.convolve!(yI, Z[I,J], xJ, csxJ, i, j_start, j_stop)
#         end
#     end
#     return y
# end

"""
    marchonintime(W0,Z,B,I)

Solve by marching-on-in-time the causal convolution problem defined by `(W0,Z,B)`
up to timestep `I`. Here, `Z` is an array of order 3 that contains a discretisation
of a time translation invariant retarded potential operator. `W0` is the inverse of
the slice `Z[:,:,1]`.
"""
function marchonintime(W0,Z,B,I)

    T = eltype(W0)
    M,N = size(W0)
    @assert M == size(B,1)

    x = zeros(T,N,I)
    y = zeros(T,N)
    csx = zeros(T,N,I)

    for i in 1:I
        # R = [ B[j][i] for j in 1:N ]
        R = B[:,i]
        # @show norm(R)
        k_start = 2
        k_stop = I

        fill!(y,0)
        ConvolutionOperators.convolve!(y,Z,x,csx,i,k_start,k_stop)
        b = R - y
        x[:,i] .+= W0 * b
        if i > 1
            csx[:,i] .= csx[:,i-1] .+ x[:,i]
        else
            csx[:,i] .= x[:,i]
        end

        (i % 10 == 0) && print(i, "[", I, "] - ")
    end
    return x
end
