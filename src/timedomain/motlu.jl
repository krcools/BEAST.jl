motsolve(eq) = td_solve(eq)

function td_solve(eq)

    V = eq.trial_space_dict[1]

    A = assemble(eq.equation.lhs, eq.test_space_dict, eq.trial_space_dict)
    T = eltype(A)
    S = zeros(T, size(A)[1:2])
    ConvolutionOperators.timeslice!(S, A, 1)

    iS = inv(S)
    b = assemble(eq.equation.rhs, eq.test_space_dict)

    nt = numfunctions(temporalbasis(V))
    marchonintime(iS, A, b, nt)
end

"""
    marchonintime(W0,Z,B,I; convhist=false)

Solve by marching-on-in-time the causal convolution problem defined by `(W0,Z,B)`
up to timestep `I`. Here, `Z` is an array of order 3 that contains a discretisation
of a time translation invariant retarded potential operator. `W0` is the inverse of
the slice `Z[:,:,1]`.

Keyword arguments:
    - 'convhist': when true, return in addition to the space-time data for the
    solution also the vector of convergence histories as returned each time step
    by the supplied solver `W0`.
"""
function marchonintime(W0,Z,B,I; convhist=false)

    T = eltype(W0)
    M,N = size(W0)
    @assert M == size(B,1)

    x = zeros(T,N,I)
    y = zeros(T,N)
    csx = zeros(T,N,I)

    ch = []
    for i in 1:I
        R = B[:,i]
        k_start = 2
        k_stop = I

        fill!(y,0)
        ConvolutionOperators.convolve!(y,Z,x,csx,i,k_start,k_stop)
        b = R - y
        xi, chi = BEAST.solve(W0, b)
        x[:,i] .+= xi
        push!(ch, chi)
        if i > 1
            csx[:,i] .= csx[:,i-1] .+ x[:,i]
        else
            csx[:,i] .= x[:,i]
        end

        (i % 10 == 0) && print(i, "[", I, "] - ")
    end

    if convhist
        return x, ch
    else
        return x
    end
end
