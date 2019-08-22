

function assemble(op::Identity,
        testnfs::AbstractTimeBasisFunction,
        trialfns::AbstractTimeBasisFunction)

    tbf = convolve(testnfs, trialfns)
    has_zero_tail = all(tbf.polys[end].data .== 0)
    @show has_zero_tail

    T = scalartype(tbf)
    if has_zero_tail
        z = zeros(T, numintervals(tbf)-1)
    else
        z = zeros(T, numfunctions(tbf))
    end

    Δt = timestep(tbf)
    #for i in eachindex(z)
    for i in 1:numintervals(tbf)-1
        p = tbf.polys[i]
        t = (i-1)*Δt
        z[i] = evaluate(p,t)
    end

    for i in numintervals(tbf):length(z)
        p = tbf.polys[end]
        t = (i-1)*Δt
        z[i] = evaluate(p,t)
    end

    return z
end

mutable struct TensorOperator <: Operator
    spatial_factor
    temporal_factor
end

⊗(A::AbstractOperator, B::AbstractOperator) = TensorOperator(A, B)
function scalartype(A::TensorOperator)
    promote_type(
        scalartype(A.spatial_factor),
        scalartype(A.temporal_factor))
end


# function assemble(operator::TensorOperator, testfns, trialfns)
#
#     test_spatial_basis  = spatialbasis(testfns)
#     trial_spatial_basis = spatialbasis(trialfns)
#
#     M = numfunctions(test_spatial_basis)
#     N = numfunctions(trial_spatial_basis)
#
#     tbf = convolve(
#         temporalbasis(testfns),
#         temporalbasis(trialfns))
#
#     kmax = numintervals(tbf) - 1
#     k0 = fill(1,    (M,N))
#     k1 = fill(kmax, (M,N))
#
#     T = promote_type(scalartype(operator), scalartype(testfns), scalartype(trialfns))
#     data = zeros(T, M, N, kmax)
#     Z = SparseND.Banded3D(k0, k1, data)
#
#     store(v, m, n, k) = (Z[m,n,k] += v)
#     assemble!(operator, testfns, trialfns, store)
#     return Z
#
# end

function allocatestorage(op::TensorOperator, test_functions, trial_functions,
    ::Type{Val{:bandedstorage}}, ::Type{LongDelays{:ignore}},)

    M = numfunctions(spatialbasis(test_functions))
    N = numfunctions(spatialbasis(trial_functions))

    time_basis_function = BEAST.convolve(
        temporalbasis(test_functions),
        temporalbasis(trial_functions))

    # has_zero_tail = all(time_basis_function.polys[end].data .== 0)
    # @show has_zero_tail

    # if has_zero_tail
        # Numintervals includes the semi-infinite interval stretching to +Inf
        # K = numintervals(time_basis_function)-1
    # else
    #     speedoflight = 1.0
    #     @warn "Assuming speed of light to be equal to 1!"
    #     Δt = timestep(time_basis_function)
    #     ct, hs = boundingbox(geometry(spatialbasis(trial_functions)).vertices)
    #     diam = 2 * sqrt(3) * hs
    #     K = ceil(Int, (numintervals(time_basis_function)-1) + diam/speedoflight/Δt)+1
    # end
    # @assert K > 0

    space_operator = op.spatial_factor
    A = assemble(space_operator, spatialbasis(test_functions), spatialbasis(trial_functions))

    K0 = Int.(A .!= 0)
    # K0 = ones(Int,M,N)
    bandwidth = numintervals(time_basis_function) - 1
    data = zeros(scalartype(op), bandwidth, M, N)
    maxk1 = bandwidth
    Z = SparseND.Banded3D(K0, data, maxk1)
    return Z, (v,m,n,k)->(Z[m,n,k] += v)
end


function allocatestorage(op::TensorOperator, test_functions, trial_functions)

    M = numfunctions(spatialbasis(test_functions))
    N = numfunctions(spatialbasis(trial_functions))

    time_basis_function = BEAST.convolve(
        temporalbasis(test_functions),
        temporalbasis(trial_functions))

    tbf = time_basis_function
    has_zero_tail = all(tbf.polys[end].data .== 0)
    @show has_zero_tail

    if has_zero_tail
        K = numintervals(time_basis_function)-1
    else
        speedoflight = 1.0
        @warn "Assuming speed of light to be equal to 1!"
        Δt = timestep(tbf)
        ct, hs = boundingbox(geometry(spatialbasis(trial_functions)).vertices)
        diam = 2 * sqrt(3) * hs
        K = ceil(Int, (numintervals(tbf)-1) + diam/speedoflight/Δt)+1
    end
    @assert K > 0

    Z = zeros(M, N, K)
    return MatrixConvolution(Z), (v,m,n,k)->(Z[m,n,k] += v)
end

function assemble!(operator::TensorOperator, testfns, trialfns, store)

    space_operator = operator.spatial_factor
    time_operator  = operator.temporal_factor

    space_testfns = spatialbasis(testfns)
    space_trialfns = spatialbasis(trialfns)

    time_testfns = temporalbasis(testfns)
    time_trialfns = temporalbasis(trialfns)

    zt = assemble(time_operator, time_testfns, time_trialfns)

    # Truncate to a reasonable size in case the tbf has a semi-infinite support
    tbf = convolve(time_testfns, time_trialfns)
    has_zero_tail = all(tbf.polys[end].data .== 0)
    if !has_zero_tail
        speedoflight = 1.0
        @warn "Assuming speed of light to be equal to 1!"
        Δt = timestep(tbf)
        ct, hs = boundingbox(geometry(space_trialfns).vertices)
        diam = 2 * sqrt(3) * hs
        kmax = ceil(Int, (numintervals(tbf)-1) + diam/speedoflight/Δt)+1
        zt = zt[1:kmax]
    end


    function store1(v,m,n)
        for (k,w) in enumerate(zt)
            store(w*v,m,n,k)
        end
    end
    assemble!(space_operator, space_testfns, space_trialfns, store1)

end


mutable struct TemporalDifferentiation
    operator
end

function assemble!(operator::TemporalDifferentiation, testfns, trialfns, store)

    trial_time_fns  = temporalbasis(trialfns)
    trial_space_fns = spatialbasis(trialfns)

    trialfns = SpaceTimeBasis(
        trial_space_fns,
        derive(trial_time_fns)
    )

    assemble!(operator.operator, testfns, trialfns, store)

end
