function assemble(op::Identity,
        testnfs::AbstractTimeBasisFunction,
        trialfns::AbstractTimeBasisFunction)
# function assemble(op,
#         testnfs::AbstractTimeBasisFunction,
#         trialfns::AbstractTimeBasisFunction)

    tbf = convolve(testnfs, trialfns)
    T = scalartype(tbf)
    z = zeros(T, numintervals(tbf)-1)
    Δt = timestep(tbf)
    for i in eachindex(z)
        p = tbf.polys[i]
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


function assemble(operator::TensorOperator, testfns, trialfns)

    test_spatial_basis  = spatialbasis(testfns)
    trial_spatial_basis = spatialbasis(trialfns)

    M = numfunctions(test_spatial_basis)
    N = numfunctions(trial_spatial_basis)

    tbf = convolve(
        temporalbasis(testfns),
        temporalbasis(trialfns))

    kmax = numintervals(tbf) - 1
    k0 = fill(1,    (M,N))
    k1 = fill(kmax, (M,N))

    T = promote_type(scalartype(operator), scalartype(testfns), scalartype(trialfns))
    data = zeros(T, M, N, kmax)
    Z = SparseND.Banded3D(k0, k1, data)

    store(v, m, n, k) = (Z[m,n,k] += v)
    assemble!(operator, testfns, trialfns, store)
    return Z

end

function assemble!(operator::TensorOperator, testfns, trialfns, store)

    space_operator = operator.spatial_factor
    time_operator  = operator.temporal_factor

    space_testfns = spatialbasis(testfns)
    space_trialfns = spatialbasis(trialfns)

    time_testfns = temporalbasis(testfns)
    time_trialfns = temporalbasis(trialfns)

    zt = assemble(time_operator, time_testfns, time_trialfns)
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
