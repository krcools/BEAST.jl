

function assemble(op::Identity,
        testnfs::AbstractTimeBasisFunction,
        trialfns::AbstractTimeBasisFunction)

    tbf = convolve(testnfs, trialfns)
    has_zero_tail = all(tbf.polys[end].data .== 0)
    # @show has_zero_tail

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

mutable struct TensorOperator <: SpaceTimeOperator
    spatial_factor
    temporal_factor
end

BEAST.defaultquadstrat(::TensorOperator, tfs, bfs) = nothing

⊗(A::AbstractOperator, B::AbstractOperator) = TensorOperator(A, B)
function scalartype(A::TensorOperator)
    promote_type(
        scalartype(A.spatial_factor),
        scalartype(A.temporal_factor))
end

function Base.:*(alpha::Number, A::TensorOperator)
    return TensorOperator(alpha*A.spatial_factor, A.temporal_factor)
end


function allocatestorage(op::TensorOperator, test_functions, trial_functions,
    ::Type{Val{:bandedstorage}}, long_delay_traits::Any)

    M = numfunctions(spatialbasis(test_functions))
    N = numfunctions(spatialbasis(trial_functions))

    time_basis_function = BEAST.convolve(
        temporalbasis(test_functions),
        temporalbasis(trial_functions))

    space_operator = op.spatial_factor
    A = assemble(space_operator, spatialbasis(test_functions), spatialbasis(trial_functions))

    K0 = ones(Int, M, N)
    bandwidth = numintervals(time_basis_function) - 1
    K1 = ones(Int,M,N) .+ (bandwidth - 1)
    T  = scalartype(op, test_functions, trial_functions)
    data = zeros(T, bandwidth, M, N)
    tail = zeros(T, M, N)

    Nt = numfunctions(temporalbasis(trial_functions))
    Z = ConvolutionOperators.ConvOp(data, K0, K1, tail, bandwidth)
    function store1(v,m,n,k)
		if Z.k0[m,n] ≤ k ≤ Z.k1[m,n]
			Z.data[k - Z.k0[m,n] + 1,m,n] += v
		elseif k == Z.k1[m,n]+1
			Z.tail[m,n] += v
		end
	end
    return ()->Z, store1
end



function assemble!(operator::TensorOperator, testfns::SpaceTimeBasis, trialfns::SpaceTimeBasis,
    store, threading::Type{Threading{:multi}};
    quadstrat=defaultquadstrat(operator, testfns, trialfns))

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
    assemble!(space_operator, space_testfns, space_trialfns,
        store1, threading)

end


mutable struct TemporalDifferentiation <: AbstractOperator
    operator
end

derive(op::AbstractOperator) = TemporalDifferentiation(op)
scalartype(op::TemporalDifferentiation) = scalartype(op.operator)
Base.:*(a::Number, op::TemporalDifferentiation) = TemporalDifferentiation(a * op.operator)

defaultquadstrat(op::TemporalDifferentiation, tfs, bfs) = defaultquadstrat(op.operator, tfs, bfs)

function allocatestorage(op::TemporalDifferentiation, testfns, trialfns,
	storage_trait, longdelays_trait)

	trial_time_fns  = temporalbasis(trialfns)
	trial_space_fns = spatialbasis(trialfns)

	trialfns = SpaceTimeBasis(
		trial_space_fns,
		derive(trial_time_fns)
	)

    return allocatestorage(op.operator, testfns, trialfns, storage_trait, longdelays_trait)
end

function assemble!(operator::TemporalDifferentiation, testfns, trialfns, store, threading = Threading{:multi};
    quadstrat=defaultquadstrat(operator, testfns, trialfns))

    trial_time_fns  = temporalbasis(trialfns)
    trial_space_fns = spatialbasis(trialfns)

    trialfns = SpaceTimeBasis(
        trial_space_fns,
        derive(trial_time_fns)
    )

    assemble!(operator.operator, testfns, trialfns, store, threading; quadstrat)

end

struct TemporalIntegration <: AbstractSpaceTimeOperator
    operator::AbstractSpaceTimeOperator
end

defaultquadstrat(op::TemporalIntegration, tfs, bfs) = defaultquadstrat(op.operator, tfs, bfs)

integrate(op::SpaceTimeOperator) = TemporalIntegration(op)
derive(op::TemporalIntegration) = op.operator
scalartype(op::TemporalIntegration) = scalartype(op.operator)
Base.:*(a::Number, op::TemporalIntegration) = TemporalIntegration(a * op.operator)

function allocatestorage(op::TemporalIntegration, testfns, trialfns,
	storage_trait::Type{Val{S}}, longdelays_trait) where {S}

	trial_time_fns  = temporalbasis(trialfns)
	trial_space_fns = spatialbasis(trialfns)

	trialfns = SpaceTimeBasis(
		trial_space_fns,
		integrate(trial_time_fns)
	)

	return allocatestorage(op.operator, testfns, trialfns, storage_trait, longdelays_trait)
end

function assemble!(operator::TemporalIntegration, testfns, trialfns, store,
    threading = Threading{:multi}; quadstrat=defaultquadstrat(operator, testfns, trialfns))

    trial_time_fns  = temporalbasis(trialfns)
    trial_space_fns = spatialbasis(trialfns)

    trialfns = SpaceTimeBasis(
        trial_space_fns,
        integrate(trial_time_fns)
    )

    assemble!(operator.operator, testfns, trialfns, store; quadstrat)

end
