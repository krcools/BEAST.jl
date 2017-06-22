export assemble

import Base: transpose, +, -, *

abstract type AbstractOperator end
abstract type Operator <: AbstractOperator end

# Support for operator transpose

type TransposedOperator <: Operator
    op::Operator
end

scalartype(op::TransposedOperator) = scalartype(op.op)

type LinearCombinationOfOperators{T} <: AbstractOperator
    coeffs::Vector{T}
    ops::Vector
end

function scalartype{T}(op::LinearCombinationOfOperators{T})
    W = promote_type(T, scalartype(op.ops[1]))
    for i in 2:length(op.ops)
        W = promote_type(W, scalartype(op.ops[i]))
    end
    W
end


function +(a::LinearCombinationOfOperators, b::Operator)
    LinearCombinationOfOperators(
        [a.coeffs;[1.0]],
        [a.ops;[b]]
    )
end

function +(a::LinearCombinationOfOperators, b::LinearCombinationOfOperators)
    LinearCombinationOfOperators(
        [a.coeffs; b.coeffs],
        [a.ops; b.ops]
    )
end

+(a::Operator, b::LinearCombinationOfOperators) = b + a
+(a::Operator, b::Number) = a + (b * Identity())
+(a::Number, b::Operator) = b + a

*(a::Number, b::Operator) = LinearCombinationOfOperators([a], [b])
*(a::Number, b::LinearCombinationOfOperators) = LinearCombinationOfOperators(a * b.coeffs, b.ops)
-(a::AbstractOperator, b::AbstractOperator) = a + (-1.0) * b

function +(a::Operator, b::Operator)
    LinearCombinationOfOperators(
        [1.0, 1.0],
        [a, b]
    )
end


transpose(op::TransposedOperator) = op.op
transpose(op::Operator) = TransposedOperator(op)


function assemble(operator::AbstractOperator, test_functions, trial_functions)
    # This is a convenience function whose only job is to allocate
    # the storage for the interaction matrix. Further dispatch on
    # operator and space types is handled by the 4-argument version
    #Z = allocatestorage(operator, test_functions, trial_functions)
    #store(v, m, n) = (Z[m,n] += v)
    Z, store = allocatestorage(operator, test_functions, trial_functions)
    assemble!(operator, test_functions, trial_functions, store)
    sdata(Z)
end

function allocatestorage(operator::AbstractOperator, test_functions, trial_functions)
    T = promote_type(
        scalartype(operator)       ,
        scalartype(test_functions) ,
        scalartype(trial_functions),
    )
    Z = SharedArray{T}(
        numfunctions(test_functions)  ,
        numfunctions(trial_functions),
    )
    fill!(Z, 0)
    store(v,m,n) = (Z[m,n] += v)
    return Z, store
end

function allocatestorage(operator::LinearCombinationOfOperators,
        test_functions::SpaceTimeBasis, trial_functions::SpaceTimeBasis)

    # TODO: remove this ugly, ugly patch
    Z, store = allocatestorage(operator.ops[end], test_functions, trial_functions)
end


function assemble!(operator::Operator, test_functions::Space, trial_functions::Space, store)

    # This method should only be called for `atomic` discrete operators, this means
    # in particular that the spaces of test and trial functions are fully conforming
    # to the Space concept and that the operator/space combinations conform the
    # concept implicitly defined by the assemble! function in itegralop.jl
    # (no more transposes or repositioning in a larger system for example)

    P = procs()
    if length(P) > 1; P = P[2:end]; end
    numchunks = length(P)
    @assert numchunks >= 1
    splits = [round(Int,s) for s in linspace(0, numfunctions(test_functions), numchunks+1)]

    T = typeof(test_functions)
    S = eltype(test_functions.fns)

    @sync begin
        for (i,p) in enumerate(P)
            I = splits[i]+1 : splits[i+1]

            fns_p = similar(test_functions.fns)
            fill!(fns_p, S[])
            for i in I
                fns_p[i] = test_functions.fns[i]
            end
            test_functions_p = T(test_functions.geo, fns_p)

            @async remotecall_wait(assemblechunk!, p, operator, test_functions_p, trial_functions, store)
        end
    end

end


function assemble!(op::TransposedOperator, tfs::Space, bfs::Space, store)

    store1(v,m,n) = store(v,n,m)
    assemble!(op.op, bfs, tfs, store1)
end


#function assemble!(op::LinearCombinationOfOperators, tfs::Space, bfs::Space, store)
function assemble!(op::LinearCombinationOfOperators, tfs::AbstractSpace, bfs::AbstractSpace, store)
#function assemble!(op::LinearCombinationOfOperators, tfs, bfs, store)
    for (a,A) in zip(op.coeffs, op.ops)
        store1(v,m,n) = store(a*v,m,n)
        assemble!(A, tfs, bfs, store1)
    end
end


# Support for direct product spaces
function assemble!(op::Operator, tfs::DirectProductSpace, bfs::Space, store)
    I = Int[0]
    for s in tfs.factors push!(I, last(I) + numfunctions(s)) end
    for (i,s) in enumerate(tfs.factors)
        store1(v,m,n) = store(v,m + I[i], n)
        assemble!(op, s, bfs, store1)
    end
end


function assemble!(op::Operator, tfs::Space, bfs::DirectProductSpace, store)
    J = Int[0]
    for s in bfs.factors push!(J, last(J) + numfunctions(s)) end
    for (j,s) in enumerate(bfs.factors)
        store1(v,m,n) = store(v,m,n + J[j])
        assemble!(op, tfs, s, store1)
    end
end

function assemble!(op::Operator, tfs::DirectProductSpace, bfs::DirectProductSpace, store)
    I = Int[0]
    for s in tfs.factors push!(I, last(I) + numfunctions(s)) end
    for (i,s) in enumerate(tfs.factors)
        store1(v,m,n) = store(v,m + I[i],n)
        assemble!(op, s, bfs, store1)
    end
end
