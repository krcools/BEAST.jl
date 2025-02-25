using .LinearSpace

struct LongDelays{T} end
struct Threading{T} end

import Base: transpose, +, -, *

abstract type AbstractOperator end
"""
*Atomic operator*: one that assemblechunk can deal with
"""
abstract type Operator <: AbstractOperator end
abstract type IntegralOperator <: Operator end

mutable struct TransposedOperator <: AbstractOperator
    op::AbstractOperator
end

scalartype(op::TransposedOperator) = scalartype(op.op)
defaultquadstrat(op::TransposedOperator, tfs::Space, bfs::Space) = defaultquadstrat(op.op, tfs, bfs)

"""
    LinearCombinationOfOperators{T} <: AbstractOperator

A linear combination of operators.
"""
mutable struct LinearCombinationOfOperators{T} <: AbstractOperator
    coeffs::Vector{T}
    ops::Vector
end

function scalartype(op::LinearCombinationOfOperators{T}) where T
    W = promote_type(T, scalartype(op.ops[1]))
    for i in 2:length(op.ops)
        W = promote_type(W, scalartype(op.ops[i]))
    end
    W
end

function derive(a::LinearCombinationOfOperators)
    coeffs = copy(a.coeffs)
    ops = [derive(op) for op in a.ops]
    return LinearCombinationOfOperators(coeffs, ops)
end


function +(a::LinearCombinationOfOperators, b::AbstractOperator)
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

+(a::AbstractOperator, b::LinearCombinationOfOperators) = b + a
+(a::AbstractOperator, b::Number) = a + (b * Identity())
+(a::Number, b::AbstractOperator) = b + a

*(a::Number, b::AbstractOperator) = LinearCombinationOfOperators([a], [b])
*(a::Number, b::LinearCombinationOfOperators) = LinearCombinationOfOperators(a * b.coeffs, b.ops)
-(a::AbstractOperator, b::AbstractOperator) = a + (-1.0) * b
-(a::AbstractOperator) = (-1.0) * a

function +(a::AbstractOperator, b::AbstractOperator)
    LinearCombinationOfOperators(
        [1.0, 1.0],
        [a, b]
    )
end


transpose(op::TransposedOperator) = op.op
transpose(op::Operator) = TransposedOperator(op)

defaultquadstrat(lc::LinearCombinationOfOperators, tfs, bfs) =
    [defaultquadstrat(op,tfs,bfs) for op in lc.ops]

"""
    assemble(operator, test_functions, trial_functions;
        storage_policy = Val{:bandedstorage},
        threading = Threading{:multi},
        quadstrat=defaultquadstrat(operator, test_functions, trial_functions))

Assemble the system matrix corresponding to the operator `operator` tested with the test functions `test_functions` and the trial functions `trial_functions`.
"""
function assemble(operator::AbstractOperator, test_functions, trial_functions;
    storage_policy = Val{:bandedstorage},
    threading = Threading{:multi},
    quadstrat=defaultquadstrat)

    Z, store = allocatestorage(operator, test_functions, trial_functions,
        storage_policy)
    # qs = quadstrat(operator, test_functions, trial_functions)
    assemble!(operator, test_functions, trial_functions,
        store, threading; quadstrat)
    return Z()
end


function assemble(A::AbstractMatrix, testfns, trialfns; kwargs...)
    @assert numfunctions(testfns) == size(A,1)
    @assert numfunctions(trialfns) == size(A,2)
    return A
end

function assemblerow(operator::AbstractOperator, test_functions, trial_functions,
    storage_policy = Val{:bandedstorage},
    long_delays_policy = LongDelays{:ignore};
    quadstrat=defaultquadstrat(operator, test_functions, trial_functions))

    Z, store = allocatestorage(operator, test_functions, trial_functions,
        storage_policy, long_delays_policy)
    assemblerow!(operator, test_functions, trial_functions, store; quadstrat)

    Z()
end

function assemblecol(operator::AbstractOperator, test_functions, trial_functions,
    storage_policy = Val{:bandestorage},
    long_delays_policy = LongDelays{:ignore};
    quadstrat=defaultquadstrat(operator, test_functions, trial_functions))

    Z, store = allocatestorage(operator, test_functions, trial_functions,
        storage_policy, long_delays_policy)
    assemblecol!(operator, test_functions, trial_functions, store; quadstrat)

    Z()
end

function allocatestorage(operator::AbstractOperator, test_functions, trial_functions,
    storage_trait=nothing, longdelays_trait=nothing)

    T = promote_type(
        scalartype(operator)       ,
        scalartype(test_functions) ,
        scalartype(trial_functions),
    )
    Z = Matrix{T}(undef,
        numfunctions(test_functions),
        numfunctions(trial_functions),
    )
    fill!(Z, 0)
    store(v,m,n) = (Z[m,n] += v)
    return ()->Z, store
end


# function allocatestorage(operator::LinearCombinationOfOperators,
#         test_functions::SpaceTimeBasis, trial_functions::SpaceTimeBasis,
#         storage_policy::Type{Val{:bandedstorage}},
#         long_delays_policy::Type{LongDelays{:ignore}})

#     # This works when are terms in the LC can share storage
#     return allocatestorage(operator.ops[end], test_functions, trial_functions,
#         storage_policy, long_delays_policy)
# end


function allocatestorage(operator::LinearCombinationOfOperators,
        test_functions::SpaceTimeBasis, trial_functions::SpaceTimeBasis,
        storage_policy::Type{Val{S}},
        long_delays_policy::Type{LongDelays{L}}) where {L,S}

    # This works when are terms in the LC can share storage
    return allocatestorage(operator.ops[end], test_functions, trial_functions,
        storage_policy, long_delays_policy)
end

struct _OffsetStore{F}
    store::F
    row_offset::Int
    col_offset::Int
end

(f::_OffsetStore)(v,m,n) = f.store(v,m + f.row_offset, n + f.col_offset)

function assemble!(operator::Operator, test_functions::Space, trial_functions::Space,
    store, threading::Type{Threading{:multi}};
    quadstrat=defaultquadstrat)

    quadstrat = quadstrat(operator, test_functions, trial_functions)

    P = Threads.nthreads()
    numchunks = P
    @assert numchunks >= 1
    splits = [round(Int,s) for s in range(0, stop=numfunctions(test_functions), length=numchunks+1)]

    Threads.@threads for i in 1:P
        lo, hi = splits[i]+1, splits[i+1]
        lo <= hi || continue
        test_functions_p = subset(test_functions, lo:hi)
        store1 = _OffsetStore(store, lo-1, 0)
        assemblechunk!(operator, test_functions_p, trial_functions, store1, quadstrat=quadstrat)
end end

function assemble!(operator::Operator, test_functions::Space, trial_functions::Space,
    store, threading::Type{Threading{:single}};
    quadstrat=defaultquadstrat)

    quadstrat = quadstrat(operator, test_functions, trial_functions)
    assemblechunk!(operator, test_functions, trial_functions, store; quadstrat)
end



function assemble!(op::TransposedOperator, tfs::Space, bfs::Space,
    store, threading = Threading{:multi};
    quadstrat=defaultquadstrat(op, tfs, bfs))

    store1(v,m,n) = store(v,n,m)
    assemble!(op.op, bfs, tfs, store1, threading; quadstrat)
end


function assemble!(op::LinearCombinationOfOperators, tfs::AbstractSpace, bfs::AbstractSpace,
    store, threading = Threading{:multi};
    quadstrat=defaultquadstrat(op, tfs, bfs))

    for (a,A) in zip(op.coeffs, op.ops)
        store1(v,m,n) = store(a*v,m,n)
        qs = quadstrat(A, tfs, bfs)
        assemble!(A, tfs, bfs, store1, threading; quadstrat=qs)
    end
end


# Support for direct product spaces
function assemble!(op::AbstractOperator, tfs::DirectProductSpace, bfs::Space,
    store, threading = Threading{:multi};
    quadstrat=defaultquadstrat(op, tfs[1], bfs))

    I = Int[0]
    for s in tfs.factors push!(I, last(I) + numfunctions(s)) end
    for (i,s) in enumerate(tfs.factors)
        store1(v,m,n) = store(v,m + I[i], n)
        assemble!(op, s, bfs, store1, threading; quadstrat)
    end
end


function assemble!(op::AbstractOperator, tfs::Space, bfs::DirectProductSpace,
    store, threading=Threading{:multi};
    quadstrat=defaultquadstrat(op, tfs, bfs[1]))

    J = Int[0]
    for s in bfs.factors push!(J, last(J) + numfunctions(s)) end
    for (j,s) in enumerate(bfs.factors)
        store1(v,m,n) = store(v,m,n + J[j])
        assemble!(op, tfs, s, store1, threading; quadstrat)
    end
end

function assemble!(op::AbstractOperator, tfs::DirectProductSpace, bfs::DirectProductSpace,
    store, threading=Threading{:multi};
    quadstrat=defaultquadstrat(op, tfs[1], bfs[1]))
    
    I = Int[0]
    for s in tfs.factors push!(I, last(I) + numfunctions(s)) end
    for (i,s) in enumerate(tfs.factors)
        store1(v,m,n) = store(v,m + I[i],n)
        assemble!(op, s, bfs, store1, threading; quadstrat)
    end
end

# TODO: Remove BlockDiagonalOperator in favour of manipulations
#       on the level of bilinear forms.

# Discretisation and assembly of these operators
# will respect the direct product structure of the
# HilbertSpace/FiniteElmeentSpace
struct BlockDiagonalOperator <: AbstractOperator
    op::AbstractOperator
end

diag(op::AbstractOperator) = BlockDiagonalOperator(op)

scalartype(op::BlockDiagonalOperator) = scalartype(op.op)
defaultquadstrat(op::BlockDiagonalOperator, U::DirectProductSpace, V::DirectProductSpace) = defaultquadstrat(op.op, U, V)
allocatestorage(op::BlockDiagonalOperator, X, Y, storage_trait, longdelays_trait) =
    allocatestorage(op.op, X, Y, storage_trait, longdelays_trait)

function assemble!(op::BlockDiagonalOperator, U::DirectProductSpace, V::DirectProductSpace,
    store, threading=Threading{:multi};
    quadstrat = defaultquadstrat(op, U, V))
    
    @assert length(U.factors) == length(V.factors)
    I = Int[0]; for u in U.factors push!(I, last(I) + numfunctions(u)) end
    J = Int[0]; for v in V.factors push!(J, last(J) + numfunctions(v)) end

    for (k,(u,v)) in enumerate(zip(U.factors, V.factors))
        store1(v,m,n) = store(v, I[k] + m, J[k] + n)
        assemble!(op.op, u, v, store1, threading; quadstrat)
    end
end

# BlockFull is default so not sure when this exists -> remove

struct BlockFullOperators <: AbstractOperator
    op::AbstractOperator
end

blocks(op::AbstractOperator) = BlockFullOperators(op)

scalartype(op::BlockFullOperators) = scalartype(op.op)
defaultquadstrat(op::BlockFullOperators, U::DirectProductSpace, V::DirectProductSpace) = defaultquadstrat(op.op, U, V)


function assemble!(op::BlockFullOperators, U::DirectProductSpace, V::DirectProductSpace,
    store, threading;
    quadstrat = defaultquadstrat(op, U, V))
    
    I = Int[0]; for u in U.factors push!(I, last(I) + numfunctions(u)) end
    J = Int[0]; for v in V.factors push!(J, last(J) + numfunctions(v)) end

    for (k,u) in enumerate(U.factors)
        for (l,v) in enumerate(V.factors)
            store1(x,m,n) = store(x, I[k]+m, J[l]+n)
            assemble!(op.op, u, v, store1, threading; quadstrat)
        end
    end
end
