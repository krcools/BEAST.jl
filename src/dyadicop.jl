struct DyadicOp{T,F,G} <: Operator
    α::T
    f::F
    g::G
end

function DyadicOp(f::F, g::G) where {F,G}
    T = scalartype(f,g)
    DyadicOp{T,F,G}(one(T), f, g)
end

Base.:*(b::Number, op::DyadicOp) = DyadicOp(b * op.α, op.f, op.g)

scalartype(op::DyadicOp) = promote_type(scalartype(op.α), scalartype(op.f), scalartype(op.g))
defaultquadstrat(op::DyadicOp, testspace, refspace) = nothing


function allocatestorage(op::DyadicOp, test_functions, trial_functions,
    storage_trait, longdelays_trait)

    T = scalartype(op, test_functions, trial_functions)
    M = length(test_functions)
    N = length(trial_functions)

    u = zeros(T,M)
    v = zeros(T,N)
    
    A = BEAST.Rank1Map{T}(u,v)
    store() = A
    freeze() = A

    return freeze, store
end


function assemble!(biop::DyadicOp, tfs::Space, bfs::Space, store,
    threading::Type{BEAST.Threading{:multi}};
    kwargs...)

    u = assemble(biop.f, tfs)
    v = assemble(biop.g, bfs)

    A = store()
    A.u .= biop.α * u
    A.v .= adjoint.(v)
end