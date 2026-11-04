struct InterpolatorOp <: Operator
end

function scalartype(op::InterpolatorOp) Float32 end
function defaultquadstrat(op::InterpolatorOp,
    tfs::SpaceWithDoFs, bfs::SpaceWithDoFs)
    nothing
end

function assemble!(operator::InterpolatorOp,
    test_functions::SpaceWithDoFs, trial_functions::SpaceWithDoFs,
    store, threading::Type{Threading{:multi}};
    quadstrat=defaultquadstrat, verbose=true, kwargs...)

    quadstrat = quadstrat(operator, test_functions, trial_functions)
    assemblechunk!(operator, test_functions, trial_functions, store; quadstrat, verbose)
end

function allocatestorage(op::InterpolatorOp,
    tfs::SpaceWithDoFs, bfs::SpaceWithDoFs,
    storage)

    A = Ref{Any}()
    function store() A end
    function freeze() return A[] end
    return freeze, store
end

function assemblechunk!(op::InterpolatorOp,
    tfs::SpaceWithDoFs, bfs::SpaceWithDoFs, store;
    quadstrat=defaultquadstrat(op, tfs, bfs), verbose=true)

    Π = tfs.interpolator
    store()[] = Π(bfs.space)

    return nothing
end