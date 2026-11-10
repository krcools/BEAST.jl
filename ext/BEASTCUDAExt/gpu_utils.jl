function launch_gpu_kernel!(gpu_kernel, args...; gpu_blocksize, problem_size)
    block_dims = gpu_blocksize isa Tuple ? gpu_blocksize : (gpu_blocksize,)
    problem_dims = problem_size isa Tuple ? problem_size : (problem_size,)

    any(iszero, problem_dims) && return
    @assert length(block_dims) == length(problem_dims)
    @assert all(>(0), block_dims)
    @assert prod(block_dims) <= 1024

    grid_dims = map(cld, problem_dims, block_dims)
    threads = length(block_dims) == 1 ? block_dims[1] : block_dims
    blocks = length(grid_dims) == 1 ? grid_dims[1] : grid_dims
    @cuda blocks=blocks threads=threads gpu_kernel(args...)

    return
end

function resolve_gpu_quadstrat(quadstrat, operator, test_space, trial_space)
    for space in (test_space, trial_space)
        mesh = geometry(space)
        element = chart(mesh, first(mesh))
        element isa CompScienceMeshes.Simplex{3,2} ||
            throw(ArgumentError(
                "GPU assembly currently supports triangular surface " *
                "elements in three dimensions"))
    end

    strategy = quadstrat(operator, test_space, trial_space)

    if CompScienceMeshes.refines(geometry(test_space), geometry(trial_space)) ||
        CompScienceMeshes.refines(geometry(trial_space), geometry(test_space))
        throw(ArgumentError(
            "GPU assembly does not support test and trial meshes in a " *
            "refinement relation"))
    end

    strategy isa Union{BEAST.DoubleNumQStrat,BEAST.DoubleNumSauterQstrat} ||
        throw(ArgumentError(
            "GPU assembly currently supports DoubleNumQStrat and " *
            "DoubleNumSauterQstrat, got $(typeof(strategy))"))

    return strategy
end
