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
