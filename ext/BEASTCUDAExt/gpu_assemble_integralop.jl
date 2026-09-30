function l2g_maps!(assembly_data, nfunctions)
    active_dofs = falses(nfunctions)
    for index in CartesianIndices(assembly_data.data)
        dof = assembly_data.data[index][1]
        dof > 0 && (active_dofs[dof] = true)
    end

    local_to_global = findall(active_dofs)
    global_to_local = zeros(Int, nfunctions)
    for (local_dof, global_dof) in enumerate(local_to_global)
        global_to_local[global_dof] = local_dof
    end

    for index in CartesianIndices(assembly_data.data)
        dof, coefficient = assembly_data.data[index]
        if dof > 0
            assembly_data.data[index] = (global_to_local[dof], coefficient)
        end
    end

    return local_to_global
end

function load_assemblydata_gpu(functions::Space, elements, assembly_data,
    ::Type{T}) where T
    element_domain = domain(first(elements))
    num_shapes = numfunctions(refspace(functions), element_domain)

    local_to_global = l2g_maps!(assembly_data, numfunctions(functions))

    rows = Int[]
    cols = Int[]
    values = T[]
    for index in CartesianIndices(assembly_data.data)
        contribution, shape, element = Tuple(index)
        dof, coefficient = assembly_data.data[contribution, shape, element]
        if dof > 0
            push!(rows, dof)
            push!(cols, num_shapes * (element - 1) + shape)
            push!(values, T(coefficient))
        end
    end

    assembly_matrix = sparse(
        rows, cols, values, length(local_to_global), num_shapes * length(elements))

    return CuArray(elements), CuSparseMatrixCSC(assembly_matrix),
        local_to_global, element_domain, coordtype(first(elements))
end

function load_assemblydata_gpu(functions::Space, ::Type{T}) where T
    elements, assembly_data, _ = BEAST.assemblydata(functions)
    return load_assemblydata_gpu(functions, elements, assembly_data, T)
end

function gpu_triangle_rule(order, ::Type{T}) where T
    points, weights = CompScienceMeshes.trgauss(order)
    rule = Vector{Tuple{SVector{2,T},T}}(undef, length(weights))
    for i in eachindex(weights)
        point = SVector{2,T}(points[1,i], points[2,i])
        rule[i] = (point, T(weights[i]))
    end
    return CuArray(rule)
end

function gpu_legendre_rule(order, ::Type{T}) where T
    rule = convert.(NTuple{2,T}, BEAST._legendre(order, zero(T), one(T)))
    return CuArray(rule)
end

function gpu_singularityflag!(singularity_map,
    test_elements::CuDeviceVector{S}, trial_elements::CuDeviceVector{S}) where S

    pair = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    num_test_elements = length(test_elements)
    num_trial_elements = length(trial_elements)

    if pair <= num_test_elements * num_trial_elements
        test_index = mod(pair - 1, num_test_elements) + 1
        trial_index = div(pair - 1, num_test_elements) + 1
        tolerance = 1.0e3 * eps(coordtype(S))

        hits = 1
        for test_vertex in vertices(test_elements[test_index])
            for trial_vertex in vertices(trial_elements[trial_index])
                hits += norm(test_vertex - trial_vertex) < tolerance
            end
        end
        singularity_map[pair, hits] = true
    end

    return nothing
end

function gpu_compact_pair_map!(pair_map, singularity_map, cumulative_counts)
    pair = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    singularity = threadIdx().y

    if pair <= size(singularity_map, 1) && singularity <= size(singularity_map, 2)
        if singularity_map[pair, singularity]
            pair_map[cumulative_counts[pair, singularity], singularity] = pair
        end
    end

    return nothing
end

function gpu_singularitydetection!(pair_map, num_pairs, test_elements, trial_elements)
    num_element_pairs = length(test_elements) * length(trial_elements)
    singularity_map = CUDA.zeros(Bool, num_element_pairs, 4)

    launch_gpu_kernel!(gpu_singularityflag!, singularity_map,
        test_elements, trial_elements;
        gpu_blocksize=256, problem_size=num_element_pairs)

    cumulative_counts = accumulate(+, singularity_map; dims=1)
    num_pairs .= Array(cumulative_counts[end, :])

    launch_gpu_kernel!(gpu_compact_pair_map!, pair_map,
        singularity_map, cumulative_counts;
        gpu_blocksize=(256, 4), problem_size=size(singularity_map))

    return nothing
end

function gpu_shapefunction_eval!(shape_values, elements, local_space, quadrule)
    element_index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    quadrature_index = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if element_index <= length(elements) && quadrature_index <= length(quadrule)
        element = elements[element_index]
        point = quadrule[quadrature_index][1]
        shape_values[element_index, quadrature_index] =
            local_space(neighborhood(element, point))
    end

    return nothing
end

@inline function gpu_doublenum_integral(integrand,
    test_element, trial_element, test_values, trial_values,
    test_local_space, trial_local_space, test_quadrule, trial_quadrule,
    ::Type{T}) where T

    num_test_shapes = numfunctions(test_local_space, domain(test_element))
    num_trial_shapes = numfunctions(trial_local_space, domain(trial_element))
    z = zeros(SMatrix{num_test_shapes,num_trial_shapes,T})

    for test_quadrature_index in eachindex(test_quadrule)
        test_point, test_weight = test_quadrule[test_quadrature_index]
        x = neighborhood(test_element, test_point)
        weighted_test_jacobian = test_weight * jacobian(x)

        for trial_quadrature_index in eachindex(trial_quadrule)
            trial_point, trial_weight = trial_quadrule[trial_quadrature_index]
            y = neighborhood(trial_element, trial_point)
            weight = weighted_test_jacobian * trial_weight * jacobian(y)
            z += weight * integrand(x, y,
                test_values[test_quadrature_index],
                trial_values[trial_quadrature_index])
        end
    end

    return z
end

function gpu_momintegral_doublenum_pair!(zlocal, operator, pair,
    test_elements, trial_elements, test_shapes, trial_shapes,
    test_local_space, trial_local_space, test_quadrule, trial_quadrule)

    num_test_elements = length(test_elements)
    test_index = mod(pair - 1, num_test_elements) + 1
    trial_index = div(pair - 1, num_test_elements) + 1

    test_element = test_elements[test_index]
    trial_element = trial_elements[trial_index]
    num_test_shapes = numfunctions(test_local_space, domain(test_element))
    num_trial_shapes = numfunctions(trial_local_space, domain(trial_element))

    integrand = BEAST.Integrand(operator, test_local_space, trial_local_space,
        test_element, trial_element)

    z = gpu_doublenum_integral(integrand,
        test_element, trial_element,
        view(test_shapes, test_index, :), view(trial_shapes, trial_index, :),
        test_local_space, trial_local_space, test_quadrule, trial_quadrule,
        eltype(zlocal))

    test_offset = num_test_shapes * (test_index - 1)
    trial_offset = num_trial_shapes * (trial_index - 1)
    for trial_shape in 1:num_trial_shapes
        for test_shape in 1:num_test_shapes
            zlocal[test_offset + test_shape, trial_offset + trial_shape] =
                z[test_shape, trial_shape]
        end
    end

    return nothing
end

function gpu_momintegral_doublenum_allpairs!(zlocal, operator,
    test_elements, trial_elements, test_shapes, trial_shapes,
    test_local_space, trial_local_space, test_quadrule, trial_quadrule)

    pair = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if pair <= length(test_elements) * length(trial_elements)
        gpu_momintegral_doublenum_pair!(zlocal, operator, pair,
            test_elements, trial_elements, test_shapes, trial_shapes,
            test_local_space, trial_local_space, test_quadrule, trial_quadrule)
    end

    return nothing
end

function gpu_momintegral_doublenum!(zlocal, operator, num_pairs, pairs,
    test_elements, trial_elements, test_shapes, trial_shapes,
    test_local_space, trial_local_space, test_quadrule, trial_quadrule)

    pair_index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if pair_index <= num_pairs
        gpu_momintegral_doublenum_pair!(zlocal, operator, pairs[pair_index],
            test_elements, trial_elements, test_shapes, trial_shapes,
            test_local_space, trial_local_space, test_quadrule, trial_quadrule)
    end

    return nothing
end

@inline function gpu_sauterschwab_reorder(
    test_vertices, trial_vertices, strategy)

    I = MVector{3,Int}(undef)
    J = MVector{3,Int}(undef)
    K = MVector{3,Int}(undef)
    L = MVector{3,Int}(undef)
    SauterSchwabQuadrature.reorder!(
        I, J, K, L, test_vertices, trial_vertices, strategy)

    return SVector(I), SVector(J)
end

@inline function gpu_sauterschwab_integral(integrand,
    test_element, trial_element, strategy,
    test_local_space, trial_local_space, ::Type{T}) where T

    num_test_shapes = numfunctions(test_local_space, domain(test_element))
    num_trial_shapes = numfunctions(trial_local_space, domain(trial_element))
    I, J = gpu_sauterschwab_reorder(
        vertices(test_element), vertices(trial_element), strategy)
    pulledback_integrand = BEAST.pulledback_integrand(
        integrand, I, test_element, J, trial_element)

    z = zeros(SMatrix{num_test_shapes,num_trial_shapes,T})
    for (eta_1, weight_1) in strategy.qps
        for (eta_2, weight_2) in strategy.qps
            for (eta_3, weight_3) in strategy.qps
                for (xi, weight_4) in strategy.qps
                    weight = weight_1 * weight_2 * weight_3 * weight_4
                    z += weight * strategy(
                        pulledback_integrand, eta_1, eta_2, eta_3, xi)
                end
            end
        end
    end

    return z
end

function gpu_momintegral_sauterschwab!(zlocal, operator, num_pairs, pairs,
    test_elements, trial_elements, test_local_space, trial_local_space, strategy)

    pair_index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if pair_index <= num_pairs
        num_test_elements = length(test_elements)
        pair = pairs[pair_index]
        test_index = mod(pair - 1, num_test_elements) + 1
        trial_index = div(pair - 1, num_test_elements) + 1

        test_element = test_elements[test_index]
        trial_element = trial_elements[trial_index]
        num_test_shapes = numfunctions(test_local_space, domain(test_element))
        num_trial_shapes = numfunctions(trial_local_space, domain(trial_element))

        integrand = BEAST.Integrand(operator, test_local_space, trial_local_space,
            test_element, trial_element)
        z = gpu_sauterschwab_integral(integrand,
            test_element, trial_element, strategy,
            test_local_space, trial_local_space, eltype(zlocal))

        test_offset = num_test_shapes * (test_index - 1)
        trial_offset = num_trial_shapes * (trial_index - 1)
        for trial_shape in 1:num_trial_shapes
            for test_shape in 1:num_test_shapes
                zlocal[test_offset + test_shape, trial_offset + trial_shape] =
                    z[test_shape, trial_shape]
            end
        end
    end

    return nothing
end

function build_matrix!(matrix, zlocal, test_assembly, trial_assembly)
    T = promote_type(eltype(zlocal), eltype(trial_assembly))
    trial_projection = CuArray{T}(
        undef, size(trial_assembly, 1), size(zlocal, 1))

    CUSPARSE.mm!('N', 'T', one(T), trial_assembly, zlocal,
        zero(T), trial_projection, 'O')
    CUSPARSE.mm!('N', 'T', one(eltype(matrix)), test_assembly, trial_projection,
        zero(eltype(matrix)), matrix, 'O')

    return nothing
end

function assemble_primer_gpu(functions::Space, quadrule, ::Type{T}) where T
    elements, assembly_data, local_to_global, element_domain, coordinate_type =
        load_assemblydata_gpu(functions, T)

    return assemble_primer_gpu(functions, elements, assembly_data,
        local_to_global, element_domain, coordinate_type, quadrule)
end

function assemble_primer_gpu(functions::Space, elements, assembly_data,
    quadrule, ::Type{T}) where T
    elements, assembly_data, local_to_global, element_domain, coordinate_type =
        load_assemblydata_gpu(functions, elements, assembly_data, T)

    return assemble_primer_gpu(functions, elements, assembly_data,
        local_to_global, element_domain, coordinate_type, quadrule)
end

function assemble_primer_gpu(functions::Space, elements, assembly_data,
    local_to_global, element_domain, coordinate_type, quadrule)
    local_space = refspace(functions)

    quadrule_d = gpu_triangle_rule(quadrule, coordinate_type)
    num_shapes = numfunctions(local_space, element_domain)
    shape_type = shapetype(local_space)
    shape_values = CuArray{SVector{num_shapes,shape_type}}(
        undef, length(elements), length(quadrule_d))

    launch_gpu_kernel!(gpu_shapefunction_eval!,
        shape_values, elements, local_space, quadrule_d;
        gpu_blocksize=(64, 4), problem_size=(length(elements), length(quadrule_d)))

    return local_to_global, (elements, assembly_data), (quadrule_d, shape_values)
end

function assemblechunk_body_gpu_device!(operator::IntegralOperator,
    test_local_space, test_elements, test_assembly, test_quadrule, test_shapes,
    trial_local_space, trial_elements, trial_assembly, trial_quadrule, trial_shapes,
    quadstrat::BEAST.DoubleNumQStrat;
    gpu_blocksize, singular_rules=nothing)

    num_test_shapes = div(size(test_assembly, 2), length(test_elements))
    num_trial_shapes = div(size(trial_assembly, 2), length(trial_elements))
    T = promote_type(scalartype(operator),
        eltype(test_assembly), eltype(trial_assembly))

    zlocal = CuArray{T}(undef,
        num_test_shapes * length(test_elements),
        num_trial_shapes * length(trial_elements))

    num_pairs = length(test_elements) * length(trial_elements)
    launch_gpu_kernel!(gpu_momintegral_doublenum_allpairs!, zlocal, operator,
        test_elements, trial_elements, test_shapes, trial_shapes,
        test_local_space, trial_local_space, test_quadrule, trial_quadrule;
        gpu_blocksize=gpu_blocksize, problem_size=num_pairs)

    matrix = CuArray{T}(
        undef, size(test_assembly, 1), size(trial_assembly, 1))
    build_matrix!(matrix, zlocal, test_assembly, trial_assembly)

    return matrix
end

function assemblechunk_body_gpu_device!(operator::IntegralOperator,
    test_local_space, test_elements, test_assembly, test_quadrule, test_shapes,
    trial_local_space, trial_elements, trial_assembly, trial_quadrule, trial_shapes,
    quadstrat::BEAST.DoubleNumSauterQstrat;
    gpu_blocksize, singular_rules=nothing)

    num_test_shapes = div(size(test_assembly, 2), length(test_elements))
    num_trial_shapes = div(size(trial_assembly, 2), length(trial_elements))
    T = promote_type(scalartype(operator),
        eltype(test_assembly), eltype(trial_assembly))

    num_element_pairs = length(test_elements) * length(trial_elements)
    zlocal = CuArray{T}(undef,
        num_test_shapes * length(test_elements),
        num_trial_shapes * length(trial_elements))
    pair_map = CuArray{Int}(undef, num_element_pairs, 4)
    num_pairs = zeros(Int, 4)
    gpu_singularitydetection!(pair_map, num_pairs, test_elements, trial_elements)

    launch_gpu_kernel!(gpu_momintegral_doublenum!, zlocal, operator,
        num_pairs[1], view(pair_map, :, 1),
        test_elements, trial_elements, test_shapes, trial_shapes,
        test_local_space, trial_local_space, test_quadrule, trial_quadrule;
        gpu_blocksize=gpu_blocksize, problem_size=num_pairs[1])

    if isnothing(singular_rules)
        coordinate_type = coordtype(eltype(test_elements))
        common_vertex = CommonVertex(gpu_legendre_rule(
            quadstrat.sauter_schwab_common_vert, coordinate_type))
        common_edge = CommonEdge(gpu_legendre_rule(
            quadstrat.sauter_schwab_common_edge, coordinate_type))
        common_face = CommonFace(gpu_legendre_rule(
            quadstrat.sauter_schwab_common_face, coordinate_type))
    else
        common_vertex = CommonVertex(singular_rules.common_vertex)
        common_edge = CommonEdge(singular_rules.common_edge)
        common_face = CommonFace(singular_rules.common_face)
    end

    for (singularity, strategy) in
        ((2, common_vertex), (3, common_edge), (4, common_face))

        launch_gpu_kernel!(gpu_momintegral_sauterschwab!, zlocal, operator,
            num_pairs[singularity], view(pair_map, :, singularity),
            test_elements, trial_elements, test_local_space, trial_local_space,
            strategy;
            gpu_blocksize=gpu_blocksize, problem_size=num_pairs[singularity])
    end

    matrix = CuArray{T}(
        undef, size(test_assembly, 1), size(trial_assembly, 1))
    build_matrix!(matrix, zlocal, test_assembly, trial_assembly)

    return matrix
end

function assemblechunk_body_gpu!(operator::IntegralOperator,
    test_local_space, test_elements, test_assembly, test_quadrule, test_shapes,
    trial_local_space, trial_elements, trial_assembly, trial_quadrule, trial_shapes,
    quadstrat::Union{BEAST.DoubleNumQStrat,BEAST.DoubleNumSauterQstrat};
    gpu_blocksize, singular_rules=nothing)

    matrix = assemblechunk_body_gpu_device!(operator,
        test_local_space, test_elements, test_assembly, test_quadrule, test_shapes,
        trial_local_space, trial_elements, trial_assembly, trial_quadrule, trial_shapes,
        quadstrat; gpu_blocksize, singular_rules)

    return Array(matrix)
end

struct GPUDenseAssemblyData{D,X,Y,R}
    device::D
    test_data::X
    trial_data::Y
    singular_rules::R
end

function gpu_dense_assembly_data(device,
    test_functions, test_elements, test_assembly, test_tiles,
    trial_functions, trial_elements, trial_assembly, trial_tiles,
    quadstrat, ::Type{T}) where T

    CUDA.device!(device)
    test_data = map(test_tiles) do tile
        elements = test_elements[tile]
        assembly_data = BEAST.AssemblyData(test_assembly.data[:, :, tile])
        assemble_primer_gpu(test_functions, elements, assembly_data,
            quadstrat.outer_rule, T)
    end
    trial_data = map(trial_tiles) do tile
        elements = trial_elements[tile]
        assembly_data = BEAST.AssemblyData(trial_assembly.data[:, :, tile])
        assemble_primer_gpu(trial_functions, elements, assembly_data,
            quadstrat.inner_rule, T)
    end

    singular_rules = nothing
    if quadstrat isa BEAST.DoubleNumSauterQstrat
        coordinate_type = coordtype(eltype(test_elements))
        singular_rules = (
            common_vertex=gpu_legendre_rule(
                quadstrat.sauter_schwab_common_vert, coordinate_type),
            common_edge=gpu_legendre_rule(
                quadstrat.sauter_schwab_common_edge, coordinate_type),
            common_face=gpu_legendre_rule(
                quadstrat.sauter_schwab_common_face, coordinate_type),
        )
    end

    CUDA.synchronize()
    return GPUDenseAssemblyData(
        CUDA.device(), test_data, trial_data, singular_rules)
end

function assemble_gpu!(operator::IntegralOperator,
    test_functions::Space, trial_functions::Space, store,
    quadstrat::Union{BEAST.DoubleNumQStrat,BEAST.DoubleNumSauterQstrat};
    gpu_blocksize,
    tilingstrat=nothing,
    nstreams::Int=1,
    devices=(CUDA.device(),))

    nstreams > 0 || throw(ArgumentError("number of CUDA streams must be positive"))
    devices = collect(devices)
    isempty(devices) && throw(ArgumentError("at least one CUDA device is required"))
    length(unique(devices)) == length(devices) ||
        throw(ArgumentError("CUDA devices must be unique"))
    if isnothing(tilingstrat)
        num_workers = length(devices) * nstreams
        tilingstrat = TilingStrategy(
            EqualTiling(num_workers), EqualTiling(1))
    end

    T = scalartype(operator, test_functions, trial_functions)

    test_geometry = geometry(test_functions)
    trial_geometry = geometry(trial_functions)
    test_tiles = tileindices(numcells(test_geometry), tilingstrat[1])
    trial_tiles = tileindices(numcells(trial_geometry), tilingstrat[2])

    test_elements, test_assembly, _ =
        BEAST.assemblydata(test_functions; onlyactives=false)
    trial_elements, trial_assembly, _ =
        BEAST.assemblydata(trial_functions; onlyactives=false)

    data_tasks = map(devices) do device
        Threads.@spawn gpu_dense_assembly_data(device,
            test_functions, test_elements, test_assembly, test_tiles,
            trial_functions, trial_elements, trial_assembly, trial_tiles,
            quadstrat, T)
    end
    device_data = fetch.(data_tasks)

    num_test_tiles = length(test_tiles)
    num_trial_tiles = length(trial_tiles)
    num_tile_pairs = num_test_tiles * num_trial_tiles
    num_workers = min(length(device_data) * nstreams, num_tile_pairs)

    result_type = Tuple{Matrix{T},Vector{Int},Vector{Int}}
    results = Channel{result_type}(2 * num_workers)
    producer = Threads.@spawn begin
        try
            @sync for worker in 1:num_workers
                Threads.@spawn let worker=worker
                    data = device_data[mod(worker - 1, length(device_data)) + 1]
                    CUDA.device!(data.device)
                    for tile_pair in worker:num_workers:num_tile_pairs
                        test_tile = mod(tile_pair - 1, num_test_tiles) + 1
                        trial_tile = div(tile_pair - 1, num_test_tiles) + 1

                        test_local_to_global,
                            (test_tile_elements, test_tile_assembly),
                            (test_quadrule, test_shapes) = data.test_data[test_tile]
                        trial_local_to_global,
                            (trial_tile_elements, trial_tile_assembly),
                            (trial_quadrule, trial_shapes) = data.trial_data[trial_tile]

                        matrix = assemblechunk_body_gpu!(operator,
                            refspace(test_functions),
                            test_tile_elements, test_tile_assembly,
                            test_quadrule, test_shapes,
                            refspace(trial_functions),
                            trial_tile_elements, trial_tile_assembly,
                            trial_quadrule, trial_shapes, quadstrat;
                            gpu_blocksize, singular_rules=data.singular_rules)

                        put!(results,
                            (matrix, test_local_to_global, trial_local_to_global))
                    end
                end
            end
        finally
            close(results)
        end
    end

    for (matrix, test_local_to_global, trial_local_to_global) in results
        for trial_dof in axes(matrix, 2)
            for test_dof in axes(matrix, 1)
                store(matrix[test_dof, trial_dof],
                    test_local_to_global[test_dof], trial_local_to_global[trial_dof])
            end
        end
    end
    fetch(producer)

    return nothing
end

function assemble_gpu!(operator::IntegralOperator,
    test_functions::Space, trial_functions::Space, store, quadstrat;
    kwargs...)

    throw(ArgumentError(
        "GPU assembly currently supports DoubleNumQStrat and DoubleNumSauterQstrat, got $(typeof(quadstrat))"))
end

function assemble!(operator::IntegralOperator,
    test_functions::Space, trial_functions::Space, store,
    threading::Type{Threading{:gpu}};
    quadstrat=BEAST.defaultquadstrat,
    gpu_blocksize=256,
    tilingstrat=nothing,
    nstreams::Int=1,
    devices=(CUDA.device(),),
    kwargs...)

    numfunctions(test_functions) == 0 && return
    numfunctions(trial_functions) == 0 && return

    strategy = resolve_gpu_quadstrat(
        quadstrat, operator, test_functions, trial_functions)
    return assemble_gpu!(operator, test_functions, trial_functions, store, strategy;
        gpu_blocksize, tilingstrat, nstreams, devices)
end
