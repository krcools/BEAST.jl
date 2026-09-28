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

function load_assemblydata_gpu(functions::Space, ::Type{T}) where T
    elements, assembly_data, _ = BEAST.assemblydata(functions)
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

function gpu_triangle_rule(order, ::Type{T}) where T
    points, weights = CompScienceMeshes.trgauss(order)
    rule = Vector{Tuple{SVector{2,T},T}}(undef, length(weights))
    for i in eachindex(weights)
        point = SVector{2,T}(points[1,i], points[2,i])
        rule[i] = (point, T(weights[i]))
    end
    return CuArray(rule)
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

function gpu_momintegral_doublenum_allpairs!(zlocal, operator,
    test_elements, trial_elements, test_shapes, trial_shapes,
    test_local_space, trial_local_space, test_quadrule, trial_quadrule)

    pair = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    num_test_elements = length(test_elements)
    num_trial_elements = length(trial_elements)

    if pair <= num_test_elements * num_trial_elements
        test_index = mod(pair - 1, num_test_elements) + 1
        trial_index = div(pair - 1, num_test_elements) + 1

        test_element = test_elements[test_index]
        trial_element = trial_elements[trial_index]
        test_domain = domain(test_element)
        trial_domain = domain(trial_element)
        num_test_shapes = numfunctions(test_local_space, test_domain)
        num_trial_shapes = numfunctions(trial_local_space, trial_domain)

        integrand = BEAST.Integrand(operator, test_local_space, trial_local_space,
            test_element, trial_element)

        T = eltype(zlocal)
        z = zeros(SMatrix{num_test_shapes,num_trial_shapes,T})
        for test_quadrature_index in eachindex(test_quadrule)
            test_point, test_weight = test_quadrule[test_quadrature_index]
            x = neighborhood(test_element, test_point)
            weighted_test_jacobian = test_weight * jacobian(x)
            test_values = test_shapes[test_index, test_quadrature_index]

            for trial_quadrature_index in eachindex(trial_quadrule)
                trial_point, trial_weight = trial_quadrule[trial_quadrature_index]
                y = neighborhood(trial_element, trial_point)
                weight = weighted_test_jacobian * trial_weight * jacobian(y)
                trial_values = trial_shapes[trial_index, trial_quadrature_index]

                z += weight * integrand(x, y, test_values, trial_values)
            end
        end

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
    trial_projection = CUDA.zeros(T, size(trial_assembly, 1), size(zlocal, 1))

    CUSPARSE.mm!('N', 'T', one(T), trial_assembly, zlocal,
        zero(T), trial_projection, 'O')
    CUSPARSE.mm!('N', 'T', one(eltype(matrix)), test_assembly, trial_projection,
        zero(eltype(matrix)), matrix, 'O')

    return nothing
end

function assemble_primer_gpu(functions::Space, quadrule, ::Type{T}) where T
    local_space = refspace(functions)
    elements, assembly_data, local_to_global, element_domain, coordinate_type =
        load_assemblydata_gpu(functions, T)

    eltype(elements) <: CompScienceMeshes.Simplex{3,2} ||
        throw(ArgumentError(
            "GPU assembly currently supports triangular surface elements in three dimensions"))

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

function assemblechunk_body_gpu!(operator::IntegralOperator,
    test_local_space, test_elements, test_assembly, test_quadrule, test_shapes,
    trial_local_space, trial_elements, trial_assembly, trial_quadrule, trial_shapes;
    gpu_blocksize)

    num_test_shapes = div(size(test_assembly, 2), length(test_elements))
    num_trial_shapes = div(size(trial_assembly, 2), length(trial_elements))
    T = promote_type(scalartype(operator),
        eltype(test_assembly), eltype(trial_assembly))

    zlocal = CUDA.zeros(T,
        num_test_shapes * length(test_elements),
        num_trial_shapes * length(trial_elements))

    num_pairs = length(test_elements) * length(trial_elements)
    launch_gpu_kernel!(gpu_momintegral_doublenum_allpairs!, zlocal, operator,
        test_elements, trial_elements, test_shapes, trial_shapes,
        test_local_space, trial_local_space, test_quadrule, trial_quadrule;
        gpu_blocksize=gpu_blocksize, problem_size=num_pairs)

    matrix = CUDA.zeros(T, size(test_assembly, 1), size(trial_assembly, 1))
    build_matrix!(matrix, zlocal, test_assembly, trial_assembly)

    return Array(matrix)
end

function assemble_gpu!(operator::IntegralOperator,
    test_functions::Space, trial_functions::Space, store,
    quadstrat::BEAST.DoubleNumQStrat; gpu_blocksize)

    T = scalartype(operator, test_functions, trial_functions)

    test_local_to_global, (test_elements, test_assembly),
        (test_quadrule, test_shapes) =
            assemble_primer_gpu(test_functions, quadstrat.outer_rule, T)

    trial_local_to_global, (trial_elements, trial_assembly),
        (trial_quadrule, trial_shapes) =
            assemble_primer_gpu(trial_functions, quadstrat.inner_rule, T)

    matrix = assemblechunk_body_gpu!(operator,
        refspace(test_functions), test_elements, test_assembly,
        test_quadrule, test_shapes,
        refspace(trial_functions), trial_elements, trial_assembly,
        trial_quadrule, trial_shapes;
        gpu_blocksize)

    for trial_dof in axes(matrix, 2)
        for test_dof in axes(matrix, 1)
            store(matrix[test_dof, trial_dof],
                test_local_to_global[test_dof], trial_local_to_global[trial_dof])
        end
    end

    return nothing
end

function assemble_gpu!(operator::IntegralOperator,
    test_functions::Space, trial_functions::Space, store, quadstrat;
    gpu_blocksize)

    throw(ArgumentError(
        "GPU assembly currently supports DoubleNumQStrat, got $(typeof(quadstrat))"))
end

function assemble!(operator::IntegralOperator,
    test_functions::Space, trial_functions::Space, store,
    threading::Type{Threading{:gpu}};
    quadstrat=BEAST.defaultquadstrat, gpu_blocksize=256, kwargs...)

    strategy = quadstrat(operator, test_functions, trial_functions)
    return assemble_gpu!(operator, test_functions, trial_functions, store, strategy;
        gpu_blocksize)
end
