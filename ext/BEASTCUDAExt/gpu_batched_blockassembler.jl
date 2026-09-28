@inline function locate_batched_block(index, offsets, num_blocks)
    lower = 1
    upper = num_blocks
    while lower < upper
        middle = (lower + upper + 1) >>> 1
        if offsets[middle] < index
            lower = middle
        else
            upper = middle - 1
        end
    end
    return lower
end

function gpu_batched_classify!(labels, total_pairs,
    pair_offsets, test_offsets, trial_offsets,
    test_elements, trial_elements, num_blocks,
    all_test_elements::CuDeviceVector{S}, all_trial_elements) where S

    pair = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if pair <= total_pairs
        block = locate_batched_block(pair, pair_offsets, num_blocks)
        local_pair = pair - pair_offsets[block] - 1
        num_test_elements = test_offsets[block + 1] - test_offsets[block]
        test_local = mod(local_pair, num_test_elements)
        trial_local = div(local_pair, num_test_elements)
        test_index = test_elements[test_offsets[block] + test_local + 1]
        trial_index = trial_elements[trial_offsets[block] + trial_local + 1]

        test_element = all_test_elements[test_index]
        trial_element = all_trial_elements[trial_index]
        tolerance = 1.0e3 * eps(coordtype(S))
        hits = 0
        for test_vertex in vertices(test_element)
            for trial_vertex in vertices(trial_element)
                hits += norm(test_vertex - trial_vertex) < tolerance
            end
        end
        labels[pair] = UInt8(hits)
    end

    return nothing
end

function gpu_batched_integrate_regular!(zstage, operator,
    pair_indices, num_pairs, pair_offsets, test_offsets, trial_offsets,
    test_elements, trial_elements, num_blocks,
    all_test_elements::CuDeviceVector{S}, all_trial_elements,
    test_shapes, trial_shapes, test_local_space, trial_local_space,
    test_quadrule, trial_quadrule) where S

    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= num_pairs
        pair = pair_indices[index]
        block = locate_batched_block(pair, pair_offsets, num_blocks)
        local_pair = pair - pair_offsets[block] - 1
        num_test_elements = test_offsets[block + 1] - test_offsets[block]
        test_local = mod(local_pair, num_test_elements)
        trial_local = div(local_pair, num_test_elements)
        test_index = test_elements[test_offsets[block] + test_local + 1]
        trial_index = trial_elements[trial_offsets[block] + trial_local + 1]

        test_element = all_test_elements[test_index]
        trial_element = all_trial_elements[trial_index]
        num_test_shapes = numfunctions(test_local_space, domain(test_element))
        num_trial_shapes = numfunctions(trial_local_space, domain(trial_element))
        integrand = BEAST.Integrand(operator,
            test_local_space, trial_local_space, test_element, trial_element)
        z = gpu_doublenum_integral(integrand,
            test_element, trial_element,
            view(test_shapes, test_index, :),
            view(trial_shapes, trial_index, :),
            test_local_space, trial_local_space, test_quadrule, trial_quadrule,
            eltype(zstage))

        offset = (pair - 1) * num_test_shapes * num_trial_shapes
        for trial_shape in 1:num_trial_shapes
            for test_shape in 1:num_test_shapes
                zstage[offset + (trial_shape - 1) * num_test_shapes + test_shape] =
                    z[test_shape, trial_shape]
            end
        end
    end

    return nothing
end

function gpu_batched_integrate_far!(zstage, operator, total_pairs,
    pair_offsets, test_offsets, trial_offsets,
    test_elements, trial_elements, num_blocks,
    all_test_elements::CuDeviceVector{S}, all_trial_elements,
    test_shapes, trial_shapes, test_local_space, trial_local_space,
    test_quadrule, trial_quadrule) where S

    pair = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if pair <= total_pairs
        block = locate_batched_block(pair, pair_offsets, num_blocks)
        local_pair = pair - pair_offsets[block] - 1
        num_test_elements = test_offsets[block + 1] - test_offsets[block]
        test_local = mod(local_pair, num_test_elements)
        trial_local = div(local_pair, num_test_elements)
        test_index = test_elements[test_offsets[block] + test_local + 1]
        trial_index = trial_elements[trial_offsets[block] + trial_local + 1]

        test_element = all_test_elements[test_index]
        trial_element = all_trial_elements[trial_index]
        num_test_shapes = numfunctions(test_local_space, domain(test_element))
        num_trial_shapes = numfunctions(trial_local_space, domain(trial_element))
        integrand = BEAST.Integrand(operator,
            test_local_space, trial_local_space, test_element, trial_element)
        z = gpu_doublenum_integral(integrand,
            test_element, trial_element,
            view(test_shapes, test_index, :),
            view(trial_shapes, trial_index, :),
            test_local_space, trial_local_space, test_quadrule, trial_quadrule,
            eltype(zstage))

        offset = (pair - 1) * num_test_shapes * num_trial_shapes
        for trial_shape in 1:num_trial_shapes
            for test_shape in 1:num_test_shapes
                zstage[offset + (trial_shape - 1) * num_test_shapes + test_shape] =
                    z[test_shape, trial_shape]
            end
        end
    end

    return nothing
end

function gpu_batched_integrate_singular!(zstage, operator,
    pair_indices, labels, num_pairs,
    pair_offsets, test_offsets, trial_offsets,
    test_elements, trial_elements, num_blocks,
    all_test_elements::CuDeviceVector{S}, all_trial_elements,
    test_local_space, trial_local_space,
    common_vertex, common_edge, common_face) where S

    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= num_pairs
        pair = pair_indices[index]
        block = locate_batched_block(pair, pair_offsets, num_blocks)
        local_pair = pair - pair_offsets[block] - 1
        num_test_elements = test_offsets[block + 1] - test_offsets[block]
        test_local = mod(local_pair, num_test_elements)
        trial_local = div(local_pair, num_test_elements)
        test_index = test_elements[test_offsets[block] + test_local + 1]
        trial_index = trial_elements[trial_offsets[block] + trial_local + 1]

        test_element = all_test_elements[test_index]
        trial_element = all_trial_elements[trial_index]
        num_test_shapes = numfunctions(test_local_space, domain(test_element))
        num_trial_shapes = numfunctions(trial_local_space, domain(trial_element))
        integrand = BEAST.Integrand(operator,
            test_local_space, trial_local_space, test_element, trial_element)

        hits = labels[pair]
        if hits == 0x01
            z = gpu_sauterschwab_integral(integrand,
                test_element, trial_element, common_vertex,
                test_local_space, trial_local_space, eltype(zstage))
        elseif hits == 0x02
            z = gpu_sauterschwab_integral(integrand,
                test_element, trial_element, common_edge,
                test_local_space, trial_local_space, eltype(zstage))
        else
            z = gpu_sauterschwab_integral(integrand,
                test_element, trial_element, common_face,
                test_local_space, trial_local_space, eltype(zstage))
        end

        offset = (pair - 1) * num_test_shapes * num_trial_shapes
        for trial_shape in 1:num_trial_shapes
            for test_shape in 1:num_test_shapes
                zstage[offset + (trial_shape - 1) * num_test_shapes + test_shape] =
                    z[test_shape, trial_shape]
            end
        end
    end

    return nothing
end

function gpu_batched_project!(output, total_cells,
    cell_offsets, row_offsets, column_offsets,
    pair_offsets, test_offsets, num_blocks, num_test_shapes, num_trial_shapes,
    row_entry_offsets, row_entry_elements, row_entry_shapes, row_entry_coefficients,
    column_entry_offsets, column_entry_elements, column_entry_shapes,
    column_entry_coefficients, zstage::CuDeviceVector{T}) where T

    cell = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if cell <= total_cells
        block = locate_batched_block(cell, cell_offsets, num_blocks)
        local_cell = cell - cell_offsets[block] - 1
        num_rows = row_offsets[block + 1] - row_offsets[block]
        local_row = mod(local_cell, num_rows)
        local_column = div(local_cell, num_rows)
        row = row_offsets[block] + local_row + 1
        column = column_offsets[block] + local_column + 1
        num_test_elements = test_offsets[block + 1] - test_offsets[block]

        value = zero(T)
        for row_entry in
            (row_entry_offsets[row] + 1):row_entry_offsets[row + 1]
            test_local = row_entry_elements[row_entry]
            test_shape = row_entry_shapes[row_entry]
            test_coefficient = row_entry_coefficients[row_entry]

            for column_entry in
                (column_entry_offsets[column] + 1):column_entry_offsets[column + 1]
                trial_local = column_entry_elements[column_entry]
                trial_shape = column_entry_shapes[column_entry]
                trial_coefficient = column_entry_coefficients[column_entry]
                pair = pair_offsets[block] + test_local - 1 +
                    (trial_local - 1) * num_test_elements + 1
                index = (pair - 1) * num_test_shapes * num_trial_shapes +
                    (trial_shape - 1) * num_test_shapes + test_shape
                value += test_coefficient * trial_coefficient * zstage[index]
            end
        end
        output[cell_offsets[block] + local_column * num_rows + local_row + 1] = value
    end

    return nothing
end

function batched_entry_table(space, ids_per_block, dof_offsets,
    elements, element_offsets, ::Type{T}) where T

    num_blocks = length(ids_per_block)
    entry_offsets = zeros(Int32, dof_offsets[end] + 1)
    for block in 1:num_blocks
        for (local_dof, dof) in enumerate(ids_per_block[block])
            entry_offsets[dof_offsets[block] + local_dof + 1] =
                length(BEAST.basisfunction(space, dof))
        end
    end
    cumsum!(entry_offsets, entry_offsets)

    num_entries = entry_offsets[end]
    entry_elements = Vector{Int32}(undef, num_entries)
    entry_shapes = Vector{Int32}(undef, num_entries)
    entry_coefficients = Vector{T}(undef, num_entries)
    fill_positions = copy(entry_offsets)

    for block in 1:num_blocks
        block_elements = view(elements,
            element_offsets[block] + 1:element_offsets[block + 1])
        for (local_dof, dof) in enumerate(ids_per_block[block])
            slot = dof_offsets[block] + local_dof
            for shape in BEAST.basisfunction(space, dof)
                local_element = searchsortedfirst(block_elements, shape.cellid)
                entry = (fill_positions[slot] += 1)
                entry_elements[entry] = Int32(local_element)
                entry_shapes[entry] = Int32(shape.refid)
                entry_coefficients[entry] = T(shape.coeff)
            end
        end
    end

    return entry_offsets, entry_elements, entry_shapes, entry_coefficients
end

function batched_worklist(assembler::GPUBlockAssembler, test_ids, trial_ids)
    num_blocks = length(test_ids)
    test_elements = Int[]
    trial_elements = Int[]
    test_offsets = zeros(Int, num_blocks + 1)
    trial_offsets = zeros(Int, num_blocks + 1)
    pair_offsets = zeros(Int, num_blocks + 1)
    row_offsets = zeros(Int, num_blocks + 1)
    column_offsets = zeros(Int, num_blocks + 1)
    cell_offsets = zeros(Int, num_blocks + 1)

    for block in 1:num_blocks
        active_test = active_element_ids(assembler.test_space, test_ids[block])
        active_trial = active_element_ids(assembler.trial_space, trial_ids[block])
        append!(test_elements, active_test)
        append!(trial_elements, active_trial)
        test_offsets[block + 1] = test_offsets[block] + length(active_test)
        trial_offsets[block + 1] = trial_offsets[block] + length(active_trial)
        pair_offsets[block + 1] = pair_offsets[block] +
            length(active_test) * length(active_trial)
        row_offsets[block + 1] = row_offsets[block] + length(test_ids[block])
        column_offsets[block + 1] =
            column_offsets[block] + length(trial_ids[block])
        cell_offsets[block + 1] = cell_offsets[block] +
            length(test_ids[block]) * length(trial_ids[block])
    end

    T = assembler.test_data.coefficient_type
    row_entries = batched_entry_table(assembler.test_space, test_ids,
        row_offsets, test_elements, test_offsets, T)
    column_entries = batched_entry_table(assembler.trial_space, trial_ids,
        column_offsets, trial_elements, trial_offsets, T)

    return (; num_blocks, test_elements, trial_elements,
        test_offsets, trial_offsets, pair_offsets,
        row_offsets, column_offsets, cell_offsets,
        row_entry_offsets=row_entries[1], row_entry_elements=row_entries[2],
        row_entry_shapes=row_entries[3], row_entry_coefficients=row_entries[4],
        column_entry_offsets=column_entries[1],
        column_entry_elements=column_entries[2],
        column_entry_shapes=column_entries[3],
        column_entry_coefficients=column_entries[4])
end

function scatter_batched_blocks!(blocks::AbstractVector{<:AbstractMatrix},
    output, worklist)

    output_host = Array(output)
    for block in 1:worklist.num_blocks
        num_rows = worklist.row_offsets[block + 1] - worklist.row_offsets[block]
        num_columns = worklist.column_offsets[block + 1] -
            worklist.column_offsets[block]
        (iszero(num_rows) || iszero(num_columns)) && continue
        data = view(output_host,
            worklist.cell_offsets[block] + 1:worklist.cell_offsets[block + 1])
        copyto!(blocks[block], reshape(data, num_rows, num_columns))
    end
    return blocks
end

function scatter_batched_blocks!(blocks::AbstractVector{<:CUDA.AnyCuMatrix},
    output, worklist)

    for block in 1:worklist.num_blocks
        num_rows = worklist.row_offsets[block + 1] - worklist.row_offsets[block]
        num_columns = worklist.column_offsets[block + 1] -
            worklist.column_offsets[block]
        (iszero(num_rows) || iszero(num_columns)) && continue
        data = view(output,
            worklist.cell_offsets[block] + 1:worklist.cell_offsets[block + 1])
        copyto!(blocks[block], reshape(data, num_rows, num_columns))
    end
    return blocks
end

function gpu_batched_blockassemble_chunk!(blocks, assembler::GPUBlockAssembler,
    test_ids, trial_ids)

    worklist = batched_worklist(assembler, test_ids, trial_ids)
    num_test_shapes = assembler.test_data.num_shapes
    num_trial_shapes = assembler.trial_data.num_shapes
    T = scalartype(assembler.operator,
        assembler.test_space, assembler.trial_space)
    total_pairs = worklist.pair_offsets[end]
    total_cells = worklist.cell_offsets[end]

    zstage = CuArray{T}(undef,
        total_pairs * num_test_shapes * num_trial_shapes)
    output = CUDA.zeros(T, total_cells)

    pair_offsets = CuArray(worklist.pair_offsets)
    test_offsets = CuArray(worklist.test_offsets)
    trial_offsets = CuArray(worklist.trial_offsets)
    test_elements = CuArray(worklist.test_elements)
    trial_elements = CuArray(worklist.trial_elements)

    if isnothing(assembler.singular_rules)
        launch_gpu_kernel!(gpu_batched_integrate_far!, zstage,
            assembler.operator, total_pairs,
            pair_offsets, test_offsets, trial_offsets,
            test_elements, trial_elements, worklist.num_blocks,
            assembler.test_data.elements, assembler.trial_data.elements,
            assembler.test_data.shapes, assembler.trial_data.shapes,
            refspace(assembler.test_space), refspace(assembler.trial_space),
            assembler.test_data.quadrule, assembler.trial_data.quadrule;
            gpu_blocksize=128, problem_size=total_pairs)
    else
        labels = CuArray{UInt8}(undef, total_pairs)
        launch_gpu_kernel!(gpu_batched_classify!, labels, total_pairs,
            pair_offsets, test_offsets, trial_offsets,
            test_elements, trial_elements, worklist.num_blocks,
            assembler.test_data.elements, assembler.trial_data.elements;
            gpu_blocksize=128, problem_size=total_pairs)

        regular_pairs = findall(iszero, labels)
        singular_pairs = findall(!iszero, labels)
        launch_gpu_kernel!(gpu_batched_integrate_regular!, zstage,
            assembler.operator, regular_pairs, length(regular_pairs),
            pair_offsets, test_offsets, trial_offsets,
            test_elements, trial_elements, worklist.num_blocks,
            assembler.test_data.elements, assembler.trial_data.elements,
            assembler.test_data.shapes, assembler.trial_data.shapes,
            refspace(assembler.test_space), refspace(assembler.trial_space),
            assembler.test_data.quadrule, assembler.trial_data.quadrule;
            gpu_blocksize=128, problem_size=length(regular_pairs))

        common_vertex = CommonVertex(assembler.singular_rules.common_vertex)
        common_edge = CommonEdge(assembler.singular_rules.common_edge)
        common_face = CommonFace(assembler.singular_rules.common_face)
        launch_gpu_kernel!(gpu_batched_integrate_singular!, zstage,
            assembler.operator, singular_pairs, labels, length(singular_pairs),
            pair_offsets, test_offsets, trial_offsets,
            test_elements, trial_elements, worklist.num_blocks,
            assembler.test_data.elements, assembler.trial_data.elements,
            refspace(assembler.test_space), refspace(assembler.trial_space),
            common_vertex, common_edge, common_face;
            gpu_blocksize=128, problem_size=length(singular_pairs))
    end

    launch_gpu_kernel!(gpu_batched_project!, output, total_cells,
        CuArray(worklist.cell_offsets),
        CuArray(worklist.row_offsets), CuArray(worklist.column_offsets),
        pair_offsets, test_offsets, worklist.num_blocks,
        num_test_shapes, num_trial_shapes,
        CuArray(worklist.row_entry_offsets),
        CuArray(worklist.row_entry_elements),
        CuArray(worklist.row_entry_shapes),
        CuArray(worklist.row_entry_coefficients),
        CuArray(worklist.column_entry_offsets),
        CuArray(worklist.column_entry_elements),
        CuArray(worklist.column_entry_shapes),
        CuArray(worklist.column_entry_coefficients), zstage;
        gpu_blocksize=256, problem_size=total_cells)

    return scatter_batched_blocks!(blocks, output, worklist)
end

function batched_block_pair_counts(assembler::GPUBlockAssembler,
    test_ids, trial_ids)

    counts = Vector{Int}(undef, length(test_ids))
    for block in eachindex(test_ids)
        num_test_elements = length(active_element_ids(
            assembler.test_space, test_ids[block]))
        num_trial_elements = length(active_element_ids(
            assembler.trial_space, trial_ids[block]))
        counts[block] = num_test_elements * num_trial_elements
    end
    return counts
end

function batched_block_chunks(assembler::GPUBlockAssembler,
    test_ids, trial_ids; budget=1 << 30)

    budget > 0 || throw(ArgumentError("GPU block budget must be positive"))
    bytes_per_pair = assembler.test_data.num_shapes *
        assembler.trial_data.num_shapes *
        sizeof(scalartype(assembler.operator,
            assembler.test_space, assembler.trial_space))
    pair_counts = batched_block_pair_counts(assembler, test_ids, trial_ids)

    chunks = UnitRange{Int}[]
    first_block = 1
    bytes = 0
    for block in eachindex(pair_counts)
        block_bytes = pair_counts[block] * bytes_per_pair
        if bytes + block_bytes > budget && block > first_block
            push!(chunks, first_block:block - 1)
            first_block = block
            bytes = 0
        end
        bytes += block_bytes
    end
    first_block <= length(pair_counts) &&
        push!(chunks, first_block:length(pair_counts))
    return chunks
end

function gpu_batched_blockassemble!(blocks, assembler::GPUBlockAssembler,
    test_ids, trial_ids; budget=1 << 30)

    length(blocks) == length(test_ids) == length(trial_ids) ||
        throw(DimensionMismatch("blocks and index collections must have equal length"))
    for block in eachindex(blocks)
        size(blocks[block]) ==
            (length(test_ids[block]), length(trial_ids[block])) ||
            throw(DimensionMismatch(
                "destination size does not match block indices"))
    end

    CUDA.device!(assembler.device)
    for chunk in batched_block_chunks(assembler, test_ids, trial_ids; budget)
        gpu_batched_blockassemble_chunk!(view(blocks, chunk), assembler,
            view(test_ids, chunk), view(trial_ids, chunk))
    end
    return blocks
end

function gpu_batched_blockassemble(assembler::GPUBlockAssembler,
    test_ids, trial_ids; kwargs...)

    length(test_ids) == length(trial_ids) ||
        throw(DimensionMismatch("index collections must have equal length"))
    T = scalartype(assembler.operator,
        assembler.test_space, assembler.trial_space)
    blocks = [zeros(T, length(test_ids[block]), length(trial_ids[block]))
        for block in eachindex(test_ids)]
    return gpu_batched_blockassemble!(
        blocks, assembler, test_ids, trial_ids; kwargs...)
end
