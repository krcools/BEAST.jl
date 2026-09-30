struct GPUBlockSpaceData{S,E,Q,V,T}
    space::S
    elements::E
    quadrule::Q
    shapes::V
    num_shapes::Int
    coefficient_type::Type{T}
end

struct GPUSauterSchwabRules{V,E,F}
    common_vertex::V
    common_edge::E
    common_face::F
end

struct GPUBlockAssembler{O,X,Y,TX,TY,Q,R,D}
    operator::O
    test_space::X
    trial_space::Y
    test_data::TX
    trial_data::TY
    quadstrat::Q
    singular_rules::R
    device::D
    gpu_blocksize::Int
end

function GPUBlockSpaceData(space::Space, quadrule, ::Type{T}) where T
    elements, assembly_data, _ = BEAST.assemblydata(space; onlyactives=false)
    _, (elements_d, _), (quadrule_d, shapes_d) =
        assemble_primer_gpu(space, elements, assembly_data, quadrule, T)
    num_shapes = numfunctions(refspace(space), domain(first(elements)))

    return GPUBlockSpaceData(
        space, elements_d, quadrule_d, shapes_d, num_shapes, T)
end

function block_assemblydata_gpu(data::GPUBlockSpaceData, ids, element_ids)
    element_map = Dict(element_id => i for (i, element_id) in enumerate(element_ids))
    rows = Int[]
    cols = Int[]
    values = data.coefficient_type[]

    for (block_dof, dof) in enumerate(ids)
        for shape in BEAST.basisfunction(data.space, dof)
            element = element_map[shape.cellid]
            push!(rows, block_dof)
            push!(cols, data.num_shapes * (element - 1) + shape.refid)
            push!(values, data.coefficient_type(shape.coeff))
        end
    end

    assembly_data = sparse(rows, cols, values,
        length(ids), data.num_shapes * length(element_ids))
    return CuSparseMatrixCSC(assembly_data)
end

function active_block_data(data::GPUBlockSpaceData, ids)
    element_ids = BEAST.active_element_ids(data.space, ids)
    element_ids_d = CuArray(element_ids)
    elements = data.elements[element_ids_d]
    shapes = data.shapes[element_ids_d, :]
    assembly_data = block_assemblydata_gpu(data, ids, element_ids)

    return elements, assembly_data, data.quadrule, shapes
end

function gpu_blockassembler(operator::IntegralOperator,
    test_space::Space, trial_space::Space;
    quadstrat=BEAST.defaultquadstrat,
    gpu_blocksize=256,
    device=CUDA.device())

    strategy = resolve_gpu_quadstrat(
        quadstrat, operator, test_space, trial_space)
    CUDA.device!(device)
    T = scalartype(operator, test_space, trial_space)
    test_data = GPUBlockSpaceData(test_space, strategy.outer_rule, T)
    trial_data = GPUBlockSpaceData(trial_space, strategy.inner_rule, T)

    singular_rules = nothing
    if strategy isa BEAST.DoubleNumSauterQstrat
        coordinate_type = coordtype(eltype(test_data.elements))
        singular_rules = GPUSauterSchwabRules(
            gpu_legendre_rule(
                strategy.sauter_schwab_common_vert, coordinate_type),
            gpu_legendre_rule(
                strategy.sauter_schwab_common_edge, coordinate_type),
            gpu_legendre_rule(
                strategy.sauter_schwab_common_face, coordinate_type))
    end

    return GPUBlockAssembler(operator, test_space, trial_space,
        test_data, trial_data, strategy, singular_rules, device, gpu_blocksize)
end

function (assembler::GPUBlockAssembler)(block, test_ids, trial_ids)
    size(block) == (length(test_ids), length(trial_ids)) ||
        throw(DimensionMismatch("destination size does not match block indices"))

    if isempty(test_ids) || isempty(trial_ids)
        fill!(block, zero(eltype(block)))
        return block
    end

    CUDA.device!(assembler.device)
    test_ids = collect(Int, test_ids)
    trial_ids = collect(Int, trial_ids)
    test_elements, test_assembly, test_quadrule, test_shapes =
        active_block_data(assembler.test_data, test_ids)
    trial_elements, trial_assembly, trial_quadrule, trial_shapes =
        active_block_data(assembler.trial_data, trial_ids)

    matrix = assemblechunk_body_gpu_device!(assembler.operator,
        refspace(assembler.test_space),
        test_elements, test_assembly, test_quadrule, test_shapes,
        refspace(assembler.trial_space),
        trial_elements, trial_assembly, trial_quadrule, trial_shapes,
        assembler.quadstrat;
        gpu_blocksize=assembler.gpu_blocksize,
        singular_rules=assembler.singular_rules)
    copyto!(block, matrix)

    return block
end
