using .Variational



mutable struct DiscreteEquation
  equation
  trial_space_dict # dictionary mapping indices into trial space to FE spaces
  test_space_dict  # dictionary mapping indices into test space to FE spaces
end

struct DiscreteBilform
    bilform
    trial_space_dict # dictionary mapping indices into trial space to FE spaces
    test_space_dict  # dictionary mapping indices into test space to FE spaces
end

struct DiscreteLinform
    linform
    test_space_dict
end

function _expand_space_mappings(sms)
    esms = []
    for sm in sms
        if first(sm) isa Vector
            j = first(sm)
            X = last(sm)
            @assert X isa BEAST.DirectProductSpace
            append!(esms, [(ji => Xi) for (ji,Xi) in zip(j,X.factors)])
        else
            append!(esms, [sm])
        end
    end
    return esms
end


function discretise(bf::BilForm, space_mappings::Pair...)

    space_mappings = _expand_space_mappings(space_mappings)

    trial_space_dict = Dict()
    test_space_dict = Dict()
    for sm in space_mappings

        found = false
        sm.first.space == bf.trial_space && (dict = trial_space_dict; found = true)
        sm.first.space == bf.test_space  && (dict = test_space_dict;  found = true)
        @assert found "Vector $(sm.first) neither in test nor in trial space"

        @assert !haskey(dict, sm.first.idx) "multiple mappings for $(sm.first)"
        dict[sm.first.idx] = sm.second
    end

    # check that all symbols where mapped
    for p in eachindex(bf.trial_space) @assert haskey(trial_space_dict,p) end
    for p in eachindex(bf.test_space)  @assert haskey(test_space_dict, p) end

    DiscreteBilform(bf, trial_space_dict, test_space_dict)
end


function discretise(lf::LinForm, space_mappings::Pair...)
    space_mappings = _expand_space_mappings(space_mappings)
    test_space_dict = Dict()
    for sm in space_mappings

        found = false
        # sm.first.space == bf.trial_space && (dict = trial_space_dict; found = true)
        sm.first.space == lf.test_space  && (dict = test_space_dict;  found = true)
        @assert found "Vector $(sm.first) not found in test space"

        @assert !haskey(dict, sm.first.idx) "multiple mappings for $(sm.first)"
        dict[sm.first.idx] = sm.second
    end

    # check that all symbols where mapped
    # for p in eachindex(bf.trial_space) @assert haskey(trial_space_dict,p) end
    for p in eachindex(lf.test_space)  @assert haskey(test_space_dict, p) end
 
    DiscreteLinform(lf, test_space_dict)
end


function discretise(eq, space_mappings::Pair...)
    space_mappings = _expand_space_mappings(space_mappings)

    trial_space_dict = Dict()
    test_space_dict = Dict()
    for sm in space_mappings

        found = false
        sm.first.space == eq.lhs.trial_space && (dict = trial_space_dict; found = true)
        sm.first.space == eq.lhs.test_space  && (dict = test_space_dict;  found = true)
        @assert found "Vector $(sm.first) neither in test nor in trial space"

        @assert !haskey(dict, sm.first.idx) "multiple mappings for $(sm.first)"
        dict[sm.first.idx] = sm.second
    end

    # check that all symbols where mapped
    for p in 1:length(eq.lhs.trial_space) @assert haskey(trial_space_dict,p) end
    for p in 1:length(eq.lhs.test_space)  @assert haskey(test_space_dict, p) end

    DiscreteEquation(eq, trial_space_dict, test_space_dict)
end


"""
    discr(eq, pairs...)

This macro provides syntactical sugar for the definition of a discretisation
of a varational formulation. Given a variational equation EQ: Find j ∈ X such
that for all k ∈ Y a(k,j) = f(k) can be discretised by stating:

    eq = @discretise EQ j∈X k∈Y
"""
macro discretise(eq, pairs...)
    r = :(BEAST.discretise($eq))
    for p in pairs
        x = p.args[2]
        X = p.args[3]
        push!(r.args, :($x=>$X))
    end
    return esc(r)
end


sysmatrix(eq::DiscreteEquation; materialize=BEAST.assemble) =
    assemble(eq.equation.lhs, eq.test_space_dict, eq.trial_space_dict, materialize=materialize)
rhs(eq::DiscreteEquation) = assemble(eq.equation.rhs, eq.test_space_dict)

assemble(dbf::DiscreteBilform; materialize=BEAST.assemble) = assemble(dbf.bilform, dbf.test_space_dict, dbf.trial_space_dict; materialize)
assemble(dlf::DiscreteLinform) = assemble(dlf.linform, dlf.test_space_dict)

function assemble(lform::LinForm, test_space_dict)
    T=ComplexF64
    terms = lform.terms
    #T = Complex{typeof(terms[1].functional.field.wavenumber)}
    #T = Complex{eltype(vertextype(test_space_dict[1].geo))}
    # I = Int[1]
    blocksizes1 = Int[]
    for p in 1:length(lform.test_space)
        X = test_space_dict[p]
        # push!(I, last(I) + numfunctions(X))
        push!(blocksizes1, numfunctions(X))
    end

    Z = zeros(T, sum(blocksizes1))
    B = PseudoBlockArray{T}(Z, blocksizes1)

    for t in terms

        α = t.coeff
        a = t.functional
        m = t.test_id
        X = test_space_dict[m]
        o = t.test_ops

        # act with the various ops on X
        for op in reverse(o)
            Y = X;
            X = op[end](op[1:end-1]..., Y)
        end

        b = assemble(a, X)
        # B[I[m] : I[m+1]-1] = α * b
        B[Block(m)] = α * b
    end

    return B
end


struct SpaceTimeData{T} <: AbstractArray{Vector{T},1}
    data::Array{T,2}
end

Base.eltype(x::SpaceTimeData{T}) where {T} = Vector{T}
Base.size(x::SpaceTimeData) = (size(x.data)[1],)
Base.getindex(x::SpaceTimeData, i::Int) = x.data[i,:]

function td_assemble(lform::LinForm, test_space_dict)

    terms = lform.terms

    T = Float32
    for term in lform.terms
        T = scalartype(T,term.coeff)
        T = scalartype(T,term.functional)
    end
    for kv in test_space_dict;  T = scalartype(T,kv[2]) end
    # @show T

    I = [numfunctions(spatialbasis(test_space_dict[i])) for i in 1:length(lform.test_space)]

    M = zeros(Int, length(test_space_dict))
    for (p,x) in test_space_dict M[p]=numfunctions(spatialbasis(x)) end

    N = [numfunctions(temporalbasis(test_space_dict[1]))]
    B = BlockArray{T}(undef, M, N)

    # row_axis = BlockArrays.blockedrange(M)
    # col_axis = BlockArrays.blockedrange(N)

    # BT = SpaceTimeData{T}
    # B = BlockArray(undef_blocks, BT, I)

    for t in terms

        α = t.coeff
        a = t.functional
        m = t.test_id
        X = test_space_dict[m]
        o = t.test_ops
        

        # act with the various ops on X
        for op in reverse(o)
            Y = X;
            X = op[end](op[1:end-1]..., Y)
        end

        b = assemble(a, X)
        B[Block(m),Block(1)] = α*b
    end

    return B
end

# function assemble_hide(bilform::BilForm, test_space_dict::Space{U}, trial_space_dict) where U

#   lhterms = bilform.terms
#   T = Complex{U} # TDOD: Fix this

#   blocksizes1 = Int[]
#   for p in 1:length(bilform.test_space)
#     X = test_space_dict[p]
#     push!(blocksizes1, numfunctions(X))
#   end

#   blocksizes2 = Int[]
#   for q in 1:length(bilform.trial_space)
#     Y = trial_space_dict[q]
#     push!(blocksizes2, numfunctions(Y))
#   end

#   # allocate the memory for the matrices
#   A = zeros(T, sum(blocksizes1), sum(blocksizes2))
#   Z = PseudoBlockArray{T}(A, blocksizes1, blocksizes2)
#   # For each block, compute the interaction matrix
#   for t in lhterms

#       α = t.coeff
#       a = t.kernel

#       m = t.test_id
#       x = test_space_dict[m]
#       for op in reverse(t.test_ops)
#           x = op[end](op[1:end-1]..., x)
#       end

#       n = t.trial_id
#       y = trial_space_dict[n]
#       for op in reverse(t.trial_ops)
#           y = op[end](op[1:end-1]..., y)
#       end

#       z = assemble(a, x, y)
#     #   @show m n norm(z)
#       Z[Block(m,n)] += α * z
#   end

#   return Z
# end


function assemble(bilform::BilForm, test_space_dict, trial_space_dict;
    materialize=BEAST.assemble)

    @assert !isempty(bilform.terms)

    M = zeros(Int, length(test_space_dict))
    N = zeros(Int, length(trial_space_dict))

    for (p,x) in test_space_dict M[p]=numfunctions(x) end
    for (p,x) in trial_space_dict N[p]=numfunctions(x) end

    U = BlockArrays.blockedrange(M)
    V = BlockArrays.blockedrange(N)

    Z = ZeroMap{Float32}(U, V)

    for term in bilform.terms
    
        x = test_space_dict[term.test_id]
        for op in reverse(term.test_ops)
            x = op[end](op[1:end-1]..., x)
        end

        y = trial_space_dict[term.trial_id]
        for op in reverse(term.trial_ops)
            y = op[end](op[1:end-1]..., y)
        end

        z = materialize(term.kernel, x, y)
        
        I = Block(term.test_id)
        J = Block(term.trial_id)
        zlifted = LiftedMap(z,I,J,U,V)

        Z = Z + term.coeff * zlifted
    end

    return Z
end


function assemble(bf::BilForm, X::DirectProductSpace, Y::DirectProductSpace)

    test_space_dict = Dict(enumerate(X.factors))
    trial_space_dict = Dict(enumerate(Y.factors))

    z = assemble(bf, test_space_dict, trial_space_dict)
end



# function td_assemble(bilform::BilForm, test_space_dict, trial_space_dict)

#     lhterms = bilform.terms

#     T = Float32
#     for term in bilform.terms
#         T = scalartype(T,term.coeff)
#         T = scalartype(T,term.kernel)
#     end
#     for kv in test_space_dict;  T = scalartype(T,kv[2]) end
#     for kv in trial_space_dict; T = scalartype(T,kv[2]) end

#     I = [numfunctions(spatialbasis(test_space_dict[i])) for i in 1:length(bilform.test_space)]
#     J = [numfunctions(spatialbasis(trial_space_dict[i])) for i in 1:length(bilform.trial_space)]

#     Z = BlockArray{Vector{T}}(zero_blocks, I, J)

#     # For each block, compute the interaction matrix
#     for t in lhterms

#         @show (t.coeff,t.kernel)
#         a = t.coeff * t.kernel

#         m = t.test_id
#         x = test_space_dict[m]
#         for op in reverse(t.test_ops)
#             x = op[end](op[1:end-1]..., x)
#         end

#         n = t.trial_id
#         y = trial_space_dict[n]
#         for op in reverse(t.trial_ops)
#             y = op[end](op[1:end-1]..., y)
#         end

#         z = assemble(a, x, y)
#         @warn "variation formulations where combinations of test and trial space recur multiple times are not supported!"
#         Z[Block(m,n)] = ConvolutionOperators.ConvOpAsMatrix(z)
 
#     end
#     return Z
# end


function td_assemble(bilform::BilForm, test_space_dict, trial_space_dict)

    lhterms = bilform.terms

    # T = Float32
    # for term in bilform.terms
    #     T = scalartype(T,term.coeff)
    #     T = scalartype(T,term.kernel)
    # end
    # for kv in test_space_dict;  T = scalartype(T,kv[2]) end
    # for kv in trial_space_dict; T = scalartype(T,kv[2]) end

    # I = [numfunctions(spatialbasis(test_space_dict[i])) for i in 1:length(bilform.test_space)]
    # J = [numfunctions(spatialbasis(trial_space_dict[i])) for i in 1:length(bilform.trial_space)]

    # rowaxis = BlockArrays.blockedrange([numfunctions(spatialbasis(test_space_dict[i])) for i in 1:length(test_space_dict)])
    # colaxis = BlockArrays.blockedrange([numfunctions(spatialbasis(trial_space_dict[i])) for i in 1:length(trial_space_dict)])

    M = zeros(Int, length(test_space_dict))
    N = zeros(Int, length(trial_space_dict))

    for (p,x) in test_space_dict M[p]=numfunctions(spatialbasis(x)) end
    for (p,x) in trial_space_dict N[p]=numfunctions(spatialbasis(x)) end

    row_axis = BlockArrays.blockedrange(M)
    col_axis = BlockArrays.blockedrange(N)

    # Z = BlockArray{Vector{T}}(zero_blocks, row_axis, col_axis)
    Z = ConvolutionOperators.ZeroConvOp(row_axis, col_axis)

    for t in lhterms

        a = t.coeff * t.kernel

        m = t.test_id
        x = test_space_dict[m]
        for op in reverse(t.test_ops)
            x = op[end](op[1:end-1]..., x)
        end

        n = t.trial_id
        y = trial_space_dict[n]
        for op in reverse(t.trial_ops)
            y = op[end](op[1:end-1]..., y)
        end

        z = assemble(a, x, y)
        Z += ConvolutionOperators.LiftedConvOp(z, row_axis, col_axis, Block(m), Block(n))
        # Z[Block(m,n)] = ConvolutionOperators.ConvOpAsMatrix(z)
 
    end
    return Z
end
