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
        sm.first.space == lf.test_space  && (dict = test_space_dict;  found = true)
        @assert found "Vector $(sm.first) not found in test space"

        @assert !haskey(dict, sm.first.idx) "multiple mappings for $(sm.first)"
        dict[sm.first.idx] = sm.second
    end

    # check that all symbols where mapped
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


function _spacedict_to_directproductspace(spacedict)
    xfactors = Vector{AbstractSpace}(undef, length(spacedict))
    for (p,x) in spacedict xfactors[p] = x end
    X = DirectProductSpace(xfactors)
end

function assemble(lform::LinForm, test_space_dict)
    X = _spacedict_to_directproductspace(test_space_dict)
    return assemble(lform, X)
end

scalartype(lf::LinForm) = scalartype(lf.terms...)
scalartype(lt::LinTerm) = scalartype(lt.coeff, lt.functional)

function assemble(lform::LinForm, X::DirectProductSpace)

    @assert !isempty(lform.terms)

    M = length.(X.factors)
    U = NestedUnitRanges.nestedrange(X, 1, numfunctions)

    T = scalartype(lform, X)
    Z = zeros(T, numfunctions(X))
    B = PseudoBlockArray{T}(Z, (U,))

    for t in lform.terms

        m = t.test_id
        x = X.factors[m]

        for op in reverse(t.test_ops) x = op[end](op[1:end-1]..., x) end

        b = assemble(t.functional, x)
        B[Block(m)] = t.coeff * b
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

    I = [numfunctions(spatialbasis(test_space_dict[i])) for i in 1:length(lform.test_space)]

    M = zeros(Int, length(test_space_dict))
    for (p,x) in test_space_dict M[p]=numfunctions(spatialbasis(x)) end

    N = [numfunctions(temporalbasis(test_space_dict[1]))]
    B = BlockArray{T}(undef, M, N)

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

function assemble(bilform::BilForm, test_space_dict, trial_space_dict;
    materialize=BEAST.assemble)

    X = _spacedict_to_directproductspace(test_space_dict)
    Y = _spacedict_to_directproductspace(trial_space_dict)

    return assemble(bilform, X, Y; materialize)
end

lift(a,I,J,U,V) = LiftedMaps.LiftedMap(a,I,J,U,V)
lift(a::ConvolutionOperators.AbstractConvOp ,I,J,U,V) =
    ConvolutionOperators.LiftedConvOp(a, U, V, I, J)

function assemble(bf::BilForm, X::DirectProductSpace, Y::DirectProductSpace;
    materialize=BEAST.assemble)

    @assert !isempty(bf.terms)

    M = numfunctions.(spatialbasis(X).factors)
    N = numfunctions.(spatialbasis(Y).factors)

    U = BlockArrays.blockedrange(M)
    V = BlockArrays.blockedrange(N)

    sum(bf.terms) do term

        x = X.factors[term.test_id]
        for op in reverse(term.test_ops)
            x = op[end](op[1:end-1]..., x)
        end

        y = Y.factors[term.trial_id]
        for op in reverse(term.trial_ops)
            y = op[end](op[1:end-1]..., y)
        end

        a = term.coeff * term.kernel
        z = materialize(a, x, y)
        lift(z, Block(term.test_id), Block(term.trial_id), U, V)
    end
end


# function assemble(bf::BilForm, X::DirectProductSpace, Y::DirectProductSpace)

#     test_space_dict = Dict(enumerate(X.factors))
#     trial_space_dict = Dict(enumerate(Y.factors))

#     z = assemble(bf, test_space_dict, trial_space_dict)
# end

# td_assemble(dbf::DiscreteBilform; materialize=BEAST.assemble) =
#     td_assemble(dbf.bilform, dbf.test_space_dict, dbf.trial_space_dict; materialize)


# function td_assemble(bilform::BilForm, test_space_dict, trial_space_dict;
#     materialize=BEAST.assemble)

#     X = _spacedict_to_directproductspace(test_space_dict)
#     Y = _spacedict_to_directproductspace(trial_space_dict)

#     return td_assemble(bilform, X, Y)
# end

# function td_assemble(bilform::BilForm,
#     X::DirectProductSpace,
#     Y::DirectProductSpace)

#     M = [length(fct.space) for fct in X.factors]
#     N = [length(fct.space) for fct in Y.factors]

#     row_axis = BlockArrays.blockedrange(M)
#     col_axis = BlockArrays.blockedrange(N)

#     sum(bilform.terms) do t

#         a = t.coeff * t.kernel

#         m = t.test_id
#         x = X.factors[m]
#         for op in reverse(t.test_ops)
#             x = op[end](op[1:end-1]..., x)
#         end

#         n = t.trial_id
#         y = Y.factors[n]
#         for op in reverse(t.trial_ops)
#             y = op[end](op[1:end-1]..., y)
#         end

#         z = assemble(a, x, y)
#         lift(z, Block(m), Block(n), row_axis, col_axis)
#     end
# end
