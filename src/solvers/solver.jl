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

function assemble(lform::LinForm, test_space_dict;kwargs...)
    X = _spacedict_to_directproductspace(test_space_dict)
    return assemble(lform, X;kwargs...)
end

scalartype(lf::LinForm) = scalartype(lf.terms...)
scalartype(lt::LinTerm) = scalartype(lt.coeff, lt.functional)

function assemble(lform::LinForm, X::DirectProductSpace;kwargs...)

    @assert !isempty(lform.terms)

    T = scalartype(lform, X)
    U = NestedUnitRanges.nestedrange(spatialbasis(X), 1, numfunctions)
    
    x = first(AbstractTrees.Leaves(X))
    N = Base.OneTo(tensordim(x,2))
    ax = _righthandside_axes(x, U, N)

    B = BlockArray{T}(undef, ax)
    fill!(B, 0)

    for t in lform.terms

        m = t.test_id
        x = X.factors[m]

        for op in reverse(t.test_ops) x = op[end](op[1:end-1]..., x) end

        b = assemble(t.functional, x;kwargs...)
        B[Block(m),Block(1)] = t.coeff * b
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
    X = _spacedict_to_directproductspace(test_space_dict)
    return td_assemble(lform, X)
end

_righthandside_axes(x::SpaceTimeBasis, U, N) = (U,N,)
_righthandside_axes(x, U, N) = (U,)

td_assemble(lform::LinForm, X::DirectProductSpace) = assemble(lform, X)

function assemble(bilform::BilForm, test_space_dict, trial_space_dict;
    kwargs...)

    X = _spacedict_to_directproductspace(test_space_dict)
    Y = _spacedict_to_directproductspace(trial_space_dict)

    return assemble(bilform, X, Y; kwargs...)
end

lift(a,I,J,U,V) = LiftedMaps.LiftedMap(a,I,J,U,V)
lift(a::ConvolutionOperators.AbstractConvOp ,I,J,U,V) =
    ConvolutionOperators.LiftedConvOp(a, U, V, I, J)

function assemble(bf::BilForm, X::DirectProductSpace, Y::DirectProductSpace;
    materialize=BEAST.assemble,kwargs...)

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
        z = materialize(a, x, y;kwargs...)
        lift(z, Block(term.test_id), Block(term.trial_id), U, V)
    end
end

function assemble(bf::BilForm, X::Space, Y::Space;kwargs...)
    @assert length(bf.terms) == 1
    assemble(bf, BEAST.DirectProductSpace([X]), BEAST.DirectProductSpace([Y]);kwargs...)
end

function assemble(bf::BilForm, pairs::Pair...;kwargs...)
    dbf = discretise(bf, pairs...)
    assemble(dbf;kwargs...)
end
# struct BilFormDirectProductSpace
# test_space
# trial_space
# terms
# end

function LinearAlgebra.dot(test::DirectProductSpace,op::Matrix{AbstractOperator},trial::DirectProductSpace)
    test_length = length(test.factors)
    trial_length = length(trial.factors)
    @assert size(op)==(test_length,trial_length)
    terms = Vector{BilTerm}()
    test_space_dict = Dict()
    trial_space_dict = Dict()
    for (i,space) in enumerate(test.factors)
        test_space_dict[i]=space
    end
    for (i,space) in enumerate(trial.factors)
        trial_space_dict[i]=space
    end

    for i in 1:test_length
        for j in 1:trial_length
            push!(terms,BilTerm(i,j,[],[],1.0,op[i,j]))
        end
    end
    BilForm([],[],terms)
end

