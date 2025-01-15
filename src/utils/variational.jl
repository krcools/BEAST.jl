module Variational

using BlockArrays
import BEAST
# import Base: start, done, next

export transposecalls!

# aux methods to hide the separation in head and args
numchds(xp) = length(xp.args) + 1
child(xp, idx) = idx == 1 ? xp.head : xp.args[idx-1]


mutable struct DepthFirst
    xp::Expr
end


"""
    depthfirst(xp)

Returns an iterator that visits all nodes in an Expr depth first. head and
args are visited before the Expr they define.
"""
depthfirst(xp) = DepthFirst(xp)


mutable struct DepthFirstState
    val
    par
    idx
end

import Base: iterate

function iterate(itr::DepthFirst)
    head = DepthFirstState(itr.xp, nothing, -1)
    state = DepthFirstState(itr.xp, head, 0)
    return iterate(itr, state)
end


function iterate(itr::DepthFirst, state)

    state.par == nothing && return nothing
    return next(itr, state)
end


function next(itr::DepthFirst, state::DepthFirstState)

    # if all children processed, move up on level
    if state.idx == numchds(state.val)
        return state.val, state.par
    end

    # move to next child
    state = DepthFirstState(state.val, state.par, state.idx+1)

    # if the next child is a leaf, return it
    chd = child(state.val, state.idx)
    if !isa(chd, Expr)
        return chd, state
    end

    # the next child is an expression; descend into it
    # and find the next valid state. There will always
    # be at least one, pointing to the child itself.
    state = DepthFirstState(chd, state, 0)
    return next(itr, state)
end


"""
    transposecall!(xp, skip=[])

Goes through the syntax tree and replace all function calls
`f(p1,p2,...,x)` with `x(f,p1,p2,...)`.
"""
function transposecalls!(xp, skip=[])
    isa(xp, Expr) || return xp
    for x in depthfirst(xp)
        if isa(x, Expr) && x.head == :call && !(x.args[1] in skip)
            @assert length(x.args) >= 2
            tmp = x.args[1]
            x.args[1] = x.args[end]
            x.args[end] = tmp
        end
    end
    return xp
end


import Base.==

export hilbertspace, @hilbertspace
export @varform, Equation
export LinForm, LinTerm, BilForm, BilTerm

import Base: +, -, *, getindex, ^, print
import LinearAlgebra: dot

mutable struct HilbertVector
    idx
    space
    opstack
end

Base.Int(hv::HilbertVector) = hv.idx

function hilbertspace(s::Symbol, numcomponents::Int)
    syms = [Symbol(s,i) for i in 1:numcomponents]
    return [HilbertVector(i,syms,[]) for i in 1:numcomponents]
end


Base.getindex(A::AbstractBlockArray, p::HilbertVector, q::HilbertVector) = A[Block(Int(p),Int(q))]
Base.getindex(u::AbstractBlockArray, p::HilbertVector) = u[Block(Int(p))]

Base.setindex!(A::AbstractBlockArray, v, p::HilbertVector, q::HilbertVector) = setindex!(A, v, Block(Int(p),Int(q)))
Base.setindex!(A::AbstractBlockArray, v, p::HilbertVector) = setindex!(A, v, Block(Int(p)))

Base.view(A::AbstractBlockArray, p::HilbertVector, q::HilbertVector) = view(A, Block(Int(p), Int(q)))
Base.view(A::AbstractBlockArray, p::HilbertVector) = view(A, Block(Int(p)))

mutable struct LinForm
  test_space
  terms
end

mutable struct LinTerm
  test_id
  test_ops
  coeff
  functional
end

mutable struct BilForm
  test_space
  trial_space
  terms
end

mutable struct BilTerm
  test_id
  trial_id
  test_ops
  trial_ops
  coeff
  kernel
end

mutable struct Equation
    lhs
    rhs
end

"""
    ==(lhs::BilForm, rhs::LinForm)

Build an equation from a left hand and right hand side
"""
==(lhs::BilForm, rhs::LinForm) = Equation(lhs, rhs)



# """
#     hilbert_space(type, g1, g2, ...)

# Returns generators defining a Hilbert space of field `type`
# """
# hilbertspace(vars::Symbol...) = [HilbertVector(i, [vars...], []) for i in 1:length(vars)]

function genspace(syms...)

    space = Vector{Symbol}()
    lengths = Int[]
    starts = Int[]
    stops = Int[]
    for sym in syms

        if sym isa Symbol
            push!(space, sym)
            push!(lengths,1)
            push!(starts,1)
            push!(stops,1)
        elseif sym isa Expr && sym.head == :ref
            base = sym.args[1]
            start = sym.args[2].args[2]
            stop = sym.args[2].args[3]
            for k in start:stop
                sym = Symbol(base,k)
                push!(space, sym)
            end
            push!(lengths,stop-start+1)
            push!(starts,start)
            push!(stops,stop)
        end
    end
    
    return space, starts, stops
end

macro hilbertspace(syms...)

    space, starts, stops = genspace(syms...)

    ex = quote end
    k = 1
    for (s, (start,stop)) in enumerate(zip(starts,stops))

        
        len = stop-start+1
        if syms[s] isa Symbol
            sym = syms[s]
            push!(ex.args, :($(esc(sym)) = HilbertVector($k,$space,[])))
            k += 1
        else
            sym = syms[s].args[1]
            push!(ex.args, :($(esc(sym)) = [HilbertVector(i,$space,[]) for i in $k:$(k+len-1)]))
            k += len
        end

    end
    return ex
end

# macro hilbertspace(syms...)

#     for sym in syms
#         @assert isa(sym, Symbol) "@hilbertspace takes a list of Symbols"
#     end

#     rhs = :(hilbertspace())
#     for sym in syms
#         push!(rhs.args, QuoteNode(sym))
#     end

#     vars = gensym()
#     xp = quote
#         $vars = $rhs
#     end
#     for (i,s) in enumerate(syms)
#         push!(xp.args, :($(esc(s)) = $vars[$i]))
#     end

#     xp
# end


"""
    call(u::HilbertVector, f, params...)
    u(f, params...)

Add another operation to the opstack of `u`.
"""
(u::HilbertVector)(f, params...) = HilbertVector(u.idx, u.space, [(f, params...); u.opstack])




"""
    getindex(f, v::HilbertVector)
    f[v]

Return a LinForm corresponding to f[v]
"""
getindex(f, v::HilbertVector) = LinForm(v.space, [LinTerm(v.idx, v.opstack, 1, f)])

function getindex(f, V::Vector{HilbertVector})
    terms = Vector{LinTerm}()
    for v in V
        term = LinTerm(v.idx, v.opstack, 1, f)
        push!(terms, term)
    end
    return LinForm(first(V).space, terms)
end

"""
    getindex(A, v::HilbertVector, u::HilbertVector)

Create a BilForm corresponding to A[v,u]
"""
function getindex(A, v::HilbertVector, u::HilbertVector)
    terms = [ BilTerm(v.idx, u.idx, v.opstack, u.opstack, 1, A) ]
    BilForm(v.space, u.space, terms)
end

function getindex(A::AbstractMatrix, v::HilbertVector, u::HilbertVector)
    terms = [ BilTerm(v.idx, u.idx, v.opstack, u.opstack, 1, A) ]
    BilForm(v.space, u.space, terms)
end


function getindex(A::BEAST.BlockDiagonalOperator, V::Vector{HilbertVector}, U::Vector{HilbertVector})
    op = A.op
    terms = Vector{BilTerm}()
    @assert length(V) == length(U)
    for (v,u) in zip(V,U)
            term = BilTerm(v.idx, u.idx, v.opstack, u.opstack, 1, op)
            push!(terms, term)
    end
    return BilForm(first(V).space, first(U).space, terms)
end

function getindex(A::BEAST.BlockFullOperators, V::Vector{HilbertVector}, U::Vector{HilbertVector})
    op = A.op
    terms = Vector{BilTerm}()
    # @assert length(V) == length(U)
    for v in V
        for u in U
            term = BilTerm(v.idx, u.idx, v.opstack, u.opstack, 1, op)
            push!(terms, term)
        end
    end
    return BilForm(first(V).space, first(U).space, terms)
end


function getindex(op::Any, V::Vector{HilbertVector}, U::Vector{HilbertVector})
    terms = Vector{BilTerm}()
    for v in V
        for u in U
            term = BilTerm(v.idx, u.idx, v.opstack, u.opstack, 1, op)
            push!(terms, term)
        end
    end
    return BilForm(first(V).space, first(U).space, terms)
end

function getindex(A::Matrix, v::HilbertVector, u::HilbertVector)
    terms = [ BilTerm(v.idx, u.idx, v.opstack, u.opstack, 1, A) ]
    BilForm(v.space, u.space, terms)
end


"Add two BilForms together"
function +(a::BilForm, b::BilForm)
    @assert a.test_space == b.test_space
    @assert a.trial_space == b.trial_space
    BilForm(a.test_space, a.trial_space, [a.terms; b.terms])
end

function+(a::LinForm, b::LinForm)
    @assert a.test_space == b.test_space
    LinForm(a.test_space, [a.terms; b.terms])
end

function *(α::Number, a::BilForm)
  b = deepcopy(a)
  for t in b.terms t.coeff *= α end
  return b
end

function *(α::Number, a::LinForm)
    b = deepcopy(a)
    for t in b.terms t.coeff *= α end
    return b
end

-(a::BilForm) = (-1 * a)
-(a::BilForm, b::BilForm) = a + (-b)

-(a::LinForm) = (-1 * a)
-(a::LinForm, b::LinForm) = a + (-b)


function print(io::IO, v::HilbertVector)
    sym = v.space[v.idx]
    ops = v.opstack
    for op in ops
        print(io, op[1], "(")
    end
    print(io, sym)
    for op in reverse(ops)
        for p in op[2:end]
            print(io, ", ", p)
        end
        print(io, ")")
    end
end


function print(io::IO, a::LinForm)
    N = length(a.terms)
    #T = typeof(a.terms[1].coeff)
    for (n,t) in enumerate(a.terms)
      u = HilbertVector(t.test_id, a.test_space, t.test_ops)
      t.coeff != 1 && print(io, t.coeff, "*")
      print(io, t.functional, "[", u, "]")
      n == N || print(io, " + ")
  end
end


function print(io::IO, f::BilForm)
    N = length(f.terms)
    #T = typeof(f.terms[1].coeff)

    for (n,t) in enumerate(f.terms)
        u = HilbertVector(t.test_id, f.test_space, t.test_ops)
        v = HilbertVector(t.trial_id, f.trial_space, t.trial_ops)
        t.coeff != 1 && print(io, t.coeff, "*")
        print(io, t.kernel, "[", u, ", ", v, "]")
        n == N || print(io, " + ")
    end
end

function print(io::IO, eq::Equation)
    print(io, eq.lhs)
    print(io, " == ")
    print(io, eq.rhs)
end


"""
    @varform <form-definition>

The Julia form compiler uses the Julia parser and meta-programming
based traversal of the AST to create a structure containing all
information required for the description of a variational formulation
from an Expr that follows closely widespread mathematical convention.

E.g:

    EFIE = @varform T[k,j] = e[k]
    MFIE = @varform 0.5*I[k,j] + K[k,j] = h[k]
    PMCH = @varform M[k,j] - η*T[k,m] + 1/η*T[l,j] + M[l,m] = e[k] + h[l]
"""
macro varform(x)
    y = transposecalls!(x, [:+, :-, :*, :^, :(==)])
    esc(y)
end


struct DirectProductKernel
    bilforms
end

function Base.getindex(A::DirectProductKernel, V::Vector{HilbertVector}, U::Vector{HilbertVector})
    terms = Vector{BilTerm}()
    @assert length(V) == length(U) == length(A.bilforms)

    for (v,u,op) in zip(V,U, A.bilforms)
        term = BilTerm(v.idx, u.idx, v.opstack, u.opstack, 1, op)
        push!(terms, term)
    end

    return BilForm(first(V).space, first(U).space, terms)
end

struct BlockDiagKernel
    bilform
end

function Base.getindex(A::BlockDiagKernel, V::Vector{HilbertVector}, U::Vector{HilbertVector})
    terms = Vector{BilTerm}()
    @assert length(V) == length(U)

    op = A.bilform
    for (v,u) in zip(V,U)
        term = BilTerm(v.idx, u.idx, v.opstack, u.opstack, 1, op)
        push!(terms, term)
    end

    return BilForm(first(V).space, first(U).space, terms)
end


struct OffDiagKernel
    bilform
end

function getindex(op::OffDiagKernel, V::Vector{HilbertVector}, U::Vector{HilbertVector})
    terms = Vector{BilTerm}()
    for (i,v) in enumerate(V)
        for (j,u) in enumerate(U)
            i == j && continue
            term = BilTerm(v.idx, u.idx, v.opstack, u.opstack, 1, op.bilform)
            push!(terms, term)
        end
    end
    return BilForm(first(V).space, first(U).space, terms)
end

end
