using BEAST
using Combinatorics
import Base: +,-,*
import Base
@hilbertspace b[1:4]
@hilbertspace t[1:4]

# nodig:    1) functies die gegeven 2 objecten(volumes) en hilbertspacede linforms teruggeven. 
#       ok  2) volumes hebben een type, PEC of HOMOGEEN
#           3) maak een interactie aan en geef deze ook telkens mee, vb: mfie, efie, pmchwt, cfie, vectorpotential
#       ok  4) multi-trace struct als wrapper voor deze types
#       ok  5) normal is outward pointing --> opslaan in object of het een outward of inward normal object is
#           6) detector for this normal 
#        ok   7) functie die matrix kan maken voor gegeven object.
#           8) functie die matrices combineert in 1 grote.
#           9) preconditioner als wrapper
#           10) post processing of a volume

@hilbertspace a
@hilbertspace b

@hilbertspace c
@hilbertspace d

t = BEAST.Identity()
q = t[a,b]

struct Zero end
+(a,b::Zero) = a
+(a::Zero,b) = b
+(a::Zero,b::Zero) = a
*(a::Zero,b) = a
*(b,a::Zero) = a

struct ObjectType end

struct Inside{T} 
    inside::T
end
struct PECEFIE
    mesh
    ω
    inside
end
PECEFIE(a,b) = PECEFIE(a,b,ObjectType())
struct PECMFIE
    mesh
    ω
    inside
end
PECMFIE(a,b) = PECMFIE(a,b,ObjectType())
struct PECCFIE
    mesh
    ω
    alpha
    inside
end
PECCFIE(a,b,c) = PECCFIE(a,b,c,ObjectType())
PEC = Union{PECEFIE,PECMFIE}
struct CavityPEC
    mesh
    ω
    inside

end


struct HOM{T}
    mesh 
    ϵ::T
    μ::T
    ω::T
    inside

end
HOM(a,b,c,d) = HOM(a,b,c,d,ObjectType())
struct FreeSpace{T} 
    ϵ::T
    μ::T
    ω::T
    inside
end
FreeSpace(a,b,c) = FreeSpace(a,b,c,ObjectType())
mutable struct Object{T} 
    index::Int
    type::T
    children::Vector{}
    parent
end
children(obj::Object) = obj.children
parent(obj::Object) = obj.parent

mutable struct World
    objecttree
    objectarray
    testhilbertspace
    trialhilbertspace
    testhilbertspacemap# given an index yields basis index, TODO write function that can be used in both directions
    trialhilbertspacemap
    testdirectproductspace
    trialdirectproductspace
end
World(a,b,c,d,e,f) = World(a,b,c,d,e,f,nothing,nothing)
Base.copy(w::World) = World(w.objecttree,w.objectarray,
w.testhilbertspace,w.trialhilbertspace,
w.testhilbertspacemap,w.trialhilbertspacemap,
w.testdirectproductspace,w.trialdirectproductspace)

function add_child!(obj,child)
    push!(obj.children,child)
end
function add_parent!(obj,parent)
    obj.parent = parent
end
"""

parent at 0 means no parent.
"""
function create_world(embeddings::Vector{Int}, objecttypes::Vector)
    @assert length(embeddings) == length(objecttypes)
    
    objectarray = [Object(i,objecttypes[i],[],nothing) for i in 1:length(objecttypes)]
    testhilbertspacemap = zeros(Int,length(embeddings))
    trialhilbertspacemap = zeros(Int,length(embeddings))
    #embeddings[embeddings.==0] = embeddings[embeddings.==0] .+ length(embeddings) .+ 1
    counter = 1
    for (ind,parent) in enumerate(embeddings[1:end])
        parent==0 && continue 
        add_child!(objectarray[parent],objectarray[ind])
        add_parent!(objectarray[ind],objectarray[parent])
        if typeof(objecttypes[ind]) != FreeSpace
            testhilbertspacemap[ind] = counter
            trialhilbertspacemap[ind] = counter
            counter += 1
        end
    end
    testhilbertspace = BEAST.hilbertspace(:t, counter-1)
    trialhilbertspace = BEAST.hilbertspace(:b, counter-1)
    return World(objectarray[end],objectarray,testhilbertspace,trialhilbertspace,testhilbertspacemap,trialhilbertspacemap)
end

function assign_basis!(w::World, strat)
    w.testdirectproductspace = BEAST.DirectProductSpace([testbasis(w.objectarray[i],strat) for (i,j) in enumerate(w.testhilbertspacemap) if j!=0])
    w.trialdirectproductspace =  BEAST.DirectProductSpace([trialbasis(w.objectarray[i],strat) for (i,j) in enumerate(w.trialhilbertspacemap) if j!=0])
end


abstract type Interaction end

struct VP <: Interaction 
    trace::Int
end # VectorPotential

VP() = VP(1)

abstract type Excitation end
struct EHExcitation <: Excitation end
struct CurrentExcitation <: Excitation end
struct VPExcitation <: Excitation
    A
    curlA
    divA
    objectids
end



# struct Outside{T} 
#     object::T
# end
parent(obj::Inside) = parent(obj.object)

epsilon(o::Union{HOM,FreeSpace}) = o.ϵ
epsilon(o::Object) = epsilon(o.type)
epsilon(o::Inside) = epsilon(o.inside)

mu(o::Union{HOM,FreeSpace}) = o.μ
mu(o::Object) = mu(o.type)
mu(o::Inside) = mu(o.inside)

omega(o::Union{HOM,FreeSpace,PEC}) = o.ω
omega(o::Object) = omega(o.type)
omega(o::Inside) = omega(o.inside)

geometry(o::Object) = geometry(o.type)
geometry(o::Union{HOM,FreeSpace,PEC}) = o.mesh
geometry(o::Inside) = geometry(o.inside)

function matrix_to_bilform(mat)
    nrows,ncols = size(mat)
    tin = BEAST.hilbertspace(:tin,nrows)
    bin = BEAST.hilbertspace(:bin,ncols)
    terms = [mat[i,j][tin[i],bin[j]] for i in 1:nrows, j in 1:ncols if typeof(mat[i,j]) != ZeroOperator]
    if length(terms) > 0
        return sum(terms)
    else
        return  ZeroOperator()[tin[1],bin[1]]
    end
end

function array_to_linform(array)
    nrows = length(array)
    tin = BEAST.hilbertspace(:tin,nrows)
    println(sum([array[i][tin[i]] for i in 1:nrows]))
    return sum([array[i][tin[i]] for i in 1:nrows])
end

# function discretise_rhs(w::World, ex::Excitation, strat::Interaction)
#     out = []
#     l = length(w.objectarray)-1 # no bases for the free space
#     test = []
#     testb = BEAST.DirectProductSpace([testbasis(obj,strat) for obj in w.objectarray])
#     for obj in w.objectarray
        
#         if !inside_interaction_matters(obj.type) 
#             push!(out,0)
#             continue
#         end
#         t = test[obj.index]
#         lins = array_to_linform(excitation(obj.type,obj.type,ex,strat))[t]
#         for c in children(obj)
#             tc = test[c.index]
#             lins += array_to_linform(excitation(c.type,obj.type,ex,strat))[tc]
#         end
#         push!(out,assemble(lins,testb))
#     end
#     return out
# end

function discretise_rhs(w::World, ex::Excitation, strat::Interaction)
    tesths = w.testhilbertspace

    testb = w.testdirectproductspace

    maptesths = w.testhilbertspacemap

    out = _discretise_rhs_matrix.(w.objectarray,Ref(tesths),Ref(testb),Ref(maptesths),Ref(ex),Ref(strat))
    return out
end
function _discretise_rhs_matrix(obj::Object{T},tesths,testb,mapt,ex,strat) where {T <: FreeSpace}
    out = []
    for child in children(obj)

        u = tesths[mapt[child.index]]

        push!(out,array_to_linform(excitation(child,obj,child.type,ex,strat))[u])

    end
    x = assemble(sum(out),testb)
    return x
end
function _discretise_rhs_matrix(obj::Object{T},tesths,testb,mapt,ex,strat) where {T <: HOM}
    out = []
    push!(out,array_to_linform(excitation(obj,obj,Inside(obj.type),ex,strat))[tesths[mapt[obj.index]]])
    for child in children(obj)
        u = tesths[mapt[child.index]]
        push!(out,array_to_linform(excitation(child,obj,child.type,ex,strat))[u])
    end
    s = sum(out)
    println(s)
    return assemble(s,testb)
end
function _discretise_rhs_matrix(obj::Object{T},tesths,testb,mapt,ex,strat) where {T <: PEC}
    return Zero()
end
function identity(w1::World,w2::World,strat)
    t = w1.trialhilbertspace
    b = w2.testhilbertspace
    @assert length(t) == length(b)
    bil = Zero()

    for (k,(tt,bb)) in enumerate(zip(t,b))
        ind = findfirst(==(k),w1.trialhilbertspacemap)
        bil += matrix_to_bilform(_identity(w1.objectarray[ind].type,strat))[tt,bb]
    end
    println(t)
    return assemble(bil,w1.trialdirectproductspace,w2.testdirectproductspace)
end

function diag_prec(w::World, strat::Interaction)
    l = length(w.objectarray)-1 # no bases for the free space
    tesths = w.testhilbertspace
    basishs = w.trialhilbertspace
    testb = w.testdirectproductspace
    trialb = w.trialdirectproductspace
    maptesths = w.testhilbertspacemap
    maptrialhs = w.trialhilbertspacemap
    out = _diag_prec.(w.objectarray,Ref(tesths),Ref(basishs),Ref(testb),Ref(trialb),Ref(maptesths),Ref(maptrialhs),Ref(strat))
    return out

end

function discretise_lhs(w::World, strat::Interaction)

    l = length(w.objectarray)-1 # no bases for the free space
    tesths = w.testhilbertspace
    basishs = w.trialhilbertspace
    testb = w.testdirectproductspace
    trialb = w.trialdirectproductspace
    maptesths = w.testhilbertspacemap
    maptrialhs = w.trialhilbertspacemap
    out = _discretise_lhs_matrix.(w.objectarray,Ref(tesths),Ref(basishs),Ref(testb),Ref(trialb),Ref(maptesths),Ref(maptrialhs),Ref(strat))
    return out
# returned the assambled matrices for each struct
end
function _discretise_lhs_matrix(obj::Object{<:HOM},tesths,basishs,testb,trialb,mapt,mapb,strat)
    t = tesths[mapt[obj.index]]
    b = basishs[mapb[obj.index]]
    # bils = BEAST.BlockDiagonalOperator(Identity())[t,b]
    bils = matrix_to_bilform(identity(obj.type,strat))[t,b]
    println(bils)
    bils -= matrix_to_bilform(interaction(obj,obj,obj,Inside(obj.type),Inside(obj.type),strat))[t,b]
    for child in children(obj)
        u = tesths[mapt[child.index]]
        v = basishs[mapb[child.index]]
        bils += matrix_to_bilform(identity(child.type,strat))[u,v]
        bils -= matrix_to_bilform(interaction(obj,child,obj,Inside(obj.type),child.type,strat))[t,v]
        bils -= matrix_to_bilform(interaction(child,obj,obj,child.type,Inside(obj.type),strat))[u,b]
        bils -= matrix_to_bilform(interaction(child,child,obj,child.type,child.typ,strat))[u,v]
    end
    if length(children(obj)) > 1
        for (c1,c2) in combinations(children(obj),2)
            t1 = tesths[mapt[c1.index]]
            t2 = tesths[mapt[c2.index]]
            b1 = basishs[mapb[c1.index]]
            b2 = basishs[mapb[c2.index]]
            bils -= matrix_to_bilform(interaction(c1,c2,obj,c1.type,c2.type,strat))[t1,b2]
            bils -= matrix_to_bilform(interaction(c2,c1,obj,c1.type,c2.type,strat))[t2,b1]
        end
    end
    println(typeof(bils.terms[1]))
    println(typeof(testb.factors[1]))
    return assemble(bils,testb,trialb)
end
function _diag_prec(obj::Object{<:HOM},tesths,basishs,testb,trialb,mapt,mapb,strat)
    t = tesths[mapt[obj.index]]
    b = basishs[mapb[obj.index]]
    # bils = BEAST.BlockDiagonalOperator(Identity())[t,b]
    bils = matrix_to_bilform(identity(obj.type,strat))[t,b]
    println(bils)
    bils -= matrix_to_bilform(interaction(obj,obj,obj,Inside(obj.type),Inside(obj.type),strat))[t,b]
    for child in children(obj)
        u = tesths[mapt[child.index]]
        v = basishs[mapb[child.index]]
        bils += matrix_to_bilform(identity(child.type,strat))[u,v]
        bils -= matrix_to_bilform(interaction(child,child,obj,child.type,child.typ,strat))[u,v]
    end
    return assemble(bils,testb,trialb)
end
function _diag_prec(obj::Object{<:PEC},tesths,basishs,testb,trialb,mapt,mapb,strat)
    return Zero()
end

function _discretise_lhs_matrix(obj::Object{<:PEC},tesths,basishs,testb,trialb,mapt,mapb,strat)
    return Zero()
end
function _discretise_lhs_matrix(obj::Object{<:FreeSpace},tesths,basishs,testb,trialb,mapt,mapb,strat)
    bilslist = []
    for child in children(obj)
        u = tesths[mapt[child.index]]
        v = basishs[mapb[child.index]]
        push!(bilslist,matrix_to_bilform(identity(child.type,strat))[u,v])
        push!(bilslist,-matrix_to_bilform(interaction(child,child,obj,child.type,child.type,strat))[u,v])
    end
    bils = sum(bilslist)
    if length(children(obj)) > 1
        for (c1,c2) in combinations(children(obj),2)
            t1 = tesths[mapt[c1.index]]
            t2 = tesths[mapt[c2.index]]
            b1 = basishs[mapb[c1.index]]
            b2 = basishs[mapb[c2.index]]
            bils -= matrix_to_bilform(interaction(c1,c2,obj,c1.type,c2.type,strat))[t1,b2]
            bils -= matrix_to_bilform(interaction(c2,c1,obj,c1.type,c2.type,strat))[t2,b1]
        end
    end
    return assemble(bils,testb,trialb)
end
function _diag_prec(obj::Object{<:FreeSpace},tesths,basishs,testb,trialb,mapt,mapb,strat)
    bilslist = []
    for child in children(obj)
        u = tesths[mapt[child.index]]
        v = basishs[mapb[child.index]]
        push!(bilslist,matrix_to_bilform(identity(child.type,strat))[u,v])
        push!(bilslist,-matrix_to_bilform(interaction(child,child,obj,child.type,child.type,strat))[u,v])
    end
    bils = sum(bilslist)
    return assemble(bils,testb,trialb)
end



# inside_interaction_matters(obj::PEC) = false
# inside_interaction_matters(obj::ObjectType) = true

#create these type of functions
# function testbasis(obj::ObjectType,::VP)
#     return
# end
# function trialbasis(obj::ObjectType,::VP)
# end

# # function interaction(testobj::ObjectType,basisobj::Inside{PEC},embedobj::HOM,::VP)
# #     return 0
# # end

# function interaction()

# end

# function excitation(testobj::ObjectType,embedobj::FreeSpace,ex::VPExcitation,::VP)

# end
# excitation(testobj::ObjectType,embedobj::ObjectType,ex::VPExcitation,::VP) = 0


