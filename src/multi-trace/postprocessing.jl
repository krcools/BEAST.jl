
#### Post processing
abstract type PostField end
struct AField <: PostField
    world::World
    result::PseudoBlockVector
    volumes::Vector{Int}
end

struct BField <: PostField
    world::World
    result::PseudoBlockVector
    volumes::Vector{Int}
end
struct Field{T}
    field::T
end
import Base
+(a::Field,b::Field) = Field(a.field + b.field)
+(a::Field{<:Nothing},b::Field) = b
+(a::Field,b::Field{<:Nothing}) = a
+(a::Field{<:Nothing},b::Field{<:Nothing}) = a
*(a::Number,b::Field) = Field(a*b.field)
*(a::Number,b::Field{<:Nothing}) = b


function calculate_field(points,f::PostField,strat)
    out = []
    for v in f.volumes
        obj = f.world.objectarray[v]
        push!(out,_calculate_field(points,obj,f,strat))
    end
    return sum(out)
end

function _calculate_field(points,obj::Object{<:FreeSpace},f::PostField,strat)
    map = f.world.trialhilbertspacemap
    out = []
    for c in children(obj)
        push!(out,sum(_calculate_field.(Ref(points),list_of_operators(obj,c,f,strat),BlockArrays.blocks(f.result[Block(map[c.index])]), f.world.trialdirectproductspace.factors[map[c.index]].factors)))
    end
    return sum(out)
end

function _calculate_field(points,obj::Object{<:HOM},f::PostField,strat)
    map = f.world.trialhilbertspacemap
    typeof(f.result)
    out = sum(_calculate_field.(Ref(points),transform(obj,strat)*list_of_operators(obj,obj,f,strat),BlockArrays.blocks(f.result[Block(map[obj.index])]), f.world.trialdirectproductspace.factors[map[obj.index]].factors))
    for c in children(obj)
        out += sum(_calculate_field.(Ref(points),list_of_operators(obj,c,f,strat),BlockArrays.blocks(f.result[Block(map[c.index])]), f.world.trialdirectproductspace.factors[map[c.index]].factors))
    end
    return out
end
function _calculate_field(points,obj::Object{<:PEC},f::PostField,strat)
    return Field(nothing)
end

function _calculate_field(points,operator::LinearCombinationOfOperators,coeffs,basis)
    out = []
    for (coeff,op) in zip(operator.coeffs,operator.ops)
        push!(out,coeff*_calculate_field(points,op,coeffs,basis))
    end
    return sum(out)
end
function _calculate_field(points,operator::ComposedOperator,coeffs,basis)
    x = potential(operator,points,coeffs,basis)

    return Field(x)
end
function _calculate_field(points,operator::ZeroOperator,coeffs,basis)
    return Field(nothing)
end

#### complement error
function complement_error(world,solution,volume::Vector{Int},strat::T) where {T}
    @warn "not tested yet"
    newstrat = T(-strat.trace)
    lhs = discretise_lhs(world,newstrat)
    out = lhs[volume].*Ref(solution)
    return norm(sum(out))/norm(solution)
end