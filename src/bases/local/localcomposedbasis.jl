abstract type _LocalBasisOperations{T} <: RefSpace{T} end
numfunctions(a::_LocalBasisOperations) = coalesce(numfunctions(a.el1) , numfunctions(a.el2))

struct _LocalBasisTimes{T,U,V} <: _LocalBasisOperations{T}
    el1::U
    el2::V
    function _LocalBasisTimes(::Type{T},el1::U,el2::V) where {T,U,V}
        new{T,U,V}(el1,el2)
    end
end

struct _LocalBasisCross{T,U,V} <: _LocalBasisOperations{T}
    el1::U
    el2::V
    function _LocalBasisCross(::Type{T},el1::U,el2::V) where {T,U,V}
        new{T,U,V}(el1,el2)
    end
end

struct _LocalBasisDot{T,U,V} <: _LocalBasisOperations{T}
    el1::U
    el2::V
    function _LocalBasisDot(::Type{T},el1::U,el2::V) where {T,U,V}
        new{T,U,V}(el1,el2)
    end
end

function (op::U where {U <: _LocalBasisOperations})(p)
    l = op.el1(p)
    r = op.el2(p)
    return _execute_operation(l,r,op)
end

function (op::NormalVector)(x::CompScienceMeshes.MeshPointNM)
    return normal(x)
end

operation(a::_LocalBasisTimes) = *
operation(a::_LocalBasisCross) = ×
operation(a::_LocalBasisDot) = (x,y) --> transpose(x)*y

_execute_operation(el1::SVector{N,<:NamedTuple},el2::SVector,op::U) where {N,U} = SVector{N}(_execute_operation_named(el1.data,el2,operation(op)))
_execute_operation(el1::SVector{N,<:NamedTuple},el2::U,op::_LocalBasisTimes) where {N,U <: Number} = SVector{N}(_execute_operation_named(el1.data,el2,operation(op)))
_execute_operation(el1::SVector,el2::SVector{N,<:NamedTuple},op::U) where {N,U} = SVector{N}(_execute_operation_named(el1,el2.data,operation(op)))
_execute_operation(el1::U,el2::SVector{N,<:NamedTuple},op::_LocalBasisTimes) where {N,U <: Number} = SVector{N}(_execute_operation_named(el1,el2.data,operation(op)))
_execute_operation(el1::SVector,el2::SVector,op::U) where {U <: Union{_LocalBasisCross,_LocalBasisDot}} = operation(op)(el1,el2)
_execute_operation(el1::U,el2::V,op::_LocalBasisTimes) where {U <: Number,V <: Number}= el1*el2
_execute_operation(el1::SVector{N,<:NamedTuple},el2::SVector{M,<:NamedTuple},op) where {N,M} = @error "multiplication of basisses not supported"

_execute_operation_named(a::NTuple{N},b,op) where {N} = ((value=op(a[1].value,b),), _execute_operation_named(Base.tail(a),b,op)...)
_execute_operation_named(a::NTuple{1},b,op) = ((value=op(a[1].value,b),),)
_execute_operation_named(b,a::NTuple{N},op) where {N} = ((value=op(b,a[1].value),), _execute_operation_named(b,Base.tail(a),op)...)
_execute_operation_named(b,a::NTuple{1},op) = ((value=op(b,a[1].value),),)