struct FunctionWrapper{T} 
    func::Function
    function FunctionWrapper(f::Function;evalpoint = @SVector [0.0,0.0,0.0])
        new{typeof(f(evalpoint))}(f)
    end
end

function (f::FunctionWrapper{T})(x)::T where {T}
    return f.func(x)
end

function (f::FunctionWrapper{T})(x::CompScienceMeshes.MeshPointNM)::T where {T}
    return f.func(cartesian(x))
end
scalartype(F::FunctionWrapper{T}) where {T} = eltype(T)
function scalartype(::NormalVector) 
    @warn "The scallartype of the NormalVector is set at Float32, if used in combination with Float64 basis or operator this is no problem, in the case of Float16 it is."
    return Float32
end

abstract type _BasisOperations{T} <: Space{T} end

struct _BasisTimes{T} <: _BasisOperations{T}
    el1
    el2
    function _BasisTimes(el1,el2)
        new{promote_type(scalartype(el1),scalartype(el2))}(el1,el2)
    end
    function _BasisTimes{T}(el1,el2) where {T}
        new{T}(el1,el2)
    end
end

struct _BasisCross{T} <: _BasisOperations{T}
    el1
    el2
    function _BasisCross(el1,el2)
        new{promote_type(scalartype(el1),scalartype(el2))}(el1,el2)
    end
    function _BasisCross{T}(el1,el2) where {T}
        new{T}(el1,el2)
    end
end

struct _BasisDot{T} <: _BasisOperations{T}
    el1
    el2
    function _BasisDot(el1,el2)
        new{promote_type(scalartype(el1),scalartype(el2))}(el1,el2)
    end
    function _BasisDot{T}(el1,el2) where {T}
        new{T}(el1,el2)
    end
end

#### wrapping of the functions
_BasisTimes(a::Function,b::Function) = _BasisTimes(FunctionWrapper(a),FunctionWrapper(b))
_BasisTimes(a::Function,b) = _BasisTimes(FunctionWrapper(a),b)
_BasisTimes(a,b::Function) = _BasisTimes(a,FunctionWrapper(b))

_BasisCross(a::Function,b::Function) = _BasisCross(FunctionWrapper(a),FunctionWrapper(b))
_BasisCross(a::Function,b) = _BasisCross(FunctionWrapper(a),b)
_BasisCross(a,b::Function) = _BasisCross(a,FunctionWrapper(b))

_BasisDot(a::Function,b::Function) = _BasisDot(FunctionWrapper(a),FunctionWrapper(b))
_BasisDot(a::Function,b) = _BasisDot(FunctionWrapper(a),b)
_BasisDot(a,b::Function) = _BasisDot(a,FunctionWrapper(b))


refspace(a::_BasisTimes{T}) where {T} = _LocalBasisTimes(T,refspace(a.el1),refspace(a.el2))
refspace(a::_BasisCross{T}) where {T} = _LocalBasisCross(T,refspace(a.el1),refspace(a.el2))
refspace(a::_BasisDot{T}) where {T} = _LocalBasisDot(T,refspace(a.el1),refspace(a.el2))
refspace(a::Function) = FunctionWrapper(a)
refspace(a::FunctionWrapper) = a 
refspace(a::NormalVector) = a

numfunctions(a::Union{NormalVector,FunctionWrapper,Function}) = missing
numfunctions(a::_BasisOperations) = coalesce(numfunctions(a.el1) , numfunctions(a.el2))

geometry(a::_BasisOperations) = coalesce(geometry(a.el1),geometry(a.el2))
basisfunction(a::_BasisOperations,i) = coalesce(basisfunction(a.el1,i),basisfunction(a.el2,i))
positions(a::_BasisOperations) = coalesce(positions(a.el1),positions(a.el2))

geometry(a::Union{NormalVector,FunctionWrapper,Function}) = missing
basisfunction(a::Union{NormalVector,FunctionWrapper,Function},i) = missing
positions(a::Union{NormalVector,FunctionWrapper,Function}) = missing


subset(a::T,I) where {T <: _BasisOperations} = T(subset(a.el1,I),subset(a.el2,I))
subset(a::Union{NormalVector,FunctionWrapper,Function},I) = a

