struct FunctionWrapper{T} 
    func::Function
end
# function FunctionWrapper(f::Function;evalpoint = @SVector [0.0,0.0,0.0])
#     FunctionWrapper{eltype(f(evalpoint))}(f)
# end
function (f::FunctionWrapper{T})(x)::T where {T}
    return f.func(x)
end

function (f::FunctionWrapper{T})(x::CompScienceMeshes.MeshPointNM)::T where {T}
    return f.func(cartesian(x))
end
scalartype(F::FunctionWrapper{T}) where {T} = eltype(T)
#scalartype(::NormalVector,s::Space) = vertextype(geometry(s))
# function scalartype(::NormalVector) 
#     @warn "The scallartype of the NormalVector is set at Float32, if used in combination with Float64 basis or operator this is no problem, in the case of Float16 it is."
#     return Float32
# end

abstract type _BasisOperations{T} <: Space{T} end

struct _BasisTimes{T} <: _BasisOperations{T}
    el1
    el2
end
_BasisTimes(el1::NormalVector,el2::Space) = _BasisTimes{promote_type(eltype(vertextype(geometry(el2))),scalartype(el2))}(el1,el2)
_BasisTimes(el2::Space,el1::NormalVector) = _BasisTimes{promote_type(eltype(vertextype(geometry(el2))),scalartype(el2))}(el2,el1)
_BasisTimes(el1::Space,el2::FunctionWrapper) = _BasisTimes{promote_type(scalartype(el1),scalartype(el2))}(el1,el2)
_BasisTimes(el1::FunctionWrapper,el2::Space) = _BasisTimes{promote_type(scalartype(el1),scalartype(el2))}(el1,el2)

struct _BasisCross{T} <: _BasisOperations{T}
    el1
    el2
end
_BasisCross(el1::NormalVector,el2::Space) = _BasisCross{promote_type(eltype(vertextype(geometry(el2))),scalartype(el2))}(el1,el2)
_BasisCross(el2::Space,el1::NormalVector) = _BasisCross{promote_type(eltype(vertextype(geometry(el2))),scalartype(el2))}(el2,el1)
_BasisCross(el1::FunctionWrapper,el2::Space) = _BasisCross{promote_type(scalartype(el1),scalartype(el2))}(el1,el2)
_BasisCross(el1::Space,el2::FunctionWrapper) = _BasisCross{promote_type(scalartype(el1),scalartype(el2))}(el1,el2)

struct _BasisDot{T} <: _BasisOperations{T}
    el1
    el2
end
_BasisDot(el1::NormalVector,el2::Space) = _BasisDot{promote_type(eltype(vertextype(geometry(el2))),scalartype(el2))}(el1,el2)
_BasisDot(el2::Space,el1::NormalVector) = _BasisDot{promote_type(eltype(vertextype(geometry(el2))),scalartype(el2))}(el2,el1)
_BasisDot(el1::Space,el2::FunctionWrapper) = _BasisDot{promote_type(scalartype(el1),scalartype(el2))}(el1,el2)
_BasisDot(el1::FunctionWrapper,el2::Space) = _BasisDot{promote_type(scalartype(el1),scalartype(el2))}(el1,el2)

# #### wrapping of the functions
# _BasisTimes(a::Function,b::Function) = _BasisTimes(FunctionWrapper(a),FunctionWrapper(b))
# _BasisTimes(a::Function,b) = _BasisTimes(FunctionWrapper(a),b)
# _BasisTimes(a,b::Function) = _BasisTimes(a,FunctionWrapper(b))

# _BasisCross(a::Function,b::Function) = _BasisCross(FunctionWrapper(a),FunctionWrapper(b))
# _BasisCross(a::Function,b) = _BasisCross(FunctionWrapper(a),b)
# _BasisCross(a,b::Function) = _BasisCross(a,FunctionWrapper(b))

# _BasisDot(a::Function,b::Function) = _BasisDot(FunctionWrapper(a),FunctionWrapper(b))
# _BasisDot(a::Function,b) = _BasisDot(FunctionWrapper(a),b)
# _BasisDot(a,b::Function) = _BasisDot(a,FunctionWrapper(b))


refspace(a::_BasisTimes{T}) where {T} = _LocalBasisTimes(T,refspace(a.el1),refspace(a.el2))
refspace(a::_BasisCross{T}) where {T} = _LocalBasisCross(T,refspace(a.el1),refspace(a.el2))
refspace(a::_BasisDot{T}) where {T} = _LocalBasisDot(T,refspace(a.el1),refspace(a.el2))
# refspace(a::Function) = FunctionWrapper(a)
refspace(a::FunctionWrapper) = a 
refspace(a::NormalVector) = a

numfunctions(a::Union{NormalVector,FunctionWrapper}) = missing
numfunctions(a::Union{NormalVector,FunctionWrapper},simp) = missing
numfunctions(a::_BasisOperations) = coalesce(numfunctions(a.el1) , numfunctions(a.el2))

geometry(a::_BasisOperations) = coalesce(geometry(a.el1),geometry(a.el2))
basisfunction(a::_BasisOperations,i) = coalesce(basisfunction(a.el1,i),basisfunction(a.el2,i))
positions(a::_BasisOperations) = coalesce(positions(a.el1),positions(a.el2))

geometry(a::Union{NormalVector,FunctionWrapper}) = missing
basisfunction(a::Union{NormalVector,FunctionWrapper},i) = missing
positions(a::Union{NormalVector,FunctionWrapper}) = missing


subset(a::T,I) where {T <: _BasisOperations} = T(subset(a.el1,I),subset(a.el2,I))
subset(a::Union{NormalVector,FunctionWrapper},I) = a

#### mathematical support will be called if there are no dedicated routines
×(n::NormalVector,s::Space) = _BasisCross(n,s)
×(s::Space,n::NormalVector) = _BasisCross(s,n)

⋅(n::NormalVector,s::Space) = _BasisDot(n,s)
⋅(s::Space,n::NormalVector) = _BasisDot(s,n)

*(n::NormalVector,s::Space) = _BasisTimes(n,s)
*(s::Space,n::NormalVector) = _BasisTimes(s,n)



# function restrict(a::_LocalBasisOperations,cell1,cell2)
#     s = sign(dot(normal(cell1),normal(cell2)))^number_of_normals(a)
#     return s*restrict(inner_refspace(a),cell1,cell2)
# end


# number_of_normals(a::NormalVector) = 1
# number_of_normals(a) = 0
# number_of_normals(a::_LocalBasisOperations) = number_of_normals(a.el1) + number_of_normals(a.el2)

# #TODO save wrapped space in type, same with refspace.
# inner_refspace(a::)