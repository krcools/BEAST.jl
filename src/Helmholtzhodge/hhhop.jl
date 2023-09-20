import LinearAlgebra: ×, ⋅
abstract type HHHOperator{T,K} <: IntegralOperator end
abstract type HHHtestOperator{T,K} <: HHHOperator{T,K} end
abstract type HHHbasisOperator{T,K} <: HHHOperator{T,K} end
abstract type HHHkernelOperator{T,K} <: HHHOperator{T,K} end
# scalartype(op::HHHOperator{T,K}) where {T, K <: Nothing} = T
# scalartype(op::HHHOperator{T,K}) where {T, K} = promote_type(T, K)
scalartype(op::HHHOperator{T,K}) where {T, K <: Nothing} = Complex
scalartype(op::HHHOperator{T,K}) where {T, K} = Complex


const i4pi = 1 / (4pi)
mutable struct HHHgreen{T,K} <: HHHkernelOperator{T,K}
    α::K
    γ::T
    op::HHHOperator
end

mutable struct HHHgradgreen{T,K} <: HHHkernelOperator{T,K}
    α::K
    γ::T
    op::HHHOperator
end

mutable struct HHHgradgreenCross{T,K} <: HHHkernelOperator{T,K}
    α::K
    γ::T
    op::HHHOperator
end

abstract type HHHdata end
struct HHHvector <:HHHdata end
struct HHHscalar <:HHHdata end


mutable struct HHHNtestCross{T,K} <: HHHtestOperator{T,K}
    op::HHHOperator{T,K}
end


mutable struct HHHNbasisCross{T,K} <: HHHbasisOperator{T,K}
    op::HHHOperator{T,K}
end

mutable struct HHHNtestDot{T,K} <: HHHtestOperator{T,K}
    op::HHHOperator{T,K}
end

mutable struct HHHNbasisnormal{T,K} <: HHHbasisOperator{T,K}
    op::HHHOperator{T,K}
end
mutable struct HHHNbasisdot{T,K} <: HHHbasisOperator{T,K}
    op::HHHOperator{T,K}
end
mutable struct HHHIdentity{T,K} <: HHHbasisOperator{T,K}
end
hhhidentity() = HHHIdentity{Complex,Complex}()
export HHH
struct BasisFunction end
basisfunction() = BasisFunction()
export basisfunction
VectorToVector = Union{HHHNtestCross,HHHNbasisCross,HHHgradgreenCross,HHHgreen}
ScalarToVector = Union{HHHNbasisnormal,HHHgradgreen}
VectorToScalar = Union{HHHNtestDot,HHHgradgreen,HHHNbasisdot}
ScalarToScalar = Union{HHHgreen}
inversemap(op::VectorToVector,::HHHvector) = inversemap(op.op,HHHvector())
inversemap(op::ScalarToVector,::HHHvector) = inversemap(op.op,HHHscalar())
inversemap(op::VectorToScalar,::HHHscalar) = inversemap(op.op,HHHvector())
inversemap(op::ScalarToScalar,::HHHscalar) = HHHscalar()
inversemap(::HHHOperator,::HHHdata) = @error "trying to convert scalar to vector or vector to scalar"
inversemap(::HHHIdentity,data::HHHdata) = data
##### Operator building

strace(op::Union{HHHkernelOperator,HHHtestOperator}) = HHHNtestCross(op)
ntrace(op::Union{HHHkernelOperator,HHHtestOperator}) = HHHNtestDot(op)
×(op::HHHgradgreen,::Nothing) = HHHgradgreenCross(op.α,op.γ,op.op)
×(op::HHHgradgreen,::NormalVector) = HHHgradgreenCross(op.α,op.γ,op.op)(HHHNbasisnormal(hhhidentity()))
×(::NormalVector,op::Union{HHHkernelOperator,HHHtestOperator}) = strace(op)
⋅(::NormalVector,op::Union{HHHkernelOperator,HHHtestOperator}) = ntrace(op)
*(::NormalVector,::BasisFunction) = HHHNbasisnormal(hhhidentity())
⋅(::NormalVector,::BasisFunction) = HHHNbasisdot(hhhidentity())
×(::NormalVector,::BasisFunction) = HHHNbasisCross(hhhidentity())
×(::NormalVector,op::HHHbasisOperator) = HHHNbasisCross(op)
⋅(::NormalVector,op::HHHbasisOperator) = HHHNbasisdot(op)


function (op::HHHkernelOperator)(Basisop::HHHbasisOperator) 
    return typeof(op)(op.α,op.γ,Basisop)
end



function (op::HHHOperator)(BasisOperator::HHHbasisOperator)
    op.op = BasisOperator
end


##### Integrand definitions
function (igd::Integrand{<:HHHOperator})(x,y,f,g)
    op = igd.operator

test = getvalue(f)
_krondot(test,op(x,y,g))
end
function (op::HHHNtestCross)(x,y,g)
    nx = normal(x)
    cross.(Ref(nx),op.op(x,y,g))
end
function (op::HHHNtestDot)(x,y,g)
    nx = normal(x)
    dot.(Ref(nx),op.op(x,y,g))
end
function (op::HHHNbasisnormal)(x,y,g)
    ny = normal(y)
    Ref(ny).*op.op(x,y,g)
end
function (op::HHHNbasisdot)(x,y,g)
    ny = normal(y)
    dot.(Ref(ny),op.op(x,y,g))
end
function (op::HHHNbasisCross)(x,y,g)
    ny = normal(y)
    cross.(Ref(ny),op.op(x,y,g))
end
function (op::HHHIdentity)(x,y,g)
    getvalue(g)
end
function (op::HHHgreen)(x,y,g)
    
    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green =  op.α * exp(-op.γ*R)*(iR*i4pi)
    green*op.op(x,y,g)
end
mydot(a::SVector{N,<:SVector},b::Base.RefValue) where {N} = dot.(a,b)
function mydot(a::SVector,b::Base.RefValue) 
    a.*b
end
function (op::HHHgradgreen)(x,y,g)
    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green =  op.α * exp(-op.γ*R)*(iR*i4pi)
    gradgreen = -(op.γ + iR) * green * (iR * r)
    mydot(op.op(x,y,g),Ref(gradgreen))
end

function (op::HHHgradgreenCross)(x,y,g)
    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = op.α * exp(-op.γ*R)*(iR*i4pi)
    gradgreen = -(op.γ + iR) * green * (iR * r)
    cross.(Ref(gradgreen),op.op(x,y,g))
end
identify(::VectorSurfaceSpace) = HHHvector()
identify(::ScalarSurfaceSpace) = HHHscalar()
kernel(op::HHHtestOperator) = kernel(op.op)
kernel(op::HHHkernelOperator) = op

function defaultquadstrat(op::HHHOperator,testspace::Space,basisspace::Space)
@assert identify(basisspace) == inversemap(op,identify(testspace))
defaultquadstrat(op,identify(testspace),identify(basisspace))
end
defaultquadstrat(op::HHHOperator,::HHHdata,::HHHdata) = DoubleNumSauterQstrat(6,7,5,5,4,3) #TODO implement strategy that is stronger for nearby triangles


sign_upon_permutation(op::HHHIdentity,I,J) = 1
sign_upon_permutation(op::HHHkernelOperator,I,J) = sign_upon_permutation(op.op,I,J) 
sign_upon_permutation(op::Union{HHHNtestCross,HHHNtestDot},I,J) = Combinatorics.levicivita(I)*sign_upon_permutation(op.op,I,J)
sign_upon_permutation(op::Union{HHHNbasisnormal,HHHNbasisCross,HHHNbasisdot},I,J) = Combinatorics.levicivita(J)*sign_upon_permutation(op.op,I,J)

