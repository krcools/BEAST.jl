import LinearAlgebra: ×, ⋅
abstract type HHHLocalOperator <: LocalOperator end
abstract type HHHOperator <: IntegralOperator end
abstract type HHHtestOperator{U} <: HHHOperator end
abstract type HHHbasisOperator <: HHHOperator end
abstract type HHHkernelOperator{T,K,U} <: HHHOperator end
# scalartype(op::HHHOperator{T,K}) where {T, K <: Nothing} = T
# scalartype(op::HHHOperator{T,K}) where {T, K} = promote_type(T, K)
# scalartype(op::HHHOperator{T,K}) where {T, K <: Nothing} = T
# scalartype(op::HHHOperator) = ComplexF64
scalartype(op::HHHLocalOperator) = Union{}

scalartype(op::HHHkernelOperator{T,K,U}) where {T, K,U} = promote_type(T, K)
scalartype(op::HHHtestOperator{U}) where {U} = scalartype(U)

scalartype(op::Type{<: HHHtestOperator{U}}) where {U} =scalartype(U)

scalartype(op::Type{<: HHHkernelOperator{T,K,U}}) where {T,K,U} = promote_type(T, K)
#const i4pi = 1 / (4pi)
struct HHHgreen{T,K,U} <: HHHkernelOperator{T,K,U}
    α::K
    γ::T
    op::U
end

struct HHHgradgreen{T,K,U} <: HHHkernelOperator{T,K,U}
    α::K
    γ::T
    op::U
end

struct HHHgradgreenCross{T,K,U} <: HHHkernelOperator{T,K,U}
    α::K
    γ::T
    op::U
end

abstract type HHHdata end
struct HHHvector <:HHHdata end
struct HHHscalar <:HHHdata end


struct HHHNtestCross{U} <: HHHtestOperator{U}
    op::U
end


struct HHHNbasisCross{U} <: HHHbasisOperator
    op::U
end

struct HHHNtestDot{U} <: HHHtestOperator{U}
    op::U
end

struct HHHNbasisnormal{U} <: HHHbasisOperator
    op::U
end
struct HHHNbasisdot{U} <: HHHbasisOperator
    op::U
end
struct HHHIdentity <: HHHbasisOperator
end
struct HHHDivergence <: HHHbasisOperator
end


mutable struct HHHNtestCrossLocal{U <: HHHLocalOperator} <: HHHLocalOperator
    op::U
end


mutable struct HHHNbasisCrossLocal{U <: HHHLocalOperator} <: HHHLocalOperator
    op::U
end

mutable struct HHHNtestDotLocal{U <: HHHLocalOperator} <: HHHLocalOperator
    op::U
end

mutable struct HHHNbasisnormalLocal{U <: HHHLocalOperator} <: HHHLocalOperator
    op::U
end
mutable struct HHHNbasisdotLocal{U <: HHHLocalOperator} <: HHHLocalOperator
    op::U
end
mutable struct HHHIdentityLocal <: HHHLocalOperator
end
mutable struct HHHDivergenceLocal <: HHHLocalOperator end

hhhidentity() = HHHIdentity()
hhhdivergence() = HHHDivergence()

export HHH
struct BasisFunction end
struct DivBasisFunction end
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
inversemap(::HHHDivergence,::HHHscalar) = HHHvector()
##### Operator building

strace(op::Union{HHHkernelOperator,HHHtestOperator}) = HHHNtestCross(op)
ntrace(op::Union{HHHkernelOperator,HHHtestOperator}) = HHHNtestDot(op)
×(op::HHHgradgreen,::Nothing) = HHHgradgreenCross(op.α,op.γ,op.op)
#×(op::HHHgradgreen,::NormalVector) = HHHgradgreenCross(op.α,op.γ,op.op)(HHHNbasisnormal(hhhidentity()))
×(::NormalVector,op::Union{HHHkernelOperator,HHHtestOperator}) = strace(op)
⋅(::NormalVector,op::Union{HHHkernelOperator,HHHtestOperator}) = ntrace(op)
*(::NormalVector,::BasisFunction) = HHHNbasisnormal(hhhidentity())
⋅(::Nabla,::BasisFunction) = hhhdivergence()
*(::NormalVector,::DivBasisFunction) = HHHNbasisnormal(hhhdivergence())
⋅(::NormalVector,::BasisFunction) = HHHNbasisdot(hhhidentity())
×(::NormalVector,::BasisFunction) = HHHNbasisCross(hhhidentity())
×(::NormalVector,op::HHHbasisOperator) = HHHNbasisCross(op)
⋅(::NormalVector,op::HHHbasisOperator) = HHHNbasisdot(op)


function (op::HHHkernelOperator)(Basisop::HHHbasisOperator) 
    return (Base.typename(typeof(op)).wrapper)(op.α,op.γ,Basisop)
end


### TODO make sure this function can be used, unpack everything and make new structure.
function (op::HHHOperator)(BasisOperator::HHHbasisOperator)
    op.op = BasisOperator
end

# ##### to generate integrand
# extract(a::Type{HHHNbasisCross{T}}) where {T} = [HHHNbasisCross;extract(T)]
# extract(a::Type{HHHNbasisdot{T}}) where {T} = [HHHNbasisdot;extract(T)]
# extract(a::Type{HHHNbasisnormal{T}}) where {T} = [HHHNbasisnormal;extract(T)]
# extract(a::Type{HHHNtestCross{T}}) where {T} = [HHHNtestCross;extract(T)]
# extract(a::Type{HHHNtestDot{T}}) where {T} = [HHHNtestDot;extract(T)]
# extract(a::Type{HHHgradgreen{T}}) where {T} = [HHHgradgreen;extract(T)]
# extract(a::Type{HHHgradgreenCross{T}}) where {T} = [HHHgradgreenCross;extract(T)]
# extract(a::Type{HHHgreen{T}}) where {T} = [HHHgreen;extract(T)]
# extract(a::Type{HHHIdentity}) = HHHIdentity
# extract(a::Type{HHHDivergence}) = HHHDivergence


# @generated function hhhintegrand(op,x,y,g)
#     println("integrand generated")
#     typelist = extract(typeof(op))



# end



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
function (op::HHHDivergence)(x,y,g)
    getdivergence(g)
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


sign_upon_permutation(op::Union{HHHIdentity,HHHDivergence},I,J) = 1
sign_upon_permutation(op::HHHkernelOperator,I,J) = sign_upon_permutation(op.op,I,J) 
sign_upon_permutation(op::Union{HHHNtestCross,HHHNtestDot},I,J) = Combinatorics.levicivita(I)*sign_upon_permutation(op.op,I,J)
sign_upon_permutation(op::Union{HHHNbasisnormal,HHHNbasisCross,HHHNbasisdot},I,J) = Combinatorics.levicivita(J)*sign_upon_permutation(op.op,I,J)

"""
(number of test normals,number of trial normals)
"""
function number_of_normals(op::HHHOperator) end 


number_of_normals(op::Union{HHHNtestCross,HHHNtestCrossLocal,HHHNtestDot,HHHNtestDotLocal}) = [1,0]+number_of_normals(op.op) 
number_of_normals(op::HHHkernelOperator) = number_of_normals(op.op)
number_of_normals(op::Union{HHHNbasisCross,HHHNbasisCrossLocal,HHHNbasisdot,HHHNbasisdotLocal,HHHNbasisnormal,HHHNbasisnormalLocal}) = [0,1]+number_of_normals(op.op)
number_of_normals(op::Union{HHHIdentity,HHHIdentityLocal,HHHDivergence,HHHDivergenceLocal}) = [0,0]



#### trace definition



alpha(op::HHHtestOperator) = alpha(op.op)
alpha(op::HHHkernelOperator) = op.α

hhhntestcrosslocal(op::ZeroOperator) = op
hhhntestcrosslocal(op::HHHLocalOperator) = HHHNtestCrossLocal(op)
hhhntestcrosslocal(op::Union{HHHNbasisnormalLocal}) = ZeroOperator()

hhhntestdotlocal(op::ZeroOperator) = op
hhhntestdotlocal(op::HHHLocalOperator) = HHHNtestDotLocal(op)
hhhntestdotlocal(op::Union{HHHNtestCrossLocal,HHHNbasisCrossLocal}) = ZeroOperator()
hhhntestdotlocal(op::HHHNbasisnormalLocal) = localoperator(op.op)

hhhnbasiscrosslocal(op::ZeroOperator) = op
hhhnbasiscrosslocal(op::HHHLocalOperator) = HHHNbasisCrossLocal(op)
hhhnbasiscrosslocal(op::Union{HHHNbasisnormalLocal}) = ZeroOperator()

hhhnbasisdotlocal(op::ZeroOperator) = op
hhhnbasisdotlocal(op::HHHLocalOperator) = HHHNbasisdotLocal(op)
hhhnbasisdotlocal(op::HHHDivergenceLocal) = HHHNbasisnormalLocal(op)
hhhnbasisdotlocal(op::Union{HHHNtestCrossLocal,HHHNbasisCrossLocal}) = ZeroOperator()
function hhhnbasisdotlocal(op::HHHIdentityLocal)
    @warn "assumed basisfunction is tangential to surface"
ZeroOperator()
end

hhhnbasisdotlocal(op::HHHNbasisnormalLocal) = op.op

hhhnbasisnormal(op::ZeroOperator) = op
hhhnbasisnormal(op::HHHLocalOperator) = HHHNbasisnormalLocal(op)

localoperator(op::HHHNtestCross) = hhhntestcrosslocal(localoperator(op.op))
localoperator(op::HHHNtestDot) = hhhntestdotlocal(localoperator(op.op))
localoperator(op::HHHNbasisCross) = hhhnbasiscrosslocal(localoperator(op.op))
localoperator(op::HHHNbasisdot) = hhhnbasisdotlocal(localoperator(op.op))
localoperator(op::HHHNbasisnormal) = hhhnbasisnormal(localoperator(op.op))
localoperator(op::HHHIdentity) = HHHIdentityLocal()
localoperator(op::HHHDivergence) = HHHDivergenceLocal()

localoperator(op::HHHgreen) = ZeroOperator()
localoperator(op::HHHgradgreen) = hhhnbasisdotlocal(localoperator(op.op))
localoperator(op::HHHgradgreenCross) = hhhnbasiscrosslocal(localoperator(op.op))

function trace(op::HHHOperator,sign)
    l = localoperator(op)
    if typeof(l) == HHHNbasisnormalLocal{HHHIdentityLocal} || typeof(l) == HHHNbasisnormalLocal{HHHDivergenceLocal}
        @warn "assumed test function is tangential to surface"
        l = ZeroOperator()
    end   

    localop = 1/2*sign*alpha(op)*localoperator(op)
    
    return op + localop
end
function normalorient(op::Union{HHHOperator,HHHLocalOperator},signtest,signtrial)
test,trial = number_of_normals(op)
signtest^test*signtrial^trial*op
end

##### defining integrand of those local operators:


function integrand(op::HHHLocalOperator,kernel,x,g,f,nt,nb)
  
dot(g.value,op(x,f,nt,nb))
end

function (op::HHHNtestCrossLocal)(x,f,nt,nb)
    nt×op.op(x,f,nt,nb)
end
function (op::HHHNtestDotLocal)(x,f,nt,nb)
    dot(nt,op.op(x,f,nt,nb))
end
function (op::HHHNbasisCrossLocal)(x,f,nt,nb)
    nb×op.op(x,f,nt,nb)
end
function (op::HHHNbasisdotLocal)(x,f,nt,nb)
    dot(nb,op.op(x,f,nt,nb))
end
function (op::HHHNbasisnormalLocal)(x,f,nt,nb)
    nb*op.op(x,f,nt,nb)
end
function (op::HHHIdentityLocal)(x,f,nt,nb)
    return f.value
end
function (op::HHHDivergenceLocal)(x,f,nt,nb)
    return f.divergence
end

function cellinteractions(biop::HHHLocalOperator, trefs::U, brefs::V, cell,tcell,bcell) where {U<:RefSpace{T},V<:RefSpace{T}} where {T}
    qr = qrtot
    num_tshs = length(qr[1][3])
    num_bshs = length(qr[1][4])

    zlocal = zeros(T, num_tshs, num_bshs)
    nt = normal(tcell)
    nb = normal(bcell)
    for q in qr

        w, mp, tvals, bvals = q[1], q[2], q[3], q[4]
        j = w * jacobian(mp)
        # kernel = kernelvals(biop, mp)
        kernel = nothing
        for m in 1 : num_tshs
            tval = tvals[m]

            for n in 1 : num_bshs
                bval = bvals[n]

                igd = integrand(biop, kernel, mp, tval, bval,nt,nb)
                zlocal[m,n] += j * igd

            end
        end
    end

    return zlocal
end
function cellinteractions_matched!(zlocal, biop::HHHLocalOperator, trefs, brefs, cell, qr,tcell=nothing,bcell=nothing)

    num_tshs = length(qr[1][3])
    num_bshs = length(qr[1][4])
    nt = normal(tcell)
    nb = normal(bcell)
    # zlocal = zeros(Float64, num_tshs, num_bshs)
    for q in qr

        w, mp, tvals, bvals = q[1], q[2], q[3], q[4]
        j = w * jacobian(mp)
        #kernel = kernelvals(biop, mp)
        kernel = nothing
        nt = normal(tcell)
        nb = normal(bcell)
        for n in 1 : num_bshs
            bval = bvals[n]
            for m in 1 : num_tshs
                tval = tvals[m]

                igd = integrand(biop, kernel, mp, tval, bval,nt,nb)
                zlocal[m,n] += j * igd
            end
        end
    end

    return zlocal
end