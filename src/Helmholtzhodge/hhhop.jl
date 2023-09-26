import LinearAlgebra: ×, ⋅
abstract type HHHLocalOperator <: LocalOperator end
abstract type HHHOperator{T,K} <: IntegralOperator end
abstract type HHHtestOperator{T,K} <: HHHOperator{T,K} end
abstract type HHHbasisOperator{T,K} <: HHHOperator{T,K} end
abstract type HHHkernelOperator{T,K} <: HHHOperator{T,K} end
# scalartype(op::HHHOperator{T,K}) where {T, K <: Nothing} = T
# scalartype(op::HHHOperator{T,K}) where {T, K} = promote_type(T, K)
scalartype(op::HHHOperator{T,K}) where {T, K <: Nothing} = Complex
scalartype(op::HHHOperator{T,K}) where {T, K} = Complex
scalartype(op::HHHLocalOperator) = Union{}

#const i4pi = 1 / (4pi)
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
mutable struct HHHDivergence{T,K} <: HHHbasisOperator{T,K}
end


mutable struct HHHNtestCrossLocal{T,K} <: HHHtestOperator{T,K}
    op::HHHOperator{T,K}
end


mutable struct HHHNbasisCrossLocal <: HHHLocalOperator
    op::HHHLocalOperator
end

mutable struct HHHNtestDotLocal <: HHHLocalOperator
    op::HHHLocalOperator
end

mutable struct HHHNbasisnormalLocal <: HHHLocalOperator
    op::HHHLocalOperator
end
mutable struct HHHNbasisdotLocal <: HHHLocalOperator
    op::HHHLocalOperator
end
mutable struct HHHIdentityLocal <: HHHLocalOperator
end
mutable struct HHHDivergenceLocal <: HHHLocalOperator end

hhhidentity() = HHHIdentity{Complex,Complex}()
hhhdivergence() = HHHDivergence{Complex,Complex}()

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
×(op::HHHgradgreen,::NormalVector) = HHHgradgreenCross(op.α,op.γ,op.op)(HHHNbasisnormal(hhhidentity()))
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
hhhntestdotlocal(op::ZeroOperator) = op
hhhntestdotlocal(op::HHHLocalOperator) = HHHNtestDotLocal(op)

hhhnbasiscrosslocal(op::ZeroOperator) = op
hhhnbasiscrosslocal(op::HHHLocalOperator) = HHHNbasisCrossLocal(op)

hhhnbasisdotlocal(op::ZeroOperator) = op
hhhnbasisdotlocal(op::HHHLocalOperator) = HHHNbasisdotLocal(op)

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
localoperator(op::HHHgradgreen) = hhhnbasisdotlocal(op.op)
localoperator(op::HHHgradgreenCross) = hhhnbasiscrosslocal(op.op)

function trace(op::HHHOperator,sign)

    localop = 1/2*sign*alpha(op)*localoperator(op)
    
    return op + localop
end
function normalorient(op::Union{HHHOperator,HHHLocalOperator},signtest,signtrial)
test,trial = number_of_normals(op)
signtest^test*signtrial^trial*op
end

##### defining integrand of those local operators:

function integrand(op::HHHLocalOperator,kernel,x,g,f,nt,nb)
dot(g[1],op(x,f[1],nt,nb))
end

function (op::HHHNtestCrossLocal)(x,f,nt,nb)
    nt×op.op(x,f,nt,nb)
end
function (op::HHHNtestDotLocal)(x,f,nt,nb)
    dot(nt,op.op(x,f,nt,nb))
end
function (op::HHHNbasisCross)(x,f,nt,nb)
    nb×op.op(x,f,nt,nb)
end
function (op::HHHNbasisdotLocal)(x,f,nt,nb)
    dot(nb,op.op(x,f,nt,nb))
end
function (op::HHHNbasisnormalLocal)(x,f,nt,nb)
    nb*op.op(x,f,nt,nb)
end
function (op::HHHIdentityLocal)(x,f,nt,nb)
    return f
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
        kernel = kernelvals(biop, mp)

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
        kernel = kernelvals(biop, mp)
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