import Base: *, div
import LinearAlgebra: ×, ⋅


#TODO add ref of cross etc to eliminate brackets in notation.

#const i4pi = 1 / (4pi)
abstract type ComposedOperatorLocal <: LocalOperator end
abstract type ComposedOperatorIntegral <: IntegralOperator end
abstract type Kernel <: ComposedOperatorIntegral end
ComposedOperator = Union{ComposedOperatorIntegral,ComposedOperatorLocal}

struct BasisFunction <: ComposedOperatorLocal end
struct DivBasisFunction <: ComposedOperatorLocal end
const B = BasisFunction()
Base.div(::BasisFunction) = DivBasisFunction()

export B
export BasisFunction


struct TestNormal <: ComposedOperatorLocal end
struct TrialNormal <: ComposedOperatorLocal end
struct TraceDirection <: ComposedOperatorLocal end


const nt = TestNormal()
const nb = TrialNormal()
export nt
export nb
struct Potential{T,M} <: AbstractOperator
    operator::T
    surface::M
end
Potential(a) = Potential(a,nothing)
surface(p::Potential) = p.surface
function Potential(linop::LinearCombinationOfOperators,surface::TraceMesh)
    out = ZeroOperator()
    for (op,coeff) in zip(linop.ops,linop.coeffs)
        out += coeff*Potential(op,surface)
    end
    return out
end

struct TraceOperator{T,M} <: AbstractOperator
    operator::T
    surface::M
end
function TraceOperator(linop::LinearCombinationOfOperators,surf::TraceMesh)
    out = ZeroOperator()
    for (op,coeff) in zip(linop.ops,linop.coeffs)
        out += coeff*TraceOperator(op,surf)
    end
    return out
end
function TraceOperator(linop::LinearCombinationOfOperators)
    out = ZeroOperator()
    for (op,coeff) in zip(linop.ops,linop.coeffs)
        out += coeff*TraceOperator(op)
    end
    return out
end
surface(p::TraceOperator) = p.surface

struct TimesLocal{U,V} <: ComposedOperatorLocal
    lhs::U
    rhs::V
end

struct DotLocal{U,V} <: ComposedOperatorLocal
    lhs::U
    rhs::V
end

struct CrossLocal{U,V} <: ComposedOperatorLocal
    lhs::U
    rhs::V
end

struct TimesIntegral{U,V} <: ComposedOperatorIntegral
    lhs::U
    rhs::V
end

struct DotIntegral{U,V} <: ComposedOperatorIntegral
    lhs::U
    rhs::V
end

struct CrossIntegral{U,V} <: ComposedOperatorIntegral
    lhs::U
    rhs::V
end

function TimesLocal(lhs::LinearCombinationOfOperators,rhs::Union{Union{ComposedOperator,NormalVector},LinearCombinationOfOperators})
    out = ZeroOperator()
    for (op,coeff) in zip(lhs.ops,lhs.coeffs)
        out += coeff*TimesLocal(op,rhs)
    end
    return out
end

function TimesLocal(lhs::Union{ComposedOperator,NormalVector},rhs::LinearCombinationOfOperators)
    out = ZeroOperator()
    for (op,coeff) in zip(rhs.ops,rhs.coeffs)
        out += coeff*TimesLocal(lhs,op)
    end
    return out
end

function DotLocal(lhs::LinearCombinationOfOperators,rhs::Union{Union{ComposedOperator,NormalVector},LinearCombinationOfOperators})
    out = ZeroOperator()
    for (op,coeff) in zip(lhs.ops,lhs.coeffs)
        out += coeff*DotLocal(op,rhs)
    end
    return out
end

function DotLocal(lhs::Union{ComposedOperator,NormalVector},rhs::LinearCombinationOfOperators)
    out = ZeroOperator()
    for (op,coeff) in zip(rhs.ops,rhs.coeffs)
        out += coeff*DotLocal(lhs,op)
    end
    return out
end

function CrossLocal(lhs::LinearCombinationOfOperators,rhs::Union{Union{ComposedOperator,NormalVector},LinearCombinationOfOperators})
    out = ZeroOperator()
    for (op,coeff) in zip(lhs.ops,lhs.coeffs)
        out += coeff*CrossLocal(op,rhs)
    end
    return out
end

function CrossLocal(lhs::Union{ComposedOperator,NormalVector},rhs::LinearCombinationOfOperators)
    out = ZeroOperator()
    for (op,coeff) in zip(rhs.ops,rhs.coeffs)
        out += coeff*CrossLocal(lhs,op)
    end
    return out
end
get_constructor(::TimesLocal) = TimesLocal
get_constructor(::TimesIntegral) = TimesLocal
get_constructor(::DotLocal) = DotLocal
get_constructor(::DotIntegral) = DotLocal
get_constructor(::CrossLocal) = CrossLocal
get_constructor(::CrossIntegral) = CrossLocal

OperationsLocal = Union{TimesLocal,DotLocal,CrossLocal}
OperationsIntegral = Union{TimesIntegral,DotIntegral,CrossIntegral}
Operations = Union{OperationsIntegral,OperationsLocal}

TimesLocal(a::ComposedOperatorIntegral,b::Union{ComposedOperatorLocal,NormalVector}) = TimesIntegral(a,b)
TimesLocal(a::Union{ComposedOperatorLocal,NormalVector},b::ComposedOperatorIntegral) = TimesIntegral(a,b)
TimesLocal(a::ComposedOperatorIntegral,b::ComposedOperatorIntegral) = TimesIntegral(a,b)

DotLocal(a::ComposedOperatorIntegral,b::Union{ComposedOperatorLocal,NormalVector}) = DotIntegral(a,b)
DotLocal(a::Union{ComposedOperatorLocal,NormalVector},b::ComposedOperatorIntegral) = DotIntegral(a,b)
DotLocal(a::ComposedOperatorIntegral,b::ComposedOperatorIntegral) = DotIntegral(a,b)

CrossLocal(a::ComposedOperatorIntegral,b::Union{ComposedOperatorLocal,NormalVector}) = CrossIntegral(a,b)
CrossLocal(a::Union{ComposedOperatorLocal,NormalVector},b::ComposedOperatorIntegral) = CrossIntegral(a,b)
CrossLocal(a::ComposedOperatorIntegral,b::ComposedOperatorIntegral) = CrossIntegral(a,b)

×(a::Union{ComposedOperator,NormalVector,LinearCombinationOfOperators},b::Union{ComposedOperator,NormalVector,LinearCombinationOfOperators}) = CrossLocal(a,b)
⋅(a::Union{ComposedOperator,NormalVector,LinearCombinationOfOperators},b::Union{ComposedOperator,NormalVector,LinearCombinationOfOperators}) = DotLocal(a,b)
*(a::Union{ComposedOperator,NormalVector,LinearCombinationOfOperators},b::Union{ComposedOperator,NormalVector,LinearCombinationOfOperators}) = TimesLocal(a,b)



scalartype(op::Union{TestNormal,TrialNormal,BasisFunction,DivBasisFunction,TraceDirection}) = Float16
scalartype(op::Operations) = promote_type(scalartype(op.lhs),scalartype(op.rhs))
scalartype(op::Potential) = scalartype(op.operator)
scalartype(op::TraceOperator) = scalartype(op.operator)

# function count_test_normals(op::Operations)
#     count_test_normals(op.lhs) + count_test_normals(op.rhs)
# end
# count_test_normals(op::Kernel) = 0
# count_test_normals(op::TestNormal) = 1
# count_test_normals(op::TrialNormal) = 0

# function count_trial_normals(op::Operations)
#     count_trial_normals(op.lhs) + count_trial_normals(op.rhs)
# end
# count_trial_normals(op::Kernel) = 0
# count_trial_normals(op::TestNormal) = 0
# count_trial_normals(op::TrialNormal) = 1

function replace_normal_by_testnormal(op::Operations)
    get_constructor(op)(replace_normal_by_testnormal(op.lhs),replace_normal_by_testnormal(op.rhs))
end
replace_normal_by_testnormal(op::Kernel) = op
replace_normal_by_testnormal(op::Union{TestNormal,TrialNormal,BasisFunction,DivBasisFunction}) = op
replace_normal_by_testnormal(op::NormalVector) = TestNormal()

function replace_normal_by_trialnormal(op::Operations)
    get_constructor(op)(replace_normal_by_trialnormal(op.lhs),replace_normal_by_trialnormal(op.rhs))
end
replace_normal_by_trialnormal(op::Kernel) = op
replace_normal_by_trialnormal(op::Union{TestNormal,TrialNormal,BasisFunction,DivBasisFunction}) = op
replace_normal_by_trialnormal(op::NormalVector) = TrialNormal()



function build_potential(op::ComposedOperator,surface::CompScienceMeshes.AbstractMesh)
    newop = replace_normal_by_trialnormal(op)
    Potential(newop,TraceMesh(surface))
end
function build_potential(op::ComposedOperator,surface::TraceMesh)
    newop = replace_normal_by_trialnormal(op)
    Potential(newop,surface)
end
function build_potential(op::ComposedOperator)
    newop = replace_normal_by_trialnormal(op)
    Potential(newop,nothing)
end


function γ(op::Operations)
    return get_constructor(op)(γ(op.lhs),γ(op.rhs))
end

γ(op::Union{TestNormal,TrialNormal,BasisFunction,DivBasisFunction}) = op

# function check_if_coincide(a,b) 
#     @warn "all meshes coincide"
#     return true
# end

function γ(op::Potential,surface::CompScienceMeshes.AbstractMesh,sign::Int)# sign + if according to normal on surface, - otherwise
#    check_if_coincide(op.surface,surface) || return op.operator
    newop = nt × (γ(op.operator) × nt)
    direction = []
    for i in 1:numcells(surface)
        c = chart(surface,i)
        push!(direction,sign*normal(c)/3)
    end
    return TraceOperator(Potential(newop,op.surface),TraceMesh(surface,SVector{length(direction),typeof(direction[1])}(direction)))
end
function γₜ(op::Potential,sign::Int)# sign + if according to normal on surface, - otherwise
#    check_if_coincide(op.surface,surface) || return op.operator
    newop = nt × (γ(op.operator) × nt)
    return TraceOperator(Potential(newop,op.surface),sign)
end
function γₛ(op::Potential,sign::Int)# sign + if according to normal on surface, - otherwise
#    check_if_coincide(op.surface,surface) || return op.operator
    newop = nt × γ(op.operator)
    return TraceOperator(Potential(newop,op.surface),sign)
end
function γₙ(op::Potential,sign::Int)# sign + if according to normal on surface, - otherwise
#    check_if_coincide(op.surface,surface) || return op.operator
    newop = nt ⋅ γ(op.operator)
    return TraceOperator(Potential(newop,op.surface),sign)
end
function τ(op::Potential,sign::Int)# sign + if according to normal on surface, - otherwise
#    check_if_coincide(op.surface,surface) || return op.operator
    newop = γ(op.operator)
    return TraceOperator(Potential(newop,op.surface),sign)
end

function γ(op::Potential,surface::TraceMesh,sign::Int)# sign + if according to normal on surface, - otherwise
#    check_if_coincide(op.surface,surface) || return op.operator
    newop = nt × (γ(op.operator) × nt)
    direction = []
    for i in 1:numcells(surface)
        c = chart(surface,i)
        push!(direction,sign*normal(c)/3)
    end
    return TraceOperator(Potential(newop,op.surface),SVector{length(direction),typeof(direction[1])}(direction)+surface)
end


function γₛ(op::Potential,surface::CompScienceMeshes.AbstractMesh,sign::Int)# sign + if according to normal on surface, - otherwise
#    check_if_coincide(op.surface,surface) || return op.operator
    newop = nt × γ(op.operator)
    # direction = []
    # for i in 1:numcells(surface)
    #     c = chart(surface,i)
    #     push!(direction,sign*normal(c)/3)
    # end
    charts = chart.(Ref(surface),1:numcells(surface))
    direction = sign/3*normal.(charts)
    return TraceOperator(Potential(newop,op.surface),TraceMesh(surface,direction))
end
function γₛ(op::Potential,surface::TraceMesh,sign::Int)# sign + if according to normal on surface, - otherwise
#    check_if_coincide(op.surface,surface) || return op.operator
    newop = nt × γ(op.operator)
    # direction = []
    # for i in 1:numcells(surface)
    #     c = chart(surface,i)
    #     push!(direction,sign*normal(c)/3)
    # end
    charts = chart.(Ref(surface),1:numcells(surface))
    direction = sign/3*normal.(charts)
    return TraceOperator(Potential(newop,op.surface),direction+surface)
end

function γₜ(op::Potential,surface::CompScienceMeshes.AbstractMesh,sign::Int)# sign + if according to normal on surface, - otherwise
    #    check_if_coincide(op.surface,surface) || return op.operator
        newop = (nt × γ(op.operator)) × nt
        charts = chart.(Ref(surface),1:numcells(surface))
        direction = sign/3*normal.(charts)
        return TraceOperator(Potential(newop,op.surface),TraceMesh(surface,direction))
    end
    function γₜ(op::Potential,surface::TraceMesh,sign::Int)# sign + if according to normal on surface, - otherwise
    #    check_if_coincide(op.surface,surface) || return op.operator
        newop =  (nt × γ(op.operator)) × nt
        charts = chart.(Ref(surface),1:numcells(surface))
        direction = sign/3*normal.(charts)
        return TraceOperator(Potential(newop,op.surface),direction+surface)
    end

function γₙ(op::Potential,surface::CompScienceMeshes.AbstractMesh,sign::Int)# sign + if according to normal on surface, - otherwise
#    check_if_coincide(op.surface,surface) || return op.operator
    newop = nt ⋅ γ(op.operator)
    charts = chart.(Ref(surface),1:numcells(surface))
    direction = sign/3*normal.(charts)
    return TraceOperator(Potential(newop,op.surface),TraceMesh(surface,direction))
end
function γₙ(op::Potential,surface::TraceMesh,sign::Int)# sign + if according to normal on surface, - otherwise
#    check_if_coincide(op.surface,surface) || return op.operator
    newop = nt ⋅ γ(op.operator)
    charts = chart.(Ref(surface),1:numcells(surface))
    direction = sign/3*normal.(charts)
    return TraceOperator(Potential(newop,op.surface),direction+surface)
end

function τ(op::Potential,surface::CompScienceMeshes.AbstractMesh,sign::Int)# sign + if according to normal on surface, - otherwise
#    check_if_coincide(op.surface,surface) || return op.operator
    newop = γ(op.operator)
    charts = chart.(Ref(surface),1:numcells(surface))
    direction = sign/3*normal.(charts)
    return TraceOperator(Potential(newop,op.surface),TraceMesh(surface,direction))
end
function τ(op::Potential,surface::TraceMesh,sign::Int)# sign + if according to normal on surface, - otherwise
#    check_if_coincide(op.surface,surface) || return op.operator
    newop = γ(op.operator)
    charts = chart.(Ref(surface),1:numcells(surface))
    direction = sign/3*normal.(charts)
    return TraceOperator(Potential(newop,op.surface),direction+surface)
end

export γₛ, γₜ, γₙ, γ, τ
#### Define new greens functions below

function greenhh3d(;
    gamma=nothing,
    wavenumber=nothing)

    gamma, wavenumber = gamma_wavenumber_handler(gamma, wavenumber)
    @assert gamma !== nothing
    GreenHH3D(gamma)
end
function gradgreenhh3d(;
    gamma=nothing,
    wavenumber=nothing)

    gamma, wavenumber = gamma_wavenumber_handler(gamma, wavenumber)
    @assert gamma !== nothing
    GradGreenHH3D(gamma)
end

struct GreenHH3D{T} <: Kernel
    gamma::T
end

struct GradGreenHH3D{T} <: Kernel
    gamma::T
end


function (op::GreenHH3D)(x,y,g)
    gamma = op.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-gamma*R)*(i4pi*iR)
    Ref(green)
end
function (op::GreenHH3D)(x::Union{SVector,Vector},y,g)
    gamma = op.gamma

    r = x - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-gamma*R)*(i4pi*iR)
    Ref(green)
end

function (op::GradGreenHH3D)(x,y,g)
    gamma = op.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-gamma*R)*(iR*i4pi)
    gradgreen = -(gamma + iR) * green * (iR * r)
    tt = Ref(gradgreen)

    if maximum(abs.(gradgreen)) === NaN
        display(R)
    end
    return tt
end
function (op::GradGreenHH3D)(x::Union{SVector,Vector},y,g)
    gamma = op.gamma

    r = x - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-gamma*R)*(iR*i4pi)
    gradgreen = -(gamma + iR) * green * (iR * r)
    tt = Ref(gradgreen)

    if maximum(abs.(gradgreen)) === NaN
        display(R)
    end
    return tt
end

γ(op::GreenHH3D) = op
γ(op::GradGreenHH3D) = op + 1/2*TraceDirection()

grad(G::GreenHH3D) = GradGreenHH3D(G.gamma)

function (::Nabla)(G::Kernel)
    grad(G)
end

scalartype(G::GreenHH3D{T}) where {T} = T
scalartype(G::GradGreenHH3D{T}) where {T} = T


struct DoubleIntegralR{T} <: Kernel end

function (op::DoubleIntegralR)(x,y,g)
    r = cartesian(x) - cartesian(y)
    Ref(r)
end
scalartype(G::DoubleIntegralR{T}) where {T} = T
γ(op::DoubleIntegralR) = op


function (op::Union{TimesIntegral,TimesLocal})(x,y,g)
    op.lhs(x,y,g).*op.rhs(x,y,g)
end
function (op::Union{DotIntegral,DotLocal})(x,y,g)
    transpose.(op.lhs(x,y,g)).*op.rhs(x,y,g)
end
function (op::Union{CrossIntegral,CrossLocal})(x,y,g)
    cross.(op.lhs(x,y,g),op.rhs(x,y,g))
end

function (op::TrialNormal)(x,y,g)
    Ref(normal(y))
end
function (op::TestNormal)(x,y,g)
    Ref(normal(x))
end
function (op::TraceDirection)(x,y,g)
    #println(direction(x)-direction(y))
    Ref(approx_sign(dot(normal(x),(direction(y)-direction(x))))*normal(x))
end
function (op::BasisFunction)(x,y,g)
    getvalue(g)
end
function (op::DivBasisFunction)(x,y,g)
    getdivergence(g)
end
function approx_sign(a::T; tol=eps(T)) where {T <: Number}
    abs(a) < tol && return zero(T)
    return sign(a)
end


function (igd::Integrand{<:ComposedOperatorIntegral})(x,y,f,g)
    _krondot(getvalue(f),igd.operator(x,y,g))
end

function integrand(op::ComposedOperatorLocal,kernel,x,y,f,g)
    _krondot(getvalue(f),op(x,y,g))
end
function integrand(op::ComposedOperator,kernel, y, f, x)
    r = op(y,x,f)
    return r
end
kernelvals(op::ComposedOperatorIntegral, y,x) = nothing
quaddata(op::ComposedOperatorIntegral,rs,els,qs::SingleNumQStrat) = quadpoints(rs,els,(qs.quad_rule,))
quadrule(op::ComposedOperatorIntegral,refspace,p,y,q,el,qdata,qs::SingleNumQStrat) = qdata[1,q]

defaultquadstrat(op::ComposedOperatorIntegral, basis) = SingleNumQStrat(6)
defaultquadstrat(op::ComposedOperatorIntegral,testspace::Space,trialspace::Space) = DoubleNumSauterQstrat(7,8,5,5,4,3) 
defaultquadstrat(op::ComposedOperatorLocal,testspace::Space,trialpsace::Space) = SingleNumQStrat(8)
#sign_upon_permutation(op::ComposedOperator,I,J) = Combinatorics.levicivita(I)^count_test_normals(op)*Combinatorics.levicivita(J)^count_trial_normals(op)


#normalorient(op::ComposedOperator,sign_test_normal,sign_trial_normal) = sign_test_normal^count_test_normals*sign_trial_normal^count_trial_normals


defaultquadstrat(op::TraceOperator,test,trial) = defaultquadstrat(op.operator,test,trial)
defaultquadstrat(op::Potential,test,trial) = defaultquadstrat(op.operator,test,trial)



# function assemble!(op::Potential, test_functions::Space, trial_functions::Space,
#     store, threading = Threading{:multi};
#     quadstrat = defaultquadstrat(op, test_functions, trial_functions))
#     # checks if surface given in basis is the same as surface given in potential
#     dsurf = surface(op)
#     surf = geometry(trial_functions)
#     @assert same_geometry(dsurf,surf)
#     nsurf = mesh(dsurf)
#     trial_functions = redefine_geometrie(trial_functions,TraceMesh(OrientedMesh(surf,nsurf),direction(dsurf)))

#     assemble!(op.operator, test_functions, trial_functions, store, threading;
#     quadstrat = quadstrat)
# end
function assemble!(op::Potential, test_functions::Space, trial_functions::Space,
    store, threading = Threading{:multi};
    quadstrat = defaultquadstrat(op, test_functions, trial_functions))
    # checks if surface given in basis is the same as surface given in potential
    dsurf = surface(op)
    surf = geometry(trial_functions)
    @assert same_geometry(dsurf,surf)
    nsurf = mesh(dsurf)
    trial_functions = redefine_geometrie(trial_functions,TraceMesh(surf,direction(dsurf)))

    assemble!(op.operator, test_functions, trial_functions, store, threading;
    quadstrat = quadstrat)
end
# function assemble!(op::Potential{T,Nothing}, test_functions::Space, trial_functions::Space,
#     store, threading = Threading{:multi};
#     quadstrat = defaultquadstrat(op, test_functions, trial_functions)) where {T}

#     surf = geometry(trial_functions)
#     trial_functions = redefine_geometrie(trial_functions,TraceMesh(OrientedMesh(surf,surf)))

#     assemble!(op.operator, test_functions, trial_functions, store, threading;
#     quadstrat = quadstrat)
# end
function assemble!(op::Potential{T,Nothing}, test_functions::Space, trial_functions::Space,
    store, threading = Threading{:multi};
    quadstrat = defaultquadstrat(op, test_functions, trial_functions)) where {T}

    surf = geometry(trial_functions)
    trial_functions = redefine_geometrie(trial_functions,TraceMesh(surf))

    assemble!(op.operator, test_functions, trial_functions, store, threading;
    quadstrat = quadstrat)
end

# function assemble!(op::TraceOperator, test_functions::Space, trial_functions::Space,
#     store, threading = Threading{:multi}; 
#     quadstrat = defaultquadstrat(op, test_functions, trial_functions))

#     dsurf = surface(op)
#     surf = geometry(test_functions)
#     println(typeof(surf))
#     println(typeof(mesh(dsurf)))
#     @assert same_geometry(dsurf,surf)
#     nsurf = mesh(dsurf)
    
#     test_functions = redefine_geometrie(test_functions,TraceMesh(OrientedMesh(surf,nsurf),direction(dsurf)))

#     assemble!(op.operator, test_functions, trial_functions, store, threading;
#     quadstrat = quadstrat)
# end

function assemble!(op::TraceOperator, test_functions::Space, trial_functions::Space,
    store, threading = Threading{:multi}; 
    quadstrat = defaultquadstrat(op, test_functions, trial_functions))

    dsurf = surface(op)
    surf = geometry(test_functions)
    println(typeof(surf))
    println(typeof(mesh(dsurf)))
    @assert same_geometry(dsurf,surf)
    nsurf = mesh(dsurf)
    
    test_functions = redefine_geometrie(test_functions,TraceMesh(surf,direction(dsurf)))

    assemble!(op.operator, test_functions, trial_functions, store, threading;
    quadstrat = quadstrat)
end
function assemble!(op::TraceOperator{T,Int}, test_functions::Space, trial_functions::Space,
    store, threading = Threading{:multi}; 
    quadstrat = defaultquadstrat(op, test_functions, trial_functions)) where {T}

    surface = geometry(test_functions)

    sign = op.surface
    direction = [sign*normal(c)/3 for c in chart.(Ref(surface),1:numcells(surface))]

    if typeof(surface) <: TraceMesh
        surf = surface + direction
    else
        surf = TraceMesh(surface,direction)
    end
    
    test_functions = redefine_geometrie(test_functions,surf)
    
    assemble!(op.operator, test_functions, trial_functions, store, threading;
    quadstrat = quadstrat)
end
