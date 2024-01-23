#abstract type MaxwellOperatorCQ3D{T,K} <: IntegralOperator end
#abstract type MaxwellOperator3DReg{T,K} <: MaxwellOperator3D{T,K} end

struct MWCQSingleLayer3D{T,U,V} <: BEAST.MaxwellOperator3D{T,U}
    Δt:: T
    α::U
    β::U
    K::V
    gamma::T
    Nt::V
    rkscheme::V
end

#scalartype(op::MaxwellOperatorCQ3D{T,K}) where {T, K <: Nothing} = T
#scalartype(op::MaxwellOperatorCQ3D{T,K}) where {T, K} = promote_type(T, K)

scalartype(op::MWCQSingleLayer3D{T,U}) where {T,U} = promote_type(T,U)
sign_upon_permutation(op::MWCQSingleLayer3D, I, J) = 1

function singlelayerCQ(dt, a, b, k, Nt)
    MWCQSingleLayer3D(dt, a, b, k, 1.0, Nt,1)
end

################################################################################
#
#  Kernel definitions
#
################################################################################

#const i4pi = 1 / (4pi)
function (igd::BEAST.Integrand{<:MWCQSingleLayer3D})(x,y,f,g)
    α = igd.operator.α
    β = igd.operator.β
    #γ = igd.operator.gamma
    K = igd.operator.K
    dt = igd.operator.Δt

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    ξ = R/dt
    #iR = 1 / R
    w_s = cq_weightBE(0,K,ξ,dt)
    w_hs = cq_weightBE(2,K,ξ,dt)

    αw_s = α * w_s/R
    βw_hs = β * w_hs/R

    _integrands(f,g) do fi,gj
        αw_s * dot(fi.value, gj.value) + βw_hs * dot(fi.divergence, gj.divergence)
    end
end

#= function (igd::Integrand{<:MWSingleLayer3DReg})(x,y,f,g)
    α = igd.operator.α
    β = igd.operator.β
    γ = igd.operator.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    γR = γ*R
    # iR = 1 / R
    green = (expm1(-γR) + γR - 0.5*γR^2) / (4pi*R)

    αG = α * green
    βG = β * green

    _integrands(f,g) do fi,gj
        αG * dot(fi.value, gj.value) + βG * dot(fi.divergence, gj.divergence)
    end
end =#

function rkstages(x) return 3 end
function allocatestorage(operator::MWCQSingleLayer3D, test_functions, trial_functions,
    storage_trait, longdelays_trait)

    Nt = operator.Nt
    p = rkstages(operator.rkscheme)

    T = promote_type(
        scalartype(operator)       ,
        scalartype(test_functions) ,
        scalartype(trial_functions),
    )
    Z = Array{T}(undef,
        numfunctions(test_functions),
        numfunctions(trial_functions),
        Nt,
        p
    )
    fill!(Z, 0)
    store(v,i,j,k,l) = (Z[i,j,k,l] += v)
    return ()->Z, store
end

function cq_weightBDF2(n,k,x,dt)
    if n>1
        return (0.5/dt)*(cq_weight(n-1,k-2,x,dt)-4*cq_weight(n-1,k-1,x,dt)+3*cq_weight(n-1,k,x,dt))
    else
        return (x/2)^(k/2)*exp(-1.5*x)*hermite(sqrt(2*x),n)/factorial(k)
    end
end

function cq_weightBE(n,k,x,dt)
    if n>0
        return (1/dt)*(cq_weightBE(n-1,k,x,dt)-cq_weightBE(n-1,k-1,x,dt))
    else
        if k<0
            return 0.0
        else
            return x^(k)*exp(-x)/factorial(k)
        end
    end
end

function hermite(x,n)
    coeffs = zeros(n+1,1)
    coeffs[end] = 1
    return Hermite(coeffs)(x)
end