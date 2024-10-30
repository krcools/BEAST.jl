abstract type MaxwellOperator3D{T,K} <: IntegralOperator end
abstract type MaxwellOperator3DReg{T,K} <: MaxwellOperator3D{T,K} end

scalartype(op::MaxwellOperator3D{T,K}) where {T, K <: Val{0}} = T
scalartype(op::MaxwellOperator3D{T,K}) where {T, K} = promote_type(T, K)

gamma(op::MaxwellOperator3D{T,Val{0}}) where {T} = zero(T)
gamma(op::MaxwellOperator3D{T,K}) where {T, K} = op.gamma


struct MWSingleLayer3D{T,U} <: MaxwellOperator3D{T,U}
  gamma::T
  α::U
  β::U
end

gamma(op::MWSingleLayer3D{Val{0}, U}) where {U} = zero(U)

scalartype(op::MWSingleLayer3D{T,U}) where {T,U} = promote_type(T,U)
# sign_upon_permutation(op::MWSingleLayer3D, I, J) = 1

MWSingleLayer3D(gamma)  = MWSingleLayer3D(gamma, -gamma, -1/(gamma))
MWWeaklySingular(gamma) = MWSingleLayer3D(gamma, 1, 0)
MWHyperSingular(gamma)  = MWSingleLayer3D(gamma, 0, 1)



export Maxwell3D

struct MWSingleLayer3DReg{T,U} <: MaxwellOperator3DReg{T,U}
    gamma::T
    α::U
    β::U
end

struct MWSingleLayer3DSng{T,U} <: MaxwellOperator3D{T,U}
    gamma::T
    α::U
    β::U
end

regularpart(op::MWSingleLayer3D) = MWSingleLayer3DReg(op.gamma, op.α, op.β)
singularpart(op::MWSingleLayer3D) = MWSingleLayer3DSng(op.gamma, op.α, op.β)


function _legendre(n,a,b)
    x, w = FastGaussQuadrature.gausslegendre(n)
    w .*= (b-a)/2
    x .= (x.+1)/2*(b-a).+a
    collect(zip(x,w))
end

# defaultquadstrat(op::MaxwellOperator3D, tfs::Space, bfs::Space) = DoubleNumWiltonSauterQStrat(2,3,6,7,5,5,4,3)

defaultquadstrat(op::MaxwellOperator3D, tfs::RTRefSpace, bfs::RTRefSpace) = DoubleNumWiltonSauterQStrat(2,3,6,7,5,5,4,3)
# defaultquadstrat(op::MaxwellOperator3D, tfs::RefSpace, bfs::RefSpace) = DoubleNumWiltonSauterQStrat(2,3,6,7,5,5,4,3)




struct MWDoubleLayer3D{T,K} <: MaxwellOperator3D{T,K}
    alpha::T
    gamma::K
end

# sign_upon_permutation(op::MWDoubleLayer3D, I, J) = 1

struct MWDoubleLayer3DSng{T,K} <: MaxwellOperator3D{T,K}
    alpha::T
    gamma::K
end

struct MWDoubleLayer3DReg{T,K} <: MaxwellOperator3DReg{T,K}
    alpha::T
    gamma::K
end

MWDoubleLayer3D(gamma) = MWDoubleLayer3D(1.0, gamma) # For legacy purposes

regularpart(op::MWDoubleLayer3D) = MWDoubleLayer3DReg(op.alpha, op.gamma)
singularpart(op::MWDoubleLayer3D) = MWDoubleLayer3DSng(op.alpha, op.gamma)



# function quadrule(op::MaxwellOperator3D, g::BDMRefSpace, f::BDMRefSpace,  i, τ, j, σ, qd,
#   qs::DoubleNumWiltonSauterQStrat)

#   hits = 0
#   dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))
#   dmin2 = floatmax(eltype(eltype(τ.vertices)))
#   for t in τ.vertices
#       for s in σ.vertices
#           d2 = LinearAlgebra.norm_sqr(t-s)
#           dmin2 = min(dmin2, d2)
#           hits += (d2 < dtol)
#       end
#   end

#   hits == 3 && return SauterSchwabQuadrature.CommonFace(qd.gausslegendre[3])
#   hits == 2 && return SauterSchwabQuadrature.CommonEdge(qd.gausslegendre[2])
#   hits == 1 && return SauterSchwabQuadrature.CommonVertex(qd.gausslegendre[1])

#   h2 = volume(σ)
#   xtol2 = 0.2 * 0.2
#   k2 = abs2(gamma(op))
#   return DoubleQuadRule(
#       qd.tpoints[1,i],
#       qd.bpoints[1,j],)
# end


# function qrdf(op::MaxwellOperator3D, g::RTRefSpace, f::RTRefSpace, i, τ, j, σ, qd)
#   # defines coincidence of points
#   dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))

#   # decides on whether to use singularity extraction
#   xtol = 0.2

#   k = norm(gamma(op))

#   hits = 0
#   xmin = xtol
#   for t in τ.vertices
#     for s in σ.vertices
#       d = norm(t-s)
#       xmin = min(xmin, k*d)
#       if d < dtol
#         hits +=1
#         break
#       end
#     end
#   end

#   xmin < xtol && return WiltonSERule(
#     qd.tpoints[1,i],
#     DoubleQuadRule(
#       qd.tpoints[2,i],
#       qd.bpoints[2,j],
#     ),
#   )
#   return DoubleQuadRule(
#     qd.tpoints[1,i],
#     qd.bpoints[1,j],
#   )

# end

################################################################################
#
#  Kernel definitions
#
################################################################################

const i4pi = 1 / (4pi)
function (igd::Integrand{<:MWSingleLayer3D})(x,y,f,g)
    α = igd.operator.α
    β = igd.operator.β
    γ = igd.operator.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1 / R
    green = exp(-γ*R)*(i4pi*iR)

    αG = α * green
    βG = β * green

    _integrands(f,g) do fi,gj
        αG * dot(fi.value, gj.value) + βG * dot(fi.divergence, gj.divergence)
    end
end



function (igd::Integrand{<:MWSingleLayer3DReg})(x,y,f,g)
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
end


function (igd::Integrand{<:MWDoubleLayer3D})(x,y,f,g)

    γ = igd.operator.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-γ*R)*(iR*i4pi)
    gradgreen = -(γ + iR) * green * (iR * r)

    fvalue = getvalue(f)
    gvalue = getvalue(g)
    G = cross.(Ref(gradgreen), gvalue)
    return _krondot(fvalue, G)
end


function (igd::Integrand{<:MWDoubleLayer3DReg})(x,y,f,g)

    γ = igd.operator.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    γR = γ*R
    iR = 1/R
    expo = exp(-γR)
    green = (expo - 1 + γR - 0.5*γR^2) * (i4pi*iR)
    gradgreen = ( -(γR + 1)*expo + (1 - 0.5*γR^2) ) * (i4pi*iR^3) * r

    fvalue = getvalue(f)
    gvalue = getvalue(g)
    G = cross.(Ref(gradgreen), gvalue)
    return _krondot(fvalue, G)
end

################################################################################
#
#  Handling of operator parameters (Helmholtz and Maxwell)
#
################################################################################

"""
    gamma_wavenumber_handler(gamma, wavenumber)

This function handles the input of `gamma` and `wavenumber`. It throws an error if both `gamma` and
`wavenumber` are provided. If neither is provided, it assumes a static problem and returns `Val(0)`
for `gamma` and `wavenumber`.

# Arguments
- `gamma`: `im` * `wavenumber` or `nothing`.
- `wavenumber`: `wavenumber` or `nothing`.

# Returns
- `gamma` and `wavenumber`: Appropriate pair `gamma` and `wavenumber`.
"""
function gamma_wavenumber_handler(gamma, wavenumber)
    if !isnothing(gamma) && !isnothing(wavenumber)
        error("Supplying both gamma and wavenumber is not supported.")

    elseif isnothing(gamma) && isnothing(wavenumber)
        # if neither gamma nor wavenumber is supplied, we are assuming a static problem
        return Val(0), Val(0)
    end

    if isnothing(gamma) && !isnothing(wavenumber)
        if iszero(real(wavenumber))
            gamma = -imag(wavenumber)
        else
            gamma = im*wavenumber
        end
    end

    return gamma, wavenumber
end

"""
    isstatic(gamma)

This function checks if the provided `gamma` value represents a static problem.
It returns true if `gamma` is of type `Val{0}` indicating a static problem.

# Arguments
- `gamma`: `gamma` value.

# Returns
- A boolean indicating whether the problem is static or not.
"""
function isstatic(gamma)
    return typeof(gamma) == Val{0}
end

function isstatic(op::MaxwellOperator3D)
    return isstatic(op.gamma)
end

function operator_parameter_handler(alpha, gamma, wavenumber)

gamma, wavenumber = gamma_wavenumber_handler(gamma, wavenumber)

    if alpha === nothing
        if isstatic(gamma) # static problem
            alpha = 1.0 # default to double precision
        else
            alpha = one(real(typeof(gamma)))
        end
    end

return alpha, gamma
end
