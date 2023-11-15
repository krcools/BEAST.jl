abstract type MWFarField <: FarField end

"""
Operator to compute the far field of a current distribution. In particular, given the current distribution ``j`` this operator allows for the computation of

```math
A j = n × ∫_Γ e^{γ ̂x ⋅ y} dy
```

where ``̂x`` is the unit vector in the direction of observation. Note that the assembly routing expects the observation directions to be normalised by the caller.
"""
struct MWFarField3D{K, U} <: MWFarField
  gamma::K
  waveimpedance::U
end
struct MWDoubleLayerFarField3D{K, U} <: MWFarField
  gamma::K
  waveimpedance::U
end

function MWFarField3D(;
  gamma=nothing,
  wavenumber=nothing,
  waveimpedance=nothing
)
    if (gamma === nothing) && (wavenumber === nothing)
        error("Supply one of (not both) gamma or wavenumber")
    end

    if (gamma !== nothing) && (wavenumber !== nothing)
        error("Supply one of (not both) gamma or wavenumber")
    end

    if gamma === nothing
        if iszero(real(wavenumber))
            gamma = -imag(wavenumber)
        else
            gamma = im*wavenumber
        end
    end

    @assert gamma !== nothing

    waveimpedance === nothing && (waveimpedance = 1.0)
    
    return MWFarField3D(gamma, waveimpedance)
end

MWFarField3D(op::MWSingleLayer3D{T,U}) where {T,U} = MWFarField3D(op.gamma, sqrt(op.α*op.β))

function MWDoubleLayerFarField3D(;
  gamma=nothing,
  wavenumber=nothing,
  waveimpedance=nothing
)
    if (gamma === nothing) && (wavenumber === nothing)
        error("Supply one of (not both) gamma or wavenumber")
    end

    if (gamma !== nothing) && (wavenumber !== nothing)
        error("Supply one of (not both) gamma or wavenumber")
    end

    if gamma === nothing
        if iszero(real(wavenumber))
            gamma = -imag(wavenumber)
        else
            gamma = im*wavenumber
        end
    end

    @assert gamma !== nothing

    waveimpedance === nothing && (waveimpedance = 1.0)
    
    return MWDoubleLayerFarField3D(gamma, waveimpedance)
end

MWDoubleLayerFarField3D(op::MWDoubleLayer3D{T}) where {T} = MWDoubleLayerFarField3D(op.gamma, 1.0)

# quaddata(op::MWFarField,rs,els) = quadpoints(rs,els,(3,))
# quadrule(op::MWFarField,refspace,p,y,q,el,qdata) = qdata[1,q]

kernelvals(op::MWFarField,y,p) = exp(op.gamma*dot(y,cartesian(p)))
function integrand(op::MWFarField3D,krn,y,f,p)
  op.waveimpedance*(y × (krn * f[1])) × y
end

function integrand(op::MWDoubleLayerFarField3D,krn,y,f,p)
  op.waveimpedance*(y × (krn * f[1]))
end

struct MWFarField3DDropConstant{K, U} <: MWFarField
  gamma::K
  coeff::U
end
kernelvals(op::MWFarField3DDropConstant,y,p) = expm1(op.gamma*dot(y,cartesian(p)))
function integrand(op::MWFarField3DDropConstant,krn,y,f,p)
  op.coeff*(y × (krn * f[1])) × y
end