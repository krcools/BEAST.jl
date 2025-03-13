abstract type MWFarField <: FarField end

struct MWFarField3D{K, U} <: MWFarField
  gamma::K
  amplitude::U
end
struct MWDoubleLayerFarField3D{K, U} <: MWFarField
  gamma::K
  amplitude::U
end

struct MWDoubleLayerRotatedFarField3D{K,U} <: MWFarField
    gamma::K
    amplitude::U
end

"""
    MWFarField3D(;gamma, amplitude)

Maxwell single layer far-field operator for 3D.
"""
function MWFarField3D(;
  gamma=nothing,
  wavenumber=nothing,
  amplitude=nothing
)
    gamma, _ = gamma_wavenumber_handler(gamma, wavenumber) 
    @assert !isstatic(gamma)

    amplitude === nothing && (amplitude = 1.0)

    return MWFarField3D(gamma, amplitude)
end

MWFarField3D(op::MWSingleLayer3D{T,U}) where {T,U} = MWFarField3D(op.gamma, sqrt(op.α*op.β))

"""
  MWDoubleLayerFarField3D(;gamma, amplitude)

Maxwell double layer far-field operator for 3D.
"""
function MWDoubleLayerFarField3D(;
  gamma=nothing,
  wavenumber=nothing,
  amplitude=nothing
)
    gamma, _ = gamma_wavenumber_handler(gamma, wavenumber)
    @assert !isstatic(gamma)

    amplitude === nothing && (amplitude = 1.0)

    return MWDoubleLayerFarField3D(gamma, amplitude)
end

MWDoubleLayerFarField3D(op::MWDoubleLayer3D{T}) where {T} = MWDoubleLayerFarField3D(op.gamma, 1.0)

# quaddata(op::MWFarField,rs,els) = quadpoints(rs,els,(3,))
# quadrule(op::MWFarField,refspace,p,y,q,el,qdata) = qdata[1,q]

kernelvals(op::MWFarField,y,p) = exp(op.gamma*dot(y,cartesian(p)))
function integrand(op::MWFarField3D,krn,y,f,p)
  op.amplitude*(y × (krn * f[1])) × y
end

function integrand(op::MWDoubleLayerFarField3D,krn,y,f,p)
  op.amplitude*(y × (krn * f[1]))
end

struct MWFarField3DDropConstant{K, U} <: MWFarField
  gamma::K
  coeff::U
end
kernelvals(op::MWFarField3DDropConstant,y,p) = expm1(op.gamma*dot(y,cartesian(p)))
function integrand(op::MWFarField3DDropConstant,krn,y,f,p)
  op.coeff*(y × (krn * f[1])) × y
end

"""
  MWDoubleLayerRotatedFarField3D

Rotated Maxwell double layer far-field operator for 3D.
"""
function MWDoubleLayerRotatedFarField3D(;
    gamma=nothing,
    wavenumber=nothing,
    amplitude=nothing
)
    gamma, _ = gamma_wavenumber_handler(gamma, wavenumber)
    @assert !isstatic(gamma)

    amplitude === nothing && (amplitude = 1.0)

    return MWDoubleLayerRotatedFarField3D(gamma, amplitude)
end

MWDoubleLayerRotatedFarField3D(op::DoubleLayerRotatedMW3D{T,U}) where {T,U} =
MWDoubleLayerRotatedFarField3D(op.gamma, T(1))

function integrand(op::MWDoubleLayerRotatedFarField3D, krn, y, f, p)
    op.amplitude * (y × ((krn * f[1]) × normal(p)))
end

LinearAlgebra.cross(::NormalVector, a::MWDoubleLayerFarField3D) = MWDoubleLayerRotatedFarField3D(a.gamma, a.amplitude)
LinearAlgebra.cross(::NormalVector, a::MWDoubleLayerRotatedFarField3D) = MWDoubleLayerFarField3D(a.gamma, -a.amplitude)
