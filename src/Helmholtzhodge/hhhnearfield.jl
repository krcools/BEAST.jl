mutable struct HHHGreenField{T,U}
    gamma::T
    op::U
end
mutable struct HHHGradGreenField{T,U}
    gamma::T
    op::U
end
mutable struct HHHGradGreenCrossField{T,U}
    gamma::T
    op::U
end
mutable struct HHHGradGreenDotField{T,U}
    gamma::T
    op::U
end

mutable struct HHHBasisNtimesField{U}
    op::U
end
mutable struct HHHIdentityField end
mutable struct HHHDivergenceField end

function HHHGreenField(;
    gamma=nothing,
    wavenumber=nothing
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

    HHHGreenField(gamma,HHHIdentityField())
end
function HHHGradGreenField(;
    gamma=nothing,
    wavenumber=nothing
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

    HHHGradGreenField(gamma,HHHIdentityField())
end
function HHHGradGreenCrossField(;
    gamma=nothing,
    wavenumber=nothing
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

    HHHGradGreenCrossField(gamma,HHHIdentityField())
end
function HHHGradGreenDotField(;
    gamma=nothing,
    wavenumber=nothing
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

    HHHGradGreenDotField(gamma,HHHIdentityField())
end

const HHHField = Union{HHHGreenField,HHHGradGreenField,HHHGradGreenCrossField}
defaultquadstrat(op::HHHField, basis) = SingleNumQStrat(2)
quaddata(op::HHHField,rs,els,qs::SingleNumQStrat) = quadpoints(rs,els,(qs.quad_rule,))
quadrule(op::HHHField,refspace,p,y,q,el,qdata,qs::SingleNumQStrat) = qdata[1,q]

function kernelvals(op::HHHField,y,p)

        γ = op.gamma
        r = y - cartesian(p)
        R = norm(r)

        γR = γ*R
        expn = exp(-γR)
        green = expn / (4pi*R)
        gradgreen = -(γ + 1/R) * green / R * r

        krn = (trans=r, dist=R, green=green, gradgreen=gradgreen)
end

function integrand(op::HHHGreenField, krn, y, fp, p)
G = krn.green
G*integrand(op.op, krn, y, fp, p)
end

function integrand(op::HHHBasisNtimesField, krn, y, fp, p)
    n = normal(p)
    n*integrand(op.op, krn, y, fp, p)
end
integrand(op::HHHIdentityField, krn,y,fp,p) = fp.value
integrand(op::HHHDivergenceField,krn,y,fp,p) = fp.divergence

function integrand(op::HHHGradGreenField, krn, y, fp, p)
    integrand(op.op, krn,y,fp,p)*krn.gradgreen
end

function integrand(op::HHHGradGreenCrossField, krn, y, fp, p)
    cross(krn.gradgreen,integrand(op.op, krn,y,fp,p))
end
function integrand(op::HHHGradGreenDotField, krn, y, fp, p)
    dot(integrand(op.op, krn,y,fp,p),krn.gradgreen)
end