import Base: show, print

Base.show(io::IO,::MIME"text/plain",op::AbstractOperator) = print(io,op)

function Base.print(io::IO,op::LinearCombinationOfOperators)
    l = length(op.coeffs)
    for (ind,coef,operator) in zip(1:l,op.coeffs,op.ops)
        if typeof(coef) <: Complex
            imag(coef) ≈ 0.0 && (coef=real(coef))
        end
        if typeof(coef) <: Real
            if ind == 1 && coef < 0
                print(io," -")
            end
            if ind > 1 && coef >=0
                print(io, " + ")
            end
            if ind > 1 && coef < 0
                print(io, " - ")
            end                

            if !(abs(coef) ≈ 1.0)
                print(io,abs(coef))
            end
            
            print(io,operator)

        else            
            if ind > 1
                print(io," + ")
            end
            print(io,coef)
            print(io,operator)


        end

    end
end

Base.print(io::IO, op::Identity) = print(io,"I")
Base.print(io::IO, op::ZeroOperator) = print(io,"Zero")
function Base.print(io::IO, op::BasisOperatorLeft) 
    print(io,"f*")
    print(io,op.operator)
end
function Base.print(io::IO, op::BasisOperatorRight) 
    print(io,op.operator)
    print(io,"*f")
end

function Base.print(io::IO, op::HHHgreen)
if !(op.α ≈ 1.0)
    @warn "alpha is not 1 and not printed"
end
print(io,"S[$(abs(op.γ))]")
print(io,op.op)
end
function Base.print(io::IO, op::HHHgradgreen)
    if !(op.α ≈ 1.0)
        @warn "alpha is not 1 and not printed"
    end
    print(io,"∇S[$(abs(op.γ))]")
    print(io,op.op)
end
function Base.print(io::IO, op::HHHgradgreenDot)
    if !(op.α ≈ 1.0)
        @warn "alpha is not 1 and not printed"
    end
    print(io,"∇S⋅[$(abs(op.γ))]")
    print(io,op.op)
end
function Base.print(io::IO, op::HHHgradgreenCross)
    if !(op.α ≈ 1.0)
        @warn "alpha is not 1 and not printed"
    end
    print(io,"∇S×[$(abs(op.γ))]×")
    print(io,op.op)
end
function Base.print(io::IO, op::Union{HHHNtestCross,HHHNtestCrossLocal})
    print(io,"nₜₑ×")
    print(io,op.op)
end
function Base.print(io::IO, op::Union{HHHNbasisCross,HHHNbasisCrossLocal})
    print(io,"nₜᵣ×")
    print(io,op.op)
end
function Base.print(io::IO, op::Union{HHHNtestDot,HHHNtestDotLocal})
    print(io,"nₜₑ⋅")
    print(io,op.op)
end
function Base.print(io::IO, op::Union{HHHNbasisnormal,HHHNbasisnormalLocal})
    print(io,"nₜᵣ")
    print(io,op.op)
end
function Base.print(io::IO, op::Union{HHHNbasisdot,HHHNbasisdotLocal})
    print(io,"nₜᵣ⋅")
    print(io,op.op)
end
function Base.print(io::IO, op::Union{HHHIdentity,HHHIdentityLocal})
    print(io,"B")
end
function Base.print(io::IO, op::Union{HHHDivergence,HHHDivergenceLocal})
    print(io,"∇⋅B")
end