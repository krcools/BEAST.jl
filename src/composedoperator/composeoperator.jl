import Base: *
import BEAST
abstract type ComposedOperatorLocal <: BEAST.LocalOperator end
abstract type ComposedOperatorIntegral <: BEAST.IntegralOperator end


struct TestNormalLocal{T,O} <: ComposedOperatorLocal 
inner::T
operation::O
end


struct TrialNormalLocal{T,O} <: ComposedOperatorLocal 
inner::T
operation::O
end

struct TestNormalIntegral{T,O} <: ComposedOperatorIntegral
    inner::T
    operation::O
end
TestNormalLocal(inner::ComposedOperatorIntegral,op) = TestNormalIntegral(inner,op)

struct TrialNormalIntegral{T,O} <: ComposedOperatorIntegral
    inner::T
    operation::O
end
TrialNormalLocal(inner::ComposedOperatorIntegral,op) = TrialNormalIntegral(inner,op)

const nt = TestNormalLocal
const nb = TrialNormalLocal

abstract type Kern end


struct Kernel{T,O,P <: Kern} <: ComposedOperatorInegral
    inner::T
    operation::O
    Kern::P
end

function (k::Kern)(operation)
return x -> Kernel(x,operation,k)
end


function (a::Type{<: ComposedOperatorLocal})(operation)
return x->a(x,operation)
end



struct HH3DGreen{T} <: Kern
    gamma::T
end
struct HH3DGradGreen{T} <: Kern
    gamma::T
end
∇(G::HH3DGreen) = HH3DGradGreen(G.wavenumber)
function (G::Kernel{T,O,HH3DGreen{Q}})(testnb,trialnb) where {T,O,Q}
    green = 
    return G.operation.(G.inner(testnb,trialnb),Ref(green))
end
function (G::Kernel{T,O,HH3DGradGreen{Q}})(testnb,trialnb) where {T,O,Q}
    gradgreen = 
    return G.operation.(G.inner(testnb,trialnb),Ref(gradgreen))
end