import Base: *
import LinearAlgebra: ×, ⋅
import BEAST
abstract type ComposedOperatorLocal <: LocalOperator end
abstract type ComposedOperatorIntegral <: IntegralOperator end



struct BasisFunction end
struct DivBasisFunction end

struct TestFunctionLocal{T <: ComposedOperatorLocal} <: ComposedOperatorLocal 
    inner::T
end
struct TestFunctionIntegral{T <: ComposedOperatorIntegral} <: ComposedOperatorIntegral 
    inner::T
end
TestFunctionLocal(inner::ComposedOperatorIntegral) = TestFunctionIntegral(inner)

struct DivTestFunctionLocal{T <: ComposedOperatorLocal} <: ComposedOperatorLocal 
    inner::T
end
struct DivTestFunctionIntegral{T <: ComposedOperatorIntegral} <: ComposedOperatorIntegral
    inner::T
end

function TestFunctionLocal(inner::LinearCombinationOfOperators)
    out = ZeroOperator()
    for op in inner.ops
        out += TestFunctionLocal(op)
    end
    return out
end

function DivTestFunctionLocal(inner::LinearCombinationOfOperators)
    out = ZeroOperator()
    for op in inner.ops
        out += DivTestFunctionLocal(op)
    end
    return out
end


struct CurlBasisFunction end

const b = BasisFunction()
const t = TestFunctionLocal
export t
export b

×(::Nabla,::BasisFunction) = CurlBasisFunction()
⋅(::Nabla,::BasisFunction) = DivBasisFunction()
⋅(::Nabla,::Type{TestFunctionLocal}) = DivTestFunctionLocal
⋅(::Nabla,::Type{TestFunctionIntegral}) = DivTestFunctionIntegral

struct TestNormalLocal{T,O,P} <: ComposedOperatorLocal 
inner::T
operation::O
factor::P
end


struct TrialNormalLocal{T,O,P} <: ComposedOperatorLocal 
inner::T
operation::O
factor::P
end


*(a::Number,b::TrialNormalLocal) = TrialNormalLocal(b.inner,b.operation,a*b.factor)
*(a::Number,b::TestNormalLocal) = TestNormalLocal(b.inner,b.operation,a*b.factor)
struct TestNormalIntegral{T,O,P} <: ComposedOperatorIntegral
    inner::T
    operation::O
    factor::P
    
end

TestNormalLocal(inner::ComposedOperatorIntegral,op) = TestNormalIntegral(inner,op)

struct TrialNormalIntegral{T,O,P} <: ComposedOperatorIntegral
    inner::T
    operation::O
    factor::P
end
TrialNormalIntegral(a,b) = TrialNormalIntegral(a,b,1.0)
TestNormalIntegral(a,b) = TestNormalIntegral(a,b,1.0)
TestNormalLocal(a,b) = TestNormalLocal(a,b,1.0)
TestNormalIntegral(a,b) = TestNormalIntegral(a,b,1.0)
TrialNormalLocal(inner::ComposedOperatorIntegral,op) = TrialNormalIntegral(inner,op)

*(a::Number,b::TrialNormalIntegral) = TrialNormalIntegral(b.inner,b.operation,a*b.factor)
*(a::Number,b::TestNormalIntegral) = TestNormalIntegral(b.inner,b.operation,a*b.factor)
const nt = TestNormalLocal
const nb = TrialNormalLocal

abstract type Kern end


struct Kernel{T,O,P <: Kern} <: ComposedOperatorIntegral
    inner::T
    operation::O
    Kern::P
end

function (k::Kern)(operation)
return x -> Kernel(x,operation,k)
end


function (a::Type{<: ComposedOperatorLocal})(operation)
return x->dist(a,x,operation)
end

function dist(a::Type{<: ComposedOperatorLocal},x::Union{ComposedOperatorIntegral,ComposedOperatorLocal},operation)
    return a(x,operation)
end
function dist(a::Type{<: ComposedOperatorLocal},x::LinearCombinationOfOperators,operation)
    out = ZeroOperator()
    for op in x.ops
        out += a(op,operation)
    end
    return out
end

struct HH3DGreen{T} <: Kern
    gamma::T
end
struct HH3DGradGreen{T} <: Kern
    gamma::T
end
function (::Nabla)(G::HH3DGreen)
    HH3DGradGreen(G.wavenumber)
end
function (G::Kernel{T,O,HH3DGreen{Q}})(testnb,trialnb,f,g) where {T,O,Q}
    green = 
    return G.operation.(Ref(green),G.inner(testnb,trialnb,f,g))
end
function (G::Kernel{T,O,HH3DGradGreen{Q}})(testnb,trialnb,f,g) where {T,O,Q}
    gradgreen = 
    return G.operation.(Ref(gradgreen),G.inner(testnb,trialnb))
end

#define with respect to from inside omega.
function traceterm(G::HH3DGreen)
ZeroOperator()
end

function traceterm(G::HH3DGradGreen)
1/2*nt(G.inner,G.operation)
end




function (igd::Integrand{<:ComposedOperatorIntegral})(x,y,f,g)
op = igd.operator
op(x,y,f,g)
end

function (op::Union{TestFunctionIntegral,TestFunctionLocal})(x,y,f,g)
    _krondot(getvalue(f),op(x,y,g))
end
function (op::Union{DivTestFunctionIntegral,DivTestFunctionLocal})(x,y,f,g)
    _krondot(getdivergence(f),op(x,y,g))
end

function integrand(op::ComposedOperatorLocal,kernel,x,y,f,g)
    op(x,y,f,g)
end

function (op::Union{TestNormalIntegral,TestNormalLocal})(x,y,g)
    op.operation(normal(x),op.inner(x,y,g))
end

function (op::Union{TrialNormalIntegral,TrialNormalLocal})(x,y,g)
    op.factor*op.operation(normal(y),op.inner(x,y,g))
end

function (op::BasisFunction)(x,y,g)
    getvalue(g)
end
function (op::DivBasisFunction)(x,y,g)
    getdivergence(g)
end
function (op::CurlBasisFunction)(x,y,g)
    getcurl(g)
end


#### support for traces to find equivalent operator.
"""
inside positive, outside negative
"""
struct Trace
direction
end

function (t::Trace)(operator)
operator+t.direction*traceterm(operator)
end





const γₜᶜ = nt((x -> -x ∘ ×))∘ nt(×)∘Trace(-1)
const γₜ = nt((x -> -x ∘ ×))∘ nt(×)∘Trace(1)

const γₛ = nt(×)∘Trace(1)
const γₛᶜ = nt(×)∘Trace(-1)
const γₙ = nt(⋅)∘Trace(1)
const γₙᶜ = nt(⋅)∘Trace(-1)

export γₙ
export γₙᶜ
export γₜ
export γₜᶜ