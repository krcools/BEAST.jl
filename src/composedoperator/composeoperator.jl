import Base: *
import LinearAlgebra: ×, ⋅
import BEAST
abstract type ComposedOperatorLocal <: LocalOperator end
abstract type ComposedOperatorIntegral <: IntegralOperator end



struct BasisFunction end
struct DivBasisFunction end

struct TestFunctionLocal{T} <: ComposedOperatorLocal
    inner::T
end
struct TestFunctionIntegral{T} <: ComposedOperatorIntegral
    inner::T
end
TestFunctionLocal(inner::ComposedOperatorIntegral) = TestFunctionIntegral(inner)

struct DivTestFunctionLocal{T} <: ComposedOperatorLocal
    inner::T
end
struct DivTestFunctionIntegral{T} <: ComposedOperatorIntegral
    inner::T
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


struct Kernel{T,O,P <: Kern} <: ComposedOperatorIntegral
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
function (::Nabla)(G::HH3DGreen)
    HH3DGradGreen(G.wavenumber)
end
function (G::Kernel{T,O,HH3DGreen{Q}})(testnb,trialnb) where {T,O,Q}
    green = 
    return G.operation.(Ref(green),G.inner(testnb,trialnb))
end
function (G::Kernel{T,O,HH3DGradGreen{Q}})(testnb,trialnb) where {T,O,Q}
    gradgreen = 
    return G.operation.(Ref(gradgreen),G.inner(testnb,trialnb))
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
    op.operation(normal(y),op.inner(x,y,g))
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

