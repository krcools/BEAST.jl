import Base.*
struct NTimesTrace{T,F} <: Functional
    field::F
end
  
NTimesTrace(f::F) where {F} = NTimesTrace{scalartype(f), F}(f)
NTimesTrace{T}(f::F) where {T,F} = NTimesTrace{T,F}(f)
scalartype(s::NTimesTrace{T}) where {T} = T

(ϕ::NTimesTrace)(p) = *(normal(p), ϕ.field(cartesian(p)))
integrand(::NTimesTrace, g, ϕ) = *(g.value, ϕ)
*(::NormalVector, f) = NTimesTrace(f)

struct Trace{T,F} <: Functional
    field::F
end
  
Trace(f::F) where {F} = Trace{scalartype(f), F}(f)
Trace{T}(f::F) where {T,F} = Trace{T,F}(f)
scalartype(s::Trace{T}) where {T} = T

(ϕ::Trace)(p) = ϕ.field(cartesian(p))
integrand(::NTimesTrace, g, ϕ) = *(g.value, ϕ)
