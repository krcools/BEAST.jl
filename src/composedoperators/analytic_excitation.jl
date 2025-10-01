struct Trace{T,F} <: Functional{T}
    field::F
end
  
Trace(f::F) where {F} = Trace{scalartype(f), F}(f)
Trace{T}(f::F) where {T,F} = Trace{T,F}(f)
scalartype(s::Trace{T}) where {T} = T

(ϕ::Trace)(p) = ϕ.field(cartesian(p))
integrand(::Trace, g, ϕ) = *(g.value, ϕ)
