# struct SourceField{T,F} <: Functional{T}
#     field::F
# end

# (s::SourceField)(p) = s.f(cartesian(p))
# integrand(f::SourceField, tval, fval) = dot(fval, tval.value)
