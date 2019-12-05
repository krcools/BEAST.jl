struct SourceField{F} <: Functional
    f::F
end

(s::SourceField)(p) = s.f(cartesian(p))
integrand(f::SourceField, tval, fval) = dot(fval, tval.value)
