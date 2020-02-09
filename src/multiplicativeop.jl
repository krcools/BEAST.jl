struct Multiplicative{f} <: LocalOperator
    field::f
end

kernelvals(biop::Multiplicative, x) = Nothing
integrand(op::Multiplicative, kernel, x, g, f) = dot(f[1], op.field(cartesian(x))*g[1])
scalartype(op::Multiplicative) = Float64
