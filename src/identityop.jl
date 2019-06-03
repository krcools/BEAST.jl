struct Identity <: LocalOperator
end

kernelvals(biop::Identity, x) = Nothing
integrand(op::Identity, kernel, x, g, f) = dot(f[1], g[1])
scalartype(op::Identity) = Union{}

struct NCross <: LocalOperator
end

kernelvals(op::NCross, mp) = Nothing()
integrand(op::NCross, kernel, x, g, f) = dot(g[1], normal(x) × f[1])
scalartype(op::NCross) = Union{}

function quaddata(op::LocalOperator, g::RTRefSpace, f::RTRefSpace, tels, bels)
    u, w = trgauss(6)
    return [(w[i],SVector(u[1,i],u[2,i])) for i in 1:length(w)]
end

function quaddata(op::LocalOperator, g::NDRefSpace, f::NDRefSpace, tels, bels)
    u, w = trgauss(6)
    return [(w[i],SVector(u[1,i],u[2,i])) for i in 1:length(w)]
end

function quaddata(op::LocalOperator, g::subReferenceSpace, f::subReferenceSpace, tels, bels)
    u, w = trgauss(6)
    return [(w[i],SVector(u[1,i],u[2,i])) for i in 1:length(w)]
end



quaddata(op::LocalOperator, g::LagrangeRefSpace, f::LagrangeRefSpace,
        tels::Vector, bels::Vector) = quaddata(op, g, f, tels, bels, Val{dimension(tels[1])})

function quaddata(op::LocalOperator, g::LagrangeRefSpace, f::LagrangeRefSpace,
        tels::Vector, bels::Vector, dim::Type{Val{1}})

    u, w = legendre(6, 0.0, 1.0)
    [(w[i],u[i]) for i in eachindex(w)]
end


function quaddata(op::LocalOperator, g::LagrangeRefSpace, f::LagrangeRefSpace,
        tels::Vector, bels::Vector, dim::Type{Val{2}})

    u, w = trgauss(6)
    [(w[i], SVector(u[1,i], u[2,i])) for i in 1:length(w)]
end


function quadrule(op::LocalOperator, ψ::RefSpace, ϕ::RefSpace, τ, qd)
    q = qd[1]
    w, p = q[1], neighborhood(τ,q[2])
    a = (w, p, ψ(p), ϕ(p))
    A = Vector{typeof(a)}(undef,length(qd))
    A[1] = a
    for i in 2:length(qd)
        q = qd[i]
        w, p = q[1], neighborhood(τ,q[2])
        A[i] = (w, p, ψ(p), ϕ(p))
    end
    return A
end
