struct Identity <: LocalOperator
end

kernelvals(biop::Identity, x) = nothing
integrand(op::Identity, kernel, x, g, f) = dot(f[1], g[1])
scalartype(op::Identity) = Union{}

struct NCross <: LocalOperator
end

kernelvals(op::NCross, mp) = nothing
integrand(op::NCross, kernel, x, g, f) = dot(g[1], normal(x) × f[1])
scalartype(op::NCross) = Union{}

function _alloc_workspace(qd, g, f, tels, bels)
    q = qd[1]
    τ = tels[1]
    w, p = q[1], neighborhood(τ,q[2])
    a = (w, p, g(p), f(p))
    A = Vector{typeof(a)}(undef,length(qd))
end

const LinearRefSpaceTriangle = Union{RTRefSpace, NDRefSpace}
function quaddata(op::LocalOperator, g::LinearRefSpaceTriangle, f::LinearRefSpaceTriangle, tels, bels)
    u, w = trgauss(6)
    qd = [(w[i],SVector(u[1,i],u[2,i])) for i in 1:length(w)]
    A = _alloc_workspace(qd, g, f, tels, bels)

    return qd, A
end

function quaddata(op::LocalOperator, g::subReferenceSpace, f::subReferenceSpace, tels, bels)
    u, w = trgauss(6)
    qd = [(w[i],SVector(u[1,i],u[2,i])) for i in 1:length(w)]
    A = _alloc_workspace(qd, g, f, tels, bels)
    return qd, A
end

const LinearRefSpaceTetr = Union{NDLCCRefSpace, NDLCDRefSpace}
function quaddata(op::LocalOperator, g::LinearRefSpaceTetr, f::LinearRefSpaceTetr, tels, bels)
    o, x, y, z = CompScienceMeshes.euclidianbasis(3)
    reftet = simplex(x,y,z,o)
    qps = quadpoints(reftet, 3)
    qd = [(w, parametric(p)) for (p,w) in qps]
    A = _alloc_workspace(qd, g, f, tels, bels)
    return qd, A
end

function quaddata(op::LocalOperator, g::LagrangeRefSpace{T,Deg,2} where {T,Deg},
    f::LagrangeRefSpace, tels::Vector, bels::Vector)
    
    u, w = legendre(6, 0.0, 1.0)
    qd = [(w[i],u[i]) for i in eachindex(w)]
    A = _alloc_workspace(qd, g, f, tels, bels)
    return qd, A
end

function quaddata(op::LocalOperator, g::LagrangeRefSpace{T,Deg,3} where {T,Deg},
    f::LagrangeRefSpace, tels::Vector, bels::Vector)

    u, w = trgauss(6)
    qd = [(w[i], SVector(u[1,i], u[2,i])) for i in 1:length(w)]
    A = _alloc_workspace(qd, g, f, tels, bels)
    return qd, A
end

function quaddata(op::LocalOperator, g::LagrangeRefSpace{T,Deg,4} where {T,Deg},
    f::LagrangeRefSpace, tels::Vector, bels::Vector)

     o, x, y, z = CompScienceMeshes.euclidianbasis(3)
     reftet = simplex(x,y,z,o)
     qps = quadpoints(reftet, 6)
     qd = [(w, parametric(p)) for (p,w) in qps]
     A = _alloc_workspace(qd, g, f, tels, bels)
     return qd, A
end


function quadrule(op::LocalOperator, ψ::RefSpace, ϕ::RefSpace, τ, (qd,A))
    for i in eachindex(qd)
        q = qd[i]
        w, p = q[1], neighborhood(τ,q[2])
        A[i] = (w, p, ψ(p), ϕ(p))
    end
    return A
end
