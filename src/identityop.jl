
"""
    Identity <: LocalOperator

The identity operator.
"""
struct Identity <: LocalOperator
end

kernelvals(biop::Identity, x) = nothing
integrand(op::Identity, kernel, x, g, f) = dot(f[1], g[1])
scalartype(op::Identity) = Union{}

"""
    NCross <: LocalOperator

The identity operator where the trial function is rotated.
"""
struct NCross <: LocalOperator
end

kernelvals(op::NCross, mp) = nothing
integrand(op::NCross, kernel, x, g, f) = dot(g[1], normal(x) × f[1])
scalartype(op::NCross) = Union{}

#= Should be no longer be needed
struct Curl <: LocalOperator
end

kernelvals(biop::Curl, x) = nothing
function integrand(op::Curl, kernel, x, g, f)
     dot(f.curl, g.value)
 end
scalartype(op::Curl) = Union{}
=#


defaultquadstrat(::LocalOperator, ::GWPDivRefSpace{<:Real,D1},::LagrangeRefSpace{T,D2,D3}) where {T,D1,D2,D3} = SingleNumQStrat(9)
function quaddata(op::LocalOperator, g::GWPDivRefSpace, f::LagrangeRefSpace{T,D2,D1} where {T,D1,D2} , tels, bels, qs::SingleNumQStrat)

    u, w = trgauss(qs.quad_rule)
    qd = [(w[i],SVector(u[1,i],u[2,i])) for i in 1:length(w)]
    A = _alloc_workspace(qd, g, f, tels, bels)

    return qd, A
end

defaultquadstrat(::LocalOperator, ::RTRefSpace,  ::LagrangeRefSpace{T,D2,D1}) where {T,D1,D2} = SingleNumQStrat(6)
# function quaddata(op::LocalOperator, g::RTRefSpace, f::LagrangeRefSpace{T,D2,D1} where {T,D1,D2} , tels, bels, qs::SingleNumQStrat)

#     u, w = trgauss(qs.quad_rule)
#     qd = [(w[i],SVector(u[1,i],u[2,i])) for i in 1:length(w)]
#     A = _alloc_workspace(qd, g, f, tels, bels)

#     return qd, A
# end

function _alloc_workspace(qd, g, f, tels, bels)
    q = qd[1]
    τ = tels[1]
    w, p = q[1], neighborhood(τ,q[2])
    a = (w, p, g(p), f(p))
    A = Vector{typeof(a)}(undef,length(qd))
end

const LinearRefSpaceTriangle = Union{RTRefSpace, RT2RefSpace, NDRefSpace, ND2RefSpace, BDMRefSpace, NCrossBDMRefSpace}
defaultquadstrat(::LocalOperator, ::LinearRefSpaceTriangle, ::LinearRefSpaceTriangle) = SingleNumQStrat(6)
# function quaddata(op::LocalOperator, g::LinearRefSpaceTriangle, f::LinearRefSpaceTriangle, tels, bels,
#         qs::SingleNumQStrat)

#     u, w = trgauss(qs.quad_rule)
#     qd = [(w[i],SVector(u[1,i],u[2,i])) for i in 1:length(w)]
#     A = _alloc_workspace(qd, g, f, tels, bels)

#     return qd, A
# end

defaultquadstrat(::LocalOperator, ::subReferenceSpace, ::subReferenceSpace) = SingleNumQStrat(6)
function quaddata(op::LocalOperator, g::subReferenceSpace, f::subReferenceSpace, tels, bels,
        qs::SingleNumQStrat)

    u, w = trgauss(qs.quad_rule)
    qd = [(w[i],SVector(u[1,i],u[2,i])) for i in 1:length(w)]
    A = _alloc_workspace(qd, g, f, tels, bels)
    return qd, A
end

const LinearRefSpaceTetr = Union{NDLCCRefSpace, NDLCDRefSpace, BDM3DRefSpace}
defaultquadstrat(::LocalOperator, ::LinearRefSpaceTetr, ::LinearRefSpaceTetr) = SingleNumQStrat(3)
function quaddata(op::LocalOperator, g::LinearRefSpaceTetr, f::LinearRefSpaceTetr, tels, bels, qs::SingleNumQStrat)

    o, x, y, z = CompScienceMeshes.euclidianbasis(3)
    reftet = simplex(x,y,z,o)
    qps = quadpoints(reftet, qs.quad_rule)
    qd = [(w, parametric(p)) for (p,w) in qps]
    A = _alloc_workspace(qd, g, f, tels, bels)
    return qd, A
end

defaultquadstrat(::LocalOperator, ::LagrangeRefSpace{T,D1,2}, ::LagrangeRefSpace{T,D2,2}) where {T,D1,D2} = SingleNumQStrat(6)
function quaddata(op::LocalOperator, g::LagrangeRefSpace{T,Deg,2} where {T,Deg},
    f::LagrangeRefSpace, tels::Vector, bels::Vector, qs::SingleNumQStrat)
    U = typeof(tels[1].volume)
    u, w = legendre(qs.quad_rule, 0.0, 1.0)
    qd = [(w[i],u[i]) for i in eachindex(w)]
    A = _alloc_workspace(Tuple{U,U}.(qd), g, f, tels, bels)
    return qd, A
end

defaultquadstrat(::LocalOperator, ::LagrangeRefSpace{T,D1,3}, ::LagrangeRefSpace{T,D2,3}) where {T,D1,D2} = SingleNumQStrat(6)
# function quaddata(op::LocalOperator, g::LagrangeRefSpace{T,Deg,3} where {T,Deg},
#     f::LagrangeRefSpace, tels::Vector, bels::Vector, qs::SingleNumQStrat)

#     u, w = trgauss(qs.quad_rule)
#     qd = [(w[i], SVector(u[1,i], u[2,i])) for i in 1:length(w)]
#     A = _alloc_workspace(qd, g, f, tels, bels)
#     return qd, A
# end

defaultquadstrat(::LocalOperator, ::LagrangeRefSpace{T,D1,4}, ::LagrangeRefSpace{T,D2,4}) where {T,D1,D2} = SingleNumQStrat(6)
function quaddata(op::LocalOperator, g::LagrangeRefSpace{T,Deg,4} where {T,Deg},
    f::LagrangeRefSpace, tels::Vector, bels::Vector, qs::SingleNumQStrat)

     o, x, y, z = CompScienceMeshes.euclidianbasis(3)
     reftet = simplex(x,y,z,o)
     qps = quadpoints(reftet, qs.quad_rule)
     qd = [(w, parametric(p)) for (p,w) in qps]
     A = _alloc_workspace(qd, g, f, tels, bels)
     return qd, A
end

# support(::LinearRefSpaceTriangle) = TriangleSupport()
# support(::LinearRefSpaceTetr) = TetraSupport()

#defaultquadstrat(::LocalOperator, ::GWPDivRefSpace{<:Real,D1},
#    ::GWPDivRefSpace{<:Real,D2}) where {D1,D2} = SingleNumQStrat(9)


# function quaddata(op::LocalOperator, g::GWPDivRefSpace, f::GWPDivRefSpace,
#     tels::Vector, bels::Vector,
#     qs::SingleNumQStrat)

#     u, w = trgauss(qs.quad_rule)
#     qd = [(w[i],SVector(u[1,i],u[2,i])) for i in 1:length(w)]
#     A = _alloc_workspace(qd, g, f, tels, bels)
#     return qd, A
# end


#defaultquadstrat(::LocalOperator, ::GWPDivRefSpace{<:Real,D1},
#    ::GWPDivRefSpace{<:Real,D2}) where {D1,D2} = SingleNumQStrat(9)


#function quaddata(op::LocalOperator, g::GWPDivRefSpace, f::GWPDivRefSpace,
#    tels::Vector, bels::Vector,
#    qs::SingleNumQStrat)

#    u, w = trgauss(qs.quad_rule)
#    qd = [(w[i],SVector(u[1,i],u[2,i])) for i in 1:length(w)]
#    A = _alloc_workspace(qd, g, f, tels, bels)
#    return qd, A
#end


function quadrule(op::LocalOperator, ψ::RefSpace, ϕ::RefSpace, τ, (qd,A), qs::SingleNumQStrat)
    for i in eachindex(qd)
        q = qd[i]
        w, p = q[1], neighborhood(τ,q[2])
        A[i] = (w, p, ψ(p), ϕ(p))
    end
    return A
end


defaultquadstrat(::LocalOperator, _::Space, _::Space) = SingleNumQStrat(6)
# defaultquadstrat(::LocalOperator, ::TetraSupport, ::TetraSupport) = SingleNumQStrat(3)
# defaultquadstrat(op::LocalOperator, a::_LocalBasisOperations, b) = defaultquadstrat(op,support(a),support(b))
# defaultquadstrat(op::LocalOperator, a, b::_LocalBasisOperations) = defaultquadstrat(op,support(a),support(b))
# defaultquadstrat(op::LocalOperator, a::_LocalBasisOperations, b::_LocalBasisOperations) = defaultquadstrat(op,support(a),support(b))

# function quaddata(op::LocalOperator, g, f, tels, bels,
#     qs::SingleNumQStrat)
#     gs = support(g)
#     fs = support(f)
#     quaddata(op,g,f,gs,fs,tels,bels,qs)
# end

# function quaddata(op::LocalOperator, g, f,supg::TriangleSupport,supf::TriangleSupport, tels, bels,
#         qs::SingleNumQStrat)

#     u, w = trgauss(qs.quad_rule)
#     qd = [(w[i],SVector(u[1,i],u[2,i])) for i in 1:length(w)]
#     A = _alloc_workspace(qd, g, f, tels, bels)

#     return qd, A
# end

function quaddata(op::LocalOperator, g, f, 
    tels::Vector{<:CompScienceMeshes.Simplex{U,2}}, bels::Vector{<:CompScienceMeshes.Simplex{U,2}},
    qs::SingleNumQStrat) where {U}

    u, w = trgauss(qs.quad_rule)
    qd = [(w[i],SVector(u[1,i],u[2,i])) for i in 1:length(w)]
    A = _alloc_workspace(qd, g, f, tels, bels)

    return qd, A
end