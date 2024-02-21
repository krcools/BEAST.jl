import Base.sign
import Base.getindex
struct TraceSimplex{T,U}
    simp::T
    direction::SVector{3,U}
end
TraceSimplex(t::TraceSimplex,d::SVector{3,U}) where {U} = TraceSimplex(simplex(t),direction(t)+d)
simplex(t::TraceSimplex) = t.simp
direction(t::TraceSimplex) = t.direction
struct TraceMeshPointNM{T,U}
    neighborhood::T
    direction::SVector{3,U}
end
direction(p::TraceMeshPointNM) = p.direction

CompScienceMeshes.neighborhood(p::TraceSimplex, bary) = TraceMeshPointNM(neighborhood(p.simp,bary),p.direction)
CompScienceMeshes.coordtype(p::TraceSimplex) = coordtype(p.simp)
CompScienceMeshes.domain(ch::TraceSimplex) = CompScienceMeshes.domain(simplex(ch))
CompScienceMeshes.dimension(::Type{TraceSimplex{CompScienceMeshes.Simplex{U,D,C,N,T},T}}) where {U,D,C,N,T} = D
CompScienceMeshes.coordtype(::Type{TraceSimplex{CompScienceMeshes.Simplex{U,D,C,N,T},T}}) where {U,D,C,N,T} = T
CompScienceMeshes.verticeslist(p::TraceSimplex) = verticeslist(simplex(p))
CompScienceMeshes.permute_vertices(simp::TraceSimplex,I) = TraceSimplex(CompScienceMeshes.permute_vertices(simplex(simp),I),direction(simp))
CompScienceMeshes.overlap(a::TraceSimplex,b::TraceSimplex) = overlap(simplex(a),simplex(b))
CompScienceMeshes.overlap(a::CompScienceMeshes.Simplex,b::TraceSimplex) = overlap(a,simplex(b))
CompScienceMeshes.overlap(a::TraceSimplex,b::CompScienceMeshes.Simplex) = overlap(simplex(a),b)
CompScienceMeshes.intersection(a::CompScienceMeshes.Simplex,b::TraceSimplex) = intersection(a,simplex(b))
CompScienceMeshes.volume(a::TraceSimplex) = volume(simplex(a))
CompScienceMeshes.dimension(a::TraceSimplex) = dimension(simplex(a))
CompScienceMeshes.tangents(a::TraceSimplex,i::Int) = tangents(simplex(a),i)
CompScienceMeshes.carttobary(a::TraceSimplex,b::SVector{T}) where {T} = carttobary(simplex(a),b)
CompScienceMeshes.center(a::TraceSimplex) = TraceMeshPointNM(center(simplex(a)),direction(a))
CompScienceMeshes.normal(a::TraceSimplex) = normal(simplex(a))
CompScienceMeshes.barytocart(a::TraceSimplex,b) = barytocart(simplex(a),b)
function CompScienceMeshes.quadpoints(chart::TraceSimplex, rule)
    PV = quadpoints(CompScienceMeshes.domain(chart), rule)
    map(PV) do pv
        q = neighborhood(chart, pv[1])
        w = jacobian(q)*pv[2]
        (q,w)
    end
end

function CompScienceMeshes.intersection2(a::TraceSimplex,b::TraceSimplex)
    int = intersection2(simplex(a),simplex(b))
    return [(TraceSimplex(i[1],direction(a)),TraceSimplex(i[2],direction(b))) for i in int]
end
function CompScienceMeshes.intersection2(a::CompScienceMeshes.Simplex,b::TraceSimplex)
    int = intersection2(a,simplex(b))
    return [(TraceSimplex(i[1],direction(b)*0),TraceSimplex(i[2],direction(b))) for i in int]
end
function CompScienceMeshes.intersection2(a::TraceSimplex,b::CompScienceMeshes.Simplex)
    int = intersection2(simplex(a),b)
    return [(TraceSimplex(i[1],direction(a)),TraceSimplex(i[2],direction(a)*0)) for i in int]
end




CompScienceMeshes.cartesian(p::TraceMeshPointNM) = cartesian(p.neighborhood)
CompScienceMeshes.parametric(p::TraceMeshPointNM) = parametric(p.neighborhood)
CompScienceMeshes.chart(p::TraceMeshPointNM) = chart(p.neighborhood)
CompScienceMeshes.normal(p::TraceMeshPointNM) = normal(p.neighborhood)
CompScienceMeshes.jacobian(p::TraceMeshPointNM) = jacobian(p.neighborhood)
CompScienceMeshes.tangents(p::TraceMeshPointNM,i) = tangents(p.neighborhood,i)
CompScienceMeshes.barycentric(a::TraceMeshPointNM) = barycentric(a.neighborhood)




struct OrientedMesh{U,D1,T} <: CompScienceMeshes.AbstractMesh{U,D1,T}
    mesh::CompScienceMeshes.AbstractMesh{U,D1,T}
    normalmesh::CompScienceMeshes.AbstractMesh{U,D1,T}
end
struct TraceMesh{U,D1,T} <: CompScienceMeshes.AbstractMesh{U,D1,T}
    mesh::CompScienceMeshes.AbstractMesh{U,D1,T}
    direction::Vector{SVector{3,T}}
end
TraceMesh(a::TraceMesh,b::Vector{SVector{3,T}}) where {T} = a+b
TraceMesh(a::CompScienceMeshes.AbstractMesh{U,D,T}) where {U,D,T} = TraceMesh(a,zeros(SVector{3,T},(numcells(a))))
+(a::TraceMesh{U,D,T},b::SVector{3,T}) where {U,D,T} = TraceMesh(mesh(a),direction(a).+Ref(b))
+(a::TraceMesh{U,D,T},b::Vector{SVector{3,T}}) where {U,D,T} = TraceMesh(mesh(a),a.direction+b)
+(a::Union{SVector,Vector},b::TraceMesh) = b+a
CompScienceMeshes.indices(t::TraceMesh,i::Int) = CompScienceMeshes.indices(mesh(t),i)
CompScienceMeshes.numvertices(t::TraceMesh) = CompScienceMeshes.numvertices(mesh(t))
CompScienceMeshes.vertextype(t::TraceMesh) = CompScienceMeshes.vertextype(mesh(t))
CompScienceMeshes.universedimension(t::TraceMesh) = CompScienceMeshes.universedimension(mesh(t))

mesh(p::OrientedMesh) = p.mesh
mesh(p::TraceMesh) = p.mesh
CompScienceMeshes.normal(p::OrientedMesh) = p.normalmesh
direction(p::TraceMesh) = p.direction

function CompScienceMeshes.chart(p::OrientedMesh,i)
    c = chart(mesh(p),i)
    n1 = normal(c)
    n2 = normal(chart(normal(p),i))
    d = dot(n1,n2)
    @assert abs(d) ≈ 1.0
    sign(d) == 1 && return c
    return mirror(c)
end
function CompScienceMeshes.chart(p::TraceMesh,i)
    c = chart(mesh(p),i)
    d = direction(p)[i]
    TraceSimplex(c,d)
end

function same_geometry(m1,m2)
    return (vertices(m1) == vertices(m2))*prod(_ispermutation.(cells(m1),cells(m2)))
end

_ispermutation(a,b) = sort(a)==sort(b)
CompScienceMeshes.cells(p::OrientedMesh) = cells(mesh(p))
CompScienceMeshes.cells(p::TraceMesh) = cells(mesh(p))
CompScienceMeshes.vertices(p::TraceMesh) = vertices(mesh(p))


# function buffachristiansen(Γ::TraceMesh,γ=mesh(coordtype(Γ),1,3);kwargs...)
#     rtspace = buffachristiansen(mesh(Γ),γ;kwargs...)
#     geo = geometry(rtspace)

#     direction = [Γ.direction[CompScienceMeshes.parent(geo,i)] for i in 1:length(geo)]
#    redefine_geometrie(rtspace,TraceMesh(geo,direction))
# end

function CompScienceMeshes.barycentric_refinement(m::TraceMesh;kwargs...)
    geo = mesh(m)
    fine = barycentric_refinement(geo)
    direction = [m.direction[CompScienceMeshes.parent(fine,i)] for i in 1:length(fine)]
    return TraceMesh(fine,direction)
end
CompScienceMeshes.vertextocellmap(m::TraceMesh) = vertextocellmap(mesh(m))
Base.getindex(m::TraceMesh,i::Vector{Int}) = mesh(m)[i]