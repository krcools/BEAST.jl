import Base.sign
struct TraceSimplex{T,U}
    simp::T
    direction::SVector{3,U}
end
simplex(t::TraceSimplex) = t.simp
direction(t::TraceSimplex) = t.direction
struct TraceMeshPointNM{T,U}
    neighborhood::T
    direction::SVector{3,U}
end
direction(p::TraceMeshPointNM) = p.direction
function CompScienceMeshes.neighborhood(p::TraceSimplex, bary)
    TraceMeshPointNM(neighborhood(p.simp,bary),p.direction)
end
CompScienceMeshes.cartesian(p::TraceMeshPointNM) = cartesian(p.neighborhood)
CompScienceMeshes.parametric(p::TraceMeshPointNM) = parametric(p.neighborhood)
chart(p::TraceMeshPointNM) = chart(p.neighborhood)
CompScienceMeshes.normal(p::TraceMeshPointNM) = normal(p.neighborhood)

CompScienceMeshes.jacobian(p::TraceMeshPointNM) = jacobian(p.neighborhood)
CompScienceMeshes.tangents(p::TraceMeshPointNM,i) = tangents(p.neighborhood,i)
CompScienceMeshes.barycentric(a::TraceMeshPointNM) = barycentric(a.neighborhood)
coordtype(p::TraceSimplex) = coordtype(p.simp)



struct OrientedMesh{U,D1,T} <: CompScienceMeshes.AbstractMesh{U,D1,T}
    mesh::CompScienceMeshes.AbstractMesh{U,D1,T}
    normalmesh::CompScienceMeshes.AbstractMesh{U,D1,T}
end
struct TraceMesh{U,D1,T} <: CompScienceMeshes.AbstractMesh{U,D1,T}
    mesh::CompScienceMeshes.AbstractMesh{U,D1,T}
    direction::SVector
end
TraceMesh(a::CompScienceMeshes.AbstractMesh{U,D,T}) where {U,D,T} = TraceMesh(a,zeros(SVector{numcells(a),SVector{3,T}}))
+(a::TraceMesh{U,D,T},b::SVector{3,T}) where {U,D,T} = TraceMesh(mesh(a),direction(a).+Ref(b))
+(a::TraceMesh,b::SVector) = TraceMesh(mesh(a),a.direction+b)
+(a::SVector,b::TraceMesh) = b+a



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
    return flip_normal(c)
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

CompScienceMeshes.domain(ch::TraceSimplex) = CompScienceMeshes.domain(simplex(ch))
CompScienceMeshes.dimension(::Type{TraceSimplex{CompScienceMeshes.Simplex{U,D,C,N,T},T}}) where {U,D,C,N,T} = D
CompScienceMeshes.coordtype(::Type{TraceSimplex{CompScienceMeshes.Simplex{U,D,C,N,T},T}}) where {U,D,C,N,T} = T


function CompScienceMeshes.quadpoints(chart::TraceSimplex, rule)
    PV = quadpoints(CompScienceMeshes.domain(chart), rule)
    map(PV) do pv
        q = neighborhood(chart, pv[1])
        w = jacobian(q)*pv[2]
        (q,w)
    end
end
CompScienceMeshes.verticeslist(p::TraceSimplex) = verticeslist(simplex(p))
CompScienceMeshes.permute_vertices(simp::TraceSimplex,I) = TraceSimplex(CompScienceMeshes.permute_vertices(simplex(simp),I),direction(simp))
CompScienceMeshes.overlap(a::TraceSimplex,b::TraceSimplex) = overlap(simplex(a),simplex(b))
CompScienceMeshes.overlap(a::CompScienceMeshes.Simplex,b::TraceSimplex) = overlap(a,simplex(b))
CompScienceMeshes.overlap(a::TraceSimplex,b::CompScienceMeshes.Simplex) = overlap(simplex(a),b)

CompScienceMeshes.intersection(a::CompScienceMeshes.Simplex,b::TraceSimplex) = intersection(a,simplex(b))

# function CompScienceMeshes.intersection(a::TraceSimplex,b::Union{TraceSimplex,CompScienceMeshes.Simplex})
#     int = intersection(simplex(a),b)
#     TraceSimplex.(int,sign(a))
# end

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
CompScienceMeshes.volume(a::TraceSimplex) = volume(simplex(a))
CompScienceMeshes.dimension(a::TraceSimplex) = dimension(simplex(a))
CompScienceMeshes.tangents(a::TraceSimplex,i::Int) = tangents(simplex(a),i)
CompScienceMeshes.carttobary(a::TraceSimplex,b::SVector{T}) where {T} = carttobary(simplex(a),b)

