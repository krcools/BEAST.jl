import Base.-
struct MinusMeshPointNM{T}
    meshpoint::T
end
# struct MinusSimplex{T}
#     simplex::T
# end

normal(a::MinusMeshPointNM) = -normal(a.meshpoint)
# normal(a::MinusSimplex) = -normal(a.MinusSimplex)

function map(b::Union{CompScienceMeshes.MeshPointNM,MinusMeshPointNM},sign::Int) 
    sign==-1 && (return -b)
    sign == 1 && (return b)
    @error "sign was not 1 or -1"
end

-(b::CompScienceMeshes.MeshPointNM) = MinusMeshPointNM(b)
-(b::MinusMeshPointNM) = b.meshpoint
# -(b::Simplex) = MinusSimple(b)
# -(b::MinusSimplex) = b.simplex

Base.length(m::MinusMeshPointNM) = length(m.meshpoint)
Base.getindex(p::MinusMeshPointNM, i::Int) = p.meshpoint[i]
cartesian(m::MinusMeshPointNM) = cartesian(m.meshpoint)
parametric(m::MinusMeshPointNM) = parametric(m.meshpoint)
chart(m::MinusMeshPointNM) = chart(m.meshpoint)
barycentric(mp::MinusMeshPointNM) = barycentric(mp.meshpoint)
jacobian(mp::MinusMeshPointNM) = jacobian(mp.meshpoint)
tangents(mp::MinusMeshPointNM,i) = tangents(mp.meshpoint,i)
