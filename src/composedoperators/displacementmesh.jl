##### displacement meshes
abstract type DisplacementMesh{T} <: CompScienceMeshes.AbstractMesh{3,3,T} end

"""
GlobalDisplacementMesh(mesh,epsilon)
    every chart is displaced with epsilon allong the normal
"""
struct GlobalDisplacementMesh{T} <: DisplacementMesh{T}
    mesh::Mesh{<:Any,<:Any,T}
    epsilon 
end
struct DisplacementChart
    chart 
    epsilon
end
CompScienceMeshes.mesh(m::GlobalDisplacementMesh) = m.mesh
displacementchart(m::GlobalDisplacementMesh,i) = DisplacementChart(chart(m.mesh,i),m.epsilon)
displacementchart(m::Mesh,i) = DisplacementChart(chart(m,i),-1.0)
CompScienceMeshes.chart(m::DisplacementMesh,i) = chart(mesh(m),i)
CompScienceMeshes.cells(m::DisplacementMesh) = cells(mesh(m))
CompScienceMeshes.celltype(m::DisplacementMesh) = CompScienceMeshes.celltype(mesh(m))
CompScienceMeshes.indextype(m::DisplacementMesh,a) = CompScienceMeshes.indextype(mesh(m),a)
CompScienceMeshes.celltype(m::DisplacementMesh,a) = CompScienceMeshes.celltype(mesh(m),a)
CompScienceMeshes.indices(m::DisplacementMesh,n) = CompScienceMeshes.indices(mesh(m),n)
CompScienceMeshes.vertices(m::DisplacementMesh) = CompScienceMeshes.vertices(mesh(m))
CompScienceMeshes.numvertices(m::DisplacementMesh) = CompScienceMeshes.numvertices(mesh(m))
CompScienceMeshes.numcells(m::DisplacementMesh) = CompScienceMeshes.numcells(mesh(m))
CompScienceMeshes.vertextype(m::DisplacementMesh) = CompScienceMeshes.vertextype(mesh(m))


function displacement(a::DisplacementChart,b::DisplacementChart)
    rel_orient = sign(dot(normal(a.chart),normal(b.chart)))
    if a.epsilon*rel_orient ≈ b.epsilon
        return 0
    end
    testdisp = a.epsilon*normal(a.chart)
    trialdisp = b.epsilon*normal(b.chart)
    
    return sign(dot(trialdisp-testdisp,normal(b.chart)))
end