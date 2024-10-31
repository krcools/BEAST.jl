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

displacementchart(m::GlobalDisplacementMesh,i) = DisplacementChart(chart(m.mesh,i),m.epsilon)
displacementchart(m::Mesh,i) = DisplacementChart(chart(m,i),-1.0)
CompScienceMeshes.chart(m::DisplacementMesh,i) = chart(m.mesh,i)
function displacement(a::DisplacementChart,b::DisplacementChart)
    rel_orient = sign(dot(normal(a.chart),normal(b.chart)))
    if a.epsilon*rel_orient ≈ b.epsilon
        return 0
    end
    testdisp = a.epsilon*normal(a.chart)
    trialdisp = b.epsilon*normal(b.chart)
    
    return sign(dot(trialdisp-testdisp,normal(a.chart)))
end