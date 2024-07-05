using BEAST
using CompScienceMeshes
using Plotly
using StaticArrays
using LinearAlgebra

"""
plot_basis(mesh, lagrangebasis)

Plots the wireframe of the mesh and the positions stored in the Lagrange basis.
The positions are marked with a dot marker along with a number that indicates the index at which it is being stored in the
lagrange basis.

"""
function plot_basis(Γ5::Mesh, Y5::BEAST.LagrangeBasis)
    plot(patch(Γ5,ones(numcells(Γ5))))
    trace1 = (wireframe(Γ5))
    trace2 = (scatter3d(x=getindex.(Y5.pos,1), y=getindex.(Y5.pos,2), z=getindex.(Y5.pos,3), mode="text",text=string.(range(start=1,step=1,stop=length(Y5.pos)))))
    trace3 = (scatter3d(x=getindex.(Y5.pos,1), y=getindex.(Y5.pos,2), z=getindex.(Y5.pos,3), mode="markers", marker_size=4))
    plot([trace1,trace2,trace3])
end

"""
plot_global_vertices(mesh)

Plots the wireframe of the mesh along with its vertices with a number as its marker. The number 'n' corresponding to a vertex 'v'
indicates the index at which v occurs in the mesh.

"""
function plot_global_vertices(Γ::Mesh)
    plot(patch(Γ,ones(numcells(Γ))))
    trace1 = (wireframe(Γ))
    trace2 = scatter3d(x=getindex.(Γ.vertices,1),y=getindex.(Γ.vertices,2),z=getindex.(Γ.vertices,3),mode="text",text=string.(range(start=1,step=1,stop=numvertices(Γ))))
    plot([trace1, trace2])
end


"""
lagrangecxd1(mesh)

Function that generates a piecewise constant Lagrange basis given a mesh.

"""
function lagrangecxd1(mesh::Mesh)
    num_cells = numcells(mesh)
    T = coordtype(mesh)

    fns = Vector{Vector{BEAST.Shape{T}}}(undef, (num_cells*3))
    pos = Vector{vertextype(mesh)}(undef, (num_cells*3))

    combs = [(0,0,1),(0,1,0),(1,0,0)]
    for n in 1:num_cells
        ref_id = 3 #local barycentric coordinates all start with the id 3 and go down to 1
        charted = chart(mesh,n)
        for (ind,(i,j,k)) in enumerate(combs)
            index = ((n-1)*3 + ind)
            fns[index] = [BEAST.Shape(n,ref_id,T(1.0))]
            pos[index] = barytocart(charted,[i,j,k])
            ref_id -= 1
        end
    end
    NF = 3
    BEAST.LagrangeBasis{1,-1,NF}(mesh,fns,pos)
end

"""
lagrangecxd2(mesh)

Function that generates a piecewise constant basis of order 2 given a mesh. 
Each cell has 6 nodes and thus 6 local shape function.

"""
function lagrangecxd2(mesh::Mesh)
    num_cells = numcells(mesh)
    T = coordtype(mesh)

    fns = Vector{Vector{BEAST.Shape{T}}}(undef, (num_cells*6))
    pos = Vector{vertextype(mesh)}(undef, (num_cells*6))

    index = 1
    for n in 1:num_cells
        ref_id_cor = 3 #ref_id for the corner vertices
        ref_id_mid = 4 #ref_id for the vertices at the midpoint of each edge
        charted = chart(mesh,n)
        
        for i in 0:2
            for j in 0:2
                if(i+j <= 2) #lower triangular structure
                    k = 2 - i - j
                    if(i==2||j==2||k==2) # local corner vertices
                        fns[index] = [BEAST.Shape(n,ref_id_cor,T(1.0))]
                        pos[index] = barytocart(charted,[i/2,j/2,k/2])
                        ref_id_cor -= 1
                    else #vertices at the mid point of each edge of the current cell
                        fns[index] = [BEAST.Shape(n,ref_id_mid,T(1.0))]
                        pos[index] = barytocart(charted,[i/2,j/2,k/2])
                        ref_id_mid += 1
                    end
                    index += 1
                end
            end
        end
    end

    NF = 6
    BEAST.LagrangeBasis{2,-1,NF}(mesh,fns,pos)
end

"""
lagrangecxdp(mesh,order)

Function that generates a piece wise constant of basis of order p. Each cell has ((p+1)*(p+2)/2)
number of nodes.

"""

function lagrangecxdp(mesh::Mesh,p::Int)
    num_cells = numcells(mesh)
    T = coordtype(mesh)

    numshapes = Int((p+1)*(p+2)/2) #number of shapes in a cell

    fns = Vector{Vector{BEAST.Shape{T}}}(undef,(num_cells*numshapes))
    pos = Vector{vertextype(mesh)}(undef,(num_cells*numshapes))

    index = 1
    for n = 1:num_cells
        ref_id_cor = 3 # ref_id to address corner vertices of each cell
        ref_id_im = 4 # ref_id to address interior and points on edges
        charted = chart(mesh,n)
        
        for i in 0:p
            for j in 0:p
                if(i+j <= p) #lower triagular structure
                    k = p - i - j
                    if(i==p||j==p||k==p) #local corner vertices
                        fns[index] = [BEAST.Shape(n,ref_id_cor,T(1.0))]
                        pos[index] = barytocart(charted,[i/p,j/p,k/p])
                        ref_id_cor -= 1
                    else #local non-corner vertices
                        fns[index] = [BEAST.Shape(n,ref_id_im,T(1.0))]
                        pos[index] = barytocart(charted,[i/p,j/p,k/p])
                        ref_id_im += 1
                    end
                    index+= 1
                end
            end
        end
    end

    NF = numshapes
    BEAST.LagrangeBasis{p,-1,NF}(mesh,fns,pos)
end