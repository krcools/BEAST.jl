export lagrangecxd0, lagrangec0d1, duallagrangec0d1
export duallagrangecxd0

export lagdimension
"""
The dimension of the space of Lagrange shape functions of degree d over a simplex of
dimension n is binom(n+d,d) == binom(n+d,n)
"""
function lagdimension end


# D: degree
# C: continuity
# M: mesh type
# T: field type
type LagrangeBasis{D,C,M,T,NF} <: Space{T}
  geo::M
  fns::Vector{Vector{Shape{T}}}
end



# Constructor that automatically deduces MeshType and ScalarType but requires specification
# of the Degree, Cont, and NumFns type parameters
(::Type{LagrangeBasis{D,C,N}}){D,C,M,T,N}(g::M, f::Vector{Vector{Shape{T}}}) = LagrangeBasis{D,C,M,T,N}(g, f)

refspace{D,C,M,T,NF}(space::LagrangeBasis{D,C,M,T,NF}) = LagrangeRefSpace{T,D,dimension(geometry(space))+1,NF}()


function lagrangecxd0(mesh)

  U = universedimension(mesh)
  D1 = dimension(mesh)+1

  geometry = mesh
  num_cells = numcells(mesh)

  # create the local shapes
  fns = Vector{Vector{Shape{Float64}}}(num_cells)
  for i in 1 : num_cells
    fns[i] = [
      Shape(i, 1, 1.0)]
  end

  NF = 1
  LagrangeBasis{0,-1,NF}(geometry, fns)
end


"""
    lagrangec0d1(mesh[, bnd])

Construct the basis of continuous, piecewise linear basis functions subordinate to mesh `mesh`. Basis functions are constructed at vertices in the interionr of the mesh and on the closure of 'bnd'. In particular, leaving out the second argument creates a finite element space subject to homogeneous Dirichlet boundary conditions.
"""
function lagrangec0d1_dirichlet(mesh)

    verts = skeleton(mesh, 0)
    detached = trues(numvertices(mesh))
    for v in cells(verts)
        detached[v] = false
    end

    bnd = boundary(mesh)
    bndverts = skeleton(bnd, 0)
    notonbnd = trues(numvertices(mesh))
    for v in cells(bndverts)
        notonbnd[v] = false
    end

    vertexlist = find(notonbnd .& .!detached)
    lagrangec0d1(mesh, vertexlist, Val{dimension(mesh)+1})
end


function interior_and_junction_vertices(mesh, jct)
    verts = skeleton(mesh, 0)
    detached = trues(numvertices(mesh))
    for v in cells(verts)
        detached[v] = false
    end

    bndfaces = boundary(mesh)
    bndverts = skeleton(bndfaces, 0)
    notonbnd = trues(numvertices(mesh))
    for v in cells(bndverts)
        notonbnd[v] = false
    end

    onjct = broadcast(!, notonbnd)
    overlap_with_junction = overlap_gpredicate(jct)
    for indices in cells(bndfaces)
        bndface = simplex(vertices(bndfaces, indices))
        if overlap_with_junction(bndface)
            continue
        else
            for v in indices
                onjct[v] = false
            end
        end
    end

    vertexlist = find(broadcast(|, onjct, notonbnd) .& broadcast(!,detached))
end

"""
    duallagrangecxd0(mesh, jct) -> basis

Build dual Lagrange piecewise constant elements. Boundary nodes are only considered if they are in the interior of `jct`.
"""
function duallagrangecxd0(mesh, jct)
    vertexlist = interior_and_junction_vertices(mesh, jct)
    duallagrangecxd0(mesh, vertexlist)
end


function duallagrangecxd0(mesh, vertexlist::Vector)

    T = coordtype(mesh)
    fns = Vector{Vector{Shape{T}}}(length(vertexlist))

    fine = barycentric_refinement(mesh)
    vtoc, vton = vertextocellmap(fine)
    for (k,v) in enumerate(vertexlist)
        n = vton[v]
        F = vtoc[v,1:n]
        fns[k] = singleduallagd0(fine, F, v)
    end

    NF = 1
    LagrangeBasis{0,-1,NF}(fine, fns)
end


"""
    singleduallagd0(fine, F, v)

Build a single dual constant Lagrange element a mesh `fine`. `F` contains the indices to cells in the support and v is the index in the vertex list of the defining vertex.
"""
function singleduallagd0(fine, F, v)

    T = coordtype(fine)
    fn = Shape{T}[]
    for cellid in F
        cell = fine.faces[cellid]
        ptch = chart(fine, cell)
        coeff = 1 / volume(ptch) / length(F)
        refid = 1
        push!(fn, Shape(cellid, refid, coeff))
    end

    return fn
end

"""
    lagrangec0d1(mesh; dirichlet=[true|false]) -> basis

Build lagrangec0d1 elements, including (dirichlet=false) or excluding (dirichlet=true) those attached to boundary vertices.
"""
function lagrangec0d1(mesh; dirichlet::Bool=true)
    if dirichlet == false
        return lagrangec0d1(mesh, boundary(mesh))
    else
        return lagrangec0d1_dirichlet(mesh)
    end
end

function lagrangec0d1(mesh, jct)
    vertexlist = interior_and_junction_vertices(mesh, jct)
    lagrangec0d1(mesh, vertexlist, Val{dimension(mesh)+1})
end

function lagrangec0d1(mesh, vertexlist::Vector, ::Type{Val{3}})

      T = coordtype(mesh)
      U = universedimension(mesh)

      cellids, ncells = vertextocellmap(mesh)

      # create the local shapes
      fns = Vector{Shape{Float64}}[]
      sizehint!(fns, length(vertexlist))
      for v in vertexlist

        numshapes = ncells[v]
        numshapes == 0 && continue

        shapes = Vector{Shape{Float64}}(numshapes)
        for s in 1: numshapes
          c = cellids[v,s]
          cell = mesh.faces[c]

          localid = findfirst(cell, v)
          @assert localid != 0

          shapes[s] = Shape(c, localid, 1.0)
        end

        push!(fns, shapes)
      end

      NF = 3
      LagrangeBasis{1,0,NF}(mesh, fns)
    end


function lagrangec0d1(mesh, vertexlist, ::Type{Val{2}})

  T = Float64
  U = universedimension(mesh)
  geometry = mesh

  cellids, ncells = vertextocellmap(mesh)

  # create the local shapes
  numverts = numvertices(mesh)
  fns = Vector{Vector{Shape{Float64}}}()
  sizehint!(fns, length(vertexlist))
  for v in vertexlist

    numshapes = ncells[v]
    numshapes == 0 && continue # skip detached vertices

    shapes = Vector{Shape{Float64}}(numshapes)
    for s in 1: numshapes
      c = cellids[v,s]
      cell = mesh.faces[c]
      if cell[1] == v
        shapes[s] = Shape(c, 1, 1.0)
      elseif cell[2] == v
        shapes[s] = Shape(c, 2, 1.0)
      else
        error("Junctions not supported")
      end
    end

    push!(fns, shapes)
  end

  NF = 2
  LagrangeBasis{1,0,NF}(geometry, fns)
end






duallagrangec0d1(mesh) = duallagrangec0d1(mesh, barycentric_refinement(mesh), x->false, Val{dimension(mesh)+1})
function duallagrangec0d1(mesh, jct)
    jct_pred = inclosure_gpredicate(jct)
    refined = barycentric_refinement(mesh)
    duallagrangec0d1(mesh, refined, jct_pred, Val{dimension(mesh)+1})
end



function duallagrangec0d1(mesh, refined, jct_pred, ::Type{Val{3}})

    T = coordtype(mesh)
    num_faces = dimension(mesh)+1
    fns = Vector{Vector{Shape{T}}}(numcells(mesh))

    # store the fine mesh's vertices in an octree for fast retrieval
    fine_vertices = Octree(refined.vertices)
    uv_ctr = ones(dimension(mesh))/(dimension(mesh)+1)

    vtoc, vton = vertextocellmap(refined)
    for (i,coarse_idcs) in enumerate(cells(mesh))
        fns[i] = Vector{Shape{T}}()

        # It is assumed the vertices of this cell have the same index
        # mesh and its refinement.
        coarse_cell = chart(mesh, coarse_idcs)

        # get the index in fine.vertices of the centroid of coarse_cell
        centroid = barytocart(coarse_cell, uv_ctr)
        I = find(fine_vertices, centroid)
        @assert length(I) == 1
        centroid_id = I[1]

        # get the indx in fine.vertices of the centroid of the faces of coarse_cell
        face_center_ids = Vector{Int}(num_faces)
        for f in 1:num_faces

            # prepare the barycentric coordinate
            uv_face_ctr = ones(dimension(mesh)+1)/(dimension(mesh))
            uv_face_ctr[f] = 0
            uv_face_ctr = uv_face_ctr[1:end-1]

            face_ctr = barytocart(coarse_cell, uv_face_ctr)
            I = find(fine_vertices, face_ctr)
            @assert length(I) == 1
            face_center_ids[f] = I[1]
        end

        n = vton[centroid_id]
        for c in vtoc[centroid_id,1:n]
            fine_idcs = refined.faces[c]
            local_id = findfirst(fine_idcs, centroid_id)
            @assert local_id != 0
            shape = Shape(c, local_id, 1.0)
            push!(fns[i], shape)
        end

        for f in 1:num_faces
            v = face_center_ids[f]
            jct_pred(refined.vertices[v]) && continue
            n = vton[v]
            for c in vtoc[v,1:n]
                fine_idcs = refined.faces[c]
                local_id = findfirst(fine_idcs, v)
                @assert local_id != 0
                shape = Shape(c, local_id, 1/n/2)
                push!(fns[i], shape)
            end
        end

        for f in 1:length(coarse_idcs)
            v = coarse_idcs[f]
            jct_pred(refined.vertices[v]) && continue
            n = vton[v]
            for c in vtoc[v,1:n]
                #fine_idcs = cells(refined,c)
                fine_idcs = refined.faces[c]
                local_id = findfirst(fine_idcs, v)
                @assert local_id != 0
                shape = Shape(c, local_id, 1/n/2)
                push!(fns[i], shape)
            end
        end
    end

    NF = 3
    return LagrangeBasis{1,0,NF}(refined, fns)
end



"""
    duallagrangec0d1(originalmesh, refinedmesh)
  It is the user responsibility to provide two meshes representing the same object.
  The second mesh needs to be obtained using "barycentric_refinement(originalmesh)".
  This basis function creats the dual Lagrange basis function and return an object that contains array of shapes [fns]
  It also return a gemoetry containing the refined mesh
"""
function duallagrangec0d1(mesh, mesh2, pred, ::Type{Val{2}})
  T = Float64
  U = universedimension(mesh)
  # get the information about number of vertices, number of faces , and the maping between vertices and faces for the original mesh
  numverts1 = numvertices(mesh)
  num_cells1 = numcells(mesh)
  cellids1, ncells1=vertextocellmap(mesh)
  # obtain the refined mesh from the original mesh
  #mesh2 = barycentric_refinement(mesh)
  # get the information about number of vertices, number of faces , and the maping between vertices and faces for the refined mesh
  num_cells2 = numcells(mesh2)
  numverts2 = numvertices(mesh2)
  geometry = mesh2
  cellids2, ncells2 = vertextocellmap(mesh2)
  fns = Vector{Vector{Shape{T}}}(num_cells1)
  # We will iterate over the coarse mesh segments to assign all the functions to it.
  for segment_coarse in 1 : num_cells1
    # For the dual Lagrange there is a 6 shapes per segment
    numshapes = (ncells1[segment_coarse]*4) -2
    shapes = Vector{Shape{T}}(numshapes)
    # Now we will get all the smaller faces within the coarse segment
    #i.e The coose segment will have two points, and these tow points are connected to two segmesnts in the finer mesh
    # This will give us a 4 smaller faces per Dual lagrange basis, we store them first in all_faces
    all_faces= Array{SVector{2,Int}}(4)                      # faces in the original segment (4)
    # the follwoing code get the verteciec for each coarse segment
    # then it looks for the two faces connected to each point in the finer mesh
    # if for example the segment would connect to more than two faces we will have
    # mesh.faces[segment_coarse][1],n] and iterate over how many segments [n] are connected
    all_faces[1]=mesh2.faces[cellids2[mesh.faces[segment_coarse][1],1]]
    all_faces[2]=mesh2.faces[cellids2[mesh.faces[segment_coarse][1],2]]
    all_faces[3]=mesh2.faces[cellids2[mesh.faces[segment_coarse][2],1]]
    all_faces[4]=mesh2.faces[cellids2[mesh.faces[segment_coarse][2],2]]

    # now we now the first point of the corse segment will have the left hand side basis
    # and we know that that point is connected to faces 1,2 in the array all_faces
    # but now we want to which one of them is inner segment and which one is at the edge
    # as both of these two faces are in left hand side
    # The inner face will have the left hand side corase point in its first place for example
    # we have original point 2,3 , and the coarse segment (2,3)
    # after refinment we will have for example  points 20,2,21,3,22 with four faces (20,2),(2,21),(21,3),(3,22)
    # so the faces to the left (20,2),(2,21) we can decide the inner one if the original point(2) is at first index
    # same for right hand side faces (21,3),(3,22) if the original point(3) is in the second index
    # then the inner faces are (2,21) and (21,3)
    # the following code does the same checking process and assign the shapes for the dual lagragne right away

    # For the Left hand side faces
    if(all_faces[1][1] == mesh.faces[segment_coarse][1])
        shapes[1]= Shape(cellids2[mesh.faces[segment_coarse][1],2],1,0.5)
        shapes[2]= Shape(cellids2[mesh.faces[segment_coarse][1],1],2,0.5)
        shapes[3]= Shape(cellids2[mesh.faces[segment_coarse][1],1],1,1.0)
    elseif(all_faces[2][1] == mesh.faces[segment_coarse][1])
        shapes[1]= Shape(cellids2[mesh.faces[segment_coarse][1],2],1,1.0)
        shapes[2]= Shape(cellids2[mesh.faces[segment_coarse][1],2],2,0.5)
        shapes[3]= Shape(cellids2[mesh.faces[segment_coarse][1],1],1,0.5)
    end
    # For the Right hand side faces
    if(all_faces[3][2] == mesh.faces[segment_coarse][2])
        shapes[4]= Shape(cellids2[mesh.faces[segment_coarse][2],1],2,1.0)
        shapes[5]= Shape(cellids2[mesh.faces[segment_coarse][2],1],1,0.5)
        shapes[6]= Shape(cellids2[mesh.faces[segment_coarse][2],2],2,0.5)
    elseif(all_faces[4][2] == mesh.faces[segment_coarse][2])
        shapes[4]= Shape(cellids2[mesh.faces[segment_coarse][2],2],2,1.0)
        shapes[5]= Shape(cellids2[mesh.faces[segment_coarse][2],2],1,0.5)
        shapes[6]= Shape(cellids2[mesh.faces[segment_coarse][2],1],2,0.5)
    end
    # Now assign all of these shapes to the relevent segment in the coarse mesh
    fns[segment_coarse]=shapes
  end

  NF = 2
  LagrangeBasis{1,0,NF}(geometry, fns)
end


curl(space::LagrangeBasis{1,0}, geo, fns) = RTBasis(geo, fns)

#
# Sclar trace for Laggrange element based spaces
#
function strace(X::LagrangeBasis{1,0}, geo, fns::Vector)
    # dimension of the boundary
    n = dimension(geo)
    # degree of the space
    d = 1
    # number of local shape functions
    NF = binomial(n+d,n)

    LagrangeBasis{1,0,NF}(geo, fns)
end
