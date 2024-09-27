
"""
The dimension of the space of Lagrange shape functions of degree d over a simplex of
dimension n is binom(n+d,d) == binom(n+d,n)
"""
function lagdimension end


# D: degree
# C: continuity
# M: mesh type
# T: field type
# NF: number of local shape functions
mutable struct LagrangeBasis{D,C,M,T,NF,P} <: Space{T}
  geo::M
  fns::Vector{Vector{Shape{T}}}
  pos::Vector{P}
end



# Constructor that automatically deduces MeshType and ScalarType but requires specification
# of the Degree, Cont, and NumFns type parameters
function LagrangeBasis{D,C,N}(mesh::M, fns::Vector{Vector{Shape{T}}}, pos::Vector{P}) where {P,T,M,N,C,D}
    LagrangeBasis{D,C,M,T,N,P}(mesh, fns, pos)
end

refspace(space::LagrangeBasis{D,C,M,T,NF}) where {D,C,M,T,NF} = LagrangeRefSpace{T,D,dimension(geometry(space))+1,NF}()
subset(s::S,I) where {S <: Space} = S(s.geo, s.fns[I], s.pos[I])

function lagrangecxd0(mesh)

    U = universedimension(mesh)
    D1 = dimension(mesh)+1
    T = coordtype(mesh)
    geometry = mesh
    num_cells = numcells(mesh)

    # create the local shapes
    fns = Vector{Vector{Shape{T}}}(undef,num_cells)
    pos = Vector{vertextype(mesh)}(undef,num_cells)
    for (i,cell) in enumerate(mesh)
        fns[i] = [Shape(i, 1, T(1.0))]
        pos[i] = cartesian(center(chart(mesh, cell)))
  end

  NF = 1
  LagrangeBasis{0,-1,NF}(geometry, fns, pos)
end

"""
    unitfunctioncxd0(mesh)

Constructs a constant function with value 1 on `mesh`.
"""
function unitfunctioncxd0(mesh)

    T = coordtype(mesh)
    geometry = mesh

    # create the local shapes
    fns = Vector{Vector{Shape{T}}}(undef, 1)
    pos = Vector{vertextype(mesh)}(undef, 1)
    fns[1] = [Shape(i, 1, T(1.0)) for (i, cell) in enumerate(mesh)]

    # Arguably, the position is fairly meaningless
    # in case of a global function. Might be replaced by something
    # more useful.
    # For now, we fill it with the average position of the shape functions
    p = vertextype(mesh)(0.0, 0.0, 0.0)
    for cell in mesh
        p += cartesian(center(chart(mesh, cell)))
    end
    pos[1] = p ./ numcells(mesh)

    NF = 1
    LagrangeBasis{0,-1,NF}(geometry, fns, pos)
end

"""
    unitfunctionc0d1(mesh)

Constructs a constant function with value 1 on `mesh` consisting of linear shapes. For dirichlet=true goes to zero on the boundary.
"""
function unitfunctionc0d1(mesh; dirichlet=true)
    if dirichlet == false
        return unitfunctionc0d1(mesh, skeleton(mesh,0))
    else
        return unitfunctionc0d1_dirichlet(mesh)
    end
end

function unitfunctionc0d1_dirichlet(mesh)

    T = coordtype(mesh)
    
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

    vertexlist = findall(notonbnd .& .!detached)

    cellids, ncells = vertextocellmap(mesh)

    Cells = cells(mesh)
    Verts = vertices(mesh)

    # create the local shapes
    fns = Vector{Vector{Shape{T}}}(undef, 1)
    pos = Vector{vertextype(mesh)}(undef, 1)

    numshapes = sum(ncells[vertexlist])
    shapes = Vector{Shape{T}}(undef,numshapes)
    n = 0
    for v in vertexlist
        nshapes = ncells[v]
        nshapes == 0 && continue

        for s in 1: nshapes
            c = cellids[v,s]

            cell = Cells[c]

            localid = something(findfirst(isequal(v), cell),0)
            @assert localid != 0

            shapes[s+n] = Shape(c, localid, T(1.0))

        end
        n += nshapes
    end
    fns[1] = shapes
    p = sum(mesh.vertices[vertexlist])/length(vertexlist)
    pos[1] = p

    NF = 3
    LagrangeBasis{1,0,NF}(mesh, fns, pos)
end

function unitfunctionc0d1(mesh, nodes::CompScienceMeshes.AbstractMesh{U,1} where {U})
    Conn = connectivity(nodes, mesh, abs)
    rows = rowvals(Conn)
    vals = nonzeros(Conn)

    T = coordtype(mesh)
    P = vertextype(mesh)
    S = Shape{T}

    fns = Vector{Vector{Shape{T}}}(undef, 1)
    pos = Vector{vertextype(mesh)}(undef, 1)
    fn = Vector{S}()
    for (i,node) in enumerate(nodes)
        for k in nzrange(Conn,i)
            cellid = rows[k]
            refid  = vals[k]
            push!(fn, Shape(cellid, refid, T(1.0)))
        end
    end
    fns[1] = fn
    p = sum(nodes.vertices)/length(nodes.vertices)
    pos[1] = p

    NF = dimension(mesh) + 1
    LagrangeBasis{1,0,NF}(mesh, fns, pos)
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

    vertexlist = findall(notonbnd .& .!detached)
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

    vertexlist = findall(broadcast(|, onjct, notonbnd) .& broadcast(!,detached))
end

"""
    duallagrangecxd0(mesh, jct) -> basis

Build dual Lagrange piecewise constant elements. Boundary nodes are only considered if
they are in the interior of `jct`.

The default dual function (`interpolatory=false`) is similar to the one depicted
in Figure 3 of  Buffa et al (doi: 10.1090/S0025-5718-07-01965-5), with the
difference that each individual shape function is normalized with respect to 
the area so that overall the integral over the dual function is one.

When `interpolatory=true` is used, the function value is one on the support, and thus,
it gives rise to a partition of unity.
"""
function duallagrangecxd0(mesh, jct=CompScienceMeshes.mesh(coordtype(mesh), dimension(mesh)-1); interpolatory=false)
    vertexlist = interior_and_junction_vertices(mesh, jct)
    duallagrangecxd0(mesh, vertexlist; interpolatory=interpolatory)
end


function duallagrangecxd0(mesh, vertexlist::Vector{Int}; interpolatory=false)

    T = coordtype(mesh)

    fns = Vector{Vector{Shape{T}}}(undef,length(vertexlist))
    pos = Vector{vertextype(mesh)}()

    fine = barycentric_refinement(mesh)
    vtoc, vton = vertextocellmap(fine)
    verts = vertices(mesh)
    for (k,v) in enumerate(vertexlist)
        n = vton[v]
        F = vtoc[v,1:n]
        fns[k] = singleduallagd0(fine, F, v, interpolatory=interpolatory)
        push!(pos, verts[v])
    end

    NF = 1
    LagrangeBasis{0,-1,NF}(fine, fns, pos)
end


function duallagrangecxd0(mesh, vertices::CompScienceMeshes.AbstractMesh{U,1}; interpolatory=false) where {U}
    # vertexlist = Int[v[1] for v in vertices]
    vertexlist =Int[CompScienceMeshes.indices(vertices, v)[1] for v in vertices]
    return duallagrangecxd0(mesh, vertexlist; interpolatory=interpolatory)
end


"""
    singleduallagd0(fine, F, v; interpolatory=false)

Build a single dual constant Lagrange element a mesh `fine`. `F` contains the indices
to cells in the support and v is the index in the vertex list of the defining vertex.

The default dual function (`interpolatory=false`) is similar to the one depicted
in Figure 3 of  Buffa et al (doi: 10.1090/S0025-5718-07-01965-5), with the
difference that each individual shape function is normalized with respect to
the area so that overall the integral over the dual function is one.

When `interpolatory=true` is used, the function value is one on the support, and thus,
it gives rise to a partition of unity.
"""
function singleduallagd0(fine, F, v; interpolatory=false)

    T = coordtype(fine)
    fn = Shape{T}[]
    for cellid in F
        # cell = cells(fine)[cellid]
        ptch = chart(fine, cellid)
        coeff = interpolatory ? T(1.0) : 1 / volume(ptch) / length(F)
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
        # return lagrangec0d1(mesh, boundary(mesh))
        return lagrangec0d1(mesh, skeleton(mesh,0))
    else
        return lagrangec0d1_dirichlet(mesh)
    end
end

function lagrangec0d1(mesh, jct)
    vertexlist = interior_and_junction_vertices(mesh, jct)
    lagrangec0d1(mesh, vertexlist, Val{dimension(mesh)+1})
end

# build continuous linear Lagrange elements on a 2D manifold
function lagrangec0d1(mesh, vertexlist::Vector, ::Type{Val{3}})

    T = coordtype(mesh)
    U = universedimension(mesh)

    cellids, ncells = vertextocellmap(mesh)

    Cells = cells(mesh)
    Verts = vertices(mesh)

    # create the local shapes
    fns = Vector{Shape{T}}[]
    pos = Vector{vertextype(mesh)}()

    sizehint!(fns, length(vertexlist))
    sizehint!(pos, length(vertexlist))
    for v in vertexlist

        numshapes = ncells[v]
        numshapes == 0 && continue

        shapes = Vector{Shape{T}}(undef,numshapes)
        for s in 1: numshapes
            c = cellids[v,s]
            # cell = mesh.faces[c]
            cell = Cells[c]

            localid = something(findfirst(isequal(v), cell),0)
            @assert localid != 0

            shapes[s] = Shape(c, localid, T(1.0))
        end

        push!(fns, shapes)
        push!(pos, Verts[v])
    end

    NF = 3
    LagrangeBasis{1,0,NF}(mesh, fns, pos)
end

# for manifolds of dimension 1
function lagrangec0d1(mesh, vertexlist, ::Type{Val{2}})

    T = coordtype(mesh)
    U = universedimension(mesh)
    P = vertextype(mesh)

    geometry = mesh

    cellids, ncells = vertextocellmap(mesh)

    # create the local shapes
    numverts = numvertices(mesh)

    fns = Vector{Vector{Shape{T}}}()
    pos = Vector{P}()

    sizehint!(fns, length(vertexlist))
    sizehint!(pos, length(vertexlist))
    for v in vertexlist

        numshapes = ncells[v]
        numshapes == 0 && continue # skip detached vertices

        shapes = Vector{Shape{T}}(undef,numshapes)
        for s in 1: numshapes
            c = cellids[v,s]
            cell = mesh.faces[c]
            if cell[1] == v
                shapes[s] = Shape(c, 1, T(1.0))
            elseif cell[2] == v
                shapes[s] = Shape(c, 2, T(1.0))
            else
                error("Junctions not supported")
            end
        end

        push!(fns, shapes)
        push!(pos, mesh.vertices[v])
    end

    NF = 2
    LagrangeBasis{1,0,NF}(geometry, fns, pos)
end


function lagrangec0d1(mesh, nodes::CompScienceMeshes.AbstractMesh{U,1} where {U})

    Conn = connectivity(nodes, mesh, abs)
    rows = rowvals(Conn)
    vals = nonzeros(Conn)

    T = coordtype(mesh)
    P = vertextype(mesh)
    S = Shape{T}

    fns = Vector{Vector{S}}()
    pos = Vector{P}()
    for (i,node) in enumerate(nodes)
        fn = Vector{S}()
        for k in nzrange(Conn,i)
            cellid = rows[k]
            refid  = vals[k]
            push!(fn, Shape(cellid, refid, T(1.0)))
        end
        push!(fns,fn)
        push!(pos,cartesian(center(chart(nodes,node))))
    end

    NF = dimension(mesh) + 1
    LagrangeBasis{1,0,NF}(mesh, fns, pos)
end


function lagrangec0d2(mesh::CompScienceMeshes.AbstractMesh{U,3},
    nodes::CompScienceMeshes.AbstractMesh{U,1},
    edges::CompScienceMeshes.AbstractMesh{U,2}) where {U}

    Conn = connectivity(nodes, mesh, abs)
    rows = rowvals(Conn)
    vals = nonzeros(Conn)

    T = coordtype(mesh)
    P = vertextype(mesh)
    S = Shape{T}

    fns = Vector{Vector{S}}()
    pos = Vector{P}()
    for (i,node) in enumerate(nodes)
        fn = Vector{S}()
        for k in nzrange(Conn,i)
            cellid = rows[k]
            refid  = vals[k]
            push!(fn, Shape(cellid, refid, T(1.0)))
        end
        push!(fns,fn)
        push!(pos,cartesian(center(chart(nodes,node))))
    end

    Conn = connectivity(edges, mesh, abs)
    rows = rowvals(Conn)
    vals = nonzeros(Conn)

    for (i,edge) in enumerate(edges)
        fn = Vector{S}()
        for k in nzrange(Conn,i)
            cellid = rows[k]
            refid  = vals[k]
            push!(fn, Shape(cellid, 3+refid, T(1.0)))
        end
        push!(fns,fn)
        push!(pos,cartesian(center(chart(edges,edge))))
    end

    NF = 6
    LagrangeBasis{2,0,NF}(mesh, fns, pos)
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

    fns = Vector{Vector{Shape{T}}}(undef,numcells(mesh))
    pos = Vector{vertextype(mesh)}()

    # store the fine mesh's vertices in an octree for fast retrieval
    fine_vertices = Octree(vertices(refined))
    uv_ctr = ones(dimension(mesh))/(dimension(mesh)+1)

    vtoc, vton = vertextocellmap(refined)
    for (i,p) in enumerate(mesh)
        coarse_idcs = CompScienceMeshes.indices(mesh, p)
        coarse_chart = chart(mesh,p)

        fns[i] = Vector{Shape{T}}()
        push!(pos, cartesian(center(coarse_chart)))

        # It is assumed the vertices of this cell have the same index
        # mesh and its refinement.
        # coarse_cell = chart(mesh, coarse_idcs)

        # get the index in fine.vertices of the centroid of coarse_cell
        centroid = barytocart(coarse_chart, uv_ctr)
        I = CollisionDetection.find(fine_vertices, centroid)
        @assert length(I) == 1
        centroid_id = I[1]

        # get the indx in fine.vertices of the centroid of the faces of coarse_cell
        face_center_ids = Vector{Int}(undef,num_faces)
        for f in 1:num_faces

            # prepare the barycentric coordinate
            uv_face_ctr = ones(dimension(mesh)+1)/(dimension(mesh))
            uv_face_ctr[f] = 0
            uv_face_ctr = uv_face_ctr[1:end-1]

            face_ctr = barytocart(coarse_chart, uv_face_ctr)
            I = CollisionDetection.find(fine_vertices, face_ctr)
            @assert length(I) == 1
            face_center_ids[f] = I[1]
        end

        n = vton[centroid_id]
        for c in vtoc[centroid_id,1:n]
            fine_idcs = cells(refined)[c]
            local_id = something(findfirst(isequal(centroid_id), fine_idcs), 0)
            @assert local_id != 0
            shape = Shape(c, local_id, 1.0)
            push!(fns[i], shape)
        end

        for f in 1:num_faces
            v = face_center_ids[f]
            jct_pred(vertices(refined)[v]) && continue
            n = vton[v]
            for c in vtoc[v,1:n]
                fine_idcs = cells(refined)[c]
                local_id = something(findfirst(isequal(v), fine_idcs),0)
                @assert local_id != 0
                shape = Shape(c, local_id, 1/n/2)
                push!(fns[i], shape)
            end
        end

        for f in 1:length(coarse_idcs)
            v = coarse_idcs[f]
            jct_pred(vertices(refined)[v]) && continue
            n = vton[v]
            for c in vtoc[v,1:n]
                fine_idcs = cells(refined)[c]
                local_id = something(findfirst(isequal(v), fine_idcs),0)
                @assert local_id != 0
                shape = Shape(c, local_id, 1/n/2)
                push!(fns[i], shape)
            end
        end
    end

    @assert length(fns) == length(pos)

    NF = 3
    return LagrangeBasis{1,0,NF}(refined, fns, pos)
end



"""
    duallagrangec0d1(originalmesh, refinedmesh)

It is the user responsibility to provide two meshes representing the same object.
The second mesh needs to be obtained using "barycentric_refinement(originalmesh)".
This basis function creats the dual Lagrange basis function and return an object that contains array of shapes [fns]
It also return a gemoetry containing the refined mesh.
"""
function duallagrangec0d1(mesh, mesh2, pred, ::Type{Val{2}})
  T = coordtype(mesh)
  U = universedimension(mesh)
  # get the information about number of vertices, number of faces , and the maping between vertices and faces for the original mesh
  numverts1 = numvertices(mesh)
  num_cells1 = numcells(mesh)
  cellids1, ncells1=vertextocellmap(mesh)
  # get the information about number of vertices, number of faces , and the maping between vertices and faces for the refined mesh
  num_cells2 = numcells(mesh2)
  numverts2 = numvertices(mesh2)
  geometry = mesh2
  cellids2, ncells2 = vertextocellmap(mesh2)

  fns = Vector{Vector{Shape{T}}}(undef,num_cells1)
  pos = Vector{vertextype(mesh)}()
  # We will iterate over the coarse mesh segments to assign all the functions to it.
  for segment_coarse in 1 : num_cells1
    # For the dual Lagrange there is a 6 shapes per segment
    numshapes = (ncells1[segment_coarse]*4) -2
    shapes = Vector{Shape{T}}(undef,numshapes)
    # Now we will get all the smaller faces within the coarse segment
    #i.e The coose segment will have two points, and these tow points are connected to two segmesnts in the finer mesh
    # This will give us a 4 smaller faces per Dual lagrange basis, we store them first in all_faces
    all_faces= Array{SVector{2,Int}}(undef,4)                      # faces in the original segment (4)
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
    push!(pos, cartesian(center(chart(mesh, segment_coarse))))
  end

  NF = 2
  LagrangeBasis{1,0,NF}(geometry, fns, pos)
end

gradient(space::LagrangeBasis{1,0}, geo, fns) = NDLCCBasis(geo, fns)
# gradient(space::LagrangeBasis{1,0}, geo::CompScienceMeshes.AbstractMesh{U,3} where {U}, fns) = NDBasis(geo, fns)

curl(space::LagrangeBasis{1,0}, geo, fns) = RTBasis(geo, fns)

curl(space::LagrangeBasis{2,0}, geo, fns) = BDMBasis(geo, fns) 

gradient(space::LagrangeBasis{1,0,<:CompScienceMeshes.AbstractMesh{<:Any,2}}, geo, fns) =
    LagrangeBasis{0,-1,1}(geo, fns, space.pos)


gradient(space::LagrangeBasis{1,0,<:CompScienceMeshes.AbstractMesh{<:Any,3}}, geo, fns) =
    NDBasis(geo, fns, space.pos)

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

    trpos = deepcopy(positions(X))
    LagrangeBasis{1,0,NF}(geo, fns, trpos)
end



function extend_0_form(supp, dirichlet, x_prt, Lg_prt)

    Id = BEAST.Identity()

    bnd_supp = boundary(supp)
    supp_nodes = CompScienceMeshes.skeleton_fast(supp, 0)

    dir_compl = submesh(!in(dirichlet), bnd_supp)
    dir_compl_nodes = CompScienceMeshes.skeleton_fast(dir_compl, 0)
    int_nodes = submesh(!in(dir_compl_nodes), supp_nodes)

    if length(int_nodes) == 0
        T = scalartype(Id, Lg_prt, Lg_prt)
        S = typeof(Lg_prt)
        P = eltype(Lg_prt.pos)
        return T[], int_nodes, S(supp, Vector{Vector{Shape{T}}}(), Vector{P}())
    end

    Lg_int = BEAST.lagrangec0d1(supp, int_nodes)

    grad_Lg_prt = gradient(Lg_prt)
    grad_Lg_int = gradient(Lg_int)

    A = assemble(Id, grad_Lg_int, grad_Lg_int, threading=Threading{:single})
    a = -assemble(Id, grad_Lg_int, grad_Lg_prt, threading=Threading{:single}) * x_prt

    x_int = A \ a
    return x_int, int_nodes, Lg_int
end


function dualforms_init(Supp, Dir)
    tetrs = barycentric_refinement(Supp)
    v2t, v2n = CompScienceMeshes.vertextocellmap(tetrs)
    bnd = boundary(tetrs)
    gpred = CompScienceMeshes.overlap_gpredicate(Dir)
    dir = submesh((m,face) -> gpred(chart(m,face)), bnd)
    return tetrs, bnd, dir, v2t, v2n
end

function duallagrangec0d1(mesh::CompScienceMeshes.AbstractMesh{<:Any,3}, nodes, dirbnd)
    refd, bnd, dir, v2t, v2n = dualforms_init(mesh, dirbnd)
    dual0forms_body(mesh, refd, bnd, dir, v2t, v2n)
end

function dual0forms_body(mesh::CompScienceMeshes.AbstractMesh{<:Any,3}, refd, bnd, dir, v2t, v2n)

    T = coordtype(refd)
    S = Shape{T}
    V = vertextype(mesh)
    bfs = Vector{Vector{S}}(undef, length(mesh))
    pos = Vector{V}(undef, length(mesh))

    # Cells = cells(mesh)
    Cells = [c for c in mesh]
    num_threads = Threads.nthreads()
    for F in 1:length(mesh)
        # Cell = Cells[F]
        Cell = CompScienceMeshes.indices(mesh, Cells[F])

        myid = Threads.threadid()
        # myid == 1 && F % 20 == 0 &&
        #     println("Constructing dual 1-forms: $(F) out of $(length(mesh)).")

        idcs1 = v2t[Cell[1],1:v2n[Cell[1]]]
        idcs2 = v2t[Cell[2],1:v2n[Cell[2]]]
        idcs3 = v2t[Cell[3],1:v2n[Cell[3]]]

        supp1 = refd[idcs1] # 2D
        supp2 = refd[idcs2]
        supp3 = refd[idcs3]

        bnd_supp1 = boundary(supp1) # 1D
        bnd_supp2 = boundary(supp2)
        bnd_supp3 = boundary(supp3)

        supp23 = submesh(in(bnd_supp2), bnd_supp3) # 1D
        supp31 = submesh(in(bnd_supp3), bnd_supp1)
        supp12 = submesh(in(bnd_supp1), bnd_supp2)

        port = boundary(supp23) # 0D
        port = submesh(in(boundary(supp31)), port)
        port = submesh(in(boundary(supp12)), port)
        @assert length(port) == 1

        in_dir = in(dir)
        dir1_edges = submesh(in_dir, bnd_supp1) # 1D
        dir2_edges = submesh(in_dir, bnd_supp2)
        dir3_edges = submesh(in_dir, bnd_supp3)

        bnd_dir1 = boundary(dir1_edges) # 0D
        bnd_dir2 = boundary(dir2_edges)
        bnd_dir3 = boundary(dir3_edges)

        dir23_nodes = submesh(in(bnd_dir2), bnd_dir3) # 0D
        dir31_nodes = submesh(in(bnd_dir3), bnd_dir1)
        dir12_nodes = submesh(in(bnd_dir1), bnd_dir2)

        # Step 1: set port flux and extend to dual faces
        x0 = ones(T,1)

        Lg12_prt = BEAST.lagrangec0d1(supp12, port)
        Lg23_prt = BEAST.lagrangec0d1(supp23, port)
        Lg31_prt = BEAST.lagrangec0d1(supp31, port)

        x23, supp23_int_nodes, _ = extend_0_form(supp23, dir23_nodes, x0, Lg23_prt)
        x31, supp31_int_nodes, _ = extend_0_form(supp31, dir31_nodes, x0, Lg31_prt)
        x12, supp12_int_nodes, _ = extend_0_form(supp12, dir12_nodes, x0, Lg12_prt)

        port1_nodes = union(port, supp31_int_nodes, supp12_int_nodes)
        port2_nodes = union(port, supp12_int_nodes, supp23_int_nodes)
        port3_nodes = union(port, supp23_int_nodes, supp31_int_nodes)

        x1_prt = [x0; x31; x12]
        x2_prt = [x0; x12; x23]
        x3_prt = [x0; x23; x31]

        Lg1_prt = BEAST.lagrangec0d1(supp1, port1_nodes)
        Lg2_prt = BEAST.lagrangec0d1(supp2, port2_nodes)
        Lg3_prt = BEAST.lagrangec0d1(supp3, port3_nodes)

        x1_int, _, Lg1_int = extend_0_form(supp1, dir1_edges, x1_prt, Lg1_prt)
        x2_int, _, Lg2_int = extend_0_form(supp2, dir2_edges, x2_prt, Lg2_prt)
        x3_int, _, Lg3_int = extend_0_form(supp3, dir3_edges, x3_prt, Lg3_prt)

        # inject in the global space
        fn = BEAST.Shape{T}[]
        addf!(fn, x1_prt, Lg1_prt, idcs1)
        addf!(fn, x1_int, Lg1_int, idcs1)

        addf!(fn, x2_prt, Lg2_prt, idcs2)
        addf!(fn, x2_int, Lg2_int, idcs2)

        addf!(fn, x3_prt, Lg3_prt, idcs3)
        addf!(fn, x3_int, Lg3_int, idcs3)

        pos[F] = cartesian(CompScienceMeshes.center(chart(mesh, Cells[F])))
        bfs[F] = fn
    end

    LagrangeBasis{1,0,3}(refd, bfs, pos)
end



function lagrangecx(mesh::CompScienceMeshes.AbstractMesh{<:Any,3}; order)

    T = coordtype(mesh)
    NF = binomial(2+order, 2)
    P = vertextype(mesh)
    S = Shape{T}

    fns = Vector{Vector{S}}(undef, length(mesh) * NF)
    pos = Vector{P}(undef, length(mesh) * NF)

    idx = 1
    u = one(T)
    for (c,cell) in enumerate(mesh)
        ch = chart(mesh,cell)
        for r in 1:NF
            fns[idx] = S[S(c,r,u)]
            pos[idx] = cartesian(center(ch))
            idx += 1
        end
    end

    return LagrangeBasis{order,-1,NF}(mesh, fns, pos)
end



"""
    localindices(dof, chart, i)

Returns a vector of indices into the vector of local shape functions that correspond to
global degrees of freedom supported on sub-entity `i`, where the type of entity (nodes,
edge, face) is encoded in the type of 'dof'.
"""
function localindices(dof::_LagrangeGlobalNodesDoFs, chart::CompScienceMeshes.Simplex,
    localspace, i)

    d = dof.order
    [div((d+1)*(d+2),2), d+1, 1][[i]]
end

function localindices(dof::_LagrangeGlobalEdgeDoFs, chart::CompScienceMeshes.Simplex,
    localspace, i)

    d = dof.order
    lids, lid = Int[], 0
    if i == 1
        for i in 0:d for j in 0:d
            k = d - i - j
            k < 0 && continue
            lid += 1
            i != 0 && continue
            j in (0,d) && continue
            push!(lids, lid)
        end end
    elseif i == 2
        for i in 0:d for j in 0:d
            k = d - i - j
            k < 0 && continue
            lid += 1
            i in (0,d) && continue
            j != 0 && continue
            push!(lids, lid)
        end end
    elseif i ==3
        for i in 0:d for j in 0:d
            k = d - i - j
            k < 0 && continue
            lid += 1
            j in (0,d) && continue
            k != 0 && continue
            push!(lids, lid)
        end end
    else
        error("wrong local edge index")
    end
    return lids
end

function localindices(dof::_LagrangeGlobalFaceDoFs, chart::CompScienceMeshes.Simplex,
    localspace, i)

    d = dof.order
    lids, lid = Int[], 0
    for i in 0:d, j in 0:d
        k = d - i - j
        k < 0 && continue
        lid += 1
        (i == 0 || j == 0 || k == 0) && continue
        # @show (i,j,k)
        push!(lids, lid)
    end
    return lids
end

function localindices(dof::_LagrangeGlobalNodesDoFs, chart::CompScienceMeshes.Simplex,
    localspace::LagrangeRefSpace{<:Real,2}, i)
    return [i]
end

function localindices(dof::_LagrangeGlobalEdgeDoFs, chart::CompScienceMeshes.Simplex,
    localspace::LagrangeRefSpace{<:Real,2}, i)
    return [3+i]
end

function localindices(dof::_LagrangeGlobalFaceDoFs, chart::CompScienceMeshes.Simplex,
    localspace::LagrangeRefSpace{<:Real,2}, i)
    return []
end


function lagrangec0(mesh::CompScienceMeshes.AbstractMesh{<:Any,3}; order)

    T = coordtype(mesh)
    atol = sqrt(eps(T))

    P = vertextype(mesh)
    S = Shape{T}
    F = Vector{S}

    verts = skeleton(mesh, 0)
    edges = skeleton(mesh, 1)

    C02 = connectivity(mesh, verts, abs); R02 = rowvals(C02); V02 = nonzeros(C02)
    C12 = connectivity(mesh, edges, abs); R12 = rowvals(C12); V12 = nonzeros(C12)

    ne = order-1
    nf = div((order-2)*(order-1), 2)

    nV = length(verts)
    nE = length(edges) * ne
    nF = length(mesh) * nf
    N = nV + nE + nF

    localspace = LagrangeRefSpace{T,order,3,binomial(2+order,2)}()
    localdim = numfunctions(localspace)
    
    d = order
    fns = [S[] for n in 1:(nV+nE+nF)]
    pos = fill(point(0,0,0), nV+nE+nF)
    for cell in mesh
        cell_ch = chart(mesh, cell)
        V = R02[nzrange(C02,cell)]
        I = V02[nzrange(C02,cell)]
        for (i,v) in zip(I, V)
            vertex_ch = chart(verts, v)
            gids = [v]
            lids = localindices(_LagrangeGlobalNodesDoFs(d), cell_ch, localspace, i)
            v = globaldofs(vertex_ch, cell_ch, localspace, _LagrangeGlobalNodesDoFs(d))
            α = v[lids, :]
            β = inv(α')
            for i in axes(β,1)
                for j in axes(β,2)
                    isapprox(β[i,j], 0; atol) && continue
                    push!(fns[gids[i]], S(cell, lids[j], β[i,j]))
                end
            end
        end

        E = R12[nzrange(C12,cell)]
        I = V12[nzrange(C12,cell)]
        for (i,e) in zip(I, E)
            edge_ch = chart(edges, e)
            gids = nV + (e-1)*ne + 1: nV + e*ne
            lids = localindices(_LagrangeGlobalEdgeDoFs(d), cell_ch, localspace, i)
            @assert length(lids) == length(gids)
            v = globaldofs(edge_ch, cell_ch, localspace, _LagrangeGlobalEdgeDoFs(d))
            α = v[lids, :]
            β = inv(α')
            for i in axes(β,1)
                for j in axes(β,2)
                    isapprox(β[i,j], 0; atol) && continue
                    push!(fns[gids[i]], S(cell, lids[j], β[i,j]))
                end
            end
        end

        order < 3 && continue
        F = [cell]
        for (q,f) in enumerate(F)
            face_ch = chart(mesh, f)
            gids = nV + nE + (f-1)*nf + 1 : nV + nE + f*nf
            lids = localindices(_LagrangeGlobalFaceDoFs(d), cell_ch, localspace, 1)
            v = globaldofs(face_ch, cell_ch, localspace, _LagrangeGlobalFaceDoFs(d))
            α = v[lids, :]
            # @show α
            β = inv(α')
            for i in axes(β,1)
                for j in axes(β,2)
                    isapprox(β[i,j], 0; atol) && continue
                    push!(fns[gids[i]], S(cell, lids[j], β[i,j]))
                end
            end
        end 
    end

    for v in verts pos[v] = cartesian(center(chart(verts, v))) end
    for e in edges pos[nV + (e-1)*ne + 1: nV + e*ne] .= Ref(cartesian(center(chart(edges, e)))) end
    for f in mesh pos[nV + nE + (f-1)*nf + 1: nV + nE + f*nf] .= Ref(cartesian(center(chart(mesh, f)))) end

    return LagrangeBasis{order,0,localdim}(mesh, fns, pos)
end