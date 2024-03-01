

mutable struct RT2Basis{T,M,P} <: Space{T}
  geo::M
  fns::Vector{Vector{Shape{T}}}
  pos::Vector{P}
end

RT2Basis(geo, fns) = RT2Basis(geo, fns, Vector{vertextype(geo)}(undef,length(fns)))

#= positions(rt) = rt.pos =#
refspace(space::RT2Basis{T}) where {T} = RT2RefSpace{T}()
subset(rt::RT2Basis,I) = RT2Basis(rt.geo, rt.fns[I], rt.pos[I])

#= mutable struct ValDiv end =#



"""
    raviartthomas2(mesh, cellpairs::Array{Int,2})

Constructs the RT2 basis on the input `mesh`. The i-th RT2 basis function will
    represent a current distribution flowing from cell `cellpairs[1,i]` to
    `cellpairs[2,i]` on the mesh.

Returns an object of type `RT2Basis`, which comprises both the mesh and pairs of
    Shape objects which corresponds to the cell pairs, containing the necsessary
    coefficients and indices to compute the exact basis functions when required
    by the solver.
"""
function raviartthomas2(mesh::CompScienceMeshes.AbstractMesh{U,D1,T}, cellpairs::Array{Int,2}) where {U,D1,T}

    # combine now the pairs of monopolar RWGs in div-conforming RWGs
    numpairs = size(cellpairs,2)
    Cells = cells(mesh)
    numcells = size(Cells,1)
    functions = Vector{Vector{Shape{T}}}(undef,2*(numpairs+numcells))
    positions = Vector{vertextype(mesh)}(undef,2*(numpairs+numcells))
    chooseref1(x) = if x<0 return 0 else return -1 end
    chooseref2(x) = if x<0 return -1 else return 0 end
    for i in 1:numpairs
        if cellpairs[2,i] > 0
            c1 = cellpairs[1,i]; # cell1 = Cells[c1] #mesh.faces[c1]
            c2 = cellpairs[2,i]; # cell2 = Cells[c2] #mesh.faces[c2]

            cell1 = CompScienceMeshes.indices(mesh, c1)
            cell2 = CompScienceMeshes.indices(mesh, c2)

            e1, e2 = getcommonedge(cell1, cell2)
            functions[2*i-1] = [
              Shape{T}(c1, 2*abs(e1)+chooseref1(sign(e1)), T(+1.0)),
              Shape{T}(c2, 2*abs(e2)+chooseref1(sign(e2)), T(-1.0))]
            functions[2*i] = [
              Shape{T}(c1, 2*abs(e1)+chooseref2(sign(e1)), T(+1.0)),
              Shape{T}(c2, 2*abs(e2)+chooseref2(sign(e2)), T(-1.0))]  
            isct = intersect(cell1, cell2)
            @assert length(isct) == 2
            @assert !(cell1[abs(e1)] in isct)
            @assert !(cell2[abs(e2)] in isct)

            ctr1 = cartesian(center(chart(mesh, c1)))
            ctr2 = cartesian(center(chart(mesh, c2)))
            positions[2*i-1] = (ctr1 + ctr2) / 2
            positions[2*i] = (ctr1 + ctr2) / 2
        else
            c1 = cellpairs[1,i]
            e1 = cellpairs[2,i]
            functions[i] = [
            Shape(c1, abs(e1), T(+1.0))]
            positions[i] = cartesian(center(chart(mesh, c1)))
        end
    end

    for (i,cell) in enumerate(Cells)
        functions[2*numpairs+2*i-1] = [
              Shape{T}(i, 7, T(+1.0))]
        
        functions[2*numpairs+2*i] = [
              Shape{T}(i, 8, T(+1.0))]

        ctr1 = cartesian(center(chart(mesh, i)))
        positions[2*numpairs+2*i-1] = ctr1
        positions[2*numpairs+2*i] = ctr1
    end

    geo = mesh
    RT2Basis(geo, functions, positions)
end


function raviartthomas2(mesh, edges::CompScienceMeshes.AbstractMesh{U,2} where {U})
    cps = CompScienceMeshes.cellpairs(mesh, edges)
    # ids = findall(x -> x>0, cps[2,:])
    raviartthomas2(mesh, cps)
end

function raviartthomas2(mesh::CompScienceMeshes.AbstractMesh{U,3} where {U})
    bnd = boundary(mesh)
    edges = submesh(!in(bnd), skeleton(mesh,1))
    return raviartthomas2(mesh, edges)
end


"""
    raviartthomas2(mesh)

Conducts pre-processing on the input `mesh` by extracting the cell edges, cell pairs
    and indices required to construct the RT2 basis on the `mesh`.

Calls raviartthomas2(mesh::Mesh, cellpairs::Array{Int,2}), which constructs
    the RT2 basis on the `mesh`, using the cell pairs identified.

Returns the RT2 basis object.
"""
function raviartthomas2(mesh; sort=:spacefillingcurve)
    edges = skeleton(mesh, 1; sort)
    cps = cellpairs(mesh, edges, dropjunctionpair=true)
    ids = findall(x -> x>0, cps[2,:])
    raviartthomas2(mesh, cps[:,ids])
end

#= divergence(X::RT2Basis, geo, fns) = LagrangeBasis{0,-1,1}(geo, fns, deepcopy(positions(X)))
ntrace(X::RT2Basis, geo, fns) = LagrangeBasis{0,-1,1}(geo, fns, deepcopy(positions(X))) =#


#= function LinearAlgebra.cross(::NormalVector, s::RT2Basis)
    @assert CompScienceMeshes.isoriented(s.geo)
    fns = similar(s.fns)
    for (i,fn) in pairs(s.fns)
        fns[i] = [Shape(sh.cellid, sh.refid, -sh.coeff) for sh in fn]
    end
    ND2Basis(s.geo, fns, s.pos)
end =#
