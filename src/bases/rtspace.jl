

mutable struct RTBasis{T,M,P} <: Space{T}
  geo::M
  fns::Vector{Vector{Shape{T}}}
  pos::Vector{P}
end

RTBasis(geo, fns) = RTBasis(geo, fns, Vector{vertextype(geo)}(undef,length(fns)))

positions(rt) = rt.pos
refspace(space::RTBasis{T}) where {T} = RTRefSpace{T}()
subset(rt::RTBasis,I) = RTBasis(rt.geo, rt.fns[I], rt.pos[I])


mutable struct ValDiv end



"""
    raviartthomas(mesh, cellpairs::Array{Int,2})

Constructs the RT basis on the input `mesh`. The i-th RT basis function will
    represent a current distribution flowing from cell `cellpairs[1,i]` to
    `cellpairs[2,i]` on the mesh.

Returns an object of type `RTBasis`, which comprises both the mesh and pairs of
    Shape objects which corresponds to the cell pairs, containing the necsessary
    coefficients and indices to compute the exact basis functions when required
    by the solver.
"""
function raviartthomas(mesh::CompScienceMeshes.AbstractMesh{U,D1,T}, cellpairs::Array{Int,2}) where {U,D1,T}

    # combine now the pairs of monopolar RWGs in div-conforming RWGs
    numpairs = size(cellpairs,2)
    functions = Vector{Vector{Shape{T}}}(undef,numpairs)
    positions = Vector{vertextype(mesh)}(undef,numpairs)
    # Cells = cells(mesh)
    for i in 1:numpairs
        if cellpairs[2,i] > 0
            c1 = cellpairs[1,i]; # cell1 = Cells[c1] #mesh.faces[c1]
            c2 = cellpairs[2,i]; # cell2 = Cells[c2] #mesh.faces[c2]

            cell1 = CompScienceMeshes.indices(mesh, c1)
            cell2 = CompScienceMeshes.indices(mesh, c2)

            e1, e2 = getcommonedge(cell1, cell2)
            functions[i] = [
              Shape{T}(c1, abs(e1), T(+1.0)),
              Shape{T}(c2, abs(e2), T(-1.0))]
            isct = intersect(cell1, cell2)
            @assert length(isct) == 2
            @assert !(cell1[abs(e1)] in isct)
            @assert !(cell2[abs(e2)] in isct)

            ctr1 = cartesian(center(chart(mesh, c1)))
            ctr2 = cartesian(center(chart(mesh, c2)))
            positions[i] = (ctr1 + ctr2) / 2
        else
            c1 = cellpairs[1,i]
            e1 = cellpairs[2,i]
            functions[i] = [
            Shape(c1, abs(e1), T(+1.0))]
            positions[i] = cartesian(center(chart(mesh, c1)))
        end
    end

    geo = mesh
    RTBasis(geo, functions, positions)
end


function raviartthomas(mesh, edges::CompScienceMeshes.AbstractMesh{U,2} where {U})
    cps = CompScienceMeshes.cellpairs(mesh, edges)
    # ids = findall(x -> x>0, cps[2,:])
    raviartthomas(mesh, cps)
end

function raviartthomas(mesh::CompScienceMeshes.AbstractMesh{U,3} where {U})
    bnd = boundary(mesh)
    edges = submesh(!in(bnd), skeleton(mesh,1))
    return raviartthomas(mesh, edges)
end


"""
    raviartthomas(mesh)

Conducts pre-processing on the input `mesh` by extracting the cell edges, cell pairs
    and indices required to construct the RT basis on the `mesh`.

Calls raviartthomas(mesh::Mesh, cellpairs::Array{Int,2}), which constructs
    the RT basis on the `mesh`, using the cell pairs identified.

Returns the RT basis object.
"""
function raviartthomas(mesh; sort=:spacefillingcurve)
    edges = skeleton(mesh, 1; sort)
    cps = cellpairs(mesh, edges, dropjunctionpair=true)
    ids = findall(x -> x>0, cps[2,:])
    raviartthomas(mesh, cps[:,ids])
end



# """
#     raviartthomas(Γ; neumann)

# Constructs the RT space relative to boundary `neumann` of an open surface, only
#     selecting cell pairs whose common edge does lie on `γ` . (This prevents
#     the calculation of physically-impossible surface currents, such as those
#     flowing 'off the edge' of a surface.)
# """
# function raviartthomas(Γ; neumann)

#     γ = neumann

#     in_interior = interior_tpredicate(Γ)
#     on_junction = overlap_gpredicate(γ)

#     pred = c -> (in_interior(c) || on_junction(chart(Γ,c)))

#     edges = skeleton(pred, Γ, 1)
#     cps = cellpairs(Γ, edges, dropjunctionpair=true)

#     raviartthomas(Γ, cps)
# end

# raowiltonglisson = raviartthomas


"""
    portcells(Γ::Mesh, γ::Mesh)
returns an array containing cell pairs of mesh Γ
around a boundary edge that overlaps with mesh γ
"""
function portcells(Γ, γ)
  in_interior = interior_tpredicate(Γ)

  overlaps = overlap_gpredicate(γ)
  on_junction = (m,c) -> overlaps(chart(m,c))

  pred = (m,x) -> (!in_interior(m,x) && on_junction(m,x)) #check only for exterior overlapping edges
  # - if γ is defined to overlap within the structure Γ, this isn't considered a port -
  #   edges = skeleton(pred, Γ, 1) #Take only exterior edges overlapping γ segment
  edges = submesh(pred, skeleton(Γ,1))
  cps = cellpairs(Γ, edges, dropjunctionpair=true)
  return cps
end

"""
    rt_cedge(cps::Array{Int,2}, weight)
Computes single basis function with equally distributed constant current
leaving or entering port defined by cellpairs cps.  weight defines the
total current over the port and its direction (+ve = out, -ve = in)
"""
function rt_cedge(cps::Array{Int,2}, weight)
    T=typeof(weight)
  numpairs = size(cps,2)
  @assert numpairs > 0
  weight = weight / numpairs #total current leaving and entering equal 1
  functions =  Vector{Shape{T}}(undef,numpairs) #note: not a Vector{Vector}
    for i in 1:numpairs
        c1 = cps[1,i]
        e1 = cps[2,i]
        functions[i] =
        Shape(c1, abs(e1), weight)
        # this assumes cellpairs will be adjacent eachother through the loop
      end
  functions
end

"""
    rt_vedge(cps::Array{Int,2}, weight)
Computes n-1 basis function with oscillating current in(leaving and entering)
pairs of half triangles defined over port specified by cellpairs cps.
weight defines the magnitude of individual current in and out the half triangles,
and it's polarity simply defines whether to start with in or out
"""
function rt_vedge(cps::Array{Int,2}, weight)
  T=typeof(weight)
  numpairs = size(cps,2)
  @assert numpairs > 0
  #adjacent cells are considered a pair, so we have one less numpairs
  functions =  Vector{Vector{Shape{T}}}(undef,numpairs - 1)
    for i in 1:numpairs
      if i < numpairs # stop when on last cellpair
        c1 = cps[1,i]; e1 = cps[2,i]
        c2 = cps[1,i+1]; e2 = cps[2,i+1]
        functions[i] = [
          Shape(c1, abs(e1), weight),
          Shape(c2, abs(e2), -weight)]
        # this assumes cellpairs will be adjacent eachother through the loop
      end
    end
 functions
end

"""
    rt_ports(Γ::Mesh, γ::Mesh ...)
Constructs the RT space on `Γ`, relative to boundary pairs in `γ`. `γ` expects
any number of pairs-of-ports as arguments and accepts tuples, arrays,
vectors etc. e.g `rt_ports(Γ, a, b ...);` where a = [γ₁ γ₂], b = (γ₃,γ₄) etc.
`rt_ports` with no pair of ports supplied i.e `rt_ports(Γ)` reduces to the
`raviartthomas(Γ)` function. The RT space ensures current continuity in each
pair of ports. i.e. current leaving mesh Γ through γ₁ is accounted for in γ₂.

Returns the RT basis object.
"""
function rt_ports(Γ, γ...)
    #TODO: can include an extra argument for polarity applied to ports

    T = coordtype(Γ)

    internal = raviartthomas(Γ)

    fns = internal.fns
    pos = internal.pos

    for i = 1:length(γ) #loop through pairs of ports

        port1 = portcells(Γ, γ[i][1])
        port2 = portcells(Γ, γ[i][2])

        ce1 = rt_cedge(port1, T(+1.0))
        ce2 = rt_cedge(port2, T(-1.0))

        ffs = Vector{Shape{T}}(undef,length(ce1) + length(ce2))
        ffs = [ce1;ce2]

        ve1 = rt_vedge(port1, T(+1.0))
        ve2 = rt_vedge(port2, T(-1.0))

        fns = [fns;[ffs];ve1;ve2]
    end

    RTBasis(Γ, fns)

end

"""
    getindex_rtg(RT::RTBasis)
Returns the indices of the global half RWGs present in `RT`.
`RT` is typically gotten from `rt_ports`
"""
function getindex_rtg(RT::RTBasis)
    #TODO: use more reasonable condition to extract indices
    # or have rt_ports return indices with RTBasis
    idx=[]
    A = RT.fns
    for i in 1:length(A)
        if size(A[i])[1] > 2
            push!(idx,i)
        end
    end
    idx
end

divergence(X::RTBasis, geo, fns) = LagrangeBasis{0,-1,1}(geo, fns, deepcopy(positions(X)))
ntrace(X::RTBasis, geo, fns) = LagrangeBasis{0,-1,1}(geo, fns, deepcopy(positions(X)))


function LinearAlgebra.cross(::NormalVector, s::RTBasis)
    @assert CompScienceMeshes.isoriented(s.geo)
    fns = similar(s.fns)
    for (i,fn) in pairs(s.fns)
        fns[i] = [Shape(sh.cellid, sh.refid, sh.coeff) for sh in fn]
    end
    NDBasis(s.geo, fns, s.pos)
end
