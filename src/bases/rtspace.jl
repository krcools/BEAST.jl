export raviartthomas, raowiltonglisson
export portcells, rt_ports

type RTRefSpace{T} <: RefSpace{T,3} end

type RTBasis{T,M} <: Space{T}
  geo::M
  fns::Vector{Vector{Shape{T}}}
end

refspace{T}(space::RTBasis{T}) = RTRefSpace{T}()
valuetype{T}(ref::RTRefSpace{T}, charttype) = SVector{3,Tuple{Pt{universedimension(charttype),T},T}}

type ValDiv end



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
function raviartthomas(mesh, cellpairs::Array{Int,2})

    # combine now the pairs of monopolar RWGs in div-conforming RWGs
    numpairs = size(cellpairs,2)
    functions = Vector{Vector{Shape{Float64}}}(numpairs)
    for i in 1:numpairs
        if cellpairs[2,i] > 0
            c1 = cellpairs[1,i]; cell1 = mesh.faces[c1]
            c2 = cellpairs[2,i]; cell2 = mesh.faces[c2]
            e1, e2 = getcommonedge(cell1, cell2)
            functions[i] = [
              Shape(c1, abs(e1), +1.0),
              Shape(c2, abs(e2), -1.0)]
            isct = intersect(cell1, cell2)
            @assert length(isct) == 2
            @assert !(cell1[abs(e1)] in isct)
            @assert !(cell2[abs(e2)] in isct)
        else
            c1 = cellpairs[1,i]
            e1 = cellpairs[2,i]
            functions[i] = [
            Shape(c1, abs(e1), +1.0)]
        end
    end

    geo = mesh
    RTBasis(geo, functions)
end



"""
    raviartthomas(mesh)

Conducts pre-processing on the input `mesh` by extracting the cell edges, cell pairs
    and indices required to construct the RT basis on the `mesh`.

Calls raviartthomas(mesh::Mesh, cellpairs::Array{Int,2}), which constructs
    the RT basis on the `mesh`, using the cell pairs identified.

Returns the RT basis object.
"""
function raviartthomas(mesh)
    edges = skeleton(mesh, 1)
    cps = cellpairs(mesh, edges, dropjunctionpair=true)
    ids = find(x -> x>0, cps[2,:])
    raviartthomas(mesh, cps[:,ids])
end



"""
    raviartthomas(Γ, γ)

Constructs the RT space relative to boundary `γ` of an open surface, only
    selecting cell pairs whose common edge does not lie on `γ` . (This prevents
    the calculation of physically-impossible surface currents, such as those
    flowing 'off the edge' of a surface.)

Calls raviartthomas(Γ, duals) which constructs
    the RT basis on the `mesh`, using the cell pairs identified.
"""
function raviartthomas(Γ, γ)

    in_interior = interior_tpredicate(Γ)

    overlaps = overlap_gpredicate(γ)
    on_junction = c -> overlaps(simplex(vertices(Γ,c)))

    pred = x -> (in_interior(x) || on_junction(x))
    edges = skeleton(pred, Γ, 1)

    duals = cellpairs(Γ, edges, dropjunctionpair=true)
    raviartthomas(Γ, duals)
end

raowiltonglisson = raviartthomas


"""
    portcells(Γ::Mesh, γ::Mesh)
returns an array containing cell pairs of mesh Γ
around a boundary edge that overlaps with mesh γ
"""
function portcells(Γ, γ)
  in_interior = interior_tpredicate(Γ)

  overlaps = overlap_gpredicate(γ)
  on_junction = c -> overlaps(simplex(vertices(Γ,c)))

  pred = x -> (!in_interior(x) && on_junction(x)) #check only for exterior overlapping edges
  # - if γ is defined to overlap within the structure Γ, this isn't considered a port -
  edges = skeleton(pred, Γ, 1) #Take only exterior edges overlapping γ segment
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
  numpairs = size(cps,2)
  @assert numpairs > 0
  weight = weight / numpairs #total current leaving and entering equal 1
  functions =  Vector{Shape{Float64}}(numpairs) #note: not a Vector{Vector}
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
  numpairs = size(cps,2)
  @assert numpairs > 0
  #adjacent cells are considered a pair, so we have one less numpairs
  functions =  Vector{Vector{Shape{Float64}}}(numpairs - 1)
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
    rt_ports(Γ::Mesh, γ₁::Mesh, γ₂::Mesh)
Construct the RT space relative to boundary γ₁ and γ₂,
ensuring current leaving mesh Γ through γ₁ is accounted for
in γ₂ and vice versa
"""
function rt_ports(Γ, γ₁, γ₂)
  internal = raviartthomas(Γ)

  port1 = portcells(Γ, γ₁); port2 = portcells(Γ, γ₂)
  ce1 = rt_cedge(port1, 1.0); ce2 = rt_cedge(port2, -1.0)
  ffs = Vector{Shape{Float64}}(length(ce1) + length(ce2))
  ffs = [ce1;ce2]

  ve1 = rt_vedge(port1, 1.0); ve2 = rt_vedge(port2, -1.0)
  fns = [[ffs];internal.fns;ve1;ve2]
  RTBasis(Γ, fns)
end

divergence(ref::RTRefSpace, sh, el) = Shape(sh.cellid, 1, sh.coeff/volume(el))
divergence(sp::RTBasis, geo, fns) = LagrangeBasis{0,-1,1}(geo, fns)


function (ϕ::RTRefSpace)(mp)

    u, v = parametric(mp)
    j = jacobian(mp)

    tu = tangents(mp,1)
    tv = tangents(mp,2)

    d = 2/j

    return SVector((
        ((tu*(u-1) + tv*v    ) / j, d),
        ((tu*u     + tv*(v-1)) / j, d),
        ((tu*u     + tv*v    ) / j, d)
    ))
end



function restrict{T}(ϕ::RTRefSpace{T}, dom1, dom2)

    K = numfunctions(ϕ)
    D = dimension(dom1)

    @assert K == 3
    @assert D == 2
    @assert D == dimension(dom2)

    Q = zeros(T,K,K)
    for i in 1:K

        # find the center of edge i of dom2
        a = dom2.vertices[mod1(i+1,D+1)]
        b = dom2.vertices[mod1(i+2,D+1)]
        c = (a + b) / 2

        # find the outer binormal there
        t = b - a
        l = norm(t)
        n = dom2.normals[1]
        m = cross(t, n) / l

        u = carttobary(dom1, c)
        x = neighborhood(dom1, u)

        y = ϕ(x)

        for j in 1:K
            Q[j,i] = dot(y[j][1], m) * l
        end
    end

    return Q
end


ntrace(X::RTBasis, geo, fns) = LagrangeBasis{0,-1,1}(geo, fns)

"""
    ntrace(refspace, element, localindex, face)

Compute the normal trace of all local shape functions on `elements` belonging to
`refspace` on `face`. This function returns a matrix expressing the traces of local
shape functions in `refspace` as linear combinations of functions in the local
trace space. Cf. `restrict`. `localindex` is the index of `face` in the enumeration
of faces of `elements`. In many special cases knowing this index allows for highly
optimised implementations.
"""
function ntrace(x::RTRefSpace, el, q, fc)
    t = zeros(scalartype(x),1,3)
    t[q] = 1 / volume(fc)
    return t
end
