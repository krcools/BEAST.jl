

abstract type AbstractSpace end
abstract type Space{T} <: AbstractSpace end

Base.length(s::AbstractSpace) = numfunctions(s)
Base.in(x, s::AbstractSpace) = (x => s)

"""
    scalartype(s)

The scalar field over which the argument to a basis function or operator's integration
kernel are defined. This is always a salar data type, even if the function or kernel
is defined over a multi-dimensional space.
"""
scalartype(s::RefSpace{T}) where {T} = T
scalartype(s::Space{T}) where {T} = T

"""
    refspace(basis)

Returns the ReferenceSpace of local shape functions on which the basis is built.
"""
function refspace end


# """
#     numfunctions(r)

# Return the number of functions in a `Space` or `RefSpace`.
# """
# numfunctions(rs::RefSpace{T,D}) where {T,D} = D


"""
    scalartype(x)

The scalar field over which the values of a global or local basis function, or an
operator are defined. This should always be a scalar type, even if the basis or
operator takes on values in a vector or tensor space. This data type is used to
determine the `eltype` of assembled discrete operators.
"""
function scalartype end
scalartype(x1, xs...) = Base.promote_type(scalartype(x1), scalartype(xs...))
scalartype(x1::Number) = typeof(x1)
scalartype(T::Type, x) = Base.promote_type(T, scalartype(x))


"""
    geometry(basis)

Returns an iterable collection of geometric elements on which the functions
in `basis` are defined. The order the elements are encountered needs correspond
to the element indices used in the data structure returned by `assemblydata`.
"""
geometry(s::Space) = s.geo

"""
    basisfunction(basis, i)

Returns a vector of the shape functions defining the `i`th function
of the `basis`.
"""
basisfunction(s::Space, i) = s.fns[i]

"""
    numfunctions(basis)

Number of functions in the basis.
"""
numfunctions(space::Space) = length(space.fns)

struct DirectProductSpace{T,S<:AbstractSpace} <: AbstractSpace
    factors::Vector{S}
end

function DirectProductSpace(factors::Vector{S}) where {S<:AbstractSpace}
    @assert !isempty(factors)
    T = scalartype(factors...)
    return DirectProductSpace{T,S}(factors)
end

Base.getindex(dps::DirectProductSpace, i) = dps.factors[i]

defaultquadstrat(op, tfs::DirectProductSpace, bfs::DirectProductSpace) = defaultquadstrat(op, tfs.factors[1], bfs.factors[1])
defaultquadstrat(op, tfs::Space, bfs::DirectProductSpace) = defaultquadstrat(op, tfs, bfs.factors[1])
defaultquadstrat(op, tfs::DirectProductSpace, bfs::Space) = defaultquadstrat(op, tfs.factors[1], bfs)

# defaultquadstrat(op, tfs::DirectProductSpace, bfs::DirectProductSpace) = defaultquadstrat(op, tfs.factors[1], bfs.factors[1])
# defaultquadstrat(op, tfs::RefSpace, bfs::DirectProductSpace) = defaultquadstrat(op, tfs, bfs.factors[1])
# defaultquadstrat(op, tfs::DirectProductSpace, bfs::RefSpace) = defaultquadstrat(op, tfs.factors[1], bfs)
# scalartype(sp::DirectProductSpace{T}) where {T} = T

# export cross, ×
export ×

function Base.:+(x::AbstractSpace...)
    T = scalartype(x...)
    return DirectProductSpace{T, AbstractSpace}([x...])
end

cross(a::Space{T}, b::Space{T}) where {T} = DirectProductSpace{T,Space{T}}(Space{T}[a,b])
cross(a::DirectProductSpace{T}, b::Space{T}) where {T} = DirectProductSpace{T,Space{T}}([a.factors; b])
numfunctions(S::DirectProductSpace) = sum([numfunctions(s) for s in S.factors])
Base.length(S::DirectProductSpace) = numfunctions(S)
scalartype(s::DirectProductSpace{T}) where {T} = T
geometry(x::DirectProductSpace) = weld(x.geo...)


AbstractTrees.children(x::AbstractSpace) = ()
AbstractTrees.children(x::DirectProductSpace) = x.factors

Base.iterate(x::DirectProductSpace) = iterate(x.factors)
Base.iterate(x::DirectProductSpace, state) = iterate(x.factors, state)

struct Shape{T}
  cellid::Int
  refid::Int
  coeff::T
end

struct AssemblyData{T}
    data::Array{Tuple{Int,T},3}
end

Base.getindex(ad::AssemblyData, c, r) = ADIterator(c,r,size(ad.data,1),ad)

struct AssemblyDataEl{T}
    element_index::Int
    ad::AssemblyData{T}
end
Base.length(adp::AssemblyDataEl) = size(adp.ad.data,2)
Base.getindex(ad::AssemblyData, c) = AssemblyDataEl(c,ad)
function Base.getindex(adp::AssemblyDataEl, r)
    ad = adp.ad
    ADIterator(adp.element_index,r,size(ad.data,1),ad)
end

struct ADIterator{T}
    c::Int
    r::Int
    I::Int
    ad::AssemblyData{T}
end

function Base.iterate(it::ADIterator, i = 1)
    (it.I < i || it.ad.data[i,it.r,it.c][1] < 1) && return nothing
    (it.ad.data[i,it.r,it.c], i+1)
end


function add!(bf::Vector{Shape{T}}, cellid, refid, coeff) where T
    for (i,sh) in pairs(bf)
        if sh.cellid == cellid && sh.refid == refid
            bf[i] = Shape(cellid, refid, sh.coeff + T(coeff))
            return nothing
        end
    end
    push!(bf, Shape(cellid, refid, T(coeff)))
    return nothing
end


"""
    charts, admap, act_to_global = assemblydata(basis; onlyactives=true)

Given a `basis` this function returns a data structure containing the information
required for matrix assemble, that is, the vector `charts` containing `Simplex` elements,
a variable `admap` of type `AssemblyData`, and a mapping from indices of actively
used simplices to global simplices.

When `onlyactives` is `true`, another layer of indices is introduced to filter out all
cells of the mesh that are not in the union of the support of the basis functions
(i.e., when the basis functions are defined only on a part of the mesh).

`admap` is, in essence, a three-dimensional array of named tuples, which,
by wrapping it in the struct `AssemblyData`, allows the definition of iterators.
The tuple consists of the two entries

```
admap[i,r,c].globalindex
admap[i,r,c].coefficient
```

Here, `c` and `r` are indices in the iterable set of (active) simplices and the
set of shape functions on each cell/simplex: `r` ranges from 1 to the number of
shape functions on a cell/simplex, `c` ranges from 1 to the number of active
simplices, and `i` ranges from 1 to the number of maximal number of basis functions,
where any of the shape functions contributes to. 

For example, for continuous piecewise linear lagrange functions (c0d1), each of the
three shape functions on a triangle are associated with exactly one Lagrange function,
and therefore `i` is limited to 1.

*Note*: When `onlyactives=false`, the indices `c` correspond to
the position of the corresponding cell/simplex whilst iterating over `geometry(basis)`.
When `onlyactives=true`, then `act_to_global(c)` correspond to the position of the
corresponding cell/simplex whilst iterating over `geometry(basis)`.

For a triplet `(i,r,c)`, `globalindex` is the index in the `basis` of the
`i`th basis function that has a contribution from shape function `r` on
(active) cell/simplex `c`. `coefficient` is the coefficient of that contribution in the
linear combination defining that basis function in terms of shape
function.
"""
function assemblydata(basis::Space; onlyactives=true)

    @assert numfunctions(basis) != 0

    T = scalartype(basis)

    geo = geometry(basis)
    num_cells = numcells(geo)

    num_bfs  = numfunctions(basis)
    
    ch = chart(geo, first(geo))
    dom = domain(ch)
    num_refs = numfunctions(refspace(basis), dom)

    # Determine the number of functions defined over a given cell
    # and per local shape function.
    celltonum = make_celltonum(num_cells, num_refs, num_bfs, basis)

    # In general, a basis function space might only be defined
    # over a small portion of the underlying mesh. To avoid
    # the inefficient iterating of cells, which are not in the support of
    # any of the basis functions, we filter out only those cells, over 
    # which at least one basis function is defined.
    if onlyactives
        active, index_among_actives, num_active_cells, act_to_global =
            index_actives(num_cells, celltonum)
    else
        active = trues(num_cells)
        num_active_cells = num_cells
        index_among_actives = collect(1:num_cells)
        act_to_global = collect(1:num_cells)
    end

    num_active_cells == 0 && return nothing

    # Generate the a vector of Simplexes associated
    # with the active cells only
    elements = instantiate_charts(geo, num_active_cells, active)

    # Determine the maximal number of functions associated with a
    # local shape function
    max_celltonum = maximum(celltonum)
    fill!(celltonum, 0)
    data = fill((0,zero(T)), max_celltonum, num_refs, num_active_cells)
    for b in 1 : num_bfs
        for shape in basisfunction(basis, b)
            c = shape.cellid
            l = index_among_actives[c]
            @assert 0 < l <= num_active_cells
            r = shape.refid
            w = shape.coeff
            k = (celltonum[c,r] += 1)
            data[k,r,l] = (b,w)
        end
    end

    return elements, AssemblyData(data), act_to_global
end

"""
    celltonum = make_celltonum(num_cells, num_refs, num_bfs, basis)

Computes the array `celltonum[c,r]` that, given the index
`c` of a cell of the mesh associated with the `basis` and the index `r`
of the shape function in the refspace associated with `basis`, provides the number
of basis functions where the `r`th shape function on cell `c` is used in their definition.
"""
function make_celltonum(num_cells, num_refs, num_bfs, basis)
    celltonum = zeros(Int, num_cells, num_refs)
    for b in 1 : num_bfs
        basisfunc = basisfunction(basis, b)
        for shape in basisfunc
            c = shape.cellid
            r = shape.refid
            celltonum[c,r] += 1
        end
    end
    return celltonum
end

"""
    active, index_among_actives, num_active_cells, act_to_global =
        index_actives(num_cells, celltonum)

Given the number of cells of the mesh `num_cells` and an array indicating, which
shape functions of each cell are used in the definition of the basis, `celltonum`,
`index_actives(num_cells, celltonum)` computes a Boolean vector `active` that indicates for
each cell index if the cell is in the support of any basis function. For this reduced set
of active cells, a new index set is introduced, where `index_among_actives` is a vector
providing the mapping from the cells of the mesh to the indices of the active cells,
`num_active_cells` provides the number of all active cells, and `act_to_global` is a
vector providing the mapping from the indices of the active cells to the indices of
cells of the underlying mesh of the basis function spce.
"""
function index_actives(num_cells, celltonum)
    @assert num_cells == size(celltonum,1)
    active = falses(num_cells)
    index_among_actives = fill(0, num_cells)
    act_to_global = fill(0, num_cells)
    num_active_cells = 0
    for i in 1:num_cells
        if maximum(@view celltonum[i,:]) > 0
            num_active_cells += 1
            active[i] = true
            index_among_actives[i] = num_active_cells
            act_to_global[num_active_cells] = i
        end
    end
    resize!(act_to_global, num_active_cells)
    return active, index_among_actives, num_active_cells, act_to_global
end

"""
    instantiate_charts(geo, num_active_cells, active)

Returns a vector of all actively used simplices.
"""
function instantiate_charts(geo, num_active_cells, active)
    @assert length(geo) != 0
    E = typeof(chart(geo, first(geo)))
    elements = Vector{E}(undef,num_active_cells)
    j = 1
    for (i,p) in enumerate(geo)
        active[i] || continue
        elements[j] = chart(geo, p)
        j += 1
    end
    return elements
end


using SparseArrays
function Base.:*(space::Space, A::SparseMatrixCSC)
    fns = similar(space.fns, size(A,2))
    pos = similar(space.pos, size(A,2))
    rows = rowvals(A)
    vals = nonzeros(A)
    for i in axes(A,2)
        fns[i] = Shape{scalartype(space)}[]
        for k in nzrange(A,i)
            m = rows[k]
            c = vals[k]
            for (s,sh) in enumerate(space.fns[m])
                add!(fns[i], sh.cellid, sh.refid, sh.coeff * c)
            end
            pos[i] += space.pos[m]
        end
        pos[i] /= length(nzrange(A,i))
    end
    return (typeof(space))(space.geo, fns, pos)
end


function addf!(fn::Vector{<:Shape}, x::Vector, space::Space, idcs::Vector{Int})
    for (m,bf) in enumerate(space.fns)
        for sh in bf
            cellid = idcs[sh.cellid]
            BEAST.add!(fn, cellid, sh.refid, sh.coeff * x[m])
        end
    end
end

function support(s::BEAST.AbstractSpace, index::Int)
    geo = geometry(s)
    s1 = subset(s,[index])
    charts, ad, activecells = BEAST.assemblydata(s1)
    return geo[activecells]
end

function functionvals(s::BEAST.Space, index::Int, n=3)

    s1 = subset(s,[index])
    charts, ad, a2g = BEAST.assemblydata(s1)
    support = geometry(s)[a2g]
    
    vals = Any[]
    ctrs = Any[]
    refs = refspace(s)
    for (p,ch) in enumerate(charts)
        for i in 1:n-1
            for j in 1:n-1
                i+j < n || continue
                val = zero(BEAST.valuetype(refs, typeof(ch)))
                ct = CompScienceMeshes.neighborhood(ch, (i/n,j/n))
                fx = refs(ct)
                for r in eachindex(fx)
                    for (m,a) in ad[p,r]
                        @assert m == 1
                        val += a * fx[r].value
                    end
                end
                push!(ctrs, cartesian(ct))
                push!(vals, val)
            end
        end
    end

    return ctrs, vals
end



function eval(s::BEAST.Space, i::Int, cellid, u)

    ch = chart(cellid, s.geo)
    p = neighborhood(ch, u)

    ϕ = refpsace(s)
    ϕp = ϕ(p)

    U = valuetype(ϕ, ch)

    r = zero(U)
    for sh in s.fns[i]
        sh.cellid == cellid || continue
        r += sh.coeff * ϕp[sh.refid].value
    end

    return r
end