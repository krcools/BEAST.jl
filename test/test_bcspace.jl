using Base.Test

using CompScienceMeshes
using BEAST

"""
    isdivconforming(space)

Returns for every basis functions in the space the total of the absolute
values of non-compensated skeleton flux, i.e. flux flowing into from one
side of an edge and not flowing out on the other side. Should be all zeros
for divergence conforming functions, at least on a closed surface.
"""
function isdivconforming(space)

    geo = geometry(space)
    mesh = geo

    edges = skeleton(mesh,1)
    D = connectivity(edges, mesh, abs)
    rows = rowvals(D)
    vals = nonzeros(D)

    Flux = zeros(Float64, numcells(mesh),3)
    TotalFlux = zeros(Float64, numcells(edges), numfunctions(space))

    for i in 1 : numfunctions(space)

        fill!(Flux, 0)
        bfs = BEAST.basisfunction(space,i)
        for bf in bfs
            c = bf.cellid
            e = bf.refid
            x = bf.coeff
            Flux[c,e] += x
        end

        for E in 1 : numcells(edges)
            for j in nzrange(D,E)
                F = rows[j]
                e = vals[j]
                @assert 1 <= e <= 3
                TotalFlux[E,i] += Flux[F,e]
            end
        end
    end

    return TotalFlux
end

function touches_predicate(mesh)

    verts = skeleton(mesh,0)
    verts = sort(verts.faces, lt=lexless)

    edges = skeleton(mesh,1)
    edges = sort(edges.faces, lt=lexless)

    function pred(s)

        num_hits = 0
        for i in 1:length(s)
            if !isempty(searchsorted(verts, SVector(s[i]), lt=lexless))
                num_hits += 1
            end
        end

        @assert num_hits < 3
        num_hits == 0 && return false
        num_hits == 1 && return true
        isempty(searchsorted(edges, s, lt=lexless)) && return true
        return false
    end

    return pred
end

function interior(mesh::Mesh)
    D = dimension(mesh)
    edges = skeleton(mesh, D-1)
    pred = interior_tpredicate(mesh)
    submesh(pred, edges)
end

meshfile = Pkg.dir("BEAST","test","sphere2.in")
mesh = readmesh(meshfile)
@test numvertices(mesh) == 160
@test numcells(mesh) == 316

rt = raviartthomas(mesh)
@test numfunctions(rt) == 316 * 3 / 2

fine = barycentric_refinement(mesh)
edges = skeleton(mesh, 1)

bc = buffachristiansen(mesh)
@test numfunctions(bc) == 316 * 3 / 2

lc = isdivconforming(rt)
@test maximum(lc) < eps(Float64) * 1000
println("RT space is div-conforming")

lc = isdivconforming(bc);
@test maximum(lc) < eps(Float64) * 1000
println("BC space is div-conforming")


# Now repeat the exercise with an open mesh
mesh = meshrectangle(1.0, 1.0, 0.2);
fine = barycentric_refinement(mesh);

rt = raviartthomas(mesh)
bc = buffachristiansen(mesh)

@test numfunctions(rt) == 65
@test numfunctions(bc) == 65

int_pred = interior_tpredicate(mesh)
bnd_pred(s) = !int_pred(s)

leaky_edges = find(sum(abs.(isdivconforming(rt)),1))
@test length(leaky_edges) == 0

bnd = boundary(mesh)
bndtch_pred = touches_predicate(bnd)
edges = interior(mesh)
bndtch_edges = find(bndtch_pred, edges.faces)
leaky_edges = find(sum(abs.(isdivconforming(bc)),1))
@test bndtch_edges == leaky_edges


## Test the charge of BC functions
meshfile = Pkg.dir("BEAST","test","sphere2.in")
mesh = readmesh(meshfile)
bc = buffachristiansen(mesh)
fine = geometry(bc)
charges = zeros(numcells(fine))

for fn in bc.fns
    abs_charge = 0.0
    net_charge = 0.0
    fill!(charges,0)
    for  sh in fn
        cellid = sh.cellid

        #cell = simplex(vertices(fine, fine.faces[cellid]))
        net_charge += sh.coeff
        charges[cellid] += sh.coeff
    end
    abs_charge = sum(abs.(charges))
    @test net_charge + 1 ≈ 1
    @test abs_charge ≈ 2
end
