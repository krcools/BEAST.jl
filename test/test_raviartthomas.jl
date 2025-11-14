using Test
using BEAST
using CompScienceMeshes

function neighbortest(X)
    Γ = X.geo
    for (e,f) in enumerate(X.fns)
        i = Γ.faces[f[1].cellid].indices
        j = Γ.faces[f[2].cellid].indices
        ij = intersect(i,j)
        @assert length(ij) == 2
        ri = f[1].refid
        rj = f[2].refid
        #@assert !(i[ri] in ij)
        if (j[rj] in ij)
            @show e
            @show f
            @show i, ri
            @show j, rj
            @assert false
        end
    end
end

for T in [Float32,Float64]
    tol = eps(T) * 10^3

    mesh = meshrectangle(T(1.0), T(1.0), T(0.5))
    @test numvertices(mesh) == 9

    idcs = mesh.faces[1]
    # @test size(idcs) == (3,)
    @test idcs isa CompScienceMeshes.SimplexGraph{3}

    # verts = vertices(mesh, idcs)
    verts = chart(mesh, idcs).vertices
    @test size(verts) == (3,)

    faces = skeleton(mesh, 2)
    idcs = faces.faces[1]
    verts = chart(mesh, idcs).vertices
    local p = simplex(verts, Val{2})
    @test volume(p) == T(1/8)

    local edges = skeleton(mesh,1)
    @test numcells(edges) == 16

    local cps = cellpairs(mesh, edges)
    @test size(cps) == (2,16)

    # select only inner edges
    I = findall(x->(x>0), cps[2,:])
    cps = cps[:,I]
    @test size(cps) == (2,8)

    # build the Raviart-Thomas elements
    rt = raviartthomas(mesh, cps)
    @test numfunctions(rt) == 8



    neighbortest(rt)
end