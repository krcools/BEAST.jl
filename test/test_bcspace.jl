using Test
using SparseArrays

using CompScienceMeshes
using BEAST
using StaticArrays

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
    T=eltype(vertextype(mesh))
    edges = skeleton(mesh,1)
    D = connectivity(edges, mesh, abs)
    rows = rowvals(D)
    vals = nonzeros(D)

    Flux = zeros(T, numcells(mesh),3)
    TotalFlux = zeros(T, numcells(edges), numfunctions(space))

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
    verts = sort(verts.faces, lt=isless)

    edges = skeleton(mesh,1)
    edges = sort(edges.faces, lt=isless)

    function pred(s)

        num_hits = 0
        for i in 1:length(s)
            if !isempty(searchsorted(verts, SVector(s[i]), lt=isless))
                num_hits += 1
            end
        end

        @assert num_hits < 3
        num_hits == 0 && return false
        num_hits == 1 && return true
        isempty(searchsorted(edges, s, lt=isless)) && return true
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

#meshfile = Pkg.dir("BEAST","test","sphere2.in")
for T in [Float32, Float64]
    meshfile = joinpath(dirname(@__FILE__),"assets","sphere316.in")
    mesh = readmesh(meshfile,T=T)
    @test numvertices(mesh) == 160
    @test numcells(mesh) == 316

    rt = raviartthomas(mesh)
    @test numfunctions(rt) == 316 * 3 / 2

    local fine = barycentric_refinement(mesh)
    local edges = skeleton(mesh, 1)

    bc = buffachristiansen(mesh)
    @test numfunctions(bc) == 316 * 3 / 2

    lc = isdivconforming(rt)
    @test maximum(lc) < eps(T) * 1000
    println("RT space is div-conforming")

    lc = isdivconforming(bc);
    @test maximum(lc) < eps(T) * 1000
    println("BC space is div-conforming")

    # Now repeat the exercise with an open mesh
    mesh = meshrectangle(T(1.0), T(1.0), T(0.2));
    fine = barycentric_refinement(mesh);

    rt = raviartthomas(mesh)
    bc = buffachristiansen(mesh)

    @test numfunctions(rt) == 65
    @test numfunctions(bc) == 65

    int_pred = interior_tpredicate(mesh)
    bnd_pred(s) = !int_pred(s)

    leaky_edges = findall(sum(abs.(isdivconforming(rt)),dims=1) .!= 0)
    @test length(leaky_edges) == 0

    bnd = boundary(mesh)
    bndtch_pred = touches_predicate(bnd)
    edges = interior(mesh)
    bndtch_edges = findall(bndtch_pred, cells(edges))
    leaky_edges = findall(vec(sum(abs.(isdivconforming(bc)),dims=1)) .!= 0)
    @test bndtch_edges == leaky_edges


    ## Test the charge of BC functions
    #meshfile = Pkg.dir("BEAST","test","sphere2.in")
    meshfile = joinpath(dirname(@__FILE__),"assets","sphere316.in")
    mesh = readmesh(meshfile,T=T)
    bc = buffachristiansen(mesh)
    fine = geometry(bc)
    charges = zeros(numcells(fine))

    for fn in bc.fns
        abs_charge = T(0.0)
        net_charge = T(0.0)
        fill!(charges,0)
        for  _sh in fn
            cellid = _sh.cellid

            #cell = simplex(vertices(fine, fine.faces[cellid]))
            net_charge += _sh.coeff
            charges[cellid] += _sh.coeff
        end
        abs_charge = sum(abs.(charges))
        @test net_charge + 1 ≈ 1
        @test abs_charge ≈ 2
    end
end

# THe BC construction function should throw for non-oriented surfaces
# width, height, h = 1.0, 0.5, 0.05
G1 = meshrectangle(1.0, 1.0, 0.25)
G2 = CompScienceMeshes.rotate(G1, 0.5π * x̂)
G = CompScienceMeshes.weld(G1,G2)
@test_throws AssertionError buffachristiansen(G)
# G3 = CompScienceMeshes.rotate(G1, 1.0π * x̂)

@testitem "Neumann BC" begin
    using CompScienceMeshes
    using LinearAlgebra

    Γ = meshrectangle(1.0, 1.0, 0.1, 3)
    all_edges = skeleton(Γ, 1)
    bnd_edges = boundary(Γ)
    int_edges = setminus(all_edges, bnd_edges)

    X = raviartthomas(Γ, all_edges)
    Y = Ydi = BEAST.buffachristiansen(Γ, bnd_edges)

    T = Maxwell3D.singlelayer(wavenumber=1.0)
    Id = BEAST.NCross()

    Txx = assemble(T, X, X)
    Tyy = assemble(T, Y, Y)

    Ixy = assemble(Id, X, Y)
    IYX = BEAST.GMRES(Ixy)
    IXY = BEAST.GMRES(Ixy')

    𝗧 = Matrix(Txx)
    𝕋 = Matrix(Tyy)
    𝗜 = Matrix(Ixy)
    𝗜⁻¹ = inv(𝗜)
    𝗜⁻ᵀ = transpose(𝗜⁻¹)
    𝗣 = 𝗜⁻ᵀ * 𝕋 * 𝗜⁻¹
    𝗣𝗧 = 𝗣*𝗧

    @test cond(𝗧) > 1e3
    @test cond(𝕋) > 1e3
    @show cond(𝗜) < 20
    @show cond(𝗣*𝗧) < 20
end