@testitem "lagrangecx order=3 - global" begin
    using CompScienceMeshes

    projectdir = joinpath(dirname(pathof(BEAST)),"..")
    m = readmesh(joinpath(projectdir, "test/assets/sphere2.in"))

    lagspace3 = BEAST.lagrangecx(m, order=3)
    @test numfunctions(lagspace3) == 10 * length(m)
    @test refspace(lagspace3) == BEAST.LagrangeRefSpace{Float64,3,3,10}()
end

@testitem "lagrangecxd2: local, ref" begin
    using CompScienceMeshes

    T = Float64
    s = simplex(
        point(1,0,0),
        point(0,1,0),
        point(0,0,0),
    )
    p = CompScienceMeshes.center(s)

    ϕ = BEAST.LagrangeRefSpace{Float64,2,3,6}()
    A = ϕ(p)

    # test the dimension
    @test length(A) == binomial(4,2)

    # test the partition of unity property
    valp = sum(a.value for a in A)
    crlp = sum(a.curl for a in A)
    @test valp ≈ 1
    @test crlp ≈ point(0,0,0) atol=sqrt(eps(T))

    u = T(0.2); du = eps(T) * 1000
    v = T(0.6); dv = eps(T) * 1000

    p00 = neighborhood(s, (u,v))
    p10 = neighborhood(s, (u+du,v))
    p01 = neighborhood(s, (u, v+dv))

    ϕ00 = ϕ(p00)
    ϕ10 = ϕ(p10)
    ϕ01 = ϕ(p01)

    # @show [x.value for x in ϕ(neighborhood(s, (0.0, 0.0)))]

    tu = tangents(p00,1)
    tv = tangents(p00,2)
    j = jacobian(p00)

    for (f00, f10, f01) in zip(ϕ00, ϕ10, ϕ01)
        dfdu = (f10.value - f00.value)/du
        dfdv = (f01.value - f00.value)/dv
        curl_num = (dfdv * tu - dfdu * tv) / j
        curl_ana = f00.curl
        @test curl_num ≈ curl_ana atol=sqrt(eps(T))*100
    end
end

@testitem "lagrangecxd3: local, ref" begin
    using CompScienceMeshes

    T = Float64
    s = simplex(
        point(1,0,0),
        point(0,1,0),
        point(0,0,0),
    )
    p = center(s)

    ϕ = BEAST.LagrangeRefSpace{Float64,3,3,10}()
    A = ϕ(p)

    # test the dimension
    @test length(A) == binomial(5,2)

    # test the partition of unity property
    valp = sum(a.value for a in A)
    crlp = sum(a.curl for a in A)
    @test valp ≈ 1
    @test crlp ≈ point(0,0,0) atol=sqrt(eps(T))

    u = T(0.2); du = eps(T) * 1000
    v = T(0.6); dv = eps(T) * 1000

    p00 = neighborhood(s, (u,v))
    p10 = neighborhood(s, (u+du,v))
    p01 = neighborhood(s, (u, v+dv))

    ϕ00 = ϕ(p00)
    ϕ10 = ϕ(p10)
    ϕ01 = ϕ(p01)

    tu = tangents(p00,1)
    tv = tangents(p00,2)
    j = jacobian(p00)

    for (f00, f10, f01) in zip(ϕ00, ϕ10, ϕ01)
        dfdu = (f10.value - f00.value)/du
        dfdv = (f01.value - f00.value)/dv
        curl_num = (-dfdv * tu + dfdu * tv) / j
        curl_ana = f00.curl
        @test curl_num ≈ curl_ana atol=sqrt(eps(T))*100
    end
end


@testitem "lagrangecxd3: local, generic simplex" begin
    using CompScienceMeshes

    s = simplex(
        point(3,0,0),
        point(0,2,1),
        point(-1,-1,-1),
    )
    p = center(s)

    ϕ = BEAST.LagrangeRefSpace{Float64,3,3,10}()
    A = ϕ(p)

    # test the dimension
    @test length(A) == binomial(5,2)

    # test the partition of unity property
    valp = sum(a.value for a in A)
    @test valp ≈ 1
end

@testitem "lagrangecxd3: self-interpolate" begin
    using CompScienceMeshes
    using LinearAlgebra

    ϕ = BEAST.LagrangeRefSpace{Float64,3,3,10}()
    ch = simplex(
        point(3,0,0),
        point(0,2,1),
        point(-1,-1,-1))
    ch1toch2 = simplex(
        point(1,0),
        point(0,1),
        point(0,0))

    Q = BEAST.interpolate(ϕ, ch, ϕ, ch, ch1toch2)
    # display(round.(Q, digits=3))
    Id = Matrix{Float64}(I, 10, 10)
    @test Q ≈ Id
end


@testitem "lagrnacexcd3: interpolate generic poly" begin

    using CompScienceMeshes
    using LinearAlgebra

    ϕ = BEAST.LagrangeRefSpace{Float64,3,3,10}()
    ch = simplex(
        point(3,0,0),
        point(0,2,1),
        point(-1,-1,-1))

    function fields(p)
        u, v = parametric(p)
        return [
            1 + u + v + u^2 + u*v + v^2 + u^3 + u^2*v + u*v^2 + v^3,
            1 + u + u^2 + u^3,
            1 + u + u^2,
            1 + u,
            1
        ]
    end

    Q = BEAST.interpolate(fields, ϕ, ch)
    nbds = [
        neighborhood(ch, (0.2, 0.2)),
        neighborhood(ch, (0.2, 0.6)),
        neighborhood(ch, (0.6, 0.2)),
        neighborhood(ch, (1/3, 1/3)),
        neighborhood(ch, (0.0, 0.0))]

    for p in nbds
        basis = ϕ(p)

        vals = fields(p)
        for j in eachindex(vals)
            val1 = vals[j]
            val2 = sum(Q[j,i] * b.value for (i,b) in enumerate(basis))
            @test val1≈val2 atol=1e-8
        end
        # println()
    end
end


@testitem "lagc0d2: support size" begin
    using CompScienceMeshes
    using SparseArrays

    order = 2
    projectdir = joinpath(dirname(pathof(BEAST)),"..")
    m = readmesh(joinpath(projectdir, "test/assets/sphere2.in"))

    verts = skeleton(m,0)
    edges = skeleton(m,1)

    println(pathof(CompScienceMeshes))

    lagspace3 = BEAST.lagrangec0(m, order=order)
    @test numfunctions(lagspace3) == 
        length(verts) +
        length(edges)*(order-1) +
        length(m)*div((order-1)*(order-2),2)
    @test refspace(lagspace3) == BEAST.LagrangeRefSpace{Float64,2,3,6}()

    conn20 = connectivity(verts, m, x -> 1)
    num_adjacent_faces = vec(sum(conn20, dims=1))

    nv = length(verts)
    ne = length(edges) * (order-1)
    nf = length(m) * div((order-1) * (order-2),2)
    @test length.(lagspace3.fns[1:length(verts)]) == num_adjacent_faces
    @test all(length.(lagspace3.fns[nv+1:nv+ne]) .== 2)
    @test all(length.(lagspace3.fns[nv+ne+1:nv+ne+nf]) .== 1)
end


@testitem "lagc0d3: support size" begin
    using CompScienceMeshes
    using SparseArrays
    order = 3

    projectdir = joinpath(dirname(pathof(BEAST)),"..")
    m = readmesh(joinpath(projectdir, "test/assets/sphere2.in"))

    verts = skeleton(m,0)
    edges = skeleton(m,1)

    println(pathof(CompScienceMeshes))

    lagspace3 = BEAST.lagrangec0(m, order=3)
    @test numfunctions(lagspace3) == 
        length(verts) +
        length(edges)*(order-1) +
        length(m)*div((order-1)*(order-2),2)
    @test refspace(lagspace3) == BEAST.LagrangeRefSpace{Float64,3,3,10}()

    conn20 = connectivity(verts, m, x -> 1)
    num_adjacent_faces = vec(sum(conn20, dims=1))

    nv = length(verts)
    ne = length(edges) * (order-1)
    nf = length(m) * div((order-1) * (order-2),2)
    @test length.(lagspace3.fns[1:nv]) == num_adjacent_faces
    @test all(length.(lagspace3.fns[nv+1:nv+ne]) .== 2)
    @test all(length.(lagspace3.fns[nv+ne+1:nv+ne+nf]) .== 1)
end


# @testitem "lagrangec0 order=3 - continuity" begin
#     using CompScienceMeshes
#     using SparseArrays
#     order = 3

#     projectdir = joinpath(dirname(pathof(BEAST)),"..")
#     m = readmesh(joinpath(projectdir, "test/assets/sphere2.in"))

#     verts = skeleton(m,0)
#     edges = skeleton(m,1)

#     println(pathof(CompScienceMeshes))

#     lagspace3 = BEAST.lagrangec0(m, order=3)
#     @test numfunctions(lagspace3) == 
#         length(verts) +
#         length(edges)*(order-1) +
#         length(m)*div((order-1)*(order-2),2)
#     @test refspace(lagspace3) == BEAST.LagrangeRefSpace{Float64,3,3,10}()

#     conn20 = connectivity(verts, m)
# end


# @testitem "Lagc0d3: alternative construction" begin
#     using CompScienceMeshes
#     using SparseArrays
#     import Main.InteractiveUtils
#     order = 3

#     projectdir = joinpath(dirname(pathof(BEAST)),"..")
#     m = readmesh(joinpath(projectdir, "test/assets/sphere2.in"))

#     verts = skeleton(m,0)
#     edges = skeleton(m,1)

#     println(pathof(CompScienceMeshes))

#     lagspace3 = BEAST.lagrangec0(m, order=3)
#     @test numfunctions(lagspace3) == 
#         length(verts) +
#         length(edges)*(order-1) +
#         length(m)*div((order-1)*(order-2),2)
#     # @show InteractiveUtils.@which refspace(lagspace3)
#     @test refspace(lagspace3) == BEAST.LagrangeRefSpace{Float64,3,3,10}()

#     conn20 = connectivity(verts, m)
# end