@testitem "Non-conf mesh: scalar, cuboid" tags=[:diagnostics] begin
    using CompScienceMeshes
    @show Threads.nthreads()

    function CompScienceMeshes.volume(m::CompScienceMeshes.AbstractMesh)
        sum(volume(chart(m,c)) for c in m)
    end

    struct ConstantIntegralOp{T} <: BEAST.IntegralOperator
        constant::T
    end

    BEAST.scalartype(op::ConstantIntegralOp{T}) where {T} = T

    function (igd::BEAST.Integrand{<:ConstantIntegralOp})(p,q,f,g)
        c = igd.operator.constant

        BEAST._integrands(f,g) do fi,gj
            dot(fi.value, gj.value) * c
        end
    end

    h1 = 0.2 * 0.75^3
    h2 = 0.2 * 0.75^4

    Γ1 = meshcuboid(1.0, 1.0, 1.0, h1)
    Γ2 = meshcuboid(1.0, 1.0, 1.0, h2)

    Γ11 = skeleton(Γ1,1)
    Γ21 = skeleton(Γ2,1)

    Γ10 = skeleton(Γ1,0)
    Γ20 = skeleton(Γ2,0)

    @show length(Γ1) - length(Γ11) + length(Γ10)
    @show length(Γ2) - length(Γ21) + length(Γ20)

    v1 = volume(Γ1)
    v2 = volume(Γ2)

    X1 = lagrangec0d1(Γ1)
    X2 = lagrangec0d1(Γ2)

    op = ConstantIntegralOp(1.0)

    qs1 = BEAST.DoubleNumQStrat(3,3)
    qs2 = BEAST.NonConformingIntegralOpQStrat(qs1)

    A1 = assemble(op, X1, X2; quadstrat=qs1)
    A2 = assemble(op, X1, X2; quadstrat=qs2)
    u1 = fill(1.0, numfunctions(X1))
    u2 = fill(1.0, numfunctions(X2))

    r1 = dot(u1, A1*u2)
    r2 = dot(u1, A2*u2)
    @show r1
    @show r2
    @show v1 * v2
    @test v1 * v2 - r1 < eps() * 100
    @test v1 * v2 - r2 < eps() * 100
end


@testitem "Non-conf mesh: scalar, smallboxbigbox" tags=[:diagnostics] begin
    using CompScienceMeshes
    @show Threads.nthreads()

    function CompScienceMeshes.volume(m::CompScienceMeshes.AbstractMesh)
        sum(volume(chart(m,c)) for c in m)
    end

    struct ConstantIntegralOp{T} <: BEAST.IntegralOperator
        constant::T
    end

    BEAST.scalartype(op::ConstantIntegralOp{T}) where {T} = T

    function (igd::BEAST.Integrand{<:ConstantIntegralOp})(p,q,f,g)
        c = igd.operator.constant

        BEAST._integrands(f,g) do fi,gj
            dot(fi.value, gj.value) * c
        end
    end

    # h1 = 0.2 * 0.75^3
    # h2 = 0.2 * 0.75^4
    hs = collect(0.25 * 0.835 .^ (0:1:14))
    href = hs[end]

    geofile = joinpath(pkgdir(BEAST), "test", "assets", "smallboxbigbox.geo")
    r = 0.5
    
    physicals = ["Gamma1", "Gamma2", "Gamma3"]
    err1 = zeros(length(hs)-1, length(physicals))
    err2 = zeros(length(hs)-1, length(physicals))
    for (i,physical) in enumerate(physicals)
        Γ2 = CompScienceMeshes.meshgeo(geofile; h = href, r=r, dim=2, physical=physical)
        Γ21 = skeleton(Γ2,1)
        Γ20 = skeleton(Γ2,0)
        @show length(Γ2)
        @show length(Γ2) - length(Γ21) + length(Γ20)
        v2 = volume(Γ2)
        X2 = lagrangec0d1(Γ2, Γ20)
        for (j,h) in enumerate(hs[1:end-1])
        
            Γ1 = CompScienceMeshes.meshgeo(geofile; h = h, r=r, dim=2, physical=physical)
            Γ11 = skeleton(Γ1,1)
            Γ10 = skeleton(Γ1,0)
            @show length(Γ1)
            @show length(Γ1) - length(Γ11) + length(Γ10)
            v1 = volume(Γ1)
            X1 = lagrangec0d1(Γ1, Γ10)

            op = ConstantIntegralOp(1.0)
            qs1 = BEAST.DoubleNumQStrat(3,3)
            qs2 = BEAST.NonConformingIntegralOpQStrat(qs1)

            A = assemble(op, X1, X2; quadstrat=qs1)
            u1 = fill(1.0, numfunctions(X1))
            r1 = dot(u1, A*u2)
            A = assemble(op, X1, X2; quadstrat=qs2)
            u2 = fill(1.0, numfunctions(X2))
            r2 = dot(u1, A*u2)

            @show v1 * v2 - r1
            @show v1 * v2 - r2

            err1[j,i] = v1 * v2 - r1
            err2[j,i] = v1 * v2 - r2

            @test v1 * v2 - r1 < eps() * 1e6
            @test v1 * v2 - r2 < eps() * 1e6
        end
    end
end



@testitem "Non-conf mesh - vectorial" tags=[:diagnostics] begin
    using CompScienceMeshes

    function CompScienceMeshes.volume(m::CompScienceMeshes.AbstractMesh)
        sum(CompScienceMeshes.volume(chart(m,c)) for c in m)
    end

    struct ConstantIntegralOp{T} <: BEAST.IntegralOperator
        constant::T
    end

    BEAST.scalartype(op::ConstantIntegralOp{T}) where {T} = T

    function (igd::BEAST.Integrand{<:ConstantIntegralOp})(p,q,f,g)
        c = igd.operator.constant

        BEAST._integrands(f,g) do fi,gj
            dot(fi.value, gj.value) * c
        end
    end

    E1(x) = point(1.0, 1.0, 0.0)
    E2(x) = point(0.0, 1.0, 1.0)
    BEAST.scalartype(::typeof(E1)) = Float64
    BEAST.scalartype(::typeof(E2)) = Float64

    h1 = 0.2 * 0.75^3
    h2 = 0.2 * 0.75^4

    Γ1 = meshcuboid(1.0, 1.0, 1.0, h1; generator=:gmsh)
    Γ2 = meshcuboid(1.0, 1.0, 1.0, h2; generator=:gmsh)

    pred = (m,c) -> begin
        x = chart(m,c) |> center |> cartesian
        !(x[3] ≈ 0.0)
    end
    Γ1 = submesh(pred, Γ1)
    Γ2 = submesh(pred, Γ2)

    v1 = CompScienceMeshes.volume(Γ1)
    v2 = CompScienceMeshes.volume(Γ2)

    X1 = raviartthomas(Γ1, skeleton(Γ1,1))
    X2 = raviartthomas(Γ2, skeleton(Γ2,1))

    Id = BEAST.Identity()
    I11 = assemble(Id, X1, X1)
    I22 = assemble(Id, X2, X2)

    u11 = I11 \ assemble(n × E1, X1)
    u12 = I22 \ assemble(n × E1, X2)
    u21 = I11 \ assemble(n × E2, X1)
    u22 = I22 \ assemble(n × E2, X2)

    op = ConstantIntegralOp(1.0)

    qs1 = BEAST.DoubleNumQStrat(3,3)
    qs2 = BEAST.NonConformingIntegralOpQStrat(qs1)

    A12 = assemble(op, X1, X2; quadstrat=qs2)
    A11 = assemble(op, X1, X1; quadstrat=qs1)
    A22 = assemble(op, X2, X2; quadstrat=qs1)

    @show r12 = dot(u11, A12*u22)
    @show r11 = dot(u11, A11*u21)
    @show r22 = dot(u12, A22*u22)

    @test abs(r12 - r11) < 1e-10
    @test abs(r12 - r22) < 1e-10 
end


@testitem "Non-conf mesh - vectorial/singular" tags=[:diagnostics] begin
    using CompScienceMeshes

    function CompScienceMeshes.volume(m::CompScienceMeshes.AbstractMesh)
        sum(CompScienceMeshes.volume(chart(m,c)) for c in m)
    end

    struct ConstantIntegralOp{T} <: BEAST.IntegralOperator
        constant::T
    end

    BEAST.scalartype(op::ConstantIntegralOp{T}) where {T} = T

    function (igd::BEAST.Integrand{<:ConstantIntegralOp})(p,q,f,g)
        c = igd.operator.constant

        BEAST._integrands(f,g) do fi,gj
            dot(fi.value, gj.value) * c
        end
    end

    E1(x) = point(1.0, 1.0, 0.0)
    E2(x) = point(0.0, 1.0, 1.0)
    BEAST.scalartype(::typeof(E1)) = Float64
    BEAST.scalartype(::typeof(E2)) = Float64

    h1 = 0.2 * 0.75^3
    h2 = 0.2 * 0.75^4

    Γ1 = meshcuboid(1.0, 1.0, 1.0, h1; generator=:gmsh)
    Γ2 = meshcuboid(1.0, 1.0, 1.0, h2; generator=:gmsh)

    pred = (m,c) -> begin
        x = chart(m,c) |> center |> cartesian
        !(x[3] ≈ 0.0)
    end
    Γ1 = submesh(pred, Γ1)
    Γ2 = submesh(pred, Γ2)

    v1 = CompScienceMeshes.volume(Γ1)
    v2 = CompScienceMeshes.volume(Γ2)

    X1 = raviartthomas(Γ1, skeleton(Γ1,1))
    X2 = raviartthomas(Γ2, skeleton(Γ2,1))

    Id = BEAST.Identity()
    I11 = assemble(Id, X1, X1)
    I22 = assemble(Id, X2, X2)

    u11 = I11 \ assemble(n × E1, X1)
    u12 = I22 \ assemble(n × E1, X2)
    u21 = I11 \ assemble(n × E2, X1)
    u22 = I22 \ assemble(n × E2, X2)

    op = Maxwell3D.singlelayer(gamma=1.0)

    q = 4
    qs1 = BEAST.DoubleNumSauterQstrat(3+q,3+q,3+q,3+q,3+q,3+q)
    # qs1 = BEAST.DoubleNumWiltonBogaertQStrat(3+q,3+q,3+q,3+1+q)
    qs2 = BEAST.NonConformingIntegralOpQStrat(qs1)

    A12 = assemble(op, X1, X2; quadstrat=qs2)
    A11 = assemble(op, X1, X1; quadstrat=qs1)
    A22 = assemble(op, X2, X2; quadstrat=qs1)

    @show r12 = dot(u11, A12*u22)
    @show r11 = dot(u11, A11*u21)
    @show r22 = dot(u12, A22*u22)

    @show abs(r12 - r11)
    @show abs(r12 - r22)
end