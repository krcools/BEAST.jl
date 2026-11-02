@testitem "globaldofs: interpolatory on edge" begin
    using CompScienceMeshes

    T = Float64
    face_ch = simplex(
        point(3,0,0),
        point(2,0,-1),
        point(0,0,1),
    )
    edge_ch = simplex(
        point(0,0,1),
        point(2,0,-1),
    )

    localspace = BEAST.LagrangeRefSpace{T,3,3,10}()
    L = BEAST._LagrangeGlobalEdgeDoFs(3)
    dofs = BEAST.globaldofs(edge_ch, face_ch, localspace, L)

    @test size(dofs) == (10,2)
    A = [
        0 0
        0 0
        0 0
        0 0
        0 0
        0 0
        1 0
        0 0
        0 1
        0 0]
    @test dofs ≈ A
end