@testitem "GWPRefSpace eval" begin
    using CompScienceMeshes

    T = Float64
    Dim, Nverts, Degree = 2, 3, 3

    dom =  CompScienceMeshes.ReferenceSimplex{Dim,T,Nverts}()
    ϕ = BEAST.GWPCurlRefSpace{T,Degree}()

    u = (one(T)/3, one(T)/3)
    vals = ϕ(dom, u)
end