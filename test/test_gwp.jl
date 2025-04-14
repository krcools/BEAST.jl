@testitem "refspace: dimension" begin
    using CompScienceMeshes

    T = Float64
    Dim, Nverts, Degree = 2, 3, 3

    dom =  CompScienceMeshes.ReferenceSimplex{Dim,T,Nverts}()
    ϕ = BEAST.GWPCurlRefSpace{T,Degree}()

    u = (one(T)/3, one(T)/3)
    vals = ϕ(dom, u)

    nf = (Degree+1)*(Degree+3)
    @test length(vals) == nf
end


@testitem "refspace: self-interpolate" begin
    using CompScienceMeshes
    using LinearAlgebra

    T, Dim, Nverts, Degree = Float64, 2, 3, 2
    dom =  CompScienceMeshes.ReferenceSimplex{Dim,T,Nverts}()
    ϕ = BEAST.GWPCurlRefSpace{T,Degree}()

    fields(p) = [v.value for v in ϕ(dom, p)]
    # fields(p) = [ϕ(dom,p)[1].value]
    coeffs = BEAST.interpolate(fields, ϕ, dom)

    nf = numfunctions(ϕ, dom)
    @test coeffs ≈ Matrix{T}(I,nf,nf) atol=sqrt(eps(T))
    # display(round.(coeffs, digits=3))

end


@testitem "confspace: interpolate lowest order" begin
    using CompScienceMeshes
    using LinearAlgebra

    T, Dim, Nverts, Degree = Float64, 2, 3, 2
    supp = simplex(
        point(T,3,0,0),
        point(T,0,2,0),
        point(T,0,0,-1))

    ϕ = BEAST.GWPCurlRefSpace{T,Degree}()
    ψ = BEAST.NDRefSpace{T}()

    fields(p) = [v.value for v in ψ(p)]
    coeffs = BEAST.interpolate(fields, ϕ, supp)

    nf1 = numfunctions(ψ, domain(supp))
    nf2 = numfunctions(ϕ, domain(supp))
    @test size(coeffs) == (nf1, nf2)
    # display(round.(coeffs, digits=3))

    pts = [
        point(T, 0.3, 0.1),
        point(T, 0.1, 0.3),
        point(T, 0.5, 0.0)]

    for (i,u) in enumerate(pts)
        p = neighborhood(supp, u)
        basisvals = ϕ(p)
        fieldvals = ψ(p)
        r = sum(c * v.value for (c,v) in zip(coeffs[i,:], basisvals))
        s = fieldvals[i].value
        @test r ≈ s atol=sqrt(eps(T))
    end
end