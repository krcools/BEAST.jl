@testitem "conformity" begin

using CompScienceMeshes
    τ = simplex(
        point(0.0624999999994547, 0.6856034446626162, 0.0),
        point(0.06944444444542179, 0.6735753140519203, 0.0),
        point(0.05555730739992443, 0.6735783483387338, 0.0))

    σ = simplex(
        point(0.07017543859225508, 0.6988976942759024, 0.0),
        point(0.06249997795123046, 0.6856034064739717, 0.0),
        point(0.06140343344926745, 0.6837043913625904, 0.0))

    for (k,λ) in pairs(faces(τ))
        for (l,μ) in pairs(faces(σ))
            if CompScienceMeshes.overlap(λ, μ)
                global i = k
                global j = l
            end
    end end

    @test i == 2
    @test j == 3

    τs, σs = BEAST._conforming_refinement_touching_triangles(τ,σ,i,j)

    # Ideally both should be true but the round-off causes CompScienceMeshes.overlap to
    # incorrectly report.
    @test BEAST._test_conformity(τs[1], σs[1])==true
    @test BEAST._test_conformity(τs[2], σs[1])==false

end