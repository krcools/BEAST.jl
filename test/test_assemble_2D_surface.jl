using Test

import Permutations
import StaticArrays
using SparseArrays
using CompScienceMeshes
using BEAST

"""
Test if the assemble methode for a 2D mesh with udim=2 works. 
Test are constructing a rectangle with udim=2 and udim=3 with langrange and nedelec interpolation. 
And assemble these two end checking that matrix from rectangle of udim=2 is equal to the matrix from the ractangle of udim=3.
"""



# overwrite the existing method of CompScienceMeshes for the positive value of the volume.
function CompScienceMeshes._normals(tangents::StaticArrays.SVector{2,StaticArrays.SVector{2,T}}, ::Type{Val{0}}) where {T}

    t = tangents[1]
    s = tangents[2]
    v = abs(t[1]*s[2] - t[2]*s[1])/2
    # n[3] = tangents[1] × tangents[2]
    # l = norm(n)

    P = StaticArrays.SVector{2,T}
    StaticArrays.SVector{0,P}(), v
end


function permutate_vector(X, Y)
    tol = sqrt(eps(eltype(X[1])))

    permut = Vector{Int32}()  
    temp = collect(1:length(X))

    for p in Y
        index = findfirst(isapprox(p;atol = tol), X)
        @assert !isnothing(index)   # vertex from Y exist in X
        @assert temp[index] != 0    # vertex from Y unique in X

        temp[index] = 0
        append!(permut, index)
    end

    for i in temp
        if i !=0
            push!(permut, i)
    end end

    return permut
end






@testset "assemble 2D langrange" begin
    h = 1/5
    Ω2D = meshrectangle(1.0, 1.0, h, 2) # udim = 2
    Ω3D = meshrectangle(1.0, 1.0, h, 3) # udim = 3

    #permutate_mesh(Ω2D, Ω3D)

    # lagrange spaces
    X2D = lagrangec0d1(Ω2D, skeleton(Ω2D,0))
    X3D = lagrangec0d1(Ω3D, skeleton(Ω3D,0))

    σ =  permutate_vector(X2D.pos, map(vec -> StaticArrays.SVector(vec[1], vec[2]), X3D.pos))

    @test map(vec -> StaticArrays.SVector(vec[1], vec[2]), X3D.pos) == X2D.pos[σ]

    #@test X2D.geo.vertices == map(x-> x[1:2], X3D.geo.vertices)

    I = BEAST.Identity()

    A2D = assemble(I, X2D, X2D)
    A3D = assemble(I, X3D, X3D)

    @test sum(abs.(permute(sparse(A2D), σ, σ) - sparse(A3D))) ≈ 0 atol=1e-8

end

@testset "assemble 2D nedelec" begin
    h = 1/5
    Ω2D = meshrectangle(1.0, 1.0, h, 2) # udim = 2
    Ω3D = meshrectangle(1.0, 1.0, h, 3) # udim = 3

    #permutate_mesh(Ω2D, Ω3D)

    # lagrange spaces
    X2D = lagrangec0d1(Ω2D, skeleton(Ω2D,0))
    X3D = lagrangec0d1(Ω3D, skeleton(Ω3D,0))

    σ =  permutate_vector(X2D.pos, map(vec -> StaticArrays.SVector(vec[1], vec[2]), X3D.pos))

    @test map(vec -> StaticArrays.SVector(vec[1], vec[2]), X3D.pos) == X2D.pos[σ]

    @hilbertspace u
    @hilbertspace v

    I = BEAST.Identity()
    IGG = @varform I[gradient(u), gradient(v)]

    A2D = assemble(IGG, X2D, X2D)
    A3D = assemble(IGG, X3D, X3D)

    @test sum(abs.(permute(sparse(A2D), σ, σ) - sparse(A3D))) ≈ 0 atol=1e-8
end