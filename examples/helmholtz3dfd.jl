using CompScienceMeshes
using BEAST

m = meshsphere(1.0, 0.3)
vertices = skeleton(m, 0)
X = lagrangec0d1(m)

Test.@test numfunctions(X) == numcells(vertices)

κ = 1.0
γ = im*κ
S = HH3DHyperSingularFDBIO(γ)
