using BEAST
using CompScienceMeshes
using LinearAlgebra
using Test

#steps for the terminal if not open properly: 1. using AnalytIntegralD0

#make a sphereAgl
Γ1 = meshsphere(radius=1.0, h=1.0)
Γ2 = meshsphere(radius=2.0, h=1.0)
# here for testing
# vertices = Γ1.vertices
# faces = Γ1.faces
# n = length(faces)

X1 = lagrangecxd0(Γ1)
X2 = lagrangecxd0(Γ2)
# t = Helmholtz3D.doublelayer(;) #hypersingular operator(cannot be tested). transposed double layer
t = Helmholtz3D.singlelayer(;)
# t = Helmholtz3D.doublelayer_transposed(;)

##

#get result 1
# @time begin
    Agl = assemble(t, X1, X2, quadstrat = BEAST.GumerovKanekoDuraiswamiStrat()) #threading=BEAST.Threading{:single}
    # Agl = assemble(t, X1, X2, BEAST.GumerovKanekoDuraiswamiStrat())  
# end

#get result 2
# @time begin
    Ab1 = assemble(t,X1,X2,quadstrat = BEAST.DoubleNumWiltonSauterQStrat(5,6,5,6,10,10,10,10))
# end
# Ab3 = assemble(BEAST.HH3DDoubleLayerTransposedFDBIO{Float64, Float64}(1.0, -0.0),X1,X2,quadstrat = BEAST.DoubleNumWiltonSauterQStrat(6,7,6,7,10,10,10,10))
#checking how close the values are
diff1 = Agl .- Ab1
diff_norm1 = norm(diff1)
perc_diff1 = norm(diff1)/norm(Ab1)

@test perc_diff1 < 3e-5
@test diff_norm1 < 3e-6

# # checking whether result is more accurate
# Ab2 = assemble(t,X1,X2,quadstrat = BEAST.DoubleNumWiltonSauterQStrat(6,7,6,7,30,30,30,30))
# diff2 = Agl .- Ab2
# diff_norm2 = norm(diff2)
# perc_diff2 = norm(diff2)/norm(Ab2)

# @test perc_diff1 < 1e-6
# @test diff_norm2 <= diff_norm1

# for i in 1:44
#     for j in 1:212
#         if abs(diff1[i,j]) > 1e-6
#         println("Case index: ", i, ", ", j)
#         # println("Case vertices x:", vertices[faces[i][1]],vertices[faces[i][2]],vertices[faces[i][3]])
#         # println("Case vertices y:", vertices[faces[j][1]],vertices[faces[j][2]],vertices[faces[j][3]])
#         e1 = Agl[i,j]
#         e2 = Ab1[i,j]
#         d1 = diff1[i,j]
#         println("Value from GL: $e1, from BEAST: $e2")
#         println("Difference: $d1")
#         println()
#         end
#     end
# end