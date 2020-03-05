using BEAST
using CompScienceMeshes

Ω4 = CompScienceMeshes.tetmeshsphere(1,0.4)
Ω2 = CompScienceMeshes.tetmeshsphere(1,0.2)

X4 = BEAST.nedelecc3d(Ω4)
X2 = BEAST.nedelecc3d(Ω2)
ttrX4 = BEAST.ttrace(X4,boundary(Ω4))
ttrX2 = BEAST.ttrace(X2,boundary(Ω2))

Id = BEAST.Identity()
@time begin
    assemble(Id,ttrX4,ttrX4)
end
@time begin
    assemble(Id,ttrX2,ttrX2)
end
