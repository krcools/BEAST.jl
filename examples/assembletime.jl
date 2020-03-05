using BEAST
using CompScienceMeshes

Ω4 = CompScienceMeshes.tetmeshsphere(1,0.4)
Ω2 = CompScienceMeshes.tetmeshsphere(1,0.2)

faces4 = [sort(c) for c in cells(skeleton(boundary(Ω4),2))]
function is_fint(face)
    (sort(face) in faces4)
end
Γ4 = submesh(is_fint, boundary(Ω4))

faces2 = [sort(c) for c in cells(skeleton(boundary(Ω2),2))]
function is_fint(face)
    (sort(face) in faces2)
end
Γ2 = submesh(is_fint, boundary(Ω2))

X4 = BEAST.nedelecc3d(Ω4)
X2 = BEAST.nedelecc3d(Ω2)
ttrX4 = BEAST.ttrace(X4,Γ4)
ttrX2 = BEAST.ttrace(X2,Γ2)

Id = BEAST.Identity()
ttrX4 = BEAST.ttrace(X4,Γ4)

@time begin
    assemble(Id,ttrX4,ttrX4)
end
@time begin
    assemble(Id,ttrX2,ttrX2)
end
