using CompScienceMeshes
using BEAST

using Test

o, x, y, z = euclidianbasis(3)
tet = simplex(x,y,z,o)
ctr = center(tet)

T = Float64
iref = BEAST.LagrangeRefSpace{T,1,4,4}()
oref = BEAST.NDLCCRefSpace{T}()
ishp = BEAST.Shape{T}(1, 3, 1.0)

output = BEAST.gradient(iref, ishp, tet)

oval = point(0,0,0)
for oshp in output
    global oval += oref(ctr)[oshp.refid].value * oshp.coeff
end

@test oval ≈ point(0,0,1)
