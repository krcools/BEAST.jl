using CompScienceMeshes
using BEAST

using Test
for T in [Float32, Float64]
    local o, x, y, z = euclidianbasis(3,T)
    tet = simplex(x,y,z,o)
    local ctr = center(tet)


    iref = BEAST.LagrangeRefSpace{T,1,4,4}()
    oref = BEAST.NDLCCRefSpace{T}()
    ishp = BEAST.Shape{T}(1, 3, 1.0)

    output = BEAST.gradient(iref, ishp, tet)

    global oval = point(T,0,0,0)
    for oshp in output
        oval += oref(ctr)[oshp.refid].value * oshp.coeff
    end

    @test oval ≈ point(T,0,0,1)
end