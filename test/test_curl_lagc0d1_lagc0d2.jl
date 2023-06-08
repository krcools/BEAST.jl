#test curl lagc0d2

using CompScienceMeshes
using BEAST

using Test
for T in [Float32, Float64]
    for j in [1,2,3]
        local o, x, y, z = euclidianbasis(3,T)
        triang = simplex(x,y,o)
        ctr = center(triang)
    


        iref1 = BEAST.LagrangeRefSpace{T,2,3,6}()
        oref1 = BEAST.BDMRefSpace{T}()
        ishp1 = [BEAST.Shape{T}(1, j, 1.0), BEAST.Shape{T}(1, mod1(j+1,3)+3, 0.5), BEAST.Shape{T}(1, mod1(j+2,3)+3, 0.5)]

        # we use the property that a nodal lagc0d2 function plus 1/2 times the two adiacent edge lagc0d2 functions is equal to the nodal lagc0d1 function at the same nodes, and the same must be true for their surfcurl
        
        outputbdm = [BEAST.curl(iref1, ishp1[1], triang), BEAST.curl(iref1, ishp1[2], triang), BEAST.curl(iref1, ishp1[3], triang)]

        global obdm = point(T,0,0,0)
        
        for i in [1,2,3]
            for oshp in outputbdm[i]
            obdm += oref1(ctr)[oshp.refid].value * oshp.coeff
            end
        end

        
        iref2 = BEAST.LagrangeRefSpace{T,1,3,3}()
        oref2 = BEAST.RTRefSpace{T}()
        ishp2 = BEAST.Shape{T}(1, j, 1.0)

        outputrt = BEAST.curl(iref2, ishp2, triang)

        global ort = point(T,0,0,0)
        
        for oshp in outputrt
            ort += oref2(ctr)[oshp.refid].value * oshp.coeff
        end
    

        @test ort ≈ obdm
    end
end
