@testitem "local basis: nxBDM and BDM consistency" begin
    using CompScienceMeshes
    #testing local value in the center of a triangle 
    for T in [Float64]
        for j in [1,2,3,4,5,6]
            o, x, y, z = euclidianbasis(3,T)
            triang = simplex(x,y,o)
            ctr = center(triang)
            n = normal(triang)
            oref = BEAST.NCrossBDMRefSpace{T}()
            oref2= BEAST.BDMRefSpace{T}()
            oshp=BEAST.Shape(1,j,1.0)
            oshp2=BEAST.Shape(1,j,1.0)

            o1 = oref(ctr)[oshp.refid].value * oshp.coeff
            o2=n × oref2(ctr)[oshp2.refid].value * oshp2.coeff
            @test o1 ≈ o2
        end
    end

    m=meshsphere(radius=1.0, h=1.0)
    Z = BEAST.brezzidouglasmarini(m)
    NZ = BEAST.ncrossbdm(m)
    for i=(1,11,36,52,111)
        for j=(1,2)
            @test  Z.fns[i][j].coeff ≈ NZ.fns[i][j].coeff
            @test  Z.fns[i][j].refid ≈ NZ.fns[i][j].refid
            @test  Z.fns[i][j].cellid ≈ NZ.fns[i][j].cellid
        end
    end

    #testing the values on a mesh through Gram matrices
    Id = BEAST.Identity()
    qs = BEAST.SingleNumQStrat(8)

    Gzz= assemble(Id,Z,Z,quadstrat=qs)
    Gnznz= assemble(Id,NZ,NZ,quadstrat=qs)
    @test Gzz ≈ Gnznz
end

