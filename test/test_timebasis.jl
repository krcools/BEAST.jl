using BEAST
using CompScienceMeshes
using Test

# function Base.length(it::BEAST.ADIterator)
#     l = 0
#     for (m,w) in it
#         l += 1
#     end
#     l
# end
#
# tbf = BEAST.timeaxisc0d1(0,10.0,50)
#
# @test numfunctions(tbf) == 49
# @test refspace(tbf) == LagrangeRefSpace{Float64,1}()
#
# ad = assemblydata(tbf)
# @test size(ad.data) == (1,2,50)
# @test length(ad[1,1]) == 0
# @test length(ad[1,2]) == 1
# @test length(ad[50,1]) == 1
# @test length(ad[50,2]) == 0
# for c in 2:49
#     for r in 1:2
#         @test length(ad[c,r]) == 1
#     end
# end


tbf = BEAST.timebasisc0d1(1.0, 10)
timeels, timead = BEAST.assemblydata(tbf)
#ad = BEAST.temporalassemblydata(tbf)
;
