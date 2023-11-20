mutable struct AcusticSingleLayerTDIO{T} <: RetardedPotential{T}
    "speed of sound in medium"
    speed_of_light::T #for compatibility with the rest of the package I left light instead of sound
    "weight"
    weight::T
    "number of temporal differentiations"
    diffs::Int
end

function Base.:*(a::Number, op::AcusticSingleLayerTDIO)
	@info "scalar product a * op (acusticSL)"
	AcusticSingleLayerTDIO(
		op.speed_of_light,
		a * op.weight,
		op.diffs)
end

AcusticSingleLayerTDIO(;speedofsound) = AcusticSingleLayerTDIO(speedofsound, one(speedofsound), 0)

module TDAcustic3D
import ...BEAST

function acusticsinglelayer(;speedofsound, numdiffs)
    @assert numdiffs >= 0
        #controllare con Kristof che tipo di schema temporale viene utilizzato e come si con integrali rispetto al tempo extra
        #return BEAST.integrate(BEAST.MWSingleLayerTDIO(speedoflight,-1/speedoflight,-speedoflight,2,0))
	return BEAST.AcusticSingleLayerTDIO(speedofsound,one(speedofsound),numdiffs)
end

end #of the module

export TDAcustic3D

defaultquadstrat(::AcusticSingleLayerTDIO, tfs, bfs) = AllAnalyticalQStrat(1)





function quaddata(op::AcusticSingleLayerTDIO, testrefs, trialrefs, timerefs,
    testels, trialels, timeels, quadstrat::AllAnalyticalQStrat)
    return 0.0
end

quadrule(op::AcusticSingleLayerTDIO, testrefs, trialrefs, timerefs,
        p, testel, q, trialel, r, timeel, qd, ::AllAnalyticalQStrat) = ZuccottiStrat(1.0)

function quaddata(op::AcusticSingleLayerTDIO, testrefs, trialrefs, timerefs,
        testels, trialels, timeels, quadstrat::OuterNumInnerAnalyticQStrat)
    
    dmax = numfunctions(timerefs)-1
    bn = binomial.((0:dmax),(0:dmax)')
    
    V = eltype(testels[1].vertices)
    ws = WiltonInts84.workspace(V)
    # quadpoints(testrefs, testels, (3,)), bn, ws
    quadpoints(testrefs, testels, (quadstrat.outer_rule,)), bn, ws
end
    
    
quadrule(op::AcusticSingleLayerTDIO, testrefs, trialrefs, timerefs,
        p, testel, q, trialel, r, timeel, qd, ::OuterNumInnerAnalyticQStrat) = WiltonInts84Strat(qd[1][1,p],qd[2],qd[3])

function momintegrals!(z, op::AcusticSingleLayerTDIO, g::LagrangeRefSpace{T,0,3}, f::LagrangeRefSpace{T,0,3}, t::MonomialBasis{T,0,1}, τ, σ, ι, qr::ZuccottiStrat) where T
        if ι[1]<0
           t1=0.0
        else
            t1=ι[1]
        end
          t2=ι[2]
        
            
          @assert t2 > t1
        
          sos = op.speed_of_light 

          a1index,a2index,a3index,b1index,b2index,b3index=0,0,0,0,0,0
    Tt = eltype(eltype(τ.vertices))
    hits = 0
    dtol = 1.0e3 * eps(Tt)
    dmin2 = floatmax(Tt)
    for t in 1:3
        for s in 1:3
            d2 = LinearAlgebra.norm_sqr(τ[t]-σ[s])
            dmin2 = min(dmin2, d2)
            if (d2 < dtol) 
                hits+=1
                if hits==1
                    a1index =t
                    b1index = s 
                elseif hits==2
                    a2index=t
                    b2index=s
                end
            end
        end
    end
    #print(τ[1]," ",τ[2]," ",τ[3]," ",σ[1]," ",σ[2]," ",σ[3]," ",t1," ",t2," ",)
          if hits==3 
             z[1,1,1]+=TimeDomainBEMInt.intcoinctriangles(τ[1],τ[2],τ[3],t1,t2)
          elseif hits==2
           
            if mod1(a1index +1,3)==a2index
                a3index=mod1(a1index-1,3)
            else
                a3index=mod1(a1index+1,3)
            end

            if mod1(b1index+1,3)==b2index
                b3index=mod1(a1index-1,3)
            else
                b3index=mod1(a1index+1,3)
            end

            z[1,1,1]+=TimeDomainBEMInt.inttriangletriangleadjacent(τ[a1index],τ[a2index],τ[a3index],σ[b1index],σ[b2index],σ[a3index],t1,t2)
          elseif hits==1
                a2index,a3index=mod1(a1index+1,3),mod1(a1index+2,3)
                b2index,b3index=mod1(b1index+1,3),mod1(b1index+2,3)

                z[1,1,1]+=TimeDomainBEMInt.intcommonvertex(τ[a1index],τ[a2index],τ[a3index],σ[b2index],σ[b3index],t1,t2)
          else
             z[1,1,1]+=inttriangletriangle(τ[1],τ[2],τ[3],σ[1],σ[2],σ[3],t1,t2)
          end
          #for the moment sos=1 but I will correct this
end