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
#nothing goes in hybrid qr, allanalytical goes in zuccottirule




function quaddata(op::AcusticSingleLayerTDIO, testrefs, trialrefs, timerefs,
    testels, trialels, timeels, quadstrat::AllAnalyticalQStrat)
    return nothing
end

quadrule(op::AcusticSingleLayerTDIO, testrefs, trialrefs, timerefs,
        p, testel, q, trialel, r, timeel, qd, ::AllAnalyticalQStrat) = ZuccottiRule(1.0)
    


function quaddata(operator::AcusticSingleLayerTDIO,
            test_local_space, trial_local_space, time_local_space,
            test_element, trial_element, time_element, quadstrat::Nothing)
    
        dmax = numfunctions(time_local_space)-1
        bn = binomial.((0:dmax),(0:dmax)')
    
        V = eltype(test_element[1].vertices)
        ws = WiltonInts84.workspace(V)
        order = 4
        @show order
        quadpoints(test_local_space, test_element, (order,)), bn, ws
    
    end
    
    
    # See: ?BEAST.quadrule for help
    function quadrule(operator::AcusticSingleLayerTDIO,
            test_local_space, trial_local_space, time_local_space,
            p, test_element, q, trial_element, r, time_element,
            quad_data, quadstrat::Nothing)
    
        # WiltonInts84Strat(quad_data[1,p])
        qd = quad_data
        HybridZuccottiWiltonStrat(qd[1][1,p],qd[2],qd[3])
    
    end
    
    
    function innerintegrals!(zlocal, operator::AcusticSingleLayerTDIO,
            test_point,
            test_local_space, trial_local_space, time_local_space,
            test_element, trial_element, time_element,
            quad_rule::HybridZuccottiWiltonStrat, quad_weight)
    
        # error("Here!!!")
    
        dx = quad_weight
        x = cartesian(test_point)
        # n = normal(test_point)
    
        # a = trial_element[1]
        # ξ = x - dot(x -a, n) * n
    
        r = time_element[1]
        R = time_element[2]
        @assert r < R
    
        N = max(degree(time_local_space), 1)
        ∫G, ∫vG, ∫∇G = WiltonInts84.wiltonints(
            trial_element[1],
            trial_element[2],
            trial_element[3],
            x, r, R, Val{2}, quad_rule.workspace)
    
        a = dx / (4*pi)
        #D = operator.num_diffs
        D=0
        @assert D == 0
        @assert numfunctions(test_local_space)  == 1
        @assert numfunctions(trial_local_space) == 1
    
        @inline function tmRoR_sl(d, iG)
            sgn = isodd(d) ? -1 : 1
            r = sgn * iG[d+2]
        end
    
        # bns = quad_rule.binomials
    
        @assert D == 0
        for k in 1 : numfunctions(time_local_space)
            d = k - 1
            d < D && continue
            q = reduce(*, d-D+1:d ,init=1)
            zlocal[1,1,k] += a * q * tmRoR_sl(d-D, ∫G)
        end # k
    end

    

qpclineline=10^8#quasiparallel case controller
qpclinetriang=10^10
qpctriangtriang=10^13  

const qpc=[qpclineline, qpclinetriang, qpctriangtriang] 

function momintegrals!(z, op::AcusticSingleLayerTDIO, g::LagrangeRefSpace{T,0,3}, f::LagrangeRefSpace{T,0,3}, t::MonomialBasis{T,0,1}, τ, σ, ι, qr::ZuccottiRule) where T
        
        t1=ι[1]
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
    
          if hits==3 
            #print("\np1=",τ[1],"\np2=",τ[2],"\np3=",τ[3],"\nt1=",t1,"\nt2=",t2)
        
             z[1,1,1]+=(TimeDomainBEMInt.intcoinctriangles(τ[1],τ[2],τ[3],t1,t2))/(4*π)
          elseif hits==2
            #pay attention with double layer and index permutation or in whatever case has nx 
            if mod1(a1index +1,3)==a2index
                a3index=mod1(a1index-1,3)
            else
                a3index=mod1(a1index+1,3)
            end

            if mod1(b1index+1,3)==b2index
                b3index=mod1(b1index-1,3)
            else
                b3index=mod1(b1index+1,3)
            end
            #print("\np1=",τ[a1index],"\np2=",τ[a2index],"\np3=",τ[a3index],"\nv1=",σ[b1index],"\nv2=",σ[b2index],"\nv3=",σ[b3index],"\nt1=",t1,"\nt2=",t2)
                

            z[1,1,1]+=(TimeDomainBEMInt.inttriangletriangleadjacent(τ[a1index],τ[a2index],τ[a3index],σ[b1index],σ[b2index],σ[b3index],t1,t2,qpc))/(4*π)
          elseif hits==1
                a2index,a3index=mod1(a1index+1,3),mod1(a1index+2,3)
                b2index,b3index=mod1(b1index+1,3),mod1(b1index+2,3)
                print("\na1=",τ[a1index],"\na2=",τ[a2index],"\na3=",τ[a3index],"\nb1=",σ[b1index],"\nb2=",σ[b2index],"\nb3=",σ[b3index],"\nt1=",t1,"\nt2=",t2)
                z[1,1,1]+=(TimeDomainBEMInt.intcommonvertex(τ[a1index],τ[a2index],τ[a3index],σ[b2index],σ[b3index],t1,t2,qpc))/(4*π)
                
          else
            #print("\np1=",τ[1],"\np2=",τ[2],"\np3=",τ[3],"\nv1=",σ[1],"\nv2=",σ[2],"\nv3=",σ[3],"\nt1=",t1,"\nt2=",t2)
             z[1,1,1]+=(inttriangletriangle(τ[1],τ[2],τ[3],σ[1],σ[2],σ[3],t1,t2,qpc))/(4*π)
          end
          #for the moment sos=1 but I will correct this
end


function momintegrals!(z, op::AcusticSingleLayerTDIO, g::LagrangeRefSpace{T,0,3}, f::LagrangeRefSpace{T,0,3}, t::MonomialBasis{T,0,1}, τ, σ, ι, qr::HybridZuccottiWiltonStrat) where T
        
    t1=ι[1]
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

        if hits==3 
            #print("\np1=",τ[1],"\np2=",τ[2],"\np3=",τ[3],"\nt1=",t1,"\nt2=",t2)
            entry=(TimeDomainBEMInt.intcoinctriangles(τ[1],τ[2],τ[3],t1,t2))/(4*π)
            if norm(entry)>10.0 || entry<-0.0001
                println("quadrature ",τ[1]," ",τ[2]," ",τ[3]," ", t1," ",t2," ",entry)
                XW = qr.outer_quad_points
                for p in 1 : length(XW)
                    x = XW[p].point
                    w = XW[p].weight
                    innerintegrals!(z, op, x, g, f, t, τ, σ, ι, qr, w)
                end 
            else
                z[1,1,1]+=max(entry,0.0)
            end
        elseif hits==2
           #= XW = qr.outer_quad_points
            for p in 1 : length(XW)
                x = XW[p].point
                w = XW[p].weight
                innerintegrals!(z, op, x, g, f, t, τ, σ, ι, qr, w)
            end=#
             #pay attention with double layer and index permutation or whatever case has n cross product
            if mod1(a1index +1,3)==a2index
                a3index=mod1(a1index-1,3)
            else
                a3index=mod1(a1index+1,3)
            end

            if mod1(b1index+1,3)==b2index
                b3index=mod1(b1index-1,3)
            else
                b3index=mod1(b1index+1,3)
            end
            #print("\np1=",τ[a1index],"\np2=",τ[a2index],"\np3=",τ[a3index],"\nv1=",σ[b1index],"\nv2=",σ[b2index],"\nv3=",σ[b3index],"\nt1=",t1,"\nt2=",t2)
                
                entry=(TimeDomainBEMInt.inttriangletriangleadjacent(τ[a1index],τ[a2index],τ[a3index],σ[b1index],σ[b2index],σ[b3index],t1,t2,qpc))/(4*π)
                if norm(entry)>10.0 || entry<-0.0001
                    println("quadrature ",τ[a1index]," ",τ[a2index]," ",τ[a3index]," ",σ[b1index]," ",σ[b2index]," ",σ[b3index], t1," ",t2," ",entry)
                    XW = qr.outer_quad_points
                    for p in 1 : length(XW)
                        x = XW[p].point
                        w = XW[p].weight
                        innerintegrals!(z, op, x, g, f, t, τ, σ, ι, qr, w)
                    end 
                else
                    z[1,1,1]+=max(0.0,entry)
                end
        elseif hits==1
            #=XW = qr.outer_quad_points
            for p in 1 : length(XW)
                x = XW[p].point
                w = XW[p].weight
                innerintegrals!(z, op, x, g, f, t, τ, σ, ι, qr, w)
            end=#
               a2index,a3index=mod1(a1index+1,3),mod1(a1index+2,3)
                b2index,b3index=mod1(b1index+1,3),mod1(b1index+2,3)
                #print("\np1=",τ[a1index],"\np2=",τ[a2index],"\np3=",τ[a3index],"\nv1=",σ[b1index],"\nv2=",σ[b2index],"\nv3=",σ[b3index],"\nt1=",t1,"\nt2=",t2)
                entry=(TimeDomainBEMInt.intcommonvertex(τ[a1index],τ[a2index],τ[a3index],σ[b2index],σ[b3index],t1,t2,qpc))/(4*π)
                if norm(entry)>10.0 || entry<-0.0001
                    println("quadrature ",τ[a1index]," ",τ[a2index]," ",τ[a3index]," ",σ[b2index]," ",σ[b3index], t1," ",t2," ",entry, " ",)
                    
                    XW = qr.outer_quad_points
                    for p in 1 : length(XW)
                        x = XW[p].point
                        w = XW[p].weight
                        innerintegrals!(z, op, x, g, f, t, τ, σ, ι, qr, w)
                    end 
                else
                    z[1,1,1]+=max(0.0,entry)
                end
        else
            entry=0.1#(inttriangletriangle(τ[1],τ[2],τ[3],σ[1],σ[2],σ[3],t1,t2,qpc))/(4*π)
            #print("\np1=",τ[1],"\np2=",τ[2],"\np3=",τ[3],"\nv1=",σ[1],"\nv2=",σ[2],"\nv3=",σ[3],"\nt1=",t1,"\nt2=",t2)
            if norm(entry)>0.0 || entry<-0.0001
                XW = qr.outer_quad_points
                for p in 1 : length(XW)
                    x = XW[p].point
                    w = XW[p].weight
                    innerintegrals!(z, op, x, g, f, t, τ, σ, ι, qr, w)
                end
                
            else
                z[1,1,1]+=max(0.0,entry)
            end
        end
      #for the moment sos=1 but I will correct this
end