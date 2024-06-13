function minmax1d(vertex,edge)
        T = eltype(τ[1])
        m = norm(vertex-edge[1])
        M = m
        s=edge[1]-edge[2]
        s/=norm(s)
        ev1=edge[1]-vertex
        x0=(edge[1]-dot(ev1,s)*s)
        a=(edge[2]-x0)*s
        b=(edge[1]-x0)*s
        if a<=0 && b>=0
           m=norm(vertex-x0)
           abs(a)<abs(b) && (M=norm(vertex-edge[2]))
            
        else
            for j in 1:length(edge)
                q = edge[j]
                d = norm(vertex-q)
                d < m && (m=d)
                d > M && (M=d)
            end
        end
        return m, M
end 

function rings1d(τ, σ, ΔR)
	m, M = minmaxdist(τ, σ)
	r0 = floor(Int, m/ΔR) + 1
	r1 = ceil(Int, M/ΔR+1)
	r0 : r1
end



function quaddata(op::AcusticSingleLayerTDIO, testrefs, trialrefs, timerefs,
    testels, trialels, timeels, quadstrat::AllAnalyticalQStrat)
    
    dimU=dimension(testels)
    dimV=dimension(trialels)
    #rigerenerare delta R

    if dimU+dimV==1
    #testelsboundary=skeleton(testels,dimU-1)
   # trialelsboundary=skeleton(trialels,dimV-1) 
    
        numnodes=length(testels)
        numedges=length(trialels)

        datavertexedge=Array{TimeDomainBEMInt.edgevertexgeo{T,P}, 2}(undef, numnodes, numedges)
        rings=Array{UnitRange{Int},2}(undef, numnodes, numedges)
        datarings=Array{Vector{Tuple{Int,Vector}},2}(undef,numnodes,numedges)#il type va bene

        #fill datarings with zeross !

        for p in 1:numnodes
            τ = chart(testels,p)#testels[p]
            for q in 1:numedges
                σ = chart(trialels,q)
                edgevertgeo=TimeDomainBEMInt.edgevertexinteraction(τ,σ[1],σ[2])
                datavertexedge[p,q]=edgevertgeo
                a,b=edgevertegeo.extint0[1],edgevertegeo.extint0[2]
                rngs=rings1d(τ,σ,ΔR)
                rings[p,q]=rngs
                datarings[p,q]=[0,[0.0,0.0]]
                for r in rngs
                    r > numfunctions(timebasisfunction) && continue #serve?
                    ι = ring(r,ΔR)

                    # compute interactions between reference shape functions
                    #fill!(z, 0)
                    rp=τ #se e un simplex ok se no va messo chart(τ,1).vertices credo
                    t2=ι[2]#needs a check
                    extint=TimeDomainBEMInt.edgevertexinteraction(t2,edgevertgeo)
                    push!(datarings[p,q],extint)
    
                    # qr = quadrule(op, U, V, W, p, τ, q, σ, r, ι, qd, quadstrat)
                    #momintegrals!(z, op, U, V, W, τ, σ, ι, qr)
                end
            end
        end

        return datavertexedge
    else
        return "devo ancora scrivere"
    end
end


nnodes=length(nodes)
edge1=chart(edges,p)
edge2=chart(edges,q)
cnnct=connectivity(edges,nodes)
vertind1=cnnct[1:nnodes,p].nzind
vertsgn1=cnnct[1:nnodes,p].nzval
vertind2=cnnct[1:nnodes,q].nzind
vertsgn2=cnnct[1:nnodes,q].nzval

if vertsgn1[1]==1
    a1,a2=edge1[1],edge1[2]
else
    a2,a1=edge1[1],edge1[2]
end

if vertsgn2[1]==1
    b1,b2=edge1[1],edge1[2]
else
    b2,b1=edge1[1],edge1[2]
end

geo1,rings1,datarings1=edgevertexgeo[vertind1[1],q],rings[vertind1[1],q],datarings[vertind1[1],q]
geo2,rings2,datarings2=edgevertexgeo[vertind1[2],q],rings[vertind1[2],q],datarings[vertind1[2],q]
geo3,rings3,datarings3=edgevertexgeo[vertind2[1],p],rings[vertind2[1],p],datarings[vertind2[1],p]
geo4,rings4,datarings4=edgevertexgeo[vertind2[2],p],rings[vertind2[2],p],datarings[vertind2[2],p]

geo=[geo1,geo2,geo3,geo4]
rings=[rings1,rings2,rings3,rings4]
datarings=[datarings1,datarings2,datarings3,datarings4]



function intlinelineglobal(a1,a2,b1,b2,geo,rings,datatimes,parcontrol,UB::Type{Val{N}}) where N
        
    #nedges=length(edges)
   
    
     #vertices=[a1,a2,a1′,a2′]
    vertices=[a1,a2,b1,b2]
    l12=norm(vertices[1]-vertices[2])
    l12′=norm(vertices[3]-vertices[4])
    

   #geo1=edgevertexinteraction(a1,a1′,a2′) 
    #geo2=edgevertexinteraction(a2,a1′,a2′)
    #geo3=edgevertexinteraction(a1′,a1,a2)
    #geo4=edgevertexinteraction(a2′,a1,a2)
    
    #datatime=Array{Tuple}(undef,2,4)?
    I = maketuple(eltype(a1), UB)
    K = maketuple(typeof(a1), UB)
    
    x=geo[3].tangent

    #z=cross(a12′,x)
    #J=norm(z)
    #if J blablabla
    #z /= J
    #h=dot(r22,z)
    hdir=cross(geo[1].tangent,x)
    n=hdir/norm(hdir)
        sgnn=[+1,-1,-1,+1]
        h=dot(a2-b2,n)
        sgnh=[+1,-1,+1,-1]
        angletot=0.0
        dminv=Vector{eltype(edge1[1])}(undef, 4)
        ξ=Vector{typeof(edge1[1])}(undef, 4)
        for j in 1:4
            dminv[j]=geo[j].dmin
            dmaxv[j]=geo[j].dmax 
            v=vertices[j]
            ξ[j]=v-n*h*sgnh[j]*sgnn[j] 
            angletot+=anglecontribution(ξ[j],sgnn[j]*n,geo[j])
        end
    if abs(angletot-2π)<100*eps(eltype(edge1[1]))
        dmin=abs(h)
    else
        dmin=min(dminv[1],dminv[2],dminv[3],dminv[4])
    end

    dmax=max(dmaxv[1],dmaxv[2],dmaxv[3],dmaxv[4])

    r0 = floor(Int, dmin/ΔR) + 1
	r1 = ceil(Int, dmax/ΔR+1) #recuperare deltaR
	ringtot = r0 : r1

    allint=Vector{typeof((I,K))}(undef,r1-r0+2)
    fill!(allint,(I,K))
    if norm(hdir) < (parcontrol[1])*eps(typeof(temp1))
        I=intparallelsegment(a1,a2,b1,b2,temp1,temp2)[1] #attenzione qui non compatibile con quello che stiamo scrivendo
    else
        n=hdir/norm(hdir)
        sgnn=[+1,-1,-1,+1]
        h=dot(a2-a2′,n)
        sgnh=[+1,-1,+1,-1]
        for j in 1:4  
            for i in ringtot[1]:(relrings[j][1]-1)
                        
                        P,Q  = arcsegcontribution(v,ξ[j],sgnn[j]*n,sgnh[j]*h,geo[j],0,[0,0],i*ΔR,UB)
                        allint[i-ringtot[1]+2][1] = add(allint[i-ringtot[1]+2][1],P)
                        allint[i-ringtot[1]+2][2] = add(allint[i-ringtot[1]+2][2],Q)
                        P,Q = arcsegcontribution(v,ξ[j],-sgnn[j]*n,sgnh[j]*h,geo[j],0,[0,0],(i-1)*ΔR,UB)
                        allint[i-ringtot[1]+2][1] = add(allint[i-ringtot[1]+2][1],P)
                        allint[i-ringtot[1]+2][2] = add(allint[i-ringtot[1]+2][2],Q)
            end
            for i in relrings[j]
                
                    #shall I put some check like i*deltaR > h
                        P, Q = arcsegcontribution(v,ξ[j],sgnn[j]*n,sgnh[j]*h,geo[j],datarings[j][i-relrings[j][1]+2][1],datarings[j][i-relrings[j][1]+2][2],i*ΔR,UB) #salviamo 
                        allint[i-ringtot[1]+2][1] = add(allint[i-ringtot[1]+2][1],P)
                        allint[i-ringtot[1]+2][2] = add(allint[i-ringtot[1]+2][2],Q)
                        P, Q = arcsegcontribution(v,ξ[j],-sgnn[j]*n,sgnh[j]*h,geo[j],datarings[j][i-relrings[j][1]+1][1],datarings[j][i-relrings[j][1]+1][2],(i-1)*ΔR,UB)                 
                        allint[i-ringtot[1]+2][1] = add(allint[i-ringtot[1]+2][1],P)
                        allint[i-ringtot[1]+2][2] = add(allint[i-ringtot[1]+2][2],Q)
            end #probabilmente va bene cosi anche con ceil(int,frac) invece di ceil(int,frac+1)
            allint[i-ringtot[1]+2][1]=multiply(allint[i-ringtot[1]+2][1],1/(l12′*l12*norm(hdir)))
            allint[i-ringtot[1]+2][2]=multiply(allint[i-ringtot[1]+2][2],1/(l12′*l12*norm(hdir))) 
            #=   for i in (relrings[j][2]+1):ringtot[j][2]
        
                    #shall I put some check like i*deltaR > h
                        saveP[i],saveQ[i]  = arcsegcontribution(v,ξ[j],sgnn[j]*n,sgnh[j]*h,geo[j],1,[a,b],i*ΔR,UB) #dadefinire a e b 
                        save,and,subtract = arcsegcontribution(v,ξ[j],sgnn[j]*n,sgnh[j]*h,geo[j],1,[a,b],(i-1)*ΔR,UB)
            end =# #sembra che non serva a causa di ceil(int,frac+1)!
        end
            #I+=(1/abs(r12′[2]*r12[1]))*(3*(temp2^2-temp1^2)*d[1]-2*(temp2^3-temp1^3)*d[2])
    end
    
    return ringtot,allint #missing buidgrad since it is not yet adapted for int line line
end

