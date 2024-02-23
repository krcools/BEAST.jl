using BEAST
using CompScienceMeshes
using LinearAlgebra
using AbstractTrees
using StaticArrays

function Z(κ;rows=[1,2,3,4],cols=[1,2,3,4],trace=true,dual=true,tr=1) # build to test with nxt in dual true
    G = BEAST.greenhh3d(wavenumber=κ)
    gradG = BEAST.∇(G)
    if trace
        ∇Gx =  BEAST.build_potential(gradG×B)
        Gn = BEAST.build_potential(G*(n*B))
        #Gnx = BEAST.build_potential(G*(n × B))
        ∇G = BEAST.build_potential(gradG*B)
        ∇Gdotn = BEAST.build_potential(nb ⋅ (gradG*B))
        ∇Gdot = BEAST.build_potential(B ⋅ gradG)
        
        Gr = BEAST.build_potential(G*B)
        ∇G∇B = BEAST.build_potential(gradG*BEAST.div(B))
        ∇Gxn = BEAST.build_potential(gradG×(n*B))
    else
        ∇Gx =  gradG×B
        Gn = G*(nb*B)
        #Gnx = G*(n × B)
        ∇G = gradG*B
        ∇Gdotn = nb ⋅ (gradG*B)
        ∇Gdot = B ⋅ gradG
        
        Gr = G*B
        ∇G∇B = gradG*BEAST.div(B)
        ∇Gxn = gradG×(nb*B)

    end
    
    if dual && trace
        int = -[γₜ(∇G∇B,tr)+κ^2*γₜ(Gr,tr) -γₜ(∇Gxn,tr) γₜ(∇Gx,tr) ZeroOperator();
        γₙ(∇Gx,tr) -γₙ(Gn,tr) γₙ(Gr,tr) -γₙ(∇G,tr);
        γₜ(∇Gx,tr) -γₜ(Gn,tr) γₜ(Gr,tr) -γₜ(∇G,tr);
        ZeroOperator() -τ(∇Gdotn,tr) τ(∇Gdot,tr)  κ^2*τ(Gr,tr)]
        
    elseif dual
        int = -[-nt×(nt×(∇G∇B))-κ^2*nt×(nt×(Gr)) nt×(nt×(∇Gxn)) -nt×(nt×(∇Gx)) ZeroOperator();
        nt ⋅(∇Gx) -nt ⋅(Gn) nt ⋅(Gr) -nt ⋅(∇G);
        -nt×(nt×(∇Gx)) nt×(nt×(Gn)) -nt×(nt×(Gr)) nt×(nt×(∇G));
        ZeroOperator() -(∇Gdotn) (∇Gdot)  κ^2*(Gr)]
        
    elseif trace
        int = -[γₛ(∇G∇B,tr)+κ^2*γₛ(Gr,tr) -γₛ(∇Gxn,tr) γₛ(∇Gx,tr) ZeroOperator();
        γₙ(∇Gx,tr) -γₙ(Gn,tr) γₙ(Gr,tr) -γₙ(∇G,tr);
        γₛ(∇Gx,tr) -γₛ(Gn,tr) γₛ(Gr,tr) -γₛ(∇G,tr);
        ZeroOperator() -τ(∇Gdotn,tr) τ(∇Gdot,tr)  κ^2*τ(Gr,tr)]
        
    else
        int = -[nt×(∇G∇B)+κ^2*nt×(Gr) -nt×(∇Gxn) nt×(∇Gx) ZeroOperator();
        nt ⋅(∇Gx) -nt ⋅(Gn) nt ⋅(Gr) -nt ⋅(∇G);
        nt×(∇Gx) -nt×(Gn) nt×(Gr) -nt×(∇G);
        ZeroOperator() -(∇Gdotn) (∇Gdot)  κ^2*(Gr)]
        
    end
    return BEAST.matrix_to_bilform(int[rows,cols])
end
function Zp(κ;rows=[1,2,3,4],cols=[1,2,3,4],trace=true,dual=true,tr=1) # build to test with nxt in dual true
    G = BEAST.greenhh3d(wavenumber=κ)
    gradG = BEAST.∇(G)
    tr *=1
    if trace
        ∇Gx =  BEAST.build_potential(gradG×B)
        Gn = BEAST.build_potential(G*(n*B))
        #Gnx = BEAST.build_potential(G*(n × B))
        ∇G = BEAST.build_potential(gradG*B)
        ∇Gdotn = BEAST.build_potential(nb ⋅ (gradG*B))
        ∇Gdot = BEAST.build_potential(B ⋅ gradG)
        
        Gr = BEAST.build_potential(G*B)
        ∇G∇B = BEAST.build_potential(gradG*BEAST.div(B))
        ∇Gxn = BEAST.build_potential(gradG×(n*B))
    else
        ∇Gx =  gradG×B
        Gn = G*(nb*B)
        #Gnx = G*(n × B)
        ∇G = gradG*B
        ∇Gdotn = nb ⋅ (gradG*B)
        ∇Gdot = B ⋅ gradG
        
        Gr = G*B
        ∇G∇B = gradG*BEAST.div(B)
        ∇Gxn = gradG×(nb*B)

    end
    
    if dual && trace
        int = -[γₜ(∇G∇B,tr)+κ^2*γₜ(Gr,tr) -γₜ(∇Gxn,tr) γₜ(∇Gx,tr) ZeroOperator();
        γₙ(∇Gx,tr) -γₙ(Gn,tr) γₙ(Gr,tr) -γₙ(∇G,tr);
        γₜ(∇Gx,tr) -γₜ(Gn,tr) γₜ(Gr,tr) -γₜ(∇G,tr);
        ZeroOperator() -τ(∇Gdotn,tr) τ(∇Gdot,tr)  κ^2*τ(Gr,tr)][[3,4,1,2],[3,4,1,2]]
        
    elseif dual
        int = -[-nt×(nt×(∇G∇B))-κ^2*nt×(nt×(Gr)) nt×(nt×(∇Gxn)) -nt×(nt×(∇Gx)) ZeroOperator();
        nt ⋅(∇Gx) -nt ⋅(Gn) nt ⋅(Gr) -nt ⋅(∇G);
        -nt×(nt×(∇Gx)) nt×(nt×(Gn)) -nt×(nt×(Gr)) nt×(nt×(∇G));
        ZeroOperator() -(∇Gdotn) (∇Gdot)  κ^2*(Gr)]
        
    elseif trace
        int = -[γₛ(∇G∇B,tr)+κ^2*γₛ(Gr,tr) -γₛ(∇Gxn,tr) γₛ(∇Gx,tr) ZeroOperator();
        γₙ(∇Gx,tr) -γₙ(Gn,tr) γₙ(Gr,tr) -γₙ(∇G,tr);
        γₛ(∇Gx,tr) -γₛ(Gn,tr) γₛ(Gr,tr) -γₛ(∇G,tr);
        ZeroOperator() -τ(∇Gdotn,tr) τ(∇Gdot,tr)  κ^2*τ(Gr,tr)][[3,4,1,2],[3,4,1,2]]
        
    else
        int = -[nt×(∇G∇B)+κ^2*nt×(Gr) -nt×(∇Gxn) nt×(∇Gx) ZeroOperator();
        nt ⋅(∇Gx) -nt ⋅(Gn) nt ⋅(Gr) -nt ⋅(∇G);
        nt×(∇Gx) -nt×(Gn) nt×(Gr) -nt×(∇G);
        ZeroOperator() -(∇Gdotn) (∇Gdot)  κ^2*(Gr)][[3,4,1,2],[3,4,1,2]]
        
    end
    return BEAST.matrix_to_bilform(int[rows,cols])
end
function excitation(A,curlA,divA;rows=[1,2,3,4],dual=true)
    if dual
    out = [((n × FunctionExcitation{ComplexF64}(curlA))×n),
    BEAST.NDotTrace{ComplexF64}(FunctionExcitation{ComplexF64}(A)),
    ((n × FunctionExcitation{ComplexF64}(A))×n),
    FunctionExcitation{ComplexF64}(divA)]
    else
        out = [((n × FunctionExcitation{ComplexF64}(curlA))),
        BEAST.NDotTrace{ComplexF64}(FunctionExcitation{ComplexF64}(A)),
        ((n × FunctionExcitation{ComplexF64}(A))),
        FunctionExcitation{ComplexF64}(divA)]
    end
    return BEAST.array_to_linform(out[rows])
end
function solve(Z,b,X;strat,kwargs...)
    if strat ==:LU
        return BEAST.solve(b,Z,X;kwargs...)
    elseif strat == :GMRES
        return BEAST.gmres_ch(b,Z,X;kwargs...)
    end
    @error "no existing strategy given, use :LU or :GMRES"
end

parent(ind,tree) = tree[ind]
children(ind,tree) = findall(==(ind),tree)
### physical constants
ϵ0 = 8.854e-12
μ0 = 4π*1e-7
ω = 10.0^8*2*pi
κ0 = ω*sqrt(ϵ0*μ0)

#### define meshes
hh = 0.3
Γ1 = meshcuboid(1.1, 1.0, 1.1, hh)
Γ2 =  (@SVector [-1.0,0.0,0.0]) + BEAST.TraceMesh(-Mesh([point(-x,y,z) for (x,y,z) in vertices(Γ1)], deepcopy(cells(Γ1))))
Γ3 = (@SVector [0.0,0.0,-1.0]) + BEAST.TraceMesh(-translate(Mesh([point(x,y,-z-1.1) for (x,y,z) in vertices(Γ1)], deepcopy(cells(Γ1))),[0.0,0.0,0.0]))

Γ = [Γ1,Γ2,Γ3]
Tree = [0,0,0] # give index in which material is

HOM = [1,2] #indices of homogeneous domains (without free space)
HomPars = Dict(0=>(ϵ0,μ0),1=>(ϵ0,μ0*3),2=>(ϵ0*2,μ0*2))#,3=>(ϵ0*3,μ0))#index --> (eps,mu)
κ = Dict(i => ω*sqrt(prod(HomPars[i])) for i in HOM)
κ[0]=κ0
EFIE = [] #indices of pec domains modeled with efie
MFIE = [] #indices of pec domains modeled with mfie
CFIE = [3] #indices of pec domains modeled with cfie
α = Dict(3=>0.2)#index --> alpha

#### Incident field
A(x) = x[1]*[0.0,0.0,1.0]*sqrt(ϵ0*μ0)*exp(-1.0im*κ0*x[3])
curlA(x) = -[0.0,1.0,0.0]*sqrt(ϵ0*μ0)*exp(-1.0im*κ0 *x[3])
divA(x) = -im*κ0 *x[1]*sqrt(ϵ0*μ0)*exp(-1.0im*κ0 *x[3])
graddivA(x) = -im*κ0 *sqrt(ϵ0*μ0)*exp(-1.0im*κ0 *x[3])*[1.0,0.0,-im*κ0*x[1]]
#### spaces
Xdi(Γ) = BEAST.DirectProductSpace([raviartthomas(Γ),lagrangec0d1(Γ)])
Xni(Γ) = BEAST.DirectProductSpace([raviartthomas(Γ),lagrangecxd0(Γ)])
Ydi(Γ) = BEAST.DirectProductSpace([BEAST.buffachristiansen(Γ),duallagrangec0d1(Γ)])
Yni(Γ) = BEAST.DirectProductSpace([BEAST.buffachristiansen(Γ),duallagrangecxd0(Γ)])
Xi(Γ) = BEAST.DirectProductSpace([raviartthomas(Γ),lagrangec0d1(Γ),raviartthomas(Γ),lagrangecxd0(Γ)])
Yi(Γ) = BEAST.DirectProductSpace([BEAST.buffachristiansen(Γ),duallagrangecxd0(Γ),BEAST.buffachristiansen(Γ),duallagrangec0d1(Γ)])# this one maybe omgedraaid

#### incident field trace
Xin = excitation(A,curlA,divA)
Xind = excitation(A,curlA,divA;rows = [3,4])
Xinn = excitation(A,curlA,divA; rows = [1,2])
#### incident field trace non dual
Xin_notdual = excitation(A,curlA,divA;dual=false)
Xind_notdual = excitation(A,curlA,divA;rows = [3,4],dual=false)
Xinn_notdual = excitation(A,curlA,divA; rows = [1,2],dual=false)
#@error "check if correct order of traces"
##### identities
id() = BEAST.matrix_to_bilform(diagm([-NCross(),Identity(),-NCross(),Identity()]))
idN() = BEAST.matrix_to_bilform(diagm([-NCross(),Identity()]))
idnotdual() = BEAST.matrix_to_bilform(diagm([Identity(),Identity(),Identity(),Identity()]))
idNnotdual() = BEAST.matrix_to_bilform(diagm([Identity(),Identity()]))
##### from here do nothing anymore.
N = length(Γ)
Q = Dict(i=>diagm([1,1,HomPars[parent(i,Tree)][2]/HomPars[i][2],HomPars[i][1]/HomPars[parent(i,Tree)][1]]) for i in HOM)
Qp = Dict(i=>diagm([HomPars[parent(i,Tree)][2]/HomPars[i][2],HomPars[i][1]/HomPars[parent(i,Tree)][1],1,1]) for i in HOM)
t = BEAST.hilbertspace(:t, length(Γ))
b = BEAST.hilbertspace(:b, length(Γ))

##### define space
perm = sortperm([HOM...,EFIE...,CFIE...,MFIE...])
X = [BEAST.DirectProductSpace[Xi(Γ[i]) for i in HOM];BEAST.DirectProductSpace[Xdi(Γ[i]) for i in [EFIE...,CFIE...,MFIE...]]][perm]
Y = [BEAST.DirectProductSpace[Yi(Γ[i]) for i in HOM];BEAST.DirectProductSpace[Yni(Γ[i]) for i in [EFIE...,CFIE...,MFIE...]]][perm]
##### define equation

eqs1 = BEAST.Equation[(Qp[i]*Z(κ[i];tr=-1)*(Q[i]^-1))[t[i],b[i]] +
        -sum(BEAST.BilForm[(Qp[i]*Z(κ[i];tr=-1))[t[i],b[j]] for j in HOM ∩ children(i,Tree)]) +
        -sum(BEAST.BilForm[(Qp[i]*Z(κ[i];cols=[3,4],tr=-1))[t[i],b[j]] for j in [EFIE...,CFIE...,MFIE...] ∩ children(i,Tree)]) ==0 for i in HOM]

# deqs1 = BEAST.discretise.(eqs1, Ref.((t .∈ X))..., Ref.((b .∈ X))...)


##### define equation
eqs2hom = begin BEAST.Equation[-sum(BEAST.BilForm[Z(κ0)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
           -sum(BEAST.BilForm[Z(κ0;cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==-Xin[t[ci]]
           for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ HOM] end
eqs2efie = begin BEAST.Equation[-sum(BEAST.BilForm[Z(κ0;rows=[3,4])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -sum(BEAST.BilForm[Z(κ0;rows=[3,4],cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==-Xind[t[ci]]
            for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ EFIE] end
eqs2mfie = begin BEAST.Equation[-sum(BEAST.BilForm[Z(κ0;rows=[1,2])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -sum(BEAST.BilForm[Z(κ0;rows=[1,2],cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) +
            -idN()[t[ci],b[ci]] ==-Xinn[t[ci]]
            for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ MFIE]end
eqs2cefie =begin  BEAST.Equation[-α[ci]*sum(BEAST.BilForm[Z(κ0;rows=[3,4])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -α[ci]*sum(BEAST.BilForm[Z(κ0;rows=[3,4],cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==
            -α[ci]*Xind[t[ci]] for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end
eqs2cmfie = begin BEAST.Equation[-(1-α[ci])* sum(BEAST.BilForm[Z(κ0;rows=[1,2])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -(1-α[ci])*sum(BEAST.BilForm[Z(κ0;rows=[1,2],cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) +
            -(1-α[ci])*idN()[t[ci],b[ci]] == -(1-α[ci])*Xinn[t[ci]]
            for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end

# deqs2hom = BEAST.discretise(eqs2hom, (t.∈X),(b.∈X))    
# deqs2efie = BEAST.discretise(eqs2efie, (t.∈X),(b.∈X))    
# deqs2mfie = BEAST.discretise(eqs2mfie, (t.∈Y),(b.∈X))    
# deqs2cefie = BEAST.discretise(eqs2cefie, (t.∈X),(b.∈X))    
# deqs2cmfie = BEAST.discretise(eqs2cmfie, (t.∈Y),(b.∈X))    


##### define equation
eqs3hom = begin BEAST.Equation[(Z(κ[i])*(Q[i]^-1))[t[ci],b[i]]+
            -sum(BEAST.BilForm[Z(κ[i])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
           -sum(BEAST.BilForm[Z(κ[i];cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) == 0
           for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ HOM] end
eqs3efie = begin BEAST.Equation[(Z(κ[i];rows=[3,4])*(Q[i]^-1))[t[ci],b[i]]+
            -sum(BEAST.BilForm[Z(κ[i];rows=[3,4])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -sum(BEAST.BilForm[Z(κ[i];rows=[3,4],cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) == 0
            for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ EFIE] end
eqs3mfie = begin BEAST.Equation[(Z(κ[i];rows=[1,2])*(Q[i]^-1))[t[ci],b[i]]+
            -sum(BEAST.BilForm[Z(κ[i];rows=[1,2])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -sum(BEAST.BilForm[Z(κ[i];rows=[1,2],cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) +
            -idN()[t[ci],b[ci]] ==0
            for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ MFIE]end
eqs3cefie = begin BEAST.Equation[α[ci]*(Z(κ[i];rows=[3,4])*(Q[i]^-1))[t[ci],b[i]]+
            -α[ci]*sum(BEAST.BilForm[Z(κ[i];rows=[3,4])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -α[ci]*sum(BEAST.BilForm[Z(κ[i];rows=[3,4],cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==0
            for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end
eqs3cmfie = begin BEAST.Equation[(1-α[ci])*(Z(κ[i];rows=[1,2])*(Q[i]^-1))[t[ci],b[i]]+
            -(1-α[ci])* sum(BEAST.BilForm[Z(κ[i];rows=[1,2])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -(1-α[ci])*sum(BEAST.BilForm[Z(κ[i];rows=[1,2],cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) +
            -(1-α[ci])*idN()[t[ci],b[ci]] == 0
            for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end

# deqs3hom = BEAST.discretise.(eqs3hom, Ref.((t.∈X))...,Ref.((b.∈X))...)    
# deqs3efie = BEAST.discretise.(eqs3efie, (t.∈X),(b.∈X))    
# deqs3mfie = BEAST.discretise.(eqs3mfie, (t.∈Y),(b.∈X))    
# deqs3cefie = BEAST.discretise.(eqs3cefie, (t.∈X),(b.∈X))    
# deqs3cmfie = BEAST.discretise.(eqs3cmfie, (t.∈Y),(b.∈X))    

## sum the equations in the two parts in a pmchwt fassion,

symeq = -sum(eqs1)+sum(eqs2cefie)+sum(eqs2efie)+sum(eqs2hom)+sum(eqs3cefie)+sum(eqs3efie)+sum(eqs3hom)
asymeq = sum(eqs2cmfie)+sum(eqs2mfie)+sum(eqs3cmfie)+sum(eqs3mfie)

symfilled = typeof(symeq) <: BEAST.Equation
asymfilled = typeof(asymeq) <: BEAST.Equation


symfilled && (Dsymeq = BEAST.discretise(symeq, (t.∈X)..., (b.∈X)...))
asymfilled && (Dasymeq = BEAST.discretise(asymeq, (t.∈Y)..., (b.∈X)...))
#assemble system
Zs,Za,bs,ba = 0,0,0,0
asymfilled && ((Za,ba,xx,yy) = assemble(Dasymeq))
symfilled && ((Zs,bs,xx,yy) = assemble(Dsymeq))

#unprec system

Zunprec = Za*I+Zs*I
bunprec = ba .+ bs 

solvestrat = :GMRES
solveargs = Dict(
    :maxiter => 20000,
    :restart => 20000,
    :tol => real(sqrt(eps())))
outnew = solve(Zunprec,bunprec,yy;strat=solvestrat,solveargs...)


#generate_plot_field



#preconditioner




##### define equation

peqs1 = BEAST.Equation[(Q[i]*Zp(κ[i];tr=-1)*(Qp[i]^-1))[t[i],b[i]] +
        -sum(BEAST.BilForm[(Q[i]*Zp(κ[i];tr=-1))[t[i],b[j]] for j in HOM ∩ children(i,Tree)]) +
        -sum(BEAST.BilForm[(Q[i]*Zp(κ[i];cols=[3,4],tr=-1))[t[i],b[j]] for j in [EFIE...,CFIE...,MFIE...] ∩ children(i,Tree)]) ==0 for i in HOM]
pdeqs1 = BEAST.Equation[(Q[i]*Zp(κ[i];tr=-1)*(Qp[i]^-1))[t[i],b[i]] ==0 for i in HOM]
pnuleqs1 = BEAST.Equation[(Q[i]*Zp(κ[i];tr=0)*(Qp[i]^-1))[t[i],b[i]] +
-sum(BEAST.BilForm[(Q[i]*Zp(κ[i];tr=0))[t[i],b[j]] for j in HOM ∩ children(i,Tree)]) +
-sum(BEAST.BilForm[(Q[i]*Zp(κ[i];cols=[3,4],tr=0))[t[i],b[j]] for j in [EFIE...,CFIE...,MFIE...] ∩ children(i,Tree)]) ==0 for i in HOM]

# deqs1 = BEAST.discretise.(eqs1, Ref.((t .∈ X))..., Ref.((b .∈ X))...)

##### define equation
peqs2hom = begin BEAST.Equation[-sum(BEAST.BilForm[Zp(κ0)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
           -sum(BEAST.BilForm[Zp(κ0;cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==-Xin[t[ci]]
           for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ HOM] end

pdeqs2hom = begin BEAST.Equation[-sum(BEAST.BilForm[Zp(κ0)[t[ci],b[ci]] for j in [HOM ∩ [ci]]])+
        -sum(BEAST.BilForm[Zp(κ0;cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ [ci]]) ==-Xin[t[ci]]
        for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ HOM] end


peqs2efie = begin BEAST.Equation[-sum(BEAST.BilForm[Zp(κ0;rows=[3,4])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
        -sum(BEAST.BilForm[Zp(κ0;rows=[3,4],cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==-Xind[t[ci]]
        for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ EFIE] end
# peqs2mfie = begin BEAST.Equation[-sum(BEAST.BilForm[Zp(κ0;rows=[1,2])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
#             -sum(BEAST.BilForm[Zp(κ0;rows=[1,2],cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) +
#             -idN()[t[ci],b[ci]] ==-Xinn[t[ci]]
#             for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ MFIE]end
peqs2cefie =begin  BEAST.Equation[-α[ci]*sum(BEAST.BilForm[Zp(κ0;rows=[3,4])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -α[ci]*sum(BEAST.BilForm[Zp(κ0;rows=[3,4],cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==
            -α[ci]*Xind[t[ci]] for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end
# peqs2cmfie = begin BEAST.Equation[-(1-α[ci])* sum(BEAST.BilForm[Zp(κ0;rows=[1,2])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
#             -(1-α[ci])*sum(BEAST.BilForm[Zp(κ0;rows=[1,2],cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) +
#             -(1-α[ci])*idN()[t[ci],b[ci]] == -(1-α[ci])*Xinn[t[ci]]
#             for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end

# deqs2hom = BEAST.discretise(eqs2hom, (t.∈X),(b.∈X))    
# deqs2efie = BEAST.discretise(eqs2efie, (t.∈X),(b.∈X))    
# deqs2mfie = BEAST.discretise(eqs2mfie, (t.∈Y),(b.∈X))    
# deqs2cefie = BEAST.discretise(eqs2cefie, (t.∈X),(b.∈X))    
# deqs2cmfie = BEAST.discretise(eqs2cmfie, (t.∈Y),(b.∈X))    


##### define equation
peqs3hom = begin BEAST.Equation[(Zp(κ[i])*(Qp[i]^-1))[t[ci],b[i]]+
            -sum(BEAST.BilForm[Zp(κ[i])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
           -sum(BEAST.BilForm[Zp(κ[i];cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) == 0
           for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ HOM] end
peqs3efie = begin BEAST.Equation[(Zp(κ[i];rows=[3,4])*(Qp[i]^-1))[t[ci],b[i]]+
            -sum(BEAST.BilForm[Zp(κ[i];rows=[3,4])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -sum(BEAST.BilForm[Zp(κ[i];rows=[3,4],cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) == 0
            for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ EFIE] end
# peqs3mfie = begin BEAST.Equation[(Zp(κ[i];rows=[1,2])*(Q[i]^-1))[t[ci],b[i]]+
#             -sum(BEAST.BilForm[Zp(κ[i];rows=[1,2])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
#             -sum(BEAST.BilForm[Zp(κ[i];rows=[1,2],cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) +
#             -idN()[t[ci],b[ci]] ==0
#             for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ MFIE]end
peqs3cefie = begin BEAST.Equation[α[ci]*(Zp(κ[i];rows=[3,4])*(Qp[i]^-1))[t[ci],b[i]]+
            -α[ci]*sum(BEAST.BilForm[Zp(κ[i];rows=[3,4])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -α[ci]*sum(BEAST.BilForm[Zp(κ[i];rows=[3,4],cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==0
            for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end
# peqs3cmfie = begin BEAST.Equation[(1-α[ci])*(Zp(κ[i];rows=[1,2])*(Q[i]^-1))[t[ci],b[i]]+
#             -(1-α[ci])* sum(BEAST.BilForm[Zp(κ[i];rows=[1,2])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
#             -(1-α[ci])*sum(BEAST.BilForm[Zp(κ[i];rows=[1,2],cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) +
#             -(1-α[ci])*idN()[t[ci],b[ci]] == 0
#             for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end


#discretise system
psymeq = -sum(peqs1)+sum(peqs2cefie)+sum(peqs2efie)+sum(peqs2hom)+sum(peqs3cefie)+sum(peqs3efie)+sum(peqs3hom)
pdsymeq = -sum(pdeqs1)+sum(pdeqs2hom)
pnulsymeq = sum(pnuleqs1)
Dsymeq = BEAST.discretise(psymeq, (t.∈Y)..., (b.∈Y)...)
# quadstrat
qs(::BEAST.LocalOperator, a, b) = BEAST.SingleNumQStrat(3)
qs(op::BEAST.ComposedOperatorIntegral,testspace,trialspace) = BEAST.DoubleNumSauterQstrat(5,5,5,5,2,2) 
qs(op::BEAST.ComposedOperatorLocal,testspace,trialpsace) = BEAST.SingleNumQStrat(3)

Pyy = assemble(Dsymeq.equation.lhs,BEAST.DirectProductSpace(Y),BEAST.DirectProductSpace(Y);quadstratfunction=qs)


# construct gramm matrix
idtest() = BEAST.matrix_to_bilform(diagm([-NCross(),Identity(),-NCross(),Identity()]))

idNtest() = BEAST.matrix_to_bilform(diagm([-NCross(),Identity()]))

gramm = sum([idtest()[t[i],b[i]] for i in HOM])+sum([idNtest()[t[i],b[i]] for i in [EFIE...,CFIE...,MFIE...]])
Gxy = assemble(gramm,BEAST.DirectProductSpace(X),BEAST.DirectProductSpace(Y))
Ginvyx = BEAST.GMRESSolver(Gxy;abstol=1e-15, maxiter=1_000, restart=1_000, verbose=false)
Ginvxy = BEAST.GMRESSolver(Gxy';abstol=1e-15, maxiter=1_000, restart=1_000, verbose=false)

Zprec = Za*I+Ginvxy*Pyy*Ginvyx*Zs*I
bprec = ba .+ Ginvxy*Pyy*Ginvyx*bs 

solvestrat = :GMRES
solveargs = Dict(
    :maxiter => 20000,
    :restart => 20000,
    :tol => real(sqrt(eps())))
outnewprec = solve(Zprec,bprec,yy;strat=solvestrat,solveargs...)

u = outnewprec[1]

###### exit dict: mesh, number of iterations, inputdata, cond number both, eigenvalue distribution both cases


##### post processing


##### nearfield
using BlockArrays
function nearfield_B(u,X,κ,points;sign=1,HOM=true,Q=diagm([1.0,1.0,1.0,1.0]))#nearfield of single volume
    G = BEAST.greenhh3d(wavenumber=κ)
    gradG = BEAST.∇(G)
    ∇Gx =  gradG×B
    Gn = G*(nb*B)
    #Gnx = G*(n × B)
    ∇G = gradG*B
    ∇Gdotn = nb ⋅ (gradG*B)
    ∇Gdot = B ⋅ gradG
    
    Gr = G*B
    ∇G∇B = gradG*BEAST.div(B)
    ∇Gxn = gradG×(nb*B)
    

    if HOM
        @hilbertspace T D R N
        B = -Q[1,1]*potential(∇G∇B, points, u[T], X.factors[1])-Q[1,1]*potential(κ^2*Gr, points, u[T], X.factors[1])
        B .+= Q[2,2]*potential(∇Gxn, points, u[D], X.factors[2])
        B .+= -Q[3,3]*potential(∇Gx, points, u[R], X.factors[3])
        
    else
        @hilbertspace R N
        B = -potential(∇Gx, points, u[R], X.factors[1])
        
    end

    return sign*B
end

function nearfield_A(u,X,κ,points;sign=1,HOM=true,Q=diagm([1.0,1.0,1.0,1.0]))#nearfield of single volume
    G = BEAST.greenhh3d(wavenumber=κ)
    gradG = BEAST.∇(G)
    ∇Gx =  gradG×B
    Gn = G*(nb*B)
    #Gnx = G*(n × B)
    ∇G = gradG*B
    ∇Gdotn = nb ⋅ (gradG*B)
    ∇Gdot = B ⋅ gradG
    
    Gr = G*B
    ∇G∇B = gradG*BEAST.div(B)
    ∇Gxn = gradG×(nb*B)
    

    if HOM
        @hilbertspace T D R N
        A = -Q[1,1]*potential(∇Gx, points, u[T], X.factors[1])
        A .+= Q[2,2]*potential(Gn, points, u[D], X.factors[2])
        A .+= -Q[3,3]*potential(Gr, points, u[R], X.factors[3])
        A .+= Q[4,4]*potential(∇G, points, u[N], X.factors[4])
    else
        @hilbertspace R N
        A = -potential(Gr, points, u[R], X.factors[1])
        A .+= potential(∇G, points, u[N], X.factors[2])
    end

    return sign*A
end
function nearfield_phi(u,X,κ,ω,points;sign=1,HOM=true,Q=diagm([1.0,1.0,1.0,1.0]))#nearfield of single volume
    G = BEAST.greenhh3d(wavenumber=κ)
    gradG = BEAST.∇(G)
    ∇Gx =  gradG×B
    Gn = G*(nb*B)
    #Gnx = G*(n × B)
    ∇G = gradG*B
    ∇Gdotn = nb ⋅ (gradG*B)
    ∇Gdot = B ⋅ gradG
    
    Gr = G*B
    ∇G∇B = gradG*BEAST.div(B)
    ∇Gxn = gradG×(nb*B)
    

    if HOM
        @hilbertspace T D R N
        A = Q[2,2]*potential(∇Gdotn, points, u[D], X.factors[2])
        A .+= -Q[3,3]*potential(∇Gdot, points, u[R], X.factors[3])
        A .+= -κ^2* Q[4,4]*potential(Gr, points, u[N], X.factors[4])
    else
        @hilbertspace R N
        A = -potential(∇Gdot, points, u[R], X.factors[1])
        A .+= -κ^2* potential(Gr, points, u[N], X.factors[2])
    end

    return sign*A/κ^2*im*ω
end
function nearfield_E(u,X,κ,ω,points;sign=1,HOM=true,Q=diagm([1.0,1.0,1.0,1.0]))#nearfield of single volume
    G = BEAST.greenhh3d(wavenumber=κ)
    gradG = BEAST.∇(G)
    graddivG = BEAST.graddiv(G)
    ∇Gx =  gradG×B
    Gn = G*(nb*B)
    #Gnx = G*(n × B)
    ∇G = gradG*B
    ∇Gdotn = nb ⋅ (gradG*B)
    graddivGdotn = graddivG*(nb*B)
    graddivGdot = graddivG*B
    ∇Gdot = B ⋅ gradG
    
    Gr = G*B
    ∇G∇B = gradG*BEAST.div(B)
    ∇Gxn = gradG×(nb*B)
    

    if HOM
        @hilbertspace T D R N
        p = Q[2,2]*potential(graddivGdotn, points, u[D], X.factors[2])
        p .+= -Q[3,3]*potential(graddivGdot, points, u[R], X.factors[3])
        p .+= -κ^2* Q[4,4]*potential(∇G, points, u[N], X.factors[4])
    else
        @hilbertspace R N
        p = -potential(graddivGdot, points, u[R], X.factors[1])
        p .+= -κ^2* potential(∇G, points, u[N], X.factors[2])
    end
    if HOM
        @hilbertspace T D R N
        A = -Q[1,1]*potential(∇Gx, points, u[T], X.factors[1])
        A .+= Q[2,2]*potential(Gn, points, u[D], X.factors[2])
        A .+= -Q[3,3]*potential(Gr, points, u[R], X.factors[3])
        A .+= Q[4,4]*potential(∇G, points, u[N], X.factors[4])
    else
        @hilbertspace R N
        A = -potential(Gr, points, u[R], X.factors[1])
        A .+= potential(∇G, points, u[N], X.factors[2])
    end
    return -sign*A/κ^2*im*ω-sign*im*ω*A
end
Xs = range(-2.0,2.0,length=150)
Zs = range(-1.5,2.5,length=100)
pts = [point(x,0.5,z) for z in Zs, x in Xs]

A1 = [nearfield_A(u[Block(i)],X[i],κ[parent(i,Tree)],pts) for i in HOM]
A2 = [nearfield_A(u[Block(i)],X[i],κ[i],pts;sign=-1,Q=Q[i]^-1) for i in HOM]
A3 = [nearfield_A(u[Block(i)],X[i],κ[parent(i,Tree)],pts;HOM=false) for i in [EFIE...,MFIE...,CFIE...]]
Atot = sum(A1)+sum(A2)+sum(A3)+A.(pts)

B1 = [nearfield_B(u[Block(i)],X[i],κ[parent(i,Tree)],pts) for i in HOM]
B2 = [nearfield_B(u[Block(i)],X[i],κ[i],pts;sign=-1,Q=Q[i]^-1) for i in HOM]
B3 = [nearfield_B(u[Block(i)],X[i],κ[parent(i,Tree)],pts;HOM=false) for i in [EFIE...,MFIE...,CFIE...]]
Btot = sum(B1)+sum(B2)+sum(B3)+curlA.(pts)

H1 = [1/HomPars[parent(i,Tree)][2]*nearfield_B(u[Block(i)],X[i],κ[parent(i,Tree)],pts) for i in HOM]
H2 = [1/HomPars[i][2]*nearfield_B(u[Block(i)],X[i],κ[i],pts;sign=-1,Q=Q[i]^-1) for i in HOM]
H3 = [1/HomPars[parent(i,Tree)][2]*nearfield_B(u[Block(i)],X[i],κ[parent(i,Tree)],pts;HOM=false) for i in [EFIE...,MFIE...,CFIE...]]
Htot = sum(H1)+sum(H2)+sum(H3)+1/μ0*curlA.(pts)

E1 = [nearfield_E(u[Block(i)],X[i],κ[parent(i,Tree)],ω,pts) for i in HOM]
E2 = [nearfield_E(u[Block(i)],X[i],κ[i],ω,pts;sign=-1,Q=Q[i]^-1) for i in HOM]
E3 = [nearfield_E(u[Block(i)],X[i],κ[parent(i,Tree)],ω,pts;HOM=false) for i in [EFIE...,MFIE...,CFIE...]]
Etot = sum(E1)+sum(E2)+sum(E3)-im*ω*A.(pts)-graddivA.(pts)/κ^2*im*ω

##### complement error, different from aps paper, traces from inside are compared, for the reset same formula as in aps paper, for pec all are treated with mfie, as efie yields zero traces in denominator.






##### define equation
ceqs2hom = begin BEAST.Equation[-sum(BEAST.BilForm[Z(κ0;tr=-1,dual=false)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
           -sum(BEAST.BilForm[Z(κ0;cols=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==-Xin_notdual[t[ci]]
           for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ HOM] end
# ceqs2efie = begin BEAST.Equation[-sum(BEAST.BilForm[Z(κ0;rows=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
#             -sum(BEAST.BilForm[Z(κ0;rows=[3,4],cols=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==-Xind_notdual[t[ci]]
#             for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ EFIE] end
ceqs2mfie = begin BEAST.Equation[-sum(BEAST.BilForm[Z(κ0;rows=[1,2],tr=-1,dual=false)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -sum(BEAST.BilForm[Z(κ0;rows=[1,2],cols=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) +
             ==-Xinn_notdual[t[ci]]
            for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ [MFIE...,EFIE...]]end
# ceqs2cefie =begin  BEAST.Equation[-α[ci]*sum(BEAST.BilForm[Z(κ0;rows=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
#             -α[ci]*sum(BEAST.BilForm[Z(κ0;rows=[3,4],cols=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==
#             -α[ci]*Xind_notdual[t[ci]] for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end
ceqs2cmfie = begin BEAST.Equation[-(1-α[ci])* sum(BEAST.BilForm[Z(κ0;rows=[1,2],tr=-1,dual=false)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -(1-α[ci])*sum(BEAST.BilForm[Z(κ0;rows=[1,2],cols=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) +
            == -(1-α[ci])*Xinn_notdual[t[ci]]
            for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end

# deqs2hom = BEAST.discretise(eqs2hom, (t.∈X),(b.∈X))    
# deqs2efie = BEAST.discretise(eqs2efie, (t.∈X),(b.∈X))    
# deqs2mfie = BEAST.discretise(eqs2mfie, (t.∈Y),(b.∈X))    
# deqs2cefie = BEAST.discretise(eqs2cefie, (t.∈X),(b.∈X))    
# deqs2cmfie = BEAST.discretise(eqs2cmfie, (t.∈Y),(b.∈X))    


##### define equation
ceqs3hom = begin BEAST.Equation[(Z(κ[i];tr=-1,dual=false)*(Q[i]^-1))[t[ci],b[i]]+
            -sum(BEAST.BilForm[Z(κ[i];tr=-1,dual=false)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
           -sum(BEAST.BilForm[Z(κ[i];cols=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) == 0
           for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ HOM] end
# ceqs3efie = begin BEAST.Equation[(Z(κ[i];rows=[3,4],tr=-1,dual=false)*(Q[i]^-1))[t[ci],b[i]]+
#             -sum(BEAST.BilForm[Z(κ[i];rows=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
#             -sum(BEAST.BilForm[Z(κ[i];rows=[3,4],cols=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) == 0
#             for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ EFIE] end
ceqs3mfie = begin BEAST.Equation[(Z(κ[i];rows=[1,2],tr=-1,dual=false)*(Q[i]^-1))[t[ci],b[i]]+
            -sum(BEAST.BilForm[Z(κ[i];rows=[1,2],tr=-1,dual=false)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -sum(BEAST.BilForm[Z(κ[i];rows=[1,2],cols=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)])  ==0
            for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ [MFIE...,EFIE...]]end
# ceqs3cefie = begin BEAST.Equation[α[ci]*(Z(κ[i];rows=[3,4],tr=-1,dual=false)*(Q[i]^-1))[t[ci],b[i]]+
#             -α[ci]*sum(BEAST.BilForm[Z(κ[i];rows=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
#             -α[ci]*sum(BEAST.BilForm[Z(κ[i];rows=[3,4],cols=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==0
#             for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end
ceqs3cmfie = begin BEAST.Equation[(1-α[ci])*(Z(κ[i];rows=[1,2],tr=-1,dual=false)*(Q[i]^-1))[t[ci],b[i]]+
            -(1-α[ci])* sum(BEAST.BilForm[Z(κ[i];rows=[1,2],tr=-1,dual=false)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -(1-α[ci])*sum(BEAST.BilForm[Z(κ[i];rows=[1,2],cols=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) == 0
            for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end



ceq = sum(ceqs2hom)+sum(ceqs3hom)
cmeq = sum(ceqs2cmfie)+sum(ceqs2mfie)+sum(ceqs3cmfie)+sum(ceqs3mfie)

Dseq = BEAST.discretise(ceq, (t.∈X)..., (b.∈X)...)

F,bsa,_,_ = assemble(Dseq)

u_error = F*u-bsa

#select in u_error en in u

function select_trace(u,trace_ind)
    usub = 0*deepcopy(u)
    N = length(blocks(u))
    for i in 1:N
        if i ∈ HOM
            usub[Block(i)][Block(trace_ind)] += u[Block(i)][Block(trace_ind)]
        elseif i∈ [EFIE...,MFIE...,CFIE...] && trace_ind ∈ [3,4]
            usub[Block(i)][Block(trace_ind)] += u[Block(i)][Block(trace_ind)]
        end
    end
    return usub
end

Gsym = assemble(diag(diag(Identity())),BEAST.DirectProductSpace(X),BEAST.DirectProductSpace(X))
Ginvsym = inv(Matrix(Gsym))

nom = sqrt.([Array(select_trace(u_error,i))'*Ginvsym*Array(select_trace(u_error,i))  for i in 1:4])
denom = sqrt.([Array(select_trace(u,i))'*Matrix(Gsym)*Array(select_trace(u,i))  for i in 1:4])
error = 1/4*sum(nom./denom)