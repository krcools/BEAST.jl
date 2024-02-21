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
    tr *=-1
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
function excitation(A,curlA,divA;rows=[1,2,3,4])
    out = [((n × FunctionExcitation{ComplexF64}(curlA))×n),
    BEAST.NDotTrace{ComplexF64}(FunctionExcitation{ComplexF64}(A)),
    ((n × FunctionExcitation{ComplexF64}(A))×n),
    FunctionExcitation{ComplexF64}(divA)]
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
Γ1 = meshcuboid(1.0, 1.0, 1.0, hh)
Γ2 =  (@SVector [-1.0,0.0,0.0]) + BEAST.TraceMesh(-Mesh([point(-x,y,z) for (x,y,z) in vertices(Γ1)], deepcopy(cells(Γ1))))
Γ3 = (@SVector [0.0,0.0,-1.0]) + BEAST.TraceMesh(-Mesh([point(x,y,-z) for (x,y,z) in vertices(Γ1)], deepcopy(cells(Γ1))))

Γ = [Γ1,Γ2,Γ3]
Tree = [0,0,0] # give index in which material is

HOM = [1,2,3] #indices of homogeneous domains (without free space)
HomPars = Dict(0=>(ϵ0,μ0),1=>(ϵ0,μ0*3),2=>(ϵ0*6,μ0*4),3=>(ϵ0*2,μ0*3))#index --> (eps,mu)
κ = Dict(i => ω*sqrt(prod(HomPars[i])) for i in HOM)
EFIE = [] #indices of pec domains modeled with efie
MFIE = [] #indices of pec domains modeled with mfie
CFIE = [] #indices of pec domains modeled with cfie
CFIEPars = Dict()#index --> alpha

#### Incident field
A(x) = x[1]*[0.0,0.0,1.0]*sqrt(ϵ0*μ0)*exp(-1.0im*κ0*x[3])
curlA(x) = -[0.0,1.0,0.0]*sqrt(ϵ0*μ0)*exp(-1.0im*κ0 *x[3])
divA(x) = -im*κ0 *x[1]*sqrt(ϵ0*μ0)*exp(-1.0im*κ0 *x[3])
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
#@error "check if correct order of traces"
##### identities
id() = BEAST.matrix_to_bilform(diagm([-NCross(),Identity(),-NCross(),Identity()]))
idN() = BEAST.matrix_to_bilform(diagm([-NCross(),Identity()]))

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

# deqs1 = BEAST.discretise.(eqs1, Ref.((t .∈ X))..., Ref.((b .∈ X))...)

##### define equation
peqs2hom = begin BEAST.Equation[-sum(BEAST.BilForm[Zp(κ0)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
           -sum(BEAST.BilForm[Zp(κ0;cols=[3,4])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==-Xin[t[ci]]
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


Dsymeq = BEAST.discretise(psymeq, (t.∈Y)..., (b.∈Y)...)

Pyy = assemble(Dsymeq.equation.lhs,BEAST.DirectProductSpace(Y),BEAST.DirectProductSpace(Y))


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
