using BEAST
using CompScienceMeshes
using LinearAlgebra
using AbstractTrees
using StaticArrays
import BEAST: FunctionExcitation, matrix_to_bilform, B, nb, nt, ZeroOperator, ∇
function T1(κ;trace=true,tr=1,kwargs...)
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
    if trace
        int = -[γₛ(∇Gx,tr) -γₛ(Gn,tr) ;
        ZeroOperator() -τ(∇Gdotn,tr) ]
        
    else
        int = -[nt×(∇Gx) -nt×(Gn) ;
        ZeroOperator() -(∇Gdotn) ]
    end
    return BEAST.matrix_to_bilform(int)
end
function T2(κ;trace=true,tr=1,kwargs...)
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

    if trace
        int = -[γₛ(∇Gx,tr) ZeroOperator();
        γₙ(Gr,tr) -γₙ(∇G,tr)]
        
    else
        int = -[nt×(∇Gx) ZeroOperator();
        nt ⋅(Gr) -nt ⋅(∇G)]
    end
    return BEAST.matrix_to_bilform(int)
end
function Aa(κ;trace=true,tr=1,kwargs...)
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
    if trace
        int = -[γₛ(Gr,tr) -γₛ(∇G,tr);
        τ(∇Gdot,tr)  κ^2*τ(Gr,tr)]
        
    else
        int = -[nt×(Gr) -nt×(∇G);
        (∇Gdot)  κ^2*(Gr)]
    end
    return BEAST.matrix_to_bilform(int)
end

function Bb(κ;trace=true,tr=1,kwargs...)
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
    if trace
        int = -[γₛ(∇G∇B,tr)+κ^2*γₛ(Gr,tr) -γₛ(∇Gxn,tr);
        γₙ(∇Gx,tr) -γₙ(Gn,tr)]
        
    else
        int = -[nt×(∇G∇B)+κ^2*nt×(Gr) -nt×(∇Gxn) ;
        nt ⋅(∇Gx) -nt ⋅(Gn)]
    end
    return BEAST.matrix_to_bilform(int)
end

Z(κ;rows = [1,2],cols = [1,2],kwargs...) = BEAST.matrix_to_bilform([Bb(κ;kwargs...) T2(κ;kwargs...);T1(κ;kwargs...) Aa(κ;kwargs...)][rows,cols];kwargs...)
Zp(κ;rows = [1,2],cols = [1,2],kwargs...) = BEAST.matrix_to_bilform([Aa(κ;kwargs...) T1(κ;kwargs...);T2(κ;kwargs...) Bb(κ;kwargs...)][rows,cols];kwargs...)


function excitation_dirichlet(A,curlA,divA)
    out = [((n × FunctionExcitation{ComplexF64}(A))),
    FunctionExcitation{ComplexF64}(divA)]
    return BEAST.array_to_linform(out)
end
function excitation_neumann(A,curlA,divA)
    out = [((n × FunctionExcitation{ComplexF64}(curlA))),
    BEAST.NDotTrace{ComplexF64}(FunctionExcitation{ComplexF64}(A))]
    return BEAST.array_to_linform(out)
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


function plot_eigenvalues(M)
    eigs = eigen(M).values
    display(scatter(real.(eigs),imag.(eigs)))

end


### physical constants
ϵ0 = 8.854e-12
μ0 = 4π*1e-7
ω = 10.0^8#*2*pi
κ0 = ω*sqrt(ϵ0*μ0)


#### define meshes
hh = 0.3
Γ1 = meshcuboid(1.0, 1.0, 1.0, hh)
Γ2 =  (@SVector [-1.0,0.0,0.0]) + BEAST.TraceMesh(-Mesh([point(-x,y,z) for (x,y,z) in vertices(Γ1)], deepcopy(cells(Γ1))))
Γ3 = (@SVector [0.0,0.0,-1.0]) + BEAST.TraceMesh(-CompScienceMeshes.translate(Mesh([point(x,y,-z) for (x,y,z) in vertices(Γ1)], deepcopy(cells(Γ1))),[0.0,0.0,0.0]))

Γ = [Γ1,Γ2,Γ3]
Tree = [0,0,0] # give index in which volume material is

HOM = [1,2] #indices of homogeneous domains (without free space)
HomPars = Dict(0=>(ϵ0,μ0),1=>(ϵ0*2,μ0*2),2=>(ϵ0*2,μ0*2))#
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
Xdb = Γ -> BEAST.DirectProductSpace([raviartthomas(Γ),lagrangec0d1(Γ)])
Xnb = Γ -> BEAST.DirectProductSpace([raviartthomas(Γ),lagrangecxd0(Γ)])
Ydb = Γ -> BEAST.DirectProductSpace([BEAST.buffachristiansen(Γ),duallagrangec0d1(Γ)])
Ynb = Γ -> BEAST.DirectProductSpace([BEAST.buffachristiansen(Γ),duallagrangecxd0(Γ)])

Xdt = Γ -> BEAST.DirectProductSpace([n×raviartthomas(Γ),lagrangecxd0(Γ)])
Xnt = Γ -> BEAST.DirectProductSpace([n×(n×(n×raviartthomas(Γ))),lagrangec0d1(Γ)])
Ydt = Γ -> BEAST.DirectProductSpace([n×BEAST.buffachristiansen(Γ),duallagrangecxd0(Γ)])
Ynt = Γ -> BEAST.DirectProductSpace([n×(n×(n×BEAST.buffachristiansen(Γ))),duallagrangec0d1(Γ)])

Xinn = BEAST.array_to_linform([excitation_neumann(A,curlA,divA)])
Xind = BEAST.array_to_linform([excitation_dirichlet(A,curlA,divA)])
Xin = BEAST.array_to_linform([excitation_neumann(A,curlA,divA),excitation_dirichlet(A,curlA,divA)])

idnd = BEAST.matrix_to_bilform(diagm([Identity(),Identity()]))
@hilbertspace t1 t2
@hilbertspace b1 b2
id = idnd[t1,b1]+idnd[t2,b2]
@hilbertspace t1
@hilbertspace b1
idnd2 = idnd[t1,b1]

N = length(Γ)
Q = Dict(i=>diagm([diagm(@SVector [1.0,1.0]),diagm(SVector{2,Float64}([HomPars[parent(i,Tree)][2]/HomPars[i][2],HomPars[i][1]/HomPars[parent(i,Tree)][1]]))]) for i in HOM)
Qp = Dict(i=>diagm([diagm(SVector{2,Float64}([HomPars[parent(i,Tree)][2]/HomPars[i][2],HomPars[i][1]/HomPars[parent(i,Tree)][1]])),diagm(@SVector [1.0,1.0])]) for i in HOM)
Qinv = Dict(i=>diagm([diagm(@SVector [1.0,1.0]),diagm(SVector{2,Float64}(1.0./[HomPars[parent(i,Tree)][2]/HomPars[i][2],HomPars[i][1]/HomPars[parent(i,Tree)][1]]))]) for i in HOM)
Qpinv = Dict(i=>diagm([diagm(SVector{2,Float64}(1.0./[HomPars[parent(i,Tree)][2]/HomPars[i][2],HomPars[i][1]/HomPars[parent(i,Tree)][1]])),diagm(@SVector [1.0,1.0])]) for i in HOM)
t = BEAST.hilbertspace(:t, length(Γ))
b = BEAST.hilbertspace(:b, length(Γ))

##### define space
perm = sortperm([HOM...,EFIE...,CFIE...,MFIE...])
Xt = [BEAST.DirectProductSpace[BEAST.DirectProductSpace([Xnt(Γ[i]),Xdt(Γ[i])]) for i in HOM];BEAST.DirectProductSpace[BEAST.DirectProductSpace([Xdt(Γ[i])]) for i in [EFIE...,CFIE...,MFIE...]]][perm]
Yt = [BEAST.DirectProductSpace[BEAST.DirectProductSpace([Ydt(Γ[i]),Ynt(Γ[i])]) for i in HOM];BEAST.DirectProductSpace[BEAST.DirectProductSpace([Ydt(Γ[i])]) for i in [EFIE...,CFIE...,MFIE...]]][perm]
Xb = [BEAST.DirectProductSpace[BEAST.DirectProductSpace([Xdb(Γ[i]),Xnb(Γ[i])]) for i in HOM];BEAST.DirectProductSpace[BEAST.DirectProductSpace([Xnb(Γ[i])]) for i in [EFIE...,CFIE...,MFIE...]]][perm]
Yb = [BEAST.DirectProductSpace[BEAST.DirectProductSpace([Ynb(Γ[i]),Ydb(Γ[i])]) for i in HOM];BEAST.DirectProductSpace[BEAST.DirectProductSpace([Ynb(Γ[i])]) for i in [EFIE...,CFIE...,MFIE...]]][perm]
Xtmfie = [BEAST.DirectProductSpace[BEAST.DirectProductSpace([Xnt(Γ[i]),Xdt(Γ[i])]) for i in HOM];BEAST.DirectProductSpace[BEAST.DirectProductSpace([Ynt(Γ[i])]) for i in [EFIE...,CFIE...,MFIE...]]][perm]
##### define equation

eqs1 = BEAST.Equation[(Qp[i]*Z(κ[i];tr=-1)*(Qinv[i]))[t[i],b[i]] +
        -sum(BEAST.BilForm[(Qp[i]*Z(κ[i];tr=-1))[t[i],b[j]] for j in HOM ∩ children(i,Tree)]) +
        -sum(BEAST.BilForm[(Qp[i]*Z(κ[i];cols=[2],tr=-1))[t[i],b[j]] for j in [EFIE...,CFIE...,MFIE...] ∩ children(i,Tree)]) ==0 for i in HOM]


eqs2hom = begin BEAST.Equation[-sum(BEAST.BilForm[Z(κ0)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
        -sum(BEAST.BilForm[Z(κ0;cols=[2])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==-Xin[t[ci]]
        for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ HOM] end
eqs2efie = begin BEAST.Equation[-sum(BEAST.BilForm[Z(κ0;rows=[2])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -sum(BEAST.BilForm[Z(κ0;rows=[2],cols=[2])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==-Xind[t[ci]]
            for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ EFIE] end
eqs2mfie = begin BEAST.Equation[-sum(BEAST.BilForm[Z(κ0;rows=[1])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -sum(BEAST.BilForm[Z(κ0;rows=[1],cols=[2])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) +
            -idnd2[t[ci],b[ci]] ==-Xinn[t[ci]]
            for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ MFIE]end
eqs2cefie =begin  BEAST.Equation[-α[ci]*sum(BEAST.BilForm[Z(κ0;rows=[2])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -α[ci]*sum(BEAST.BilForm[Z(κ0;rows=[2],cols=[2])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==
            -α[ci]*Xind[t[ci]] for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end
eqs2cmfie = begin BEAST.Equation[-(1-α[ci])* sum(BEAST.BilForm[Z(κ0;rows=[1])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -(1-α[ci])*sum(BEAST.BilForm[Z(κ0;rows=[1],cols=[2])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) +
            -(1-α[ci])*idnd2[t[ci],b[ci]] == -(1-α[ci])*Xinn[t[ci]]
            for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end


eqs3hom = begin BEAST.Equation[(Z(κ[i])*(Qinv[i]))[t[ci],b[i]]+
            -sum(BEAST.BilForm[Z(κ[i])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
        -sum(BEAST.BilForm[Z(κ[i];cols=[2])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) == 0
        for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ HOM] end
eqs3efie = begin BEAST.Equation[(Z(κ[i];rows=[2])*(Qinv[i]))[t[ci],b[i]]+
            -sum(BEAST.BilForm[Z(κ[i];rows=[2])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -sum(BEAST.BilForm[Z(κ[i];rows=[2],cols=[2])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) == 0
            for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ EFIE] end
eqs3mfie = begin BEAST.Equation[(Z(κ[i];rows=[1])*(Qinv[i]))[t[ci],b[i]]+
            -sum(BEAST.BilForm[Z(κ[i];rows=[1])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -sum(BEAST.BilForm[Z(κ[i];rows=[1],cols=[2])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) +
            -idnd2[t[ci],b[ci]] ==0
            for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ MFIE]end
eqs3cefie = begin BEAST.Equation[α[ci]*(Z(κ[i];rows=[2])*(Qinv[i]))[t[ci],b[i]]+
            -α[ci]*sum(BEAST.BilForm[Z(κ[i];rows=[2])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -α[ci]*sum(BEAST.BilForm[Z(κ[i];rows=[2],cols=[2])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==0
            for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end
eqs3cmfie = begin BEAST.Equation[(1-α[ci])*(Z(κ[i];rows=[1])*(Qinv[i]))[t[ci],b[i]]+
            -(1-α[ci])* sum(BEAST.BilForm[Z(κ[i];rows=[1])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -(1-α[ci])*sum(BEAST.BilForm[Z(κ[i];rows=[1],cols=[2])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) +
            -(1-α[ci])*idnd2[t[ci],b[ci]] == 0
            for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end

## sum the equations in the two parts in a pmchwt fassion,

symeq = -sum(eqs1)+sum(eqs2cefie)+sum(eqs2efie)+sum(eqs2hom)+sum(eqs3cefie)+sum(eqs3efie)+sum(eqs3hom)
asymeq = sum(eqs2cmfie)+sum(eqs2mfie)+sum(eqs3cmfie)+sum(eqs3mfie)

symfilled = typeof(symeq) <: BEAST.Equation
asymfilled = typeof(asymeq) <: BEAST.Equation


symfilled && (Dsymeq = BEAST.discretise(symeq, (t.∈Xt)..., (b.∈Xb)...))
asymfilled && (Dasymeq = BEAST.discretise(asymeq, (t.∈Xtmfie)..., (b.∈Xb)...))
#assemble system


qsZ = [8,8,7,8,5,5,4,3]
qs(::BEAST.LocalOperator, a, b) = BEAST.SingleNumQStrat(qsZ[1])
qs(op::BEAST.ComposedOperatorLocal,testspace,trialpsace) = BEAST.SingleNumQStrat(qsZ[2])
qs(op::BEAST.ComposedOperatorIntegral,testspace,trialspace) = BEAST.DoubleNumSauterQstrat(qsZ[3:end]...) 
qs(fn::BEAST.Functional, basis) = BEAST.SingleNumQStrat(8)
asymfilled && ((Za,ba,xx,yy) = assemble(Dasymeq;quadstratfunction = qs))
symfilled && ((Zs,bs,xx,yy) = assemble(Dsymeq;quadstratfunction = qs))

Zunprec = Za*I+Zs*I
#unprec system


bunprec = ba .+ bs 

uz = solve(Zunprec,bunprec,BEAST.DirectProductSpace(Xb);strat=:GMRES,tol=real(sqrt(eps())), maxiter=20000, restart=20000)




