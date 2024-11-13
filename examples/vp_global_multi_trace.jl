using BEAST
using CompScienceMeshes
using LinearAlgebra
using SparseArrays
using PlotlyJS
using StaticArrays

ω = 10.0^8*2*pi
ϵ0 = 8.854e-12
μ0 = 4π*1e-7
κ0 = ω*√(μ0*ϵ0)

ϵr = [1,2,2]
μr = [2,1,2]

h = 0.5
fn = joinpath(@__DIR__, "assets/rectangular_torus.geo")
Γ11 = CompScienceMeshes.meshgeo(fn; dim=2, h)
Γ11 = Mesh([point(-z,y,x) for (x,y,z) in vertices(Γ11)], deepcopy(cells(Γ11)))
function onlys(a)
    c = sort.(a)
    b = [i for i in a if length(findall(==(sort(i)),c))==1]
    return b
end

Γ12 = -Mesh([point(x,y,-z-4) for (x,y,z) in vertices(Γ11)], deepcopy(cells(Γ11)))
Γ1 = weld(Γ11,Γ12;glueop=onlys)
Γ2 = Mesh([point(x,-y-4,-z-4) for (x,y,z) in vertices(Γ11)], deepcopy(cells(Γ11)))
Γ3 = -Mesh([point(x,-y-4,z) for (x,y,z) in vertices(Γ11)], deepcopy(cells(Γ11)))


Γ = [Γ1,Γ2,Γ3] # the meshes

κ = [ω*√(μ0*ϵ0*ϵri*μri) for (ϵri,μri) in zip(ϵr,μr)]

########RHS definition (lorentz gauge:)
x0 = @SVector [0.0,0.0,0.0] # used to check the type of the incidnet fields
A(x) = x[1]*sqrt(ϵ0*μ0)exp(-1im*κ0*x[3])* (@SVector [0.0,0.0,1.0])
curlA(x) = -(@SVector [0.0,1.0,0.0])*sqrt(ϵ0*μ0)*exp(-1.0im*κ0 *x[3])
phi(x) =  x[1]*exp(-1.0im*κ0 *x[3])
gradphi(x) = (@SVector [1.0,0.0,-1im*κ0*x[1]])*exp(-1.0im*κ0 *x[3])

####################################################
#
# From here on do nothing
#
####################################################
@hilbertspace a1 a2
@hilbertspace b1 b2

#### potentials + operators
G0 = BEAST.HH3DGreen(1im*κ0)
G = [BEAST.HH3DGreen(1im*k) for k in κ]
dG0 = BEAST.HH3DGradGreen(1im*κ0)
dG = [BEAST.HH3DGradGreen(1im*k) for k in κ]

    U11 = -BEAST.strace(BEAST.PotentialIntegralOperator{2}(dG0,×,b->b),-1.0)
    U12 = BEAST.strace(BEAST.PotentialIntegralOperator{2}(G0,*,b->n*b),-1.0)
    U22 = BEAST.trace(BEAST.PotentialIntegralOperator{2}(dG0,(x,y)->transpose(x)*y,b->n*b),-1.0)

DDL0 = U11[a1,b1] + U12[a1,b2] + U22[a2,b2]

    U11 = [-BEAST.strace(BEAST.PotentialIntegralOperator{2}(g,×,b->b),1.0) for g in dG]
    U12 = [BEAST.strace(BEAST.PotentialIntegralOperator{2}(g,*,b->n*b),1.0) for g in G]
    U22 = [BEAST.trace(BEAST.PotentialIntegralOperator{2}(g,(x,y)->transpose(x)*y,b->n*b),1.0) for g in dG]

DDL = [U11i[a1,b1] + er*mr*U12i[a1,b2] + U22i[a2,b2] for (U11i,U12i,U22i,er,mr) in zip(U11,U12,U22,ϵr,μr)]

    U11 = -BEAST.strace(BEAST.PotentialIntegralOperator{2}(dG0,×,b->b),-1.0)
    U21 = -BEAST.ntrace(BEAST.PotentialIntegralOperator{2}(G0,*,b->b),-1.0)
    U22 = BEAST.ntrace(BEAST.PotentialIntegralOperator{2}(dG0,*,b->b),-1.0)

NDL0 = U11[a1,b1] + U21[a2,b1] + U22[a2,b2]

    U11 = [-BEAST.strace(BEAST.PotentialIntegralOperator{2}(g,×,b->b),1.0) for g in dG]
    U21 = [-BEAST.ntrace(BEAST.PotentialIntegralOperator{2}(g,*,b->b),1.0) for g in G]
    U22 = [BEAST.ntrace(BEAST.PotentialIntegralOperator{2}(g,*,b->b),1.0) for g in dG]

NDL = [U11i[a1,b1] + er*mr*U21i[a2,b1] + U22i[a2,b2] for (U11i,U21i,U22i,er,mr) in zip(U11,U21,U22,ϵr,μr)]

    U11 = -BEAST.strace(BEAST.PotentialIntegralOperator{2}(G0,*,b->b),-1.0)
    U12 = BEAST.strace(BEAST.PotentialIntegralOperator{2}(dG0,*,b->b),-1.0)
    U21 = -BEAST.trace(BEAST.PotentialIntegralOperator{2}(dG0,(x,y)->transpose(x)*y,b->b),-1.0)
    U22 = -κ0^2* BEAST.trace(BEAST.PotentialIntegralOperator{2}(G0,*,b->b),-1.0)

NDSL0 = U11[a1,b1] + U12[a1,b2] + U21[a2,b1] + U22[a2,b2]

    U11 = [-BEAST.strace(BEAST.PotentialIntegralOperator{2}(g,*,b->b),1.0) for g in G]
    U12 = [BEAST.strace(BEAST.PotentialIntegralOperator{2}(g,*,b->b),1.0) for g in dG]
    U21 = [-BEAST.trace(BEAST.PotentialIntegralOperator{2}(g,(x,y)->transpose(x)*y,b->b),1.0) for g in dG]
    U22 = [-k^2*BEAST.trace(BEAST.PotentialIntegralOperator{2}(g,*,b->b),1.0) for (g,k) in zip(G,κ)]

NDSL = [mr*U11i[a1,b1] + 1/er*U12i[a1,b2] + 1/er*U21i[a2,b1] + 1/(er^2*mr)*U22i[a2,b2] for (U11i,U12i,U21i,U22i,er,mr) in zip(U11,U12,U21,U22,ϵr,μr)]

    U11 = -κ0^2 * BEAST.pvstrace(BEAST.PotentialIntegralOperator{2}(G0,*,b->b),-1.0)-BEAST.CompDoubleInt(B->divergence(n×B),*,G0,*,B->divergence(B))
    U12 = BEAST.strace(BEAST.PotentialIntegralOperator{2}(dG0,×,b->n*b),-1.0)
    U21 = -BEAST.ntrace(BEAST.PotentialIntegralOperator{2}(dG0,×,b->b),-1.0)
    U22 =  BEAST.ntrace(BEAST.PotentialIntegralOperator{2}(G0,*,b->n*b),-1.0)

DNSL0 = U11[a1,b1] + U12[a1,b2] + U21[a2,b1] + U22[a2,b2]

    U11 = [-k^2 * BEAST.pvstrace(BEAST.PotentialIntegralOperator{2}(g,*,b->b),1.0)-BEAST.CompDoubleInt(B->divergence(n×B),*,g,*,B->divergence(B)) for (g,k) in zip(G,κ)]
    U12 = [BEAST.strace(BEAST.PotentialIntegralOperator{2}(g,×,b->n*b),1.0) for g in dG]
    U21 = [-BEAST.ntrace(BEAST.PotentialIntegralOperator{2}(g,×,b->b),1.0) for g in dG]
    U22 = [BEAST.ntrace(BEAST.PotentialIntegralOperator{2}(g,*,b->n*b),1.0) for (g,k) in zip(G,κ)]

DNSL = [1/mr*U11i[a1,b1] + er*U12i[a1,b2] + er*U21i[a2,b1] + er^2*mr*U22i[a2,b2] for (U11i,U12i,U21i,U22i,er,mr) in zip(U11,U12,U21,U22,ϵr,μr)]

### definition per domain bilinear forms

A0 = DNSL0[a1,b1] + NDSL0[a2,b2] + DDL0[a2,b1] + NDL0[a1,b2]
Ad = [DNSLi[a1,b1] + NDSLi[a2,b2] + DDLi[a2,b1] + NDLi[a1,b2] for (DNSLi,NDSLi,DDLi,NDLi) in zip(DNSL,NDSL,DDL,NDL)]

p = BEAST.hilbertspace(:p, length(κ))
q = BEAST.hilbertspace(:q, length(κ))

LHSA = A0[p,q] + BEAST.Variational.DirectProductKernel(Ad)[p,q]

####################################################
#
# RHS definition
#
####################################################

### wrapping incident field

w_A = BEAST.FunctionWrapper{typeof(A(x0))}(A)
w_curlA = BEAST.FunctionWrapper{typeof(curlA(x0))}(curlA)
w_phi = BEAST.FunctionWrapper{typeof(phi(x0))}(phi)
w_gradphi = BEAST.FunctionWrapper{typeof(gradphi(x0))}(gradphi)

### traces RHS

t_nxA = BEAST.CrossTraceMW(w_A)
t_nxcurlA = BEAST.CrossTraceMW(w_curlA)
t_ndA = BEAST.NDotTrace(w_A)
t_div_rescA =  BEAST.Trace(w_phi)

t_d = t_nxA[a1] + -1im*κ0^2/ω*t_div_rescA[a2] 
t_n = t_nxcurlA[a1] + t_ndA[a2]
RHSAi = t_d[a2] + t_n[a1]

t_phi = BEAST.Trace(w_phi)
t_ngradphi = BEAST.NDotTrace(w_gradphi)

RHSA = RHSAi[p]

### Spaces
RT = raviartthomas.(Γ)
ND = Ref(n) .× RT
MND = Ref(n) .× (Ref(n) .× (Ref(n) .× RT))
L0 = lagrangecxd0.(Γ)
L1 = lagrangec0d1.(Γ)

Xd = [RTi×L1i for (RTi,L1i) in zip(RT,L1)]
Xn = [RTi×L0i for (RTi,L0i) in zip(RT,L0)]

Yn = [NDi×L1i for (NDi,L1i) in zip(MND,L1)]
Yd = [NDi×L0i for (NDi,L0i) in zip(ND,L0)]

Qₕ = [BEAST.DirectProductSpace([Xdi,Xni]) for (Xdi,Xni) in zip(Xd,Xn)]
Pₕ = [BEAST.DirectProductSpace([Yni,Ydi]) for (Ydi,Yni) in zip(Yd,Yn)]

deq = BEAST.discretise(LHSA==RHSA, (p .∈ Pₕ)..., (q .∈ Qₕ)...)
# Z = assemble(LHSA,BEAST.DirectProductSpace(Pₕ),BEAST.DirectProductSpace(Qₕ))
# B = assemble(RHSA,BEAST.DirectProductSpace(Pₕ))
u = solve(deq)

