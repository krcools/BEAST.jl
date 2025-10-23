abstract type QHComponent end

#Local loops + global loops
struct Loops <: QHComponent end
#Stars
struct Stars <: QHComponent end
#Local loops
struct DualLoops <: QHComponent end
#Stars + global loops
struct DualStars <: QHComponent end

abstract type ComputeStrat end
struct Direct <: ComputeStrat end
struct Iterative <: ComputeStrat end

struct QHProjector{C<:QHComponent, S<:ComputeStrat} end

function PΣ(;compStrat = Iterative)
    return QHProjector{Stars,compStrat}()
end

function PΛ(;compStrat = Iterative)
    return QHProjector{Loops,compStrat}()
end

function ℙΣ(;compStrat = Iterative)
    return QHProjector{DualStars,compStrat}()
end

function ℙΛ(;compStrat = Iterative)
    return QHProjector{DualLoops,compStrat}()
end

"""
saddlepoint(A,B,P1,P2)
Create preconditioned saddlepoint matrix from A and B with block diagonal preconditoners [P1,P2]1 and P2. 
"""
function saddlepoint(A::SparseMatrixCSC,B::SparseMatrixCSC,P1::SparseMatrixCSC,P2::SparseMatrixCSC) 
    T= eltype(A)
    nP = size(B,2)
    SP = [A    B
          B'  spzeros(T,nP,nP)]
    Pdiv = blockdiag(CholeskyFactorization(P1),CholeskyFactorization(P2))
    return GMRESSolver(SP, left_preconditioner=Pdiv, maxiter=200, restart=50, reltol=1e-8, verbose=false) 
end

function assemble(::QHProjector, X::Space; quadstrat=defaultquadstrat)
    error("Not implemented")
end

#RT Basis
function assemble(::QHProjector{Stars,Direct}, X::RTBasis; quadstrat=defaultquadstrat)
    edges = setminus(skeleton(X.geo,1), boundary(X.geo))
    Σ = Matrix(connectivity(X.geo, edges, sign))
    return Σ * pinv(Σ'*Σ) * Σ'
end

function assemble(::QHProjector{Loops,Direct}, X::RTBasis; quadstrat=defaultquadstrat)
   edges = setminus(skeleton(X.geo,1), boundary(X.geo))
   Σ = Matrix(connectivity(X.geo, edges, sign))
   return LinearAlgebra.I - Σ*pinv(Σ'*Σ)*Σ'
end

function assemble(::QHProjector{DualLoops,Direct}, X::RTBasis; quadstrat=defaultquadstrat)
    edges = setminus(skeleton(X.geo,1), boundary(X.geo))
    verts = setminus(skeleton(X.geo,0), skeleton(boundary(X.geo),0))
    Λ = Matrix(connectivity(verts, edges, sign))
    return  Λ * pinv(Λ'*Λ) * Λ'
end

function assemble(::QHProjector{DualStars,Direct}, X::RTBasis; quadstrat=defaultquadstrat)
    edges = setminus(skeleton(X.geo,1), boundary(X.geo))
    verts = setminus(skeleton(X.geo,0), skeleton(boundary(X.geo),0))
    Λ = Matrix(connectivity(verts, edges, sign))
    return  LinearAlgebra.I - Λ * pinv(Λ'*Λ) * Λ'
end

function assemble(::QHProjector{Stars,Iterative}, X::RTBasis; quadstrat=defaultquadstrat)
    #create auxilarry basis functions
    P = lagrangecxd0(X.geo)
    nX = numfunctions(X)
    nP = numfunctions(P)
    #assemble auxilarry matrices
    Gxx = assemble(Identity(),X,X;quadstrat)
    Dxx = assemble(Identity(),divergence(X),divergence(X);quadstrat)
    Gpp = assemble(Identity(),P,P;quadstrat)
    Σp = assemble(Identity(),divergence(X),P;quadstrat)
    Px0 = [Gxx
           spzeros(nP,nX)]
    P0Σ = [spzeros(nX,nX) Σp]
    SP = saddlepoint(Gxx,Σp,Gxx+Dxx,Gpp)
    return P0Σ*SP*Px0
end

function assemble(::QHProjector{Loops,Iterative}, X::RTBasis; quadstrat=defaultquadstrat)
    #create auxilarry basis functions
    P = lagrangecxd0(X.geo)
    nX = numfunctions(X)
    nP = numfunctions(P)
    #assemble auxilary matrices
    Gxx = assemble(Identity(),X,X;quadstrat)
    Dxx = assemble(Identity(),divergence(X),divergence(X);quadstrat)
    Gpp = assemble(Identity(),P,P;quadstrat)
    Σp = assemble(Identity(),divergence(X),P;quadstrat)
    Px0 = [Gxx
           spzeros(nP,nX)]
    SP = saddlepoint(Gxx,Σp,Gxx+Dxx,Gpp)
    return Px0'*SP*Px0
end

function assemble(::QHProjector{DualStars,Iterative}, X::RTBasis; quadstrat=defaultquadstrat)
    #create auxilarry basis functions
    L = lagrangec0d1(X.geo)
    nX = numfunctions(X)
    nL = numfunctions(L)
    #assemble auxilarry matrices
    Gxx = assemble(Identity(),X,X;quadstrat)
    Cll = assemble(Identity(),curl(L),curl(L);quadstrat)
    Gll = assemble(Identity(),L,L;quadstrat)
    Λp = assemble(Identity(),X,curl(L);quadstrat)
    Px0 = [Gxx
           spzeros(nL,nX)]
    SP = saddlepoint(Gxx,Λp,Gxx,Gll+Cll)
    return Px0'*SP*Px0
end

function assemble(::QHProjector{DualLoops,Iterative}, X::RTBasis; quadstrat=defaultquadstrat)
    #create auxilarry basis functions
    L = lagrangec0d1(X.geo)
    nX = numfunctions(X)
    nP = numfunctions(L)
    #assemble auxilarry matrices
    Gxx = assemble(Identity(),X,X;quadstrat)
    Cll = assemble(Identity(),curl(L),curl(L);quadstrat)
    Gll = assemble(Identity(),L,L;quadstrat)
    Λp = assemble(Identity(),X,curl(L);quadstrat)
    Px0 = [Gxx
           spzeros(nP,nX)]
    P0Λ = [spzeros(nX,nX) Λp]
    SP = saddlepoint(Gxx,Λp,Gxx,Gll+Cll)
    return P0Λ*SP*Px0
end

#GWP basis
function assemble(::QHProjector{Stars,Direct}, X::GWPDivSpace; quadstrat=defaultquadstrat)
    p = X.degree
    #create auxilarry basis functions
    P = lagrangecx(X.geo,order=p)
    #assemble auxilarry matrices
    Gxx = assemble(Identity(),X,X;quadstrat)
    Σp = assemble(Identity(),divergence(X),P;quadstrat)
    return Σp * pinv(Matrix(Σp'*inv(Matrix(Gxx))*Σp)) * Σp'
end

function assemble(::QHProjector{Loops,Direct}, X::GWPDivSpace; quadstrat=defaultquadstrat)
    p = X.degree
    #create auxilarry basis functions
    P = lagrangecx(X.geo,order=p)
    #assemble auxilarry matrices
    Gxx = assemble(Identity(),X,X;quadstrat)
    Σp = assemble(Identity(),divergence(X),P;quadstrat)
    return Gxx - Σp * pinv(Matrix(Σp'*inv(Matrix(Gxx))*Σp)) * Σp'
end

function assemble(::QHProjector{DualStars,Direct}, X::GWPDivSpace; quadstrat=defaultquadstrat)
    p = X.degree
    #create auxilarry basis functions
    L = lagrangec0(X.geo,order=p+1)
    #assemble auxilarry matrices
    Gxx = assemble(Identity(),X,X;quadstrat) 
    Λp = assemble(Identity(),X,curl(L);quadstrat)
    return Gxx - Λp * pinv(Matrix(Λp'*inv(Matrix(Gxx))*Λp)) * Λp'
end

function assemble(::QHProjector{DualLoops,Direct}, X::GWPDivSpace; quadstrat=defaultquadstrat)
    p = X.degree
    #create auxilarry basis functions
    L = lagrangec0(X.geo,order=p+1)
    #assemble auxilarry matrices
    Gxx = assemble(Identity(),X,X;quadstrat) 
    Λp = assemble(Identity(),X,curl(L);quadstrat)
    return Λp * pinv(Matrix(Λp'*inv(Matrix(Gxx))*Λp)) * Λp'
end

function assemble(::QHProjector{Stars,Iterative}, X::GWPDivSpace; quadstrat=defaultquadstrat)
    p = X.degree
    #create auxilarry basis functions
    P = lagrangecx(X.geo,order=p)
    nX = numfunctions(X)
    nP = numfunctions(P)
    #assemble auxilarry matrices
    Gxx = assemble(Identity(),X,X;quadstrat)
    Dxx = assemble(Identity(),divergence(X),divergence(X);quadstrat)
    Gpp = assemble(Identity(),P,P;quadstrat)
    Σp = assemble(Identity(),divergence(X),P;quadstrat)
    Px0 = [Gxx
           spzeros(nP,nX)]
    P0Σ = [spzeros(nX,nX) Σp]
    SP = saddlepoint(Gxx,Σp,Gxx+Dxx,Gpp)
    return P0Σ*SP*Px0
end

function assemble(::QHProjector{Loops,Iterative}, X::GWPDivSpace; quadstrat=defaultquadstrat)
    p = X.degree
    #create auxilary basis functions
    P = lagrangecx(X.geo,order=p)
    nX = numfunctions(X)
    nP = numfunctions(P)
    #assemble auxilary matrices
    Gxx = assemble(Identity(),X,X;quadstrat)
    Dxx = assemble(Identity(),divergence(X),divergence(X);quadstrat)
    Gpp = assemble(Identity(),P,P;quadstrat)
    Σp = assemble(Identity(),divergence(X),P;quadstrat)
    Px0 = [Gxx
           spzeros(nP,nX)]
    SP = saddlepoint(Gxx,Σp,Gxx+Dxx,Gpp)
    return Px0'*SP*Px0
end

function assemble(::QHProjector{DualStars,Iterative}, X::GWPDivSpace; quadstrat=defaultquadstrat)
    p = X.degree
    #create auxilarry basis functions
    L = lagrangec0(X.geo,order=p+1)
    nX = numfunctions(X)
    nL = numfunctions(L)
    #assemble auxilarry matrices
    Gxx = assemble(Identity(),X,X;quadstrat)
    Cll = assemble(Identity(),curl(L),curl(L);quadstrat)
    Gll = assemble(Identity(),L,L;quadstrat)
    #Σp = assemble(Identity(),divergence(X),P)
    Λp = assemble(Identity(),X,curl(L);quadstrat)
    Px0 = [Gxx
           spzeros(nL,nX)]
    SP = saddlepoint(Gxx,Λp,Gxx,Gll+Cll)
    return Px0'*SP*Px0
end

function assemble(::QHProjector{DualLoops,Iterative}, X::GWPDivSpace; quadstrat=defaultquadstrat)
    p = X.degree
    #create auxilarry basis functions
    L = lagrangec0(X.geo,order=p+1)
    nX = numfunctions(X)
    nP = numfunctions(L)
    #assemble auxilarry matrices
    Gxx = assemble(Identity(),X,X;quadstrat)
    Cll = assemble(Identity(),curl(L),curl(L);quadstrat)
    Gll = assemble(Identity(),L,L;quadstrat)
    #Σp = assemble(Identity(),divergence(X),P)
    Λp = assemble(Identity(),X,curl(L);quadstrat)
    Px0 = [Gxx
           spzeros(nP,nX)]
    P0Λ = [spzeros(nX,nX) Λp]
    SP = saddlepoint(Gxx,Λp,Gxx,Gll+Cll)
    return P0Λ*SP*Px0
end   


@testitem "QH-Projectors" begin 

    using CompScienceMeshes, BEAST
    using LinearAlgebra

    h = 0.5
    M = meshrectangle(1.0,1.0,h)

    X = raviartthomas(M)
    Y = BEAST.gwpdiv(M;order=1)

    G = assemble(BEAST.Identity(),Y,Y)
    iG = BEAST.cholesky(G)

    b = rand(numfunctions(Y))
    d = rand(numfunctions(X))

    PΣd = assemble(BEAST.PΣ(;compStrat =BEAST.Direct),X)
    PΛd = assemble(BEAST.PΛ(;compStrat =BEAST.Direct),X)

    ℙΣd = assemble(BEAST.ℙΣ(;compStrat =BEAST.Direct),X)
    ℙΛd = assemble(BEAST.ℙΛ(;compStrat =BEAST.Direct),X)

    @test norm(PΣd*PΛd*d)/norm(d) < sqrt(eps())
    @test norm(ℙΣd*ℙΛd*d)/norm(d) < sqrt(eps())

    PΣd = assemble(BEAST.PΣ(;compStrat =BEAST.Direct),Y)
    PΛd = assemble(BEAST.PΛ(;compStrat =BEAST.Direct),Y)

    ℙΣd = assemble(BEAST.ℙΣ(;compStrat =BEAST.Direct),Y)
    ℙΛd = assemble(BEAST.ℙΛ(;compStrat =BEAST.Direct),Y)

    @test norm(PΣd*iG*PΛd*b)/norm(b) <  sqrt(eps())
    @test norm(ℙΣd*iG*ℙΛd*b)/norm(b) <  sqrt(eps())

    PΣd = assemble(BEAST.PΣ(;compStrat =BEAST.Iterative),Y)
    PΛd = assemble(BEAST.PΛ(;compStrat =BEAST.Iterative),Y)

    ℙΣd = assemble(BEAST.ℙΣ(;compStrat =BEAST.Iterative),Y)
    ℙΛd = assemble(BEAST.ℙΛ(;compStrat =BEAST.Iterative),Y)

    @test norm(PΣd*iG*PΛd*b)/norm(b) <  sqrt(eps())
    @test norm(ℙΣd*iG*ℙΛd*b)/norm(b) <  sqrt(eps())

end 