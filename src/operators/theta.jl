abstract type Theta{S<:ComputeStrat} end
struct NedelecTheta{C<:QHComponent,S} <: Theta{S} end
struct BourhisTheta{C<:QHComponent,D<:QHComponent,S} <: Theta{S} end

function nedelecTheta(::C; compStrat = Iterative) where C<:QHComponent
    return NedelecTheta{C,compStrat}()
end

function bourhisTheta(::C,::D; compStrat = Iterative) where {C<:QHComponent,D<:QHComponent}
    return BourhisTheta{C,D,compStrat}()
end

ThetaLoops(;compStrat = Iterative) = bourhisTheta(Loops(),Stars();compStrat = compStrat)
ThetaStars(;compStrat = Iterative) = bourhisTheta(Stars(),Loops();compStrat = compStrat)

function assemble(theta::BourhisTheta{C,D,Iterative},Y::GWPDivSpace,X::GWPDivSpace;quadstrat=defaultquadstrat) where {C,D}
    Gxx   = assemble(Identity(),X,X;quadstrat)
    Gxx_r = assemble(Identity(),Y,Y;quadstrat)
    Gmix  = assemble(NCross(),Y,X;quadstrat)
    iGxx = cholesky(Gxx)
    iGxx_r = cholesky(Gxx_r)
    
    Pp = assemble(QHProjector{D,Iterative}(),X;quadstrat) 
    Pr = assemble(QHProjector{C,Iterative}(),Y;quadstrat)  

    return iGxx_r* Pr*iGxx_r*Gmix*iGxx*Pp*iGxx
end

function assemble(theta::BourhisTheta{C,D,Direct},Y::GWPDivSpace,X::GWPDivSpace;quadstrat=defaultquadstrat) where {C,D}
    Gxx   = assemble(Identity(),X,X;quadstrat)
    Gxx_r = assemble(Identity(),Y,Y;quadstrat)
    Gmix  = assemble(NCross(),Y,X;quadstrat)
    iGxx = inv(Matrix(Gxx))
    iGxx_r = inv(Matrix(Gxx_r))
    
    Pp = assemble(QHProjector{D,Direct}(),X;quadstrat) 
    Pr = assemble(QHProjector{C,Direct}(),Y;quadstrat)  

    return iGxx_r* Pr*iGxx_r*Gmix*iGxx*Pp*iGxx
end

function assemble(theta::NedelecTheta{Stars,Direct},Y::Space,X::Space;quadstrat=defaultquadstrat)
    P = lagrangecx(X.geo,order=X.degree)
    L = lagrangec0(X.geo,order=X.degree+1)
    
    Gxx =  assemble(Identity(),X,X;quadstrat)
    Gdiv = assemble(Identity(),divergence(X),P;quadstrat)
    Gll = assemble(Identity(),L,L;quadstrat)
    Gmix  = assemble(Identity(),L,P;quadstrat)
    Gcurl = assemble(Identity(),Y,curl(L);quadstrat)
    iGll = inv(Matrix(Gll))
    Gyy =  assemble(Identity(),Y,Y;quadstrat)
    iGyy = inv(Matrix(Gyy))
    
    nP = size(Gpp,2)
    nX = size(Gxx,2)
    T = eltype(Gxx)
    SP = inv(Matrix([Gxx   Gdiv
             Gdiv'  spzeros(T,nP,nP)]))  
    Px0 = [sparse(I,nX,nX)
        spzeros(T,nP,nX)]
    P0p = [spzeros(T,nP,nX)  sparse(I,nP,nP)]
    
    return -iGyy*Gcurl*iGll*Gmix*P0p*SP*Px0
end

function assemble(theta::NedelecTheta{Loops,Direct},Y::Space,X::Space;quadstrat=defaultquadstrat)
    P = lagrangecx(X.geo,order=X.degree)
    
    Gxx =  assemble(Identity(),X,X;quadstrat)
    Gdiv = assemble(Identity(),divergence(X),P;quadstrat)
    Nyx =  assemble(NCross(),Y,X;quadstrat)
    Gyy =  assemble(Identity(),Y,Y;quadstrat)
    iGyy = inv(Matrix(Gyy))

    nP = size(Gpp,2)
    nX = size(Gxx,2)
    T = eltype(Gxx)
    SP = inv(Matrix([Gxx   Gdiv
             Gdiv'  spzeros(T,nP,nP)]))     
    Px0 = [sparse(I,nX,nX)
        spzeros(T,nP,nX)]

    return iGyy*Nyx*Px0'*SP*Px0
end

function assemble(theta::NedelecTheta{Stars,Iterative},Y::Space,X::Space;quadstrat=defaultquadstrat)
    P = lagrangecx(X.geo,order=X.degree)
    L = lagrangec0(X.geo,order=X.degree+1)
    
    Gxx =  assemble(Identity(),X,X;quadstrat)
    Gdiv = assemble(Identity(),divergence(X),P;quadstrat)
    Gll = assemble(Identity(),L,L;quadstrat)
    Gmix  = assemble(Identity(),L,P;quadstrat)
    Gcurl = assemble(Identity(),Y,curl(L);quadstrat)
    D = assemble(Identity(),divergence(X),divergence(X);quadstrat)
    Gpp = assemble(Identity(),P,P;quadstrat)
    iGll = cholesky(Gll)
    Gyy =  assemble(Identity(),Y,Y;quadstrat)
    iGyy = cholesky(Gyy)
    
    SP = saddlepoint(Gxx,Gdiv,Gxx+D,Gpp)
    nP = size(Gpp,2)
    nX = size(Gxx,2)
    T = eltype(Gxx)
    Px0 = [sparse(I,nX,nX)
        spzeros(T,nP,nX)]
    P0p = [spzeros(T,nP,nX)  sparse(I,nP,nP)]
    
    return -iGyy*Gcurl*iGll*Gmix*P0p*SP*Px0
end

function assemble(theta::NedelecTheta{Loops,Iterative},Y::Space,X::Space;quadstrat=defaultquadstrat)
    P = lagrangecx(X.geo,order=X.degree)

    Gxx =  assemble(Identity(),X,X;quadstrat)
    Gdiv = assemble(Identity(),divergence(X),P;quadstrat)
    D = assemble(Identity(),divergence(X),divergence(X);quadstrat)
    Gpp = assemble(Identity(),P,P;quadstrat)
    Nyx =  assemble(NCross(),Y,X;quadstrat)
    Gyy =  assemble(Identity(),Y,Y;quadstrat)
    iGyy = cholesky(Gyy)

    SP = saddlepoint(Gxx,Gdiv,Gxx+D,Gpp)
    nP = size(Gpp,2)
    nX = size(Gxx,2)
    T = eltype(Gxx)
    Px0 = [sparse(I,nX,nX)
        spzeros(T,nP,nX)]

    return iGyy*Nyx*Px0'*SP*Px0
end