# TODO
# Make it possible to measure the complement error
# Cleanup the code, different file for post processing
# other branch Make connection with single-layer operator in different file called: connections with maxwell3D
# later Make PMCHWT, EFIE, MFIE description in the same way
# Include PEC
# Write test files
# Niet plat ding in sheets steken, halve sfeer. --> werkt (hier valt verschil efie-cfie mss nog te halen)
struct PRECVP <: Interaction 
    trace::Int
end
struct PRECVPPMCHWT2 <: Interaction 
    trace::Int
    rows
    cols
end
struct VPPMCHWT <: Interaction
    trace::Int
end
VPPMCHWT() = VPPMCHWT(1)
struct VPPMCHWT2 <: Interaction
    trace::Int
end
PRECVPPMCHWT2() = PRECVPPMCHWT2(1,[1,2,3,4],[1,2,3,4])
VPPMCHWT2() = VPPMCHWT2(1)
PRECVP() = PRECVP(1)

function testbasis(obj::Object{T}, ::VP) where {T <: HOM}
    Γ =  obj.type.mesh
    x1 = raviartthomas(Γ)
    x2 = BEAST.lagrangecxd0(Γ)
    x3 = raviartthomas(Γ)
    x4 = BEAST.lagrangec0d1(Γ)
    return BEAST.DirectProductSpace([x1,x2,x3,x4])
end

function trialbasis(obj::Object{T}, ::VP) where {T <: HOM}
    Γ =  obj.type.mesh
    x1 = raviartthomas(Γ)
    x2 = BEAST.lagrangec0d1(Γ)
    x3 = raviartthomas(Γ)
    x4 = BEAST.lagrangecxd0(Γ)
    return BEAST.DirectProductSpace([x1,x2,x3,x4])
end
function testbasis(obj::Object{T}, ::PRECVP) where {T <: HOM}
    Γ =  obj.type.mesh
    x1 = buffachristiansen(Γ)
    x2 = BEAST.duallagrangecxd0(Γ)
    x3 = buffachristiansen(Γ)
    x4 = BEAST.duallagrangec0d1(Γ)
    return BEAST.DirectProductSpace([x1,x2,x3,x4])
end

function trialbasis(obj::Object{T}, ::PRECVP) where {T <: HOM}
    Γ =  obj.type.mesh
    x1 = buffachristiansen(Γ)
    x2 = BEAST.duallagrangec0d1(Γ)
    x3 = buffachristiansen(Γ)
    x4 = BEAST.duallagrangecxd0(Γ)
    return BEAST.DirectProductSpace([x1,x2,x3,x4])
end
function testbasis(obj::Object{T}, ::PRECVPPMCHWT2) where {T <: HOM}
    Γ =  obj.type.mesh
    x1 = buffachristiansen(Γ)
    x2 = BEAST.duallagrangec0d1(Γ)
    x3 = buffachristiansen(Γ)
    x4 = BEAST.duallagrangecxd0(Γ)
    return BEAST.DirectProductSpace([x1,x2,x3,x4])
end

function trialbasis(obj::Object{T}, ::PRECVPPMCHWT2) where {T <: HOM}
    Γ =  obj.type.mesh
    x1 = buffachristiansen(Γ)
    x2 = BEAST.duallagrangec0d1(Γ)
    x3 = buffachristiansen(Γ)
    x4 = BEAST.duallagrangecxd0(Γ)
    return BEAST.DirectProductSpace([x1,x2,x3,x4])
end
function testbasis(obj::Object{T}, ::VPPMCHWT) where {T <: HOM}
    Γ =  obj.type.mesh
    x1 = raviartthomas(Γ)
    x2 = BEAST.lagrangec0d1(Γ)
    x3 = raviartthomas(Γ)
    x4 = BEAST.lagrangecxd0(Γ)
    return BEAST.DirectProductSpace([x1,x2,x3,x4])
end

function trialbasis(obj::Object{T}, ::VPPMCHWT) where {T <: HOM}
    Γ =  obj.type.mesh
    x1 = raviartthomas(Γ)
    x2 = BEAST.lagrangec0d1(Γ)
    x3 = raviartthomas(Γ)
    x4 = BEAST.lagrangecxd0(Γ)
    return BEAST.DirectProductSpace([x1,x2,x3,x4])
end
function testbasis(obj::Object{T}, ::VPPMCHWT2) where {T <: HOM}
    Γ =  obj.type.mesh
    x1 = raviartthomas(Γ)
    x2 = BEAST.lagrangec0d1(Γ)
    x3 = raviartthomas(Γ)
    x4 = BEAST.lagrangecxd0(Γ)
    return BEAST.DirectProductSpace([x1,x2,x3,x4])
end

function trialbasis(obj::Object{T}, ::VPPMCHWT2) where {T <: HOM}
    Γ =  obj.type.mesh
    x1 = raviartthomas(Γ)
    x2 = BEAST.lagrangec0d1(Γ)
    x3 = raviartthomas(Γ)
    x4 = BEAST.lagrangecxd0(Γ)
    return BEAST.DirectProductSpace([x1,x2,x3,x4])
end
# function testbasis(obj::Object{T}, ::PRECVP) where {T <: HOM}
#     Γ =  obj.type.mesh
#     x1 = buffachristiansen(Γ)
#     x2 = BEAST.lagrangec0d1(Γ)
#     x3 = buffachristiansen(Γ)
#     x4 = BEAST.lagrangec0d1(Γ)
#     return BEAST.DirectProductSpace([x1,x2,x3,x4])
# end

# function trialbasis(obj::Object{T}, ::PRECVP) where {T <: HOM}
#     Γ =  obj.type.mesh
#     x1 = buffachristiansen(Γ)
#     x2 = BEAST.lagrangec0d1(Γ)
#     x3 = buffachristiansen(Γ)
#     x4 = BEAST.lagrangec0d1(Γ)
#     return BEAST.DirectProductSpace([x1,x2,x3,x4])
# end

function testbasis(obj::Object{T}, ::VP) where {T <: PECEFIE}
    Γ =  obj.type.mesh
    x3 = raviartthomas(Γ)
    x4 = BEAST.lagrangec0d1(Γ)
    return BEAST.DirectProductSpace([x3,x4])
end

function trialbasis(obj::Object{T}, ::VP) where {T <: PECEFIE}
    Γ =  obj.type.mesh
    x3 = raviartthomas(Γ)
    x4 = BEAST.lagrangec0d1(Γ)
    return BEAST.DirectProductSpace([x3,x4])
end
function testbasis(obj::Object{T}, ::VP) where {T <: PECMFIE}
    Γ =  obj.type.mesh
    x3 = buffachristiansen(Γ)
    x4 = BEAST.lagrangec0d1(Γ)
    return BEAST.DirectProductSpace([x3,x4])
end

function trialbasis(obj::Object{T}, ::VP) where {T <: PECMFIE}
    Γ =  obj.type.mesh
    x3 = raviartthomas(Γ)
    x4 = BEAST.lagrangec0d1(Γ)
    return BEAST.DirectProductSpace([x3,x4])
end

function testbasis(obj::Object{T}, ::PRECVP) where {T <: PECEFIE}
    Γ =  obj.type.mesh
    x3 = buffachristiansen(Γ)
    x4 = BEAST.lagrangec0d1(Γ)
    return BEAST.DirectProductSpace([x3,x4])
end

function trialbasis(obj::Object{T}, ::PRECVP) where {T <: PECEFIE}
    Γ =  obj.type.mesh
    x3 = buffachristiansen(Γ)
    x4 = BEAST.lagrangec0d1(Γ)
    return BEAST.DirectProductSpace([x3,x4])
end

function testbasis(obj::Object{T}, ::PRECVP) where {T <: PECMFIE}
    Γ =  obj.type.mesh
    x3 = buffachristiansen(Γ)
    x4 = BEAST.lagrangec0d1(Γ)
    return BEAST.DirectProductSpace([x3,x4])
end

function trialbasis(obj::Object{T}, ::PRECVP) where {T <: PECMFIE}
    Γ =  obj.type.mesh
    x3 = raviartthomas(Γ)
    x4 = BEAST.lagrangec0d1(Γ)
    return BEAST.DirectProductSpace([x3,x4])
end

# function interaction(testobj::Union{Object,Inside}, 
#     trialobj::Inside{Object{<:Union{FreeSpace,HOM}}}, embobj::Object, strat::VP)

#     a = -[1 0 0 0;
#         0 1 0 0;
#         0 0 mu(trialobj)/mu(parent(trialobj)) 0;
#         0 0 0 epsilon(parent(trialobj.object))/epsilon(trialobj)]
#     return interaction(testobj,trialobj.object,embobj,strat)*a
# end


function interaction(testobj::Object,trialobj::Object,embobj::Object,
    testtype::Inside,trialtype,strat::Union{PRECVP,VP,PRECVPPMCHWT2}; sign=1,dual=true,trace=true)
    a = [1 0 0 0;
    0 1 0 0;
    0 0 mu(testobj)/mu(parent(testobj)) 0;
    0 0 0 epsilon(parent(testobj))/epsilon(testobj)]
    a^-1 * interaction(testobj,trialobj,embobj,testtype.inside,trialtype,strat;sign=-sign,dual=dual,trace=trace)
end
function interaction(testobj::Object,trialobj::Object,embobj::Object,
    testtype::Inside,trialtype,strat::Union{VPPMCHWT}; sign=1,dual=true,trace=true)
    a = [1 0 0 0;
    0 epsilon(parent(testobj))/epsilon(testobj) 0 0;
    0 0 mu(testobj)/mu(parent(testobj)) 0;
    0 0 0 1]
    a^-1 * interaction(testobj,trialobj,embobj,testtype.inside,trialtype,strat;sign=-sign,dual=dual,trace=trace)
end
function interaction(testobj::Object,trialobj::Object,embobj::Object,
    testtype::Inside,trialtype,strat::Union{VPPMCHWT2}; sign=1,dual=true,trace=true)
    a = [mu(testobj)/mu(parent(testobj)) 0 0 0;
    0 epsilon(parent(testobj))/epsilon(testobj) 0 0;
    0 0 1 0;
    0 0 0 1]
    a^-1 * interaction(testobj,trialobj,embobj,testtype.inside,trialtype,strat;sign=-sign,dual=dual,trace=trace)
end

function interaction(testobj::Object,trialobj::Object,embobj::Object,
    testtype::PECCFIE,trialtype,strat::Union{PRECVP,VP}; sign=1,dual=true,trace=true)

    testtype.alpha*interaction(testobj,trialobj,embobj,testtype.inside,trialtype,strat;sign=sign)[1:2,:]+
    (1-testtype.alpha)*interaction(testobj,trialobj,embobj,testtype.inside,trialtype,strat;sign=sign,dual=dual,trace=trace)[3:4,:]
end
function interaction(testobj::Object,trialobj::Object,embobj::Object,
    testtype::PECEFIE,trialtype,strat::Union{PRECVP,VP}; sign=1,dual=true,trace=true)

    interaction(testobj,trialobj,embobj,testtype.inside,trialtype,strat;sign=sign,dual=dual,trace=trace)[1:2,:]
end
function interaction(testobj::Object,trialobj::Object,embobj::Object,
    testtype::PECMFIE,trialtype,strat::Union{PRECVP,VP}; sign=1,dual=true,trace=true)
    interaction(testobj,trialobj,embobj,testtype.inside,trialtype,strat;sign=sign,dual=dual,trace=trace)[3:4,:]
end
function interaction(testobj::Object,trialobj::Object,embobj::Object,
    testtype::ObjectType,trialtype::Inside,strat::Union{PRECVP,VP,VPPMCHWT,VPPMCHWT2}; sign=1,dual=true,trace=true)
    a = -[1 0 0 0;
    0 1 0 0;
    0 0 mu(trialobj)/mu(parent(trialobj)) 0;
    0 0 0 epsilon(parent(trialobj))/epsilon(trialobj)]
    interaction(testobj,trialobj,embobj,testtype,trialtype.inside,strat;sign=sign,dual=dual,trace=trace)*a
end
function interaction(testobj::Object,trialobj::Object,embobj::Object,
    testtype::ObjectType,trialtype::Inside,strat::Union{PRECVPPMCHWT2}; sign=1,dual=true,trace=true)
    a = -[mu(trialobj)/mu(parent(trialobj)) 0 0 0;
    0 epsilon(parent(trialobj))/epsilon(trialobj) 0 0;
    0 0 1 0;
    0 0 0 1]
    interaction(testobj,trialobj,embobj,testtype,trialtype.inside,strat;sign=sign,dual=dual,trace=trace)*a
end
function interaction(testobj::Object,trialobj::Object,embobj::Object,
    testtype::ObjectType,trialtype::PEC,strat::Union{PRECVP,VP}; sign=1,dual=true,trace=true)
    interaction(testobj,trialobj,embobj,testtype,trialtype.inside,strat;sign=sign,dual=dual,trace=trace)[:,3:4]
end
function interaction(testobj::Object,trialobj::Object,embobj::Object,
    testtype::ObjectType,trialtype::HOM,strat::Union{PRECVP,VP,VPPMCHWT,VPPMCHWT2,PRECVPPMCHWT2}; sign=1,dual=true,trace=true)
    interaction(testobj,trialobj,embobj,testtype,trialtype.inside,strat;sign=sign,dual=dual,trace=trace)
end
function interaction(testobj::Object,trialobj::Object,embobj::Object,
    testtype::HOM,trialtype,strat::Union{PRECVP,VP,VPPMCHWT,VPPMCHWT2,PRECVPPMCHWT2}; sign=1,dual=true,trace=true)
    interaction(testobj,trialobj,embobj,testtype.inside,trialtype,strat;sign=sign,dual=dual,trace=trace)
end
function interaction(testobj::Object, trialobj::Object, embobj::Object, 
    testtype::ObjectType,trialtype::ObjectType,strat::Union{PRECVP,VP}; sign=1,dual=true,trace=true)
    tr = sign*strat.trace
    k = sqrt(epsilon(embobj)*mu(embobj))*omega(embobj)
    G = BEAST.greenhh3d(wavenumber=k)
    gradG = BEAST.∇(G)

    bs = geometry(trialobj)
    ts = geometry(testobj)
    if trace
        ∇Gx =  BEAST.build_potential(gradG×B,bs)
        Gn = BEAST.build_potential(G*(n*B),bs)
        #Gnx = BEAST.build_potential(G*(n × B),bs)
        ∇G = BEAST.build_potential(gradG*B,bs)
        ∇Gdotn = BEAST.build_potential(nb ⋅ (gradG*B),bs)
        ∇Gdot = BEAST.build_potential(B ⋅ gradG,bs)
        
        Gr = BEAST.build_potential(G*B,bs)
        ∇G∇B = BEAST.build_potential(gradG*div(B),bs)
        ∇Gxn = BEAST.build_potential(gradG×(n*B),bs)
    else
        ∇Gx =  gradG×B
        Gn = G*(nb*B)
        #Gnx = G*(n × B)
        ∇G = gradG*B
        ∇Gdotn = nb ⋅ (gradG*B)
        ∇Gdot = B ⋅ gradG
        
        Gr = G*B
        ∇G∇B = gradG*div(B)
        ∇Gxn = gradG×(nb*B)

    end
    #gebruik gamma^c --> sign = +
    # #without extra ncross
    # int = [γₛ(∇Gx,ts,trace) -γₛ(Gn,ts,trace) γₛ(Gr,ts,trace) -γₛ(∇G,ts,trace);
    # ZeroOperator() -τ(∇Gdotn,ts,trace) τ(∇Gdot,ts,trace)  k^2*τ(Gr,ts,trace);
    # γₛ(∇G∇B,ts,trace)+k^2*γₛ(Gr,ts,trace) -γₛ(∇Gxn,ts,trace) γₛ(∇Gx,ts,trace) ZeroOperator();
    # γₙ(∇Gx,ts,trace) -γₙ(Gn,ts,trace) γₙ(Gr,ts,trace) -γₙ(∇G,ts,trace)]
    # #with extra ncross
    if dual && trace
        int = [-γₜ(∇Gx,ts,tr) γₜ(Gn,ts,tr) -γₜ(Gr,ts,tr) γₜ(∇G,ts,tr);
        ZeroOperator() -τ(∇Gdotn,ts,tr) τ(∇Gdot,ts,tr)  k^2*τ(Gr,ts,tr);
        -γₜ(∇G∇B,ts,tr)-k^2*γₜ(Gr,ts,tr) γₜ(∇Gxn,ts,tr) -γₜ(∇Gx,ts,tr) ZeroOperator();
        γₙ(∇Gx,ts,tr) -γₙ(Gn,ts,tr) γₙ(Gr,ts,tr) -γₙ(∇G,ts,tr)]
        return int
    elseif dual
        int = [nt×(nt×(∇Gx)) -nt×(nt×(Gn)) nt×(nt×(Gr)) nt×(nt×(∇G));
        ZeroOperator() -(∇Gdotn) (∇Gdot)  k^2*(Gr);
        nt×(nt×(∇G∇B))+k^2*nt×(nt×(Gr)) -nt×(nt×(∇Gxn)) nt×(nt×(∇Gx)) ZeroOperator();
        nt ⋅(∇Gx) -nt ⋅(Gn) nt ⋅(Gr) -nt ⋅(∇G)]
        return int
    elseif trace
        int = [γₛ(∇Gx,ts,tr) -γₛ(Gn,ts,tr) γₛ(Gr,ts,tr) -γₛ(∇G,ts,tr);
        ZeroOperator() -τ(∇Gdotn,ts,tr) τ(∇Gdot,ts,tr)  k^2*τ(Gr,ts,tr);
        γₛ(∇G∇B,ts,tr)+k^2*γₛ(Gr,ts,tr) -γₛ(∇Gxn,ts,tr) γₛ(∇Gx,ts,tr) ZeroOperator();
        γₙ(∇Gx,ts,tr) -γₙ(Gn,ts,tr) γₙ(Gr,ts,tr) -γₙ(∇G,ts,tr)]
        return int
    else
        int = [nt×(∇Gx) -nt×(Gn) nt×(Gr) -nt×(∇G);
        ZeroOperator() -(∇Gdotn) (∇Gdot)  k^2*(Gr);
        nt×(∇G∇B)+k^2*nt×(Gr) -nt×(∇Gxn) nt×(∇Gx) ZeroOperator();
        nt ⋅(∇Gx) -nt ⋅(Gn) nt ⋅(Gr) -nt ⋅(∇G)]
        return int
    end
end
function interaction(testobj::Object, trialobj::Object, embobj::Object, 
    testtype::ObjectType,trialtype::ObjectType,strat::Union{PRECVPPMCHWT2}; sign=1,dual=true,trace=true)
    tr = sign*strat.trace
    k = sqrt(epsilon(embobj)*mu(embobj))*omega(embobj)
    G = BEAST.greenhh3d(wavenumber=k)
    gradG = BEAST.∇(G)
    switch_collums(mat) = [mat[:,3:4] mat[:,1:2]]

    zer = Matrix{AbstractOperator}(undef,4,4)
    fill!(zer,ZeroOperator())
    rows = strat.rows
    cols = strat.cols
    
    bs = geometry(trialobj)
    ts = geometry(testobj)
    if trace
        ∇Gx =  BEAST.build_potential(gradG×B,bs)
        Gn = BEAST.build_potential(G*(n*B),bs)
        #Gnx = BEAST.build_potential(G*(n × B),bs)
        ∇G = BEAST.build_potential(gradG*B,bs)
        ∇Gdotn = BEAST.build_potential(nb ⋅ (gradG*B),bs)
        ∇Gdot = BEAST.build_potential(B ⋅ gradG,bs)
        
        Gr = BEAST.build_potential(G*B,bs)
        ∇G∇B = BEAST.build_potential(gradG*div(B),bs)
        ∇Gxn = BEAST.build_potential(gradG×(n*B),bs)
    else
        ∇Gx =  gradG×B
        Gn = G*(nb*B)
        #Gnx = G*(n × B)
        ∇G = gradG*B
        ∇Gdotn = nb ⋅ (gradG*B)
        ∇Gdot = B ⋅ gradG
        
        Gr = G*B
        ∇G∇B = gradG*div(B)
        ∇Gxn = gradG×(nb*B)

    end
    #gebruik gamma^c --> sign = +
    # #without extra ncross
    # int = [γₛ(∇Gx,ts,trace) -γₛ(Gn,ts,trace) γₛ(Gr,ts,trace) -γₛ(∇G,ts,trace);
    # ZeroOperator() -τ(∇Gdotn,ts,trace) τ(∇Gdot,ts,trace)  k^2*τ(Gr,ts,trace);
    # γₛ(∇G∇B,ts,trace)+k^2*γₛ(Gr,ts,trace) -γₛ(∇Gxn,ts,trace) γₛ(∇Gx,ts,trace) ZeroOperator();
    # γₙ(∇Gx,ts,trace) -γₙ(Gn,ts,trace) γₙ(Gr,ts,trace) -γₙ(∇G,ts,trace)]
    # #with extra ncross
    if dual && trace
        int = [-γₜ(∇Gx,ts,tr) γₜ(Gn,ts,tr) -γₜ(Gr,ts,tr) γₜ(∇G,ts,tr);
        ZeroOperator() -τ(∇Gdotn,ts,tr) τ(∇Gdot,ts,tr)  k^2*τ(Gr,ts,tr);
        -γₜ(∇G∇B,ts,tr)-k^2*γₜ(Gr,ts,tr) γₜ(∇Gxn,ts,tr) -γₜ(∇Gx,ts,tr) ZeroOperator();
        γₙ(∇Gx,ts,tr) -γₙ(Gn,ts,tr) γₙ(Gr,ts,tr) -γₙ(∇G,ts,tr)]
        out = switch_collums(int)
    elseif dual
        int = [nt×(nt×(∇Gx)) -nt×(nt×(Gn)) nt×(nt×(Gr)) nt×(nt×(∇G));
        ZeroOperator() -(∇Gdotn) (∇Gdot)  k^2*(Gr);
        nt×(nt×(∇G∇B))+k^2*nt×(nt×(Gr)) -nt×(nt×(∇Gxn)) nt×(nt×(∇Gx)) ZeroOperator();
        nt ⋅(∇Gx) -nt ⋅(Gn) nt ⋅(Gr) -nt ⋅(∇G)]
        out = switch_collums(int)
    elseif trace
        int = [γₛ(∇Gx,ts,tr) -γₛ(Gn,ts,tr) γₛ(Gr,ts,tr) -γₛ(∇G,ts,tr);
        ZeroOperator() -τ(∇Gdotn,ts,tr) τ(∇Gdot,ts,tr)  k^2*τ(Gr,ts,tr);
        γₛ(∇G∇B,ts,tr)+k^2*γₛ(Gr,ts,tr) -γₛ(∇Gxn,ts,tr) γₛ(∇Gx,ts,tr) ZeroOperator();
        γₙ(∇Gx,ts,tr) -γₙ(Gn,ts,tr) γₙ(Gr,ts,tr) -γₙ(∇G,ts,tr)]
        out = switch_collums(int)
    else
        int = [nt×(∇Gx) -nt×(Gn) nt×(Gr) -nt×(∇G);
        ZeroOperator() -(∇Gdotn) (∇Gdot)  k^2*(Gr);
        nt×(∇G∇B)+k^2*nt×(Gr) -nt×(∇Gxn) nt×(∇Gx) ZeroOperator();
        nt ⋅(∇Gx) -nt ⋅(Gn) nt ⋅(Gr) -nt ⋅(∇G)]
        out = switch_collums(int)
    end
    zer[rows,cols] .+= out[rows,cols]
    return zer
end
function interaction(testobj::Object, trialobj::Object, embobj::Object, 
    testtype::ObjectType,trialtype::ObjectType,strat::VPPMCHWT; sign=1,dual=true,trace=true)
    tr = sign*strat.trace
    k = sqrt(epsilon(embobj)*mu(embobj))*omega(embobj)
    G = BEAST.greenhh3d(wavenumber=k)
    gradG = BEAST.∇(G)

    bs = geometry(trialobj)
    ts = geometry(testobj)
    if trace
        ∇Gx =  BEAST.build_potential(gradG×B,bs)
        Gn = BEAST.build_potential(G*(n*B),bs)
        #Gnx = BEAST.build_potential(G*(n × B),bs)
        ∇G = BEAST.build_potential(gradG*B,bs)
        ∇Gdotn = BEAST.build_potential(nb ⋅ (gradG*B),bs)
        ∇Gdot = BEAST.build_potential(B ⋅ gradG,bs)
        
        Gr = BEAST.build_potential(G*B,bs)
        ∇G∇B = BEAST.build_potential(gradG*div(B),bs)
        ∇Gxn = BEAST.build_potential(gradG×(n*B),bs)
    else
        ∇Gx =  gradG×B
        Gn = G*(nb*B)
        #Gnx = G*(n × B)
        ∇G = gradG*B
        ∇Gdotn = nb ⋅ (gradG*B)
        ∇Gdot = B ⋅ gradG
        
        Gr = G*B
        ∇G∇B = gradG*div(B)
        ∇Gxn = gradG×(nb*B)

    end
    #gebruik gamma^c --> sign = +
    # #without extra ncross
    # int = [γₛ(∇Gx,ts,trace) -γₛ(Gn,ts,trace) γₛ(Gr,ts,trace) -γₛ(∇G,ts,trace);
    # ZeroOperator() -τ(∇Gdotn,ts,trace) τ(∇Gdot,ts,trace)  k^2*τ(Gr,ts,trace);
    # γₛ(∇G∇B,ts,trace)+k^2*γₛ(Gr,ts,trace) -γₛ(∇Gxn,ts,trace) γₛ(∇Gx,ts,trace) ZeroOperator();
    # γₙ(∇Gx,ts,trace) -γₙ(Gn,ts,trace) γₙ(Gr,ts,trace) -γₙ(∇G,ts,trace)]
    # #with extra ncross
    if dual && trace
        int = [-γₜ(∇Gx,ts,tr) γₜ(Gn,ts,tr) -γₜ(Gr,ts,tr) γₜ(∇G,ts,tr);
        γₙ(∇Gx,ts,tr) -γₙ(Gn,ts,tr) γₙ(Gr,ts,tr) -γₙ(∇G,ts,tr);
        -γₜ(∇G∇B,ts,tr)-k^2*γₜ(Gr,ts,tr) γₜ(∇Gxn,ts,tr) -γₜ(∇Gx,ts,tr) ZeroOperator();
        ZeroOperator() -τ(∇Gdotn,ts,tr) τ(∇Gdot,ts,tr)  k^2*τ(Gr,ts,tr)]
        return int
    elseif dual
        int = [nt×(nt×(∇Gx)) -nt×(nt×(Gn)) nt×(nt×(Gr)) nt×(nt×(∇G));
        nt ⋅(∇Gx) -nt ⋅(Gn) nt ⋅(Gr) -nt ⋅(∇G);
        nt×(nt×(∇G∇B))+k^2*nt×(nt×(Gr)) -nt×(nt×(∇Gxn)) nt×(nt×(∇Gx)) ZeroOperator();
        ZeroOperator() -(∇Gdotn) (∇Gdot)  k^2*(Gr)]
        return int
    elseif trace
        int = [γₛ(∇Gx,ts,tr) -γₛ(Gn,ts,tr) γₛ(Gr,ts,tr) -γₛ(∇G,ts,tr);
        γₙ(∇Gx,ts,tr) -γₙ(Gn,ts,tr) γₙ(Gr,ts,tr) -γₙ(∇G,ts,tr);
        γₛ(∇G∇B,ts,tr)+k^2*γₛ(Gr,ts,tr) -γₛ(∇Gxn,ts,tr) γₛ(∇Gx,ts,tr) ZeroOperator();
        ZeroOperator() -τ(∇Gdotn,ts,tr) τ(∇Gdot,ts,tr)  k^2*τ(Gr,ts,tr)]
        return int
    else
        int = [nt×(∇Gx) -nt×(Gn) nt×(Gr) -nt×(∇G);
        nt ⋅(∇Gx) -nt ⋅(Gn) nt ⋅(Gr) -nt ⋅(∇G);
        nt×(∇G∇B)+k^2*nt×(Gr) -nt×(∇Gxn) nt×(∇Gx) ZeroOperator();
        ZeroOperator() -(∇Gdotn) (∇Gdot)  k^2*(Gr)]
        return int
    end
end
function interaction(testobj::Object, trialobj::Object, embobj::Object, 
    testtype::ObjectType,trialtype::ObjectType,strat::VPPMCHWT2; sign=1,dual=true,trace=true)
    tr = sign*strat.trace
    k = sqrt(epsilon(embobj)*mu(embobj))*omega(embobj)
    G = BEAST.greenhh3d(wavenumber=k)
    gradG = BEAST.∇(G)

    bs = geometry(trialobj)
    ts = geometry(testobj)
    if trace
        ∇Gx =  BEAST.build_potential(gradG×B,bs)
        Gn = BEAST.build_potential(G*(n*B),bs)
        #Gnx = BEAST.build_potential(G*(n × B),bs)
        ∇G = BEAST.build_potential(gradG*B,bs)
        ∇Gdotn = BEAST.build_potential(nb ⋅ (gradG*B),bs)
        ∇Gdot = BEAST.build_potential(B ⋅ gradG,bs)
        
        Gr = BEAST.build_potential(G*B,bs)
        ∇G∇B = BEAST.build_potential(gradG*div(B),bs)
        ∇Gxn = BEAST.build_potential(gradG×(n*B),bs)
    else
        ∇Gx =  gradG×B
        Gn = G*(nb*B)
        #Gnx = G*(n × B)
        ∇G = gradG*B
        ∇Gdotn = nb ⋅ (gradG*B)
        ∇Gdot = B ⋅ gradG
        
        Gr = G*B
        ∇G∇B = gradG*div(B)
        ∇Gxn = gradG×(nb*B)

    end
    #gebruik gamma^c --> sign = +
    # #without extra ncross
    # int = [γₛ(∇Gx,ts,trace) -γₛ(Gn,ts,trace) γₛ(Gr,ts,trace) -γₛ(∇G,ts,trace);
    # ZeroOperator() -τ(∇Gdotn,ts,trace) τ(∇Gdot,ts,trace)  k^2*τ(Gr,ts,trace);
    # γₛ(∇G∇B,ts,trace)+k^2*γₛ(Gr,ts,trace) -γₛ(∇Gxn,ts,trace) γₛ(∇Gx,ts,trace) ZeroOperator();
    # γₙ(∇Gx,ts,trace) -γₙ(Gn,ts,trace) γₙ(Gr,ts,trace) -γₙ(∇G,ts,trace)]
    # #with extra ncross
    if dual && trace
        int = [-γₜ(∇G∇B,ts,tr)-k^2*γₜ(Gr,ts,tr) γₜ(∇Gxn,ts,tr) -γₜ(∇Gx,ts,tr) ZeroOperator();
        γₙ(∇Gx,ts,tr) -γₙ(Gn,ts,tr) γₙ(Gr,ts,tr) -γₙ(∇G,ts,tr);
        -γₜ(∇Gx,ts,tr) γₜ(Gn,ts,tr) -γₜ(Gr,ts,tr) γₜ(∇G,ts,tr);
        ZeroOperator() -τ(∇Gdotn,ts,tr) τ(∇Gdot,ts,tr)  k^2*τ(Gr,ts,tr)]
        return int
    elseif dual
        int = [nt×(nt×(∇G∇B))+k^2*nt×(nt×(Gr)) -nt×(nt×(∇Gxn)) nt×(nt×(∇Gx)) ZeroOperator();
        nt ⋅(∇Gx) -nt ⋅(Gn) nt ⋅(Gr) -nt ⋅(∇G);
        nt×(nt×(∇Gx)) -nt×(nt×(Gn)) nt×(nt×(Gr)) nt×(nt×(∇G));
        ZeroOperator() -(∇Gdotn) (∇Gdot)  k^2*(Gr)]
        return int
    elseif trace
        int = [γₛ(∇G∇B,ts,tr)+k^2*γₛ(Gr,ts,tr) -γₛ(∇Gxn,ts,tr) γₛ(∇Gx,ts,tr) ZeroOperator();
        γₙ(∇Gx,ts,tr) -γₙ(Gn,ts,tr) γₙ(Gr,ts,tr) -γₙ(∇G,ts,tr);
        γₛ(∇Gx,ts,tr) -γₛ(Gn,ts,tr) γₛ(Gr,ts,tr) -γₛ(∇G,ts,tr);
        ZeroOperator() -τ(∇Gdotn,ts,tr) τ(∇Gdot,ts,tr)  k^2*τ(Gr,ts,tr)]
        return int
    else
        int = [nt×(∇G∇B)+k^2*nt×(Gr) -nt×(∇Gxn) nt×(∇Gx) ZeroOperator();
        nt ⋅(∇Gx) -nt ⋅(Gn) nt ⋅(Gr) -nt ⋅(∇G);
        nt×(∇Gx) -nt×(Gn) nt×(Gr) -nt×(∇G);
        ZeroOperator() -(∇Gdotn) (∇Gdot)  k^2*(Gr)]
        return int
    end
end
# function interaction(testobj::Object, trialobj::Object, embobj::Object, 
#     testtype::ObjectType,trialtype::ObjectType,strat::VP; sign=1)
#     trace = sign*strat.trace
#     k = sqrt(epsilon(embobj)*mu(embobj))*omega(embobj)
#     G = BEAST.greenhh3d(wavenumber=k)
#     gradG = BEAST.∇(G)


#     ∇Gx =  BEAST.build_potential(gradG×B)
#     Gn = BEAST.build_potential(G*(n*B))
#     #Gnx = BEAST.build_potential(G*(n × B),bs)
#     ∇G = BEAST.build_potential(gradG*B)
#     ∇Gdotn = BEAST.build_potential(nb ⋅ (gradG*B))
#     ∇Gdot = BEAST.build_potential(B ⋅ gradG)
    
#     Gr = BEAST.build_potential(G*B)
#     ∇G∇B = BEAST.build_potential(gradG*div(B))
#     ∇Gxn = BEAST.build_potential(gradG×(n*B))

#     #gebruik gamma^c --> sign = +
#     # #without extra ncross
#     # int = [γₛ(∇Gx,ts,trace) -γₛ(Gn,ts,trace) γₛ(Gr,ts,trace) -γₛ(∇G,ts,trace);
#     # ZeroOperator() -τ(∇Gdotn,ts,trace) τ(∇Gdot,ts,trace)  k^2*τ(Gr,ts,trace);
#     # γₛ(∇G∇B,ts,trace)+k^2*γₛ(Gr,ts,trace) -γₛ(∇Gxn,ts,trace) γₛ(∇Gx,ts,trace) ZeroOperator();
#     # γₙ(∇Gx,ts,trace) -γₙ(Gn,ts,trace) γₙ(Gr,ts,trace) -γₙ(∇G,ts,trace)]
#     # #with extra ncross
#     int = [-γₜ(∇Gx,trace) γₜ(Gn,trace) -γₜ(Gr,trace) γₜ(∇G,trace);
#     ZeroOperator() -τ(∇Gdotn,trace) τ(∇Gdot,trace)  k^2*τ(Gr,trace);
#     -γₜ(∇G∇B,trace)-k^2*γₜ(Gr,trace) γₜ(∇Gxn,trace) -γₜ(∇Gx,trace) ZeroOperator();
#     γₙ(∇Gx,trace) -γₙ(Gn,trace) γₙ(Gr,trace) -γₙ(∇G,trace)]
#     return int
# end

# function interaction(testobj::Object, trialobj::Object, embobj::Object, 
#     testtype::ObjectType,trialtype::ObjectType,strat::PRECVP; sign=1,dual=true)
#     trace = sign*strat.trace
#     k = sqrt(epsilon(embobj)*mu(embobj))*omega(embobj)
#     G = BEAST.greenhh3d(wavenumber=k)
#     gradG = BEAST.∇(G)

#     bs = barycentric_refinement(geometry(trialobj); sort = :spacefillingcurve)
#     ts = barycentric_refinement(geometry(testobj); sort = :spacefillingcurve)

#     bss = geometry(trialobj)
#     tss = geometry(testobj)

#     ∇Gx =  BEAST.build_potential(gradG×B,bs)#
#     Gn = BEAST.build_potential(G*(n*B),bss)#
#     #Gnx = BEAST.build_potential(G*(n × B),bs)
#     ∇G = BEAST.build_potential(gradG*B,bss)#
#     ∇Gdotn = BEAST.build_potential(nb ⋅ (gradG*B),bss)#
#     ∇Gdot = BEAST.build_potential(B ⋅ gradG,bs)#
    
#     Gr = BEAST.build_potential(G*B,bs)#
#     Grs = BEAST.build_potential(G*B,bss)#
#     ∇G∇B = BEAST.build_potential(gradG*div(B),bs)#
#     ∇Gxn = BEAST.build_potential(gradG×(n*B),bss)#

#     int = [-γₜ(∇Gx,ts,trace) γₜ(Gn,ts,trace) -γₜ(Gr,ts,trace) γₜ(∇G,ts,trace);
#     ZeroOperator() -τ(∇Gdotn,tss,trace) τ(∇Gdot,tss,trace)  k^2*τ(Grs,tss,trace);
#     -γₜ(∇G∇B,ts,trace)-k^2*γₜ(Gr,ts,trace) γₜ(∇Gxn,ts,trace) -γₜ(∇Gx,ts,trace) ZeroOperator();
#     γₙ(∇Gx,tss,trace) -γₙ(Gn,tss,trace) γₙ(Gr,tss,trace) -γₙ(∇G,tss,trace)]
#     return int
# end
function _identity(::ObjectType,strat::PRECVP)
    [NCross() ZeroOperator() ZeroOperator() ZeroOperator();
    ZeroOperator() Identity() ZeroOperator() ZeroOperator();
    ZeroOperator() ZeroOperator() NCross() ZeroOperator();
    ZeroOperator() ZeroOperator() ZeroOperator() Identity()]
end
_identity(p::PEC,strat::PRECVP) = _identity(p.inside,strat)[1:2,1:2]

function identity(objtype::ObjectType,strat::Union{VP,VPPMCHWT,VPPMCHWT2,PRECVPPMCHWT2})
    [NCross() ZeroOperator() ZeroOperator() ZeroOperator();
    ZeroOperator() Identity() ZeroOperator() ZeroOperator();
    ZeroOperator() ZeroOperator() NCross() ZeroOperator();
    ZeroOperator() ZeroOperator() ZeroOperator() Identity()]
end
function identity(objtype::ObjectType,strat::PRECVP)
    [NCross() ZeroOperator() ZeroOperator() ZeroOperator();
    ZeroOperator() Identity() ZeroOperator() ZeroOperator();
    ZeroOperator() ZeroOperator() NCross() ZeroOperator();
    ZeroOperator() ZeroOperator() ZeroOperator() Identity()]
end
function identity(objtype::PECCFIE,strat::Union{PRECVP,VP})
    (1-objtype.alpha)*identity(objtype.inside,strat)[3:4,3:4]
end
function identity(objtype::PECMFIE,strat::Union{PRECVP,VP})
    identity(objtype.inside,strat)[3:4,3:4]
end
function identity(objtype::PECEFIE,strat::Union{PRECVP,VP})
    [ZeroOperator() ZeroOperator();
    ZeroOperator() ZeroOperator()]
end
identity(o::HOM,strat::Union{PRECVP,VP,VPPMCHWT,VPPMCHWT2,PRECVPPMCHWT2}) = identity(o.inside,strat)
# function interaction(testobj::Inside{Object{<:Union{FreeSpace,HOM}}}, 
#     trialobj::Object, embobj::Object, 
#     strat::VP)

#     trace = -strat.trace
#     a = [1 0 0 0;
#         0 1 0 0;
#         0 0 mu(testobj)/mu(parent(testobj)) 0;
#         0 0 0 epsilon(parent(testobj.object))/epsilon(testobj)]

#         k = sqrt(epsilon(embobj)*mu(embobj))*omega(embobj)
#         G = BEAST.greenhh3d(wavenumber=k)
#         gradG = BEAST.∇(G)
    
#         bs = geometry(trialobj)
#         ts = geometry(testobj)
    
#         ∇Gx =  BEAST.build_potential(gradG×B,bs)
#         Gn = BEAST.build_potential(G*(n*B),bs)
#         #Gnx = BEAST.build_potential(G*(n × B),bs)
#         ∇G = BEAST.build_potential(gradG*B,bs)
#         ∇Gdotn = BEAST.build_potential(nb ⋅ (gradG*B),bs)
#         ∇Gdot = BEAST.build_potential(B ⋅ gradG,bs)
        
#         Gr = BEAST.build_potential(G*B,bs)
#         ∇G∇B = BEAST.build_potential(gradG*div(B),bs)
#         ∇Gxn = BEAST.build_potential(gradG×(n*B),bs)
    
#         #gebruik gamma^c --> sign = -
#         int = [γₛ(∇Gx,ts,trace) -γₛ(Gn,ts,trace) γₛ(Gr,ts,trace) -γₛ(∇G,ts,trace);
#         ZeroOperator() -τ(∇Gdotn,ts,trace) τ(∇Gdot,ts,trace)  k^2*τ(Gr,ts,trace);
#         γₛ(∇G∇B,ts,trace)+k^2*γₛ(Gr,ts,trace) -γₛ(∇Gxn,ts,trace) γₛ(∇Gx,ts,trace) ZeroOperator();
#         γₙ(∇Gx,ts,trace) -γₙ(Gn,ts,trace) γₙ(Gr,ts,trace) -γₙ(∇G,ts,trace)]
#         return (a^-1)*int
# end

# function interaction(testobj::Object, trialobj::Object, embobj::Object, strat::VP)
#     trace = strat.trace
#     k = sqrt(epsilon(embobj)*mu(embobj))*omega(embobj)
#     G = BEAST.greenhh3d(wavenumber=k)
#     gradG = BEAST.∇(G)

#     bs = geometry(trialobj)
#     ts = geometry(testobj)

#     ∇Gx =  BEAST.build_potential(gradG×B,bs)
#     Gn = BEAST.build_potential(G*(n*B),bs)
#     #Gnx = BEAST.build_potential(G*(n × B),bs)
#     ∇G = BEAST.build_potential(gradG*B,bs)
#     ∇Gdotn = BEAST.build_potential(nb ⋅ (gradG*B),bs)
#     ∇Gdot = BEAST.build_potential(B ⋅ gradG,bs)
    
#     Gr = BEAST.build_potential(G*B,bs)
#     ∇G∇B = BEAST.build_potential(gradG*div(B),bs)
#     ∇Gxn = BEAST.build_potential(gradG×(n*B),bs)

#     #gebruik gamma^c --> sign = +
#     int = [γₛ(∇Gx,ts,trace) -γₛ(Gn,ts,trace) γₛ(Gr,ts,trace) -γₛ(∇G,ts,trace);
#     ZeroOperator() -τ(∇Gdotn,ts,trace) τ(∇Gdot,ts,trace)  k^2*τ(Gr,ts,trace);
#     γₛ(∇G∇B,ts,trace)+k^2*γₛ(Gr,ts,trace) -γₛ(∇Gxn,ts,trace) γₛ(∇Gx,ts,trace) ZeroOperator();
#     γₙ(∇Gx,ts,trace) -γₙ(Gn,ts,trace) γₙ(Gr,ts,trace) -γₙ(∇G,ts,trace)]
#     return int

# end

function excitation(obj::Object,emb::Object,objtype::ObjectType,ex::VPExcitation,strat::VP)
    emb.index ∉ ex.objectids && return [0,0,0,0]
    a = [-1* ((n × FunctionExcitation{ComplexF64}(ex.A))×n),
    FunctionExcitation{ComplexF64}(ex.divA),
    -1* ((n × FunctionExcitation{ComplexF64}(ex.curlA))×n),
    NDotTrace{ComplexF64}(FunctionExcitation{ComplexF64}(ex.A))]
    @warn "fix this type"
    return a
end
function excitation(obj::Object,emb::Object,objtype::ObjectType,ex::VPExcitation,strat::VPPMCHWT)
    emb.index ∉ ex.objectids && return [0,0,0,0]
    a = [-1* ((n × FunctionExcitation{ComplexF64}(ex.A))×n),
    NDotTrace{ComplexF64}(FunctionExcitation{ComplexF64}(ex.A)),
    -1* ((n × FunctionExcitation{ComplexF64}(ex.curlA))×n),
    FunctionExcitation{ComplexF64}(ex.divA)]
    @warn "fix this type"
    return a
end
function excitation(obj::Object,emb::Object,objtype::ObjectType,ex::VPExcitation,strat::VPPMCHWT2)
    emb.index ∉ ex.objectids && return [0,0,0,0]
    a = [-1* ((n × FunctionExcitation{ComplexF64}(ex.curlA))×n),
    NDotTrace{ComplexF64}(FunctionExcitation{ComplexF64}(ex.A)),
    -1* ((n × FunctionExcitation{ComplexF64}(ex.A))×n),
    FunctionExcitation{ComplexF64}(ex.divA)]
    @warn "fix this type"
    return a
end
function excitation(testobj::Object,emb::Object,objtype::Inside,ex::VPExcitation,strat::VP)
    a = [1 0 0 0;
    0 1 0 0;
    0 0 mu(testobj)/mu(parent(testobj)) 0;
    0 0 0 epsilon(parent(testobj))/epsilon(testobj)]
    
    return a^-1*excitation(testobj,emb,objtype.inside,ex,strat)
end
function excitation(testobj::Object,emb::Object,objtype::Inside,ex::VPExcitation,strat::VPPMCHWT)
    a = [1 0 0 0;
    0 epsilon(parent(testobj))/epsilon(testobj) 0 0;
    0 0 mu(testobj)/mu(parent(testobj)) 0;
    0 0 0 1]
    
    return a^-1*excitation(testobj,emb,objtype.inside,ex,strat)
end
function excitation(testobj::Object,emb::Object,objtype::Inside,ex::VPExcitation,strat::VPPMCHWT2)
    a = [mu(testobj)/mu(parent(testobj)) 0 0 0;
    0 epsilon(parent(testobj))/epsilon(testobj) 0 0;
    0 0 1 0;
    0 0 0 1]
    
    return a^-1*excitation(testobj,emb,objtype.inside,ex,strat)
end
function excitation(testobj::Object,emb::Object,objtype::PECCFIE,ex::VPExcitation,strat::VP)

    objtype.alpha*excitation(testobj,emb,objtype.inside,ex,strat)[1:2]+
    (1-objtype.alpha)*excitation(testobj,emb,objtype.inside,ex,strat)[3:4]
end
function excitation(testobj::Object,emb::Object,objtype::PECEFIE,ex::VPExcitation,strat::VP)

    excitation(testobj,emb,objtype.inside,ex,strat)[1:2]
end
function excitation(testobj::Object,emb::Object,objtype::PECMFIE,ex::VPExcitation,strat::VP)

    excitation(testobj,emb,objtype.inside,ex,strat)[3:4]
end
excitation(testobj::Object,emb::Object,objtype::HOM,ex::VPExcitation,strat::Union{VP,VPPMCHWT,VPPMCHWT2}) = excitation(testobj,emb,objtype.inside,ex,strat)


### Post Processing
function list_of_operators(embobj::Object,basisobj::Object{<:HOM},::AField,::VP)
    k = sqrt(epsilon(embobj)*mu(embobj))*omega(embobj)
    G = BEAST.greenhh3d(wavenumber=k)
    gradG = BEAST.∇(G)

    ∇Gx =  gradG×B
    Gn = G*(nb*B)
    ∇G = gradG*B
    Gr = G*B

    return [∇Gx, -Gn, Gr, -∇G]
end

function list_of_operators(embobj::Object,basisobj::Object{<:HOM}, ::BField, ::VP)
    k = sqrt(epsilon(embobj)*mu(embobj))*omega(embobj)
    G = BEAST.greenhh3d(wavenumber=k)
    gradG = BEAST.∇(G)

    ∇Gx =  gradG×B
    Gr = G*B
    ∇G∇B = gradG*div(B)
    ∇Gxn = gradG×(nb*B)
    return [∇G∇B+k^2*Gr, -∇Gxn, ∇Gx, ZeroOperator()]
end
function list_of_operators(embobj::Object,basisobj::Object{<:PEC},::AField,::VP)
    k = sqrt(epsilon(embobj)*mu(embobj))*omega(embobj)
    G = BEAST.greenhh3d(wavenumber=k)
    gradG = BEAST.∇(G)

    ∇Gx =  gradG×B
    Gn = G*(nb*B)
    ∇G = gradG*B
    Gr = G*B

    return [Gr, -∇G]
end

function list_of_operators(embobj::Object,basisobj::Object{<:PEC}, ::BField, ::VP)
    k = sqrt(epsilon(embobj)*mu(embobj))*omega(embobj)
    G = BEAST.greenhh3d(wavenumber=k)
    gradG = BEAST.∇(G)

    ∇Gx =  gradG×B
    Gr = G*B
    ∇G∇B = gradG*div(B)
    ∇Gxn = gradG×(nb*B)
    return [∇Gx, ZeroOperator()]
end
function transform(obj,strat::VP)
    a = -[1 0 0 0;
    0 1 0 0;
    0 0 mu(obj)/mu(parent(obj)) 0;
    0 0 0 epsilon(parent(obj))/epsilon(obj)]
    return a
end


