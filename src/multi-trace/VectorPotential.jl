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
PRECVP() = PRECVP(1)

function testbasis(obj::Object{T}, ::VP) where {T <: HOM}
    Γ =  obj.type.mesh
    x1 = raviartthomas(Γ)
    x2 = BEAST.lagrangec0d1(Γ)
    x3 = raviartthomas(Γ)
    x4 = BEAST.lagrangec0d1(Γ)
    return BEAST.DirectProductSpace([x1,x2,x3,x4])
end

function trialbasis(obj::Object{T}, ::VP) where {T <: HOM}
    Γ =  obj.type.mesh
    x1 = raviartthomas(Γ)
    x2 = BEAST.lagrangec0d1(Γ)
    x3 = raviartthomas(Γ)
    x4 = BEAST.lagrangec0d1(Γ)
    return BEAST.DirectProductSpace([x1,x2,x3,x4])
end
function testbasis(obj::Object{T}, ::PRECVP) where {T <: HOM}
    Γ =  obj.type.mesh
    x1 = buffachristiansen(Γ)
    x2 = BEAST.lagrangec0d1(Γ)
    x3 = buffachristiansen(Γ)
    x4 = BEAST.lagrangec0d1(Γ)
    return BEAST.DirectProductSpace([x1,x2,x3,x4])
end

function trialbasis(obj::Object{T}, ::PRECVP) where {T <: HOM}
    Γ =  obj.type.mesh
    x1 = buffachristiansen(Γ)
    x2 = BEAST.lagrangec0d1(Γ)
    x3 = buffachristiansen(Γ)
    x4 = BEAST.lagrangec0d1(Γ)
    return BEAST.DirectProductSpace([x1,x2,x3,x4])
end

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
    testtype::Inside,trialtype,strat::Union{PRECVP,VP}; sign=1)
    a = [1 0 0 0;
    0 1 0 0;
    0 0 mu(testobj)/mu(parent(testobj)) 0;
    0 0 0 epsilon(parent(testobj))/epsilon(testobj)]
    a^-1 * interaction(testobj,trialobj,embobj,testtype.inside,trialtype,strat;sign=-sign)
end

function interaction(testobj::Object,trialobj::Object,embobj::Object,
    testtype::PECCFIE,trialtype,strat::Union{PRECVP,VP}; sign=1)

    testtype.alpha*interaction(testobj,trialobj,embobj,testtype.inside,trialtype,strat;sign=sign)[1:2,:]+
    (1-testtype.alpha)*interaction(testobj,trialobj,embobj,testtype.inside,trialtype,strat;sign=sign)[3:4,:]
end
function interaction(testobj::Object,trialobj::Object,embobj::Object,
    testtype::PECEFIE,trialtype,strat::Union{PRECVP,VP}; sign=1)

    interaction(testobj,trialobj,embobj,testtype.inside,trialtype,strat;sign=sign)[1:2,:]
end
function interaction(testobj::Object,trialobj::Object,embobj::Object,
    testtype::PECMFIE,trialtype,strat::Union{PRECVP,VP}; sign=1)
    interaction(testobj,trialobj,embobj,testtype.inside,trialtype,strat;sign=sign)[3:4,:]
end
function interaction(testobj::Object,trialobj::Object,embobj::Object,
    testtype::ObjectType,trialtype::Inside,strat::Union{PRECVP,VP}; sign=1)
    a = -[1 0 0 0;
    0 1 0 0;
    0 0 mu(trialobj)/mu(parent(trialobj)) 0;
    0 0 0 epsilon(parent(trialobj))/epsilon(trialobj)]
    interaction(testobj,trialobj,embobj,testtype,trialtype.inside,strat;sign=sign)*a
end
function interaction(testobj::Object,trialobj::Object,embobj::Object,
    testtype::ObjectType,trialtype::PEC,strat::Union{PRECVP,VP}; sign=1)
    interaction(testobj,trialobj,embobj,testtype,trialtype.inside,strat;sign=sign)[:,3:4]
end
function interaction(testobj::Object,trialobj::Object,embobj::Object,
    testtype::ObjectType,trialtype::HOM,strat::Union{PRECVP,VP}; sign=1)
    interaction(testobj,trialobj,embobj,testtype,trialtype.inside,strat;sign=sign)
end
function interaction(testobj::Object,trialobj::Object,embobj::Object,
    testtype::HOM,trialtype,strat::Union{PRECVP,VP}; sign=1)
    interaction(testobj,trialobj,embobj,testtype.inside,trialtype,strat;sign=sign)
end
# function interaction(testobj::Object, trialobj::Object, embobj::Object, 
#     testtype::ObjectType,trialtype::ObjectType,strat::VP; sign=1)
#     trace = sign*strat.trace
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
#     # #without extra ncross
#     # int = [γₛ(∇Gx,ts,trace) -γₛ(Gn,ts,trace) γₛ(Gr,ts,trace) -γₛ(∇G,ts,trace);
#     # ZeroOperator() -τ(∇Gdotn,ts,trace) τ(∇Gdot,ts,trace)  k^2*τ(Gr,ts,trace);
#     # γₛ(∇G∇B,ts,trace)+k^2*γₛ(Gr,ts,trace) -γₛ(∇Gxn,ts,trace) γₛ(∇Gx,ts,trace) ZeroOperator();
#     # γₙ(∇Gx,ts,trace) -γₙ(Gn,ts,trace) γₙ(Gr,ts,trace) -γₙ(∇G,ts,trace)]
#     # #with extra ncross
#     int = [-γₜ(∇Gx,ts,trace) γₜ(Gn,ts,trace) -γₜ(Gr,ts,trace) γₜ(∇G,ts,trace);
#     ZeroOperator() -τ(∇Gdotn,ts,trace) τ(∇Gdot,ts,trace)  k^2*τ(Gr,ts,trace);
#     -γₜ(∇G∇B,ts,trace)-k^2*γₜ(Gr,ts,trace) γₜ(∇Gxn,ts,trace) -γₜ(∇Gx,ts,trace) ZeroOperator();
#     γₙ(∇Gx,ts,trace) -γₙ(Gn,ts,trace) γₙ(Gr,ts,trace) -γₙ(∇G,ts,trace)]
#     return int
# end
function interaction(testobj::Object, trialobj::Object, embobj::Object, 
    testtype::ObjectType,trialtype::ObjectType,strat::VP; sign=1)
    trace = sign*strat.trace
    k = sqrt(epsilon(embobj)*mu(embobj))*omega(embobj)
    G = BEAST.greenhh3d(wavenumber=k)
    gradG = BEAST.∇(G)


    ∇Gx =  BEAST.build_potential(gradG×B)
    Gn = BEAST.build_potential(G*(n*B))
    #Gnx = BEAST.build_potential(G*(n × B),bs)
    ∇G = BEAST.build_potential(gradG*B)
    ∇Gdotn = BEAST.build_potential(nb ⋅ (gradG*B))
    ∇Gdot = BEAST.build_potential(B ⋅ gradG)
    
    Gr = BEAST.build_potential(G*B)
    ∇G∇B = BEAST.build_potential(gradG*div(B))
    ∇Gxn = BEAST.build_potential(gradG×(n*B))

    #gebruik gamma^c --> sign = +
    # #without extra ncross
    # int = [γₛ(∇Gx,ts,trace) -γₛ(Gn,ts,trace) γₛ(Gr,ts,trace) -γₛ(∇G,ts,trace);
    # ZeroOperator() -τ(∇Gdotn,ts,trace) τ(∇Gdot,ts,trace)  k^2*τ(Gr,ts,trace);
    # γₛ(∇G∇B,ts,trace)+k^2*γₛ(Gr,ts,trace) -γₛ(∇Gxn,ts,trace) γₛ(∇Gx,ts,trace) ZeroOperator();
    # γₙ(∇Gx,ts,trace) -γₙ(Gn,ts,trace) γₙ(Gr,ts,trace) -γₙ(∇G,ts,trace)]
    # #with extra ncross
    int = [-γₜ(∇Gx,trace) γₜ(Gn,trace) -γₜ(Gr,trace) γₜ(∇G,trace);
    ZeroOperator() -τ(∇Gdotn,trace) τ(∇Gdot,trace)  k^2*τ(Gr,trace);
    -γₜ(∇G∇B,trace)-k^2*γₜ(Gr,trace) γₜ(∇Gxn,trace) -γₜ(∇Gx,trace) ZeroOperator();
    γₙ(∇Gx,trace) -γₙ(Gn,trace) γₙ(Gr,trace) -γₙ(∇G,trace)]
    return int
end

function interaction(testobj::Object, trialobj::Object, embobj::Object, 
    testtype::ObjectType,trialtype::ObjectType,strat::PRECVP; sign=1)
    trace = sign*strat.trace
    k = sqrt(epsilon(embobj)*mu(embobj))*omega(embobj)
    G = BEAST.greenhh3d(wavenumber=k)
    gradG = BEAST.∇(G)

    bs = barycentric_refinement(geometry(trialobj); sort = :spacefillingcurve)
    ts = barycentric_refinement(geometry(testobj); sort = :spacefillingcurve)

    bss = geometry(trialobj)
    tss = geometry(testobj)
    

    ∇Gx =  BEAST.build_potential(gradG×B,bs)#
    Gn = BEAST.build_potential(G*(n*B),bss)#
    #Gnx = BEAST.build_potential(G*(n × B),bs)
    ∇G = BEAST.build_potential(gradG*B,bss)#
    ∇Gdotn = BEAST.build_potential(nb ⋅ (gradG*B),bss)#
    ∇Gdot = BEAST.build_potential(B ⋅ gradG,bs)#
    
    Gr = BEAST.build_potential(G*B,bs)#
    Grs = BEAST.build_potential(G*B,bss)#
    ∇G∇B = BEAST.build_potential(gradG*div(B),bs)#
    ∇Gxn = BEAST.build_potential(gradG×(n*B),bss)#

    int = [-γₜ(∇Gx,ts,trace) γₜ(Gn,ts,trace) -γₜ(Gr,ts,trace) γₜ(∇G,ts,trace);
    ZeroOperator() -τ(∇Gdotn,tss,trace) τ(∇Gdot,tss,trace)  k^2*τ(Grs,tss,trace);
    -γₜ(∇G∇B,ts,trace)-k^2*γₜ(Gr,ts,trace) γₜ(∇Gxn,ts,trace) -γₜ(∇Gx,ts,trace) ZeroOperator();
    γₙ(∇Gx,tss,trace) -γₙ(Gn,tss,trace) γₙ(Gr,tss,trace) -γₙ(∇G,tss,trace)]
    return int
end
function _identity(::ObjectType,strat::PRECVP)
    [NCross() ZeroOperator() ZeroOperator() ZeroOperator();
    ZeroOperator() Identity() ZeroOperator() ZeroOperator();
    ZeroOperator() ZeroOperator() NCross() ZeroOperator();
    ZeroOperator() ZeroOperator() ZeroOperator() Identity()]
end
_identity(p::PEC,strat::PRECVP) = _identity(p.inside,strat)[1:2,1:2]

function identity(objtype::ObjectType,strat::VP)
    [NCross() ZeroOperator() ZeroOperator() ZeroOperator();
    ZeroOperator() Identity() ZeroOperator() ZeroOperator();
    ZeroOperator() ZeroOperator() NCross() ZeroOperator();
    ZeroOperator() ZeroOperator() ZeroOperator() Identity()]
end
function identity(objtype::ObjectType,strat::PRECVP)
    -[NCross() ZeroOperator() ZeroOperator() ZeroOperator();
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
identity(o::HOM,strat::Union{PRECVP,VP}) = identity(o.inside,strat)
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
function excitation(testobj::Object,emb::Object,objtype::Inside,ex::VPExcitation,strat::VP)
    a = [1 0 0 0;
    0 1 0 0;
    0 0 mu(testobj)/mu(parent(testobj)) 0;
    0 0 0 epsilon(parent(testobj))/epsilon(testobj)]
    return a*excitation(testobj,emb,objtype.inside,ex,strat)
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
excitation(testobj::Object,emb::Object,objtype::HOM,ex::VPExcitation,strat::VP) = excitation(testobj,emb,objtype.inside,ex,strat)


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


