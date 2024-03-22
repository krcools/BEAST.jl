


struct SpaceTimeBasis{S,T,U} <: Space{U}
    space::S
    time::T
end

function SpaceTimeBasis(space, time)
    S = typeof(space)
    T = typeof(time)

    U = scalartype(space,time)
    return SpaceTimeBasis{S,T,U}(space,time)
end

spatialbasis(s) = s
spatialbasis(s::SpaceTimeBasis) = s.space
function spatialbasis(s::DirectProductSpace)
    r = [spatialbasis(ch) for ch in s.factors]
    return DirectProductSpace(r)
end

temporalbasis(s::SpaceTimeBasis) = s.time
function numstages(s::DirectProductSpace)
    r = [numstages(temporalbasis(ch)) for ch in s.factors]
    return r
end

⊗(a, b) = SpaceTimeBasis(a,b)
numfunctions(S::SpaceTimeBasis) = numfunctions(S.space) * numfunctions(S.time)
scalartype(st::SpaceTimeBasis) = promote_type(scalartype(st.space), scalartype(st.time))

function tensordim(s::SpaceTimeBasis,i)
    i == 1 && return numfunctions(s.space)
    i == 2 && return numfunctions(s.time)
    return 1
end

function tensordim(s, i)
    i == 1 && return numfunctions(s)
    return 1
end

