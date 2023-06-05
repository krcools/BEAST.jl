


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

spatialbasis(s::SpaceTimeBasis) = s.space
temporalbasis(s::SpaceTimeBasis) = s.time

⊗(a, b) = SpaceTimeBasis(a,b)
numfunctions(S::SpaceTimeBasis) = numfunctions(S.space) * numfunctions(S.time)
scalartype(st::SpaceTimeBasis) = promote_type(scalartype(st.space), scalartype(st.time))
