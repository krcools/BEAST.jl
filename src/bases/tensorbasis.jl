


mutable struct SpaceTimeBasis{S,T}
    space::S
    time::T
end

spatialbasis(s::SpaceTimeBasis) = s.space
temporalbasis(s::SpaceTimeBasis) = s.time

⊗(a, b) = SpaceTimeBasis(a,b)
numfunctions(S::SpaceTimeBasis) = numfunctions(S.space) * numfunctions(S.time)
scalartype(st::SpaceTimeBasis) = promote_type(scalartype(st.space), scalartype(st.time))
