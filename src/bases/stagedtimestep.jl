


"""
    StagedTimeStep{T,N}

T: the value type of the basis function.
N: the number of stages.
It corresponds to a time-space basis function where each time step has
intermediary stages given by the vertor c in a Butcher tableau (A,b,c)
"""
struct StagedTimeStep{SB, T, N}
	spatialBasis :: SB
	c  :: SVector{N,T}
	Δt :: T
	Nt :: Int
end

scalartype(sts :: StagedTimeStep{SB, T, N}) where {SB, T, N} = T;
temporalbasis(sts :: StagedTimeStep{SB, T, N}) where {SB, T, N} = timebasisdelta(sts.Δt, sts.Nt)
