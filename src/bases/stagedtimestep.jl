
export StagedTimeStep

"""
    StagedTimeStep{T,N}

T: the value type of the basis function.
N: the number of stages.
It corresponds to a time-space basis function where each time step has
intermediary stages given by the vertor c in a Butcher tableau (A,b,c)
"""
struct StagedTimeStep{T,N}
	spatialBasis
	c  :: SVector{N,T}
	Δt :: T
	Nt :: Int
end

scalartype{T,N}(sts :: StagedTimeStep{T,N}) = T;
temporalbasis{T,N}(sts :: StagedTimeStep{T,N}) = timebasisdelta(sts.Δt, sts.Nt)
