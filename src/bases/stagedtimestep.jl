


"""
    StagedTimeStep{T,N,NN}

T: the value type of the basis function.
N: the number of stages.
NN: the number of stages squared
Each time step has intermediary stages given by the vertor c in a Butcher tableau (A,b,c)
"""
struct StagedTimeStep{T, N, NN, I}
	Δt :: T
	Nt :: I
	c  :: SVector{N,T}
	A  :: SArray{Tuple{N,N},T,2,NN}
	b  :: SVector{N,T}
	zTransformedTermCount :: I
	contourRadius         :: T
end

scalartype(sts :: StagedTimeStep{T, N, NN, I}) where {T, N, NN, I} = T
temporalbasis(sts :: StagedTimeStep{T, N, NN, I}) where {T, N, NN, I} = timebasisdelta(sts.Δt, sts.Nt)

numfunctions(s::StagedTimeStep) = s.Nt

numstages(s) = 1
numstages(s::StagedTimeStep) = size(s.c,1)
