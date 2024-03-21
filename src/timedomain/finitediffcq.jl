struct FiniteDiffConvolutionQuadrature{LK, T, N}
	laplaceKernel         :: LK # function of s that returns an IntegralOperator
	method                :: N
	Δt                    :: T
	zTransformedTermCount :: Int
	contourRadius         :: T
end

abstract  type FiniteDiffMethod end

struct BDF2{T} <:FiniteDiffMethod where T
	p::BEAST.Polynomial
	dt::T
end

struct BE{T} <: FiniteDiffMethod where T
	p::BEAST.Polynomial
	dt::T
end

BDF2(dt) = BDF2(BEAST.Polynomial(1.5/dt, -2.0/dt, 0.5/dt), dt) #2-step backward difference formula

BE(dt) = BE(BEAST.Polynomial(1.0/dt,-1.0/dt), dt) #Backward Euler or 1-step backward difference formula

BDF1(dt) = BE(dt) #1-step backward difference formula

scalartype(basis :: FiniteDiffTimeStep{SB, T}) where {SB, T} = T
temporalbasis(basis :: FiniteDiffTimeStep{SB, T}) where {SB, T} = timebasisdelta(basis.Δt, basis.Nt)


scalartype(etcq::FiniteDiffConvolutionQuadrature{LK, T, N}) where {LK, T, N} = Complex{T};

function assemble(cqop :: FiniteDiffConvolutionQuadrature,
                  testfns :: FiniteDiffTimeStep,
                  trialfns :: FiniteDiffTimeStep)

	laplaceKernel = cqop.laplaceKernel;
	method = cqop.method
	Δt = cqop.Δt;
	Q = cqop.zTransformedTermCount;
	rho = cqop.contourRadius;

	test_spatial_basis  = testfns.spatialBasis;
	trial_spatial_basis = trialfns.spatialBasis;

	# Compute the Z transformed sequence.
	# Assume that the operator applied on the conjugate of s is the same as the
	# conjugate of the operator applied on s,
	# so that only half of the values are computed
	Qmax = Q>>1+1;
	M = numfunctions(test_spatial_basis);
	N = numfunctions(trial_spatial_basis);
	Tz = promote_type(scalartype(cqop), scalartype(testfns), scalartype(trialfns));
	Zz = Vector{Array{Tz,2}}(undef,Qmax);

	for q = 0:Qmax-1
		# Build a temporary matrix for each eigenvalue
		s = laplace_to_z(rho, q, Q, Δt, method);
		Zz[q+1] = assemble(laplaceKernel(s), test_spatial_basis, trial_spatial_basis);
	end

	# return the inverse Z transform
	kmax = Q;
	T = real(Tz);
	Z = zeros(T, M, N, kmax)
	#return Zz
	print("Inverse z transform dots out of 10:")
	for q = 0:kmax-1
		Z[:,:,q+1] = real_inverse_z_transform(q, rho, Q, Zz);
		isinteger((q+1)*10/kmax) ? print(".") : nothing
	end
	return ConvolutionOperators.DenseConvOp(Z)

end