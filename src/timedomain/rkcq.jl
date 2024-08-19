

"""
    RungeKuttaConvolutionQuadrature{T,N,NN}

T: the value type of the basis function.
N: the number of stages.
NN: N*N.

Performs a convolution quadrature on a laplaceKernel to represent an operator
in time domain using an implicit Runge-Kutta method.

laplaceKernel: function of the Laplace variable s that returns an IntegralOperator.
A, b: Coefficient matrix and vectors from the Butcher tableau.
Δt: time step.
zTransformedTermCount: Number of terms in the inverse Z-transform.
contourRadius: radius of circle used as integration contour for the inverse Z-transform.
"""
struct RungeKuttaConvolutionQuadrature{}
	timedomainKernel :: AbstractSpaceTimeOperator # function of s that returns an IntegralOperator
end
#scalartype(rkcq::RungeKuttaConvolutionQuadrature{LK}) where {LK} = Complex


# M = H*diagm(D)*invH
struct DiagonalizedMatrix{T,N,NN}
	H    :: SArray{Tuple{N,N},Complex{T},2,NN}
	invH :: SArray{Tuple{N,N},Complex{T},2,NN}
	D    :: SVector{N,Complex{T}}
end

# M = H*diagm(D)*invH
function diagonalizedmatrix(M :: SArray{Tuple{N,N},Complex{T},2,NN}) where {T,N,NN}
	ef = eigen(Array{Complex{T},2}(M));

	efValues = SVector{N,Complex{T}}(ef.values)  :: SVector{N,Complex{T}};
	efVectors = SArray{Tuple{N,N},Complex{T},2,NN}(ef.vectors) :: SArray{Tuple{N,N},Complex{T},2,NN};
	return DiagonalizedMatrix(efVectors, inv(efVectors), efValues);
end

function assemble(rkcq :: RungeKuttaConvolutionQuadrature,
                  testfns :: SpaceTimeBasis,
                  trialfns :: SpaceTimeBasis)

	@warn "staged assemble of the left-hand side"
	sol = rkcq.timedomainKernel.speed_of_light
	numdiffweak = rkcq.timedomainKernel.ws_diffs
	numdiffhyper = rkcq.timedomainKernel.hs_diffs

	@info "converting time-domain kernel to Laplace-domain kernel"
	LaplaceEFIO(s::T) where {T}= MWSingleLayer3D(s/sol, -s^(numdiffweak)/sol, s^(numdiffhyper)*T(sol))
	laplaceKernel = LaplaceEFIO
	
	A = testfns.time.A
	b = testfns.time.b
	Δt = testfns.time.Δt
	Q = testfns.time.zTransformedTermCount
	rho = testfns.time.contourRadius
	p = length(b) # stage count

	test_spatial_basis  = testfns.space
	trial_spatial_basis = trialfns.space

	# Compute the Z transformed sequence.
	# Assume that the operator applied on the conjugate of s is the same as the
	# conjugate of the operator applied on s,
	# so that only half of the values are computed
	Qmax = Q>>1+1
	M = numfunctions(test_spatial_basis)
	N = numfunctions(trial_spatial_basis)
	#Tz = promote_type(scalartype(rkcq), scalartype(testfns), scalartype(trialfns))
	#Tz = promote_type(scalartype(testfns), scalartype(trialfns))
	Tz = ComplexF64
	Zz = Vector{Array{Tz,2}}(undef,Qmax)
	blocksEigenvalues = Vector{Array{Tz,2}}(undef,p)
	tmpDiag = Vector{Tz}(undef,p)
	for q = 0:Qmax-1
		# Build a temporary matrix for each eigenvalue
		s = laplace_to_z(rho, q, Q, Δt, A, b)
		sFactorized = diagonalizedmatrix(s)
		for (i,sD) in enumerate(sFactorized.D)
			blocksEigenvalues[i] = assemble(laplaceKernel(sD), test_spatial_basis, trial_spatial_basis)
		end

		# Compute the Z transformed matrix by block
		Zz[q+1] = zeros(Tz, M*p, N*p)
		for m = 1:M
			for n = 1:N
				for i = 1:p
					tmpDiag[i] = blocksEigenvalues[i][m,n]
				end
				D = SVector{p,Tz}(tmpDiag);
				Zz[q+1][(m-1)*p.+(1:p),(n-1)*p.+(1:p)] = sFactorized.H * diagm(D) * sFactorized.invH
			end
		end
	end

	# return the inverse Z transform
	kmax = Q
	T = real(Tz)
	Z = zeros(T, M*p, N*p, kmax)
	for q = 0:kmax-1
		Z[:,:,q+1] = real_inverse_z_transform(q, rho, Q, Zz)
	end
	ZC = ConvolutionOperators.DenseConvOp(Z)
	return ZC

end


struct LaplaceDomainOperator
	coeff::ComplexF64
	kernel::Any
end

scalartype(op::LaplaceDomainOperator) = ComplexF64

*(a::Number, op::LaplaceDomainOperator) = LaplaceDomainOperator(a*op.coeff, op.kernel)

LaplaceDomainOperator(kernel) = LaplaceDomainOperator(Complex(1.0), kernel)

defaultquadstrat(op::LaplaceDomainOperator, tfs::SpaceTimeBasis, bfs::SpaceTimeBasis) = DoubleNumWiltonSauterQStrat(2,3,6,7,5,5,4,3)

function assemble(op::LaplaceDomainOperator, testfns::SpaceTimeBasis, trialfns::SpaceTimeBasis; quadstrat=defaultquadstrat(op, testfns, trialfns))
	laplaceKernel = op.kernel

	A = testfns.time.A
	b = testfns.time.b
	Δt = testfns.time.Δt
	Q = testfns.time.zTransformedTermCount
	rho = testfns.time.contourRadius
	p = length(b) # stage count

	test_spatial_basis  = testfns.space
	trial_spatial_basis = trialfns.space

	# Compute the Z transformed sequence.
	# Assume that the operator applied on the conjugate of s is the same as the
	# conjugate of the operator applied on s,
	# so that only half of the values are computed
	Qmax = Q>>1+1
	M = numfunctions(test_spatial_basis)
	N = numfunctions(trial_spatial_basis)
	Tz = ComplexF64
	Zz = Vector{Array{Tz,2}}(undef,Qmax)
	blocksEigenvalues = Vector{Array{Tz,2}}(undef,p)
	tmpDiag = Vector{Tz}(undef,p)

	for q = 0:Qmax-1
		# Build a temporary matrix for each eigenvalue
		s = laplace_to_z(rho, q, Q, Δt, A, b)
		sFactorized = diagonalizedmatrix(s)
		for (i,sD) in enumerate(sFactorized.D)
			blocksEigenvalues[i] = op.coeff * assemble(laplaceKernel(sD), test_spatial_basis, trial_spatial_basis, quadstrat=quadstrat)
		end

		# Compute the Z transformed matrix by block
		Zz[q+1] = zeros(Tz, M*p, N*p)
		for m = 1:M
			for n = 1:N
				for i = 1:p
					tmpDiag[i] = blocksEigenvalues[i][m,n]
				end
				D = SVector{p,Tz}(tmpDiag);
				Zz[q+1][(m-1)*p.+(1:p),(n-1)*p.+(1:p)] = sFactorized.H * diagm(D) * sFactorized.invH
			end
		end
	end

	# return the inverse Z transform
	kmax = Q
	T = real(Tz)
	Z = zeros(T, M*p, N*p, kmax)
	for q = 0:kmax-1
		Z[:,:,q+1] = real_inverse_z_transform(q, rho, Q, Zz)
	end
	ZC = ConvolutionOperators.DenseConvOp(Z)

	return ZC
end