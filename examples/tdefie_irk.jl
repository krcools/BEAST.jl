include("../src/BEAST.jl")
using CompScienceMeshes, BEAST, StaticArrays

sol = 5.0;

struct RungeKuttaConvolutionQuadrature{T,N,NN}
	laplaceKernel         # function of s that returns an IntegralOperator
	A                     :: SArray{Tuple{N,N},T,2,NN}
	b                     :: SVector{N,T}
	Δt                    :: T
	zTransformedTermCount :: Int
	contourRadius         :: T
end
scalartype{T,N,NN}(rkcq::RungeKuttaConvolutionQuadrature{T,N,NN}) = Complex{T};

LaplaceEFIO{T}(s::T) = MWSingleLayer3D(-s/sol, s*s/sol, T(sol));

# M = H*diagm(D)*invH
struct DiagonalizedMatrix{T,N,NN}
	H    :: SArray{Tuple{N,N},Complex{T},2,NN}
	invH :: SArray{Tuple{N,N},Complex{T},2,NN}
	D    :: SVector{N,Complex{T}}
end

# M = H*diagm(D)*invH
function DiagonalizedMatrix{T,N,NN}(M :: SArray{Tuple{N,N},Complex{T},2,NN})
	ef = eigfact(Array{Complex{T},2}(M));

	efValues = SVector{N,Complex{T}}(ef.values);
	efVectors = SArray{Tuple{N,N},Complex{T},2,NN}(ef.vectors);
	return DiagonalizedMatrix(efVectors, inv(efVectors), efValues);
end

	#########################################
	#                  LHS                  #
	#########################################

function BEAST.assemble(rkcq::RungeKuttaConvolutionQuadrature, testfns::BEAST.StagedTimeStep, trialfns::BEAST.StagedTimeStep)

	laplaceKernel = rkcq.laplaceKernel;
	A = rkcq.A;
	b = rkcq.b;
	Δt = rkcq.Δt;
	Q = rkcq.zTransformedTermCount;
	rho = rkcq.contourRadius;
	p = length(b);

    test_spatial_basis  = testfns.spatialBasis
    trial_spatial_basis = trialfns.spatialBasis
	Tz = promote_type(scalartype(rkcq), BEAST.scalartype(testfns), BEAST.scalartype(trialfns))

	# Compute the Z transformed sequence
    M = numfunctions(test_spatial_basis)
    N = numfunctions(trial_spatial_basis)
	Zz = Vector{Array{Tz,2}}(Q);
	for q = 0:Q-1
		# Build a temporary matrix for each eigenvalue
		s = BEAST.laplace_to_z(rho, q, Q, Δt, A, b);
		sFactorized = DiagonalizedMatrix(s);
		blocksEigenvalues = [BEAST.assemble(laplaceKernel(sD), test_spatial_basis, trial_spatial_basis) for sD in sFactorized.D];
		Zz[q+1] = zeros(Tz, M*p, N*p)

		# Compute the Z transformed matrix by block
		for m = 1:M
			for n = 1:N
				D = SVector{p,Tz}([blocksEigenvalues[i][m,n] for i in 1:p])
				Zz[q+1][(m-1)*p+(1:p),(n-1)*p+(1:p)] = sFactorized.H * diagm(D) * sFactorized.invH;
			end
		end
	end

	# return the inverse Z tranformed
	kmax = Q;
	T = real(Tz);
	Z = zeros(T, M*p, N*p, kmax)
	for q = 0:kmax-1
		Z[:,:,q+1] = real(BEAST.inverse_z_transform(q, rho, Q, Zz));
	end
    return Z

end

#####################################################################

o, x, y, z = euclidianbasis(3)

#Δt, Nt = 0.08/sol,400
Δt, Nt = 100.0/sol,200

D, Δx = 1.0, 0.45
Γ = meshsphere(D, Δx)
X = raviartthomas(Γ)

(A, b, c) = butcher_tableau_radau_2stages();
V = StagedTimeStep(X, c, Δt, Nt);

duration, delay, amplitude = 2000.0/sol, 2000.0/sol, 1.0
gaussian = creategaussian(duration, delay, duration)

direction, polarisation = z, x
E = planewave(polarisation, direction, derive(gaussian), sol)
T = RungeKuttaConvolutionQuadrature(LaplaceEFIO, A, b, Δt, 15, 1.0001);

@hilbertspace j
@hilbertspace j′
tdefie = @discretise T[j′,j] == -1E[j′]   j∈V  j′∈V
xefie_irk = solve(tdefie)

using PyPlot
semilogy(abs.(xefie_irk[2,:]))
