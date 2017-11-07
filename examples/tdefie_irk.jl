include("../src/BEAST.jl")
using CompScienceMeshes, BEAST, StaticArrays

sol = 5.0;

struct StagedTimeStep{T,N}
	spatialBasis
	c  :: SVector{N,T}
	Δt :: T
	Nt :: Int
end
scalartype{T,N}(sts::StagedTimeStep{T,N}) = T;
BEAST.temporalbasis{T,N}(sts :: StagedTimeStep{T,N}) = BEAST.timebasisdelta(sts.Δt, sts.Nt)

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

function LaplaceToZ(rho, k, N, dt, A, b)
	z = rho * exp(2*im*pi*k/N);
	s = inv(dt * (A + ones(b) * b' / (z-1)));
	return s;
end

function InverseZTransform{T}(k, rho, N, X::AbstractArray{T,1})
	result = zero(X[1]);
	for n = 0:(N-1)
		result += X[n+1] * exp(2*im*pi*k*n/N)
	end
	result *= (rho^k)/N;
	return result;
end

	#########################################
	#                  RHS                  #
	#########################################

struct TimeBasisDeltaShifted{T} <: BEAST.AbstractTimeBasisFunction
	tbf   :: BEAST.TimeBasisDelta{T}
	shift :: T
end
BEAST.scalartype(x::TimeBasisDeltaShifted) = BEAST.scalartype(x.tbf)
BEAST.numfunctions(x::TimeBasisDeltaShifted) = BEAST.numfunctions(x.tbf)
BEAST.refspace(x::TimeBasisDeltaShifted) = BEAST.refspace(x.tbf)
BEAST.timestep(x::TimeBasisDeltaShifted) = BEAST.timestep(x.tbf)
BEAST.numintervals(x::TimeBasisDeltaShifted) = BEAST.numintervals(x.tbf)

function BEAST.assemblydata(tbds::TimeBasisDeltaShifted)
	tbf = tbds.tbf

    T = BEAST.scalartype(tbf)
    Δt = BEAST.timestep(tbf)

    z = zero(BEAST.scalartype(tbf))
    w = one(BEAST.scalartype(tbf))

    num_cells = BEAST.numfunctions(tbf)
    num_refs  = 1

    max_num_funcs = 1
    num_funcs = zeros(Int, num_cells, num_refs)
    data = fill((0,z), max_num_funcs, num_refs, num_cells)

    els = [ simplex(point((i-0+tbds.shift)*Δt),point((i+1+tbds.shift)*Δt)) for i in 1:num_cells ]

    for k in 1 : BEAST.numfunctions(tbf)
        data[1,1,k] = (k,w)
    end

    return els, BEAST.AssemblyData(data)
end

function BEAST.assemble(exc::BEAST.TDFunctional, testST::StagedTimeStep)
	stageCount = length(testST.c);
    spatialBasis = testST.spatialBasis;
    Nt = testST.Nt;
	Δt = testST.Δt;
    Z = zeros(eltype(exc), BEAST.numfunctions(spatialBasis) * stageCount, Nt)
	for i = 1:stageCount
		store(v,m,k) = (Z[(m-1)*stageCount+i,k] += v);
		tbsd = TimeBasisDeltaShifted(BEAST.timebasisdelta(Δt, Nt), testST.c[i]);
		BEAST.assemble!(exc, spatialBasis ⊗ tbsd, store);
	end
    return Z
end

	#########################################
	#                  LHS                  #
	#########################################

function BEAST.assemble(rkcq::RungeKuttaConvolutionQuadrature, testfns::StagedTimeStep, trialfns::StagedTimeStep)

	laplaceKernel = rkcq.laplaceKernel;
	A = rkcq.A;
	b = rkcq.b;
	Δt = rkcq.Δt;
	Q = rkcq.zTransformedTermCount;
	rho = rkcq.contourRadius;
	p = length(b);

    test_spatial_basis  = testfns.spatialBasis
    trial_spatial_basis = trialfns.spatialBasis
	Tz = promote_type(scalartype(rkcq), scalartype(testfns), scalartype(trialfns))

	# Compute the Z transformed sequence
    M = numfunctions(test_spatial_basis)
    N = numfunctions(trial_spatial_basis)
	Zz = Vector{Array{Tz,2}}(Q);
	for q = 0:Q-1
		# Build a temporary matrix for each eigenvalue
		s = LaplaceToZ(rho, q, Q, Δt, A, b);
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
		Z[:,:,q+1] = real(InverseZTransform(q, rho, Q, Zz));
	end
    return Z

end

function BEAST.solve(eq)

    time_domain = isa(first(eq.trial_space_dict).second, BEAST.SpaceTimeBasis)
    time_domain |= isa(first(eq.trial_space_dict).second, StagedTimeStep)
    if time_domain
        return BEAST.td_solve(eq)
    end

    test_space_dict  = eq.test_space_dict
    trial_space_dict = eq.trial_space_dict

    lhs = eq.equation.lhs
    rhs = eq.equation.rhs

    b = BEAST.assemble(rhs, test_space_dict)
    Z = BEAST.assemble(lhs, test_space_dict, trial_space_dict)

    u = Z \ b

    return u
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
