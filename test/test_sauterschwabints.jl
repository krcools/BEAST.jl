using Test

using BEAST, CompScienceMeshes, SauterSchwabQuadrature, StaticArrays

t1 = simplex(
    @SVector[0.180878, -0.941848, -0.283207],
    @SVector[0.0, -0.980785, -0.19509],
    @SVector[0.0, -0.92388, -0.382683])


rt = BEAST.RTRefSpace{Float64}()
tqd = BEAST.quadpoints(rt, [t1], (12,))
bqd = BEAST.quadpoints(rt, [t1], (13,))

op1 = BEAST.MWSingleLayer3D(1.0, 250.0, 1.0)
op2 = BEAST.MWSingleLayer3D(1.0)

SE_strategy = BEAST.WiltonSEStrategy(
  tqd[1,1],
  BEAST.DoubleQuadStrategy(
	tqd[1,1],
	bqd[1,1],
  ),
)
SS_strategy = SauterSchwabQuadrature.CommonFace(8)

z_se = zeros(3,3)
z_ss = zeros(3,3)

BEAST.momintegrals!(op1, rt, rt, t1, t1, z_se, SE_strategy)
BEAST.momintegrals!(op1, rt, rt, t1, t1, z_ss, SS_strategy)

@test z_se ≈ z_ss atol=1e-4

z_se2 = zeros(3,3)
z_ss2 = zeros(3,3)

BEAST.momintegrals!(op2, rt, rt, t1, t1, z_se2, SE_strategy)
BEAST.momintegrals!(op2, rt, rt, t1, t1, z_ss2, SS_strategy)

@test z_se2 ≈ z_ss2 atol=1e-4
