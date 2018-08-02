using Test
using LinearAlgebra

using BEAST, CompScienceMeshes, SauterSchwabQuadrature, StaticArrays

t1 = simplex(
    @SVector[0.180878, -0.941848, -0.283207],
    @SVector[0.0, -0.980785, -0.19509],
    @SVector[0.0, -0.92388, -0.382683])
t2 = simplex(
    @SVector[0.373086, -0.881524, -0.289348],
    @SVector[0.180878, -0.941848, -0.283207],
    @SVector[0.294908, -0.944921, -0.141962])


# s1 = simplex(
#     @SVector[0.180878, -0.941848, -0.283207],
#     @SVector[0.0, -0.980785, -0.19509],
#     @SVector[0.0, -0.92388, -0.382683])
# s2 = simplex(
#     @SVector[0.180878, -0.941848, -0.283207],
#     @SVector[0.373086, -0.881524, -0.289348],
#     @SVector[0.294908, -0.944921, -0.141962])

# Common face case:
rt = BEAST.RTRefSpace{Float64}()
tqd = BEAST.quadpoints(rt, [t1,t2], (12,))
bqd = BEAST.quadpoints(rt, [t1,t2], (13,))

op1 = BEAST.MWSingleLayer3D(1.0, 250.0, 1.0)
op2 = BEAST.MWSingleLayer3D(1.0)
op3 = BEAST.MWDoubleLayer3D(0.0)

SE_strategy = BEAST.WiltonSEStrategy(
  tqd[1,1],
  BEAST.DoubleQuadStrategy(
	tqd[1,1],
	bqd[1,1],
  ),
)
SS_strategy = SauterSchwabQuadrature.CommonFace(BEAST._legendre(8,0.0,1.0))

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

## Common Vertex Case
SE_strategy = BEAST.WiltonSEStrategy(
  tqd[1,1],
  BEAST.DoubleQuadStrategy(
	tqd[1,1],
	bqd[1,2]))
SS_strategy = SauterSchwabQuadrature.CommonVertex(BEAST._legendre(12,0.0,1.0))

z_cv_se_1 = zeros(3,3)
z_cv_ss_1 = zeros(3,3)
z_cv_se_2 = zeros(3,3)
z_cv_ss_2 = zeros(3,3)
z_cv_se_3 = zeros(3,3)
z_cv_ss_3 = zeros(3,3)

BEAST.momintegrals!(op1, rt, rt, t1, t2, z_cv_se_1, SE_strategy)
BEAST.momintegrals!(op1, rt, rt, t1, t2, z_cv_ss_1, SS_strategy)

@test z_cv_se_1 ≈ z_cv_ss_1 atol=1e-7

BEAST.momintegrals!(op2, rt, rt, t1, t2, z_cv_se_2, SE_strategy)
BEAST.momintegrals!(op2, rt, rt, t1, t2, z_cv_ss_2, SS_strategy)

@test z_cv_se_2 ≈ z_cv_ss_2 atol=1e-7

BEAST.momintegrals!(op3, rt, rt, t1, t2, z_cv_se_3, SE_strategy)
BEAST.momintegrals!(op3, rt, rt, t1, t2, z_cv_ss_3, SS_strategy)

@test z_cv_se_3 ≈ z_cv_ss_3 atol=1e-7

## Common Edge Case:
t1 = simplex(
    @SVector[0.180878, -0.941848, -0.283207],
    @SVector[0.0, -0.980785, -0.19509],
    @SVector[0.0, -0.92388, -0.382683],
    )
t2 = simplex(
    @SVector[0.180878, -0.941848, -0.283207],
    @SVector[0.158174, -0.881178, -0.44554],
    @SVector[0.0, -0.92388, -0.382683],
    )

tqd = BEAST.quadpoints(rt, [t1,t2], (12,))
bqd = BEAST.quadpoints(rt, [t1,t2], (13,))

SE_strategy = BEAST.WiltonSEStrategy(
    tqd[1,1],
    BEAST.DoubleQuadStrategy(
        tqd[1,1],
        bqd[1,2]))
SS_strategy = SauterSchwabQuadrature.CommonEdge(BEAST._legendre(12,0.0,1.0))

z_ce_se_1 = zeros(3,3)
z_ce_ss_1 = zeros(3,3)
z_ce_se_2 = zeros(3,3)
z_ce_ss_2 = zeros(3,3)
z_ce_se_3 = zeros(3,3)
z_ce_ss_3 = zeros(3,3)

BEAST.momintegrals!(op1, rt, rt, t1, t2, z_ce_se_1, SE_strategy)
BEAST.momintegrals!(op1, rt, rt, t1, t2, z_ce_ss_1, SS_strategy)

@show norm(z_ce_se_1 - z_ce_ss_1)
@test z_ce_se_1 ≈ z_ce_ss_1 atol=1e-5

BEAST.momintegrals!(op2, rt, rt, t1, t2, z_ce_se_2, SE_strategy)
BEAST.momintegrals!(op2, rt, rt, t1, t2, z_ce_ss_2, SS_strategy)

@show norm(z_ce_se_2 - z_ce_ss_2)
@test z_ce_se_2 ≈ z_ce_ss_2 atol=1e-5

BEAST.momintegrals!(op3, rt, rt, t1, t2, z_ce_se_3, SE_strategy)
BEAST.momintegrals!(op3, rt, rt, t1, t2, z_ce_ss_3, SS_strategy)
@show norm(z_ce_se_3 - z_ce_ss_3)
@test z_ce_se_3 ≈ z_ce_ss_3 atol=1e-5

SS_strategy = SauterSchwabQuadrature.CommonEdge(BEAST._legendre(18,0.0,1.0))
z_ce_ss_3_18 = zeros(3,3)
BEAST.momintegrals!(op3, rt, rt, t1, t2, z_ce_ss_3_18, SS_strategy)
@show norm(z_ce_ss_3 - z_ce_ss_3_18)
@test z_ce_ss_3 ≈ z_ce_ss_3_18 atol=1e-14
