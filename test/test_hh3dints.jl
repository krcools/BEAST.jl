using Test
using LinearAlgebra

using BEAST, CompScienceMeshes, SauterSchwabQuadrature, StaticArrays

T=Float64


# operators
op1 = Helmholtz3D.singlelayer(gamma=T(1.0)im+T(1.0))
op2 = Helmholtz3D.doublelayer(gamma=T(1.0)im+T(1.0))
op3 = Helmholtz3D.doublelayer_transposed(gamma=T(1.0)im+T(1.0))
op4 = Helmholtz3D.hypersingular(gamma=T(1.0)im+T(1.0))
## Common face case:
@testset "Common Face" begin
  # triangles
t1 = simplex(
  T.(@SVector [0.180878, -0.941848, -0.283207]),
  T.(@SVector [0.0, -0.980785, -0.19509]),
  T.(@SVector [0.0, -0.92388, -0.382683]))

  lag0 = BEAST.LagrangeRefSpace{T,0,3,1}() # patch basis
  lag1 = BEAST.LagrangeRefSpace{T,1,3,3}() #pyramid basis

  tqd0 = BEAST.quadpoints(lag0, [t1], (13,)) #test quadrature data
  tqd1 = BEAST.quadpoints(lag1, [t1], (13,))

  bqd0 = BEAST.quadpoints(lag0, [t1], (12,)) #basis quadrature data
  bqd1 = BEAST.quadpoints(lag1, [t1], (12,))

  SE_strategy = BEAST.WiltonSERule( #wilton
    tqd0[1,1],
    BEAST.DoubleQuadRule(
    tqd0[1,1],
    bqd0[1,1],
    ),
  )

  SS_strategy = SauterSchwabQuadrature.CommonFace(BEAST._legendre(12,T(0.0),T(1.0))) #sauter
#single layer with patch-patch
  Double_strategy = BEAST.DoubleQuadRule( #doublequadstrat
    tqd0[1,1],
    bqd0[1,1]
  )

  z_se = [0.0im]
  z_ss = [0.0im]
  z_double = [0.0im]

  BEAST.momintegrals!(op1, lag0, lag0, t1, t1, z_se, SE_strategy)
  BEAST.momintegrals!(op1, lag0, lag0, t1, t1, z_ss, SS_strategy)
  BEAST.momintegrals!(op1, lag0, lag0, t1, t1, z_double, Double_strategy)

  @test z_se ≈ z_ss rtol=1e-5
  @test z_double ≈ z_ss rtol = 1e-1
  @test norm(z_se-z_ss) < norm(z_double-z_ss)

  #single layer with patch-pyramid
  SE_strategy = BEAST.WiltonSERule( #wilton
  tqd1[1,1],
  BEAST.DoubleQuadRule(
  tqd1[1,1],
  bqd0[1,1],
  ),)
  Double_strategy = BEAST.DoubleQuadRule(
    tqd1[1,1],
    bqd0[1,1],
    )

  z_se = zeros(complex(T),3,1)
  z_ss = zeros(complex(T),3,1)
  z_double = zeros(complex(T),3,1)
  BEAST.momintegrals!(op1, lag1, lag0, t1, t1, z_ss, SS_strategy)
  BEAST.momintegrals!(op1, lag1, lag0, t1, t1, z_se, SE_strategy)
  BEAST.momintegrals!(op1, lag1, lag0, t1, t1, z_double, Double_strategy)

  @test z_ss ≈ z_se rtol = 1e-5
  @test z_ss ≈ z_double rtol = 1e-1
  @test norm(z_ss-z_se) < norm(z_ss-z_double)

#single layer with pyramid-patch
  SE_strategy = BEAST.WiltonSERule( #wilton
  tqd0[1,1],
  BEAST.DoubleQuadRule(
  tqd0[1,1],
  bqd1[1,1],
  ),)
  Double_strategy = BEAST.DoubleQuadRule(
    tqd0[1,1],
    bqd1[1,1],
    )

  z_se = zeros(complex(T),1,3)
  z_ss = zeros(complex(T),1,3)
  z_double = zeros(complex(T),1,3)
  BEAST.momintegrals!(op1, lag0, lag1, t1, t1, z_ss, SS_strategy)
  BEAST.momintegrals!(op1, lag0, lag1, t1, t1, z_se, SE_strategy)
  BEAST.momintegrals!(op1, lag0, lag1, t1, t1, z_double, Double_strategy)

  @test z_ss ≈ z_se rtol = 1e-5
  @test z_ss ≈ z_double rtol = 1e-1
  @test norm(z_ss-z_se) < norm(z_ss-z_double)

  #single layer with pyramid pyramid
  SE_strategy = BEAST.WiltonSERule( #wilton
  tqd1[1,1],
  BEAST.DoubleQuadRule(
  tqd1[1,1],
  bqd1[1,1],
  ),)
  Double_strategy = BEAST.DoubleQuadRule(
    tqd1[1,1],
    bqd1[1,1],
    )

  z_se = zeros(complex(T),3,3)
  z_ss = zeros(complex(T),3,3)
  z_double = zeros(complex(T),3,3)
  BEAST.momintegrals!(op1, lag1, lag1, t1, t1, z_ss, SS_strategy)
  BEAST.momintegrals!(op1, lag1, lag1, t1, t1, z_se, SE_strategy)
  BEAST.momintegrals!(op1, lag1, lag1, t1, t1, z_double, Double_strategy)

  @test z_ss ≈ z_se rtol = 1e-5
  @test z_ss ≈ z_double rtol = 1e-1
  @test norm(z_ss-z_se) < norm(z_ss-z_double)

end
# common vertex case

# In the common vertex case, for the single layer and hypersingular,
# the double quadrature apparently is better than the wilton method.
# It is not clear to me why this is. If someone knows why this happens
# I would be grateful for an explanation.

@testset "Common vertex" begin
  t1 = simplex(
    T.(@SVector [0.180878, -0.941848, -0.283207]),
    T.(@SVector [0.0, -0.980785, -0.19509]),
    T.(@SVector [0.0, -0.92388, -0.382683]))
  t2 = simplex(
    T.(@SVector [0.373086, -0.881524, -0.289348]),
    T.(@SVector [0.180878, -0.941848, -0.283207]),
    T.(@SVector [0.294908, -0.944921, -0.141962]))
  lag0 = BEAST.LagrangeRefSpace{T,0,3,1}() # patch basis
  lag1 = BEAST.LagrangeRefSpace{T,1,3,3}() #pyramid basis

  tqd0 = BEAST.quadpoints(lag0, [t1,t2], (12,)) #test quadrature data
  tqd1 = BEAST.quadpoints(lag1, [t1,t2], (12,))

  bqd0 = BEAST.quadpoints(lag0, [t1,t2], (13,)) #basis quadrature data
  bqd1 = BEAST.quadpoints(lag1, [t1,t2], (13,))
#patch patch
  SE_strategy = BEAST.WiltonSERule(
    tqd0[1,1],
    BEAST.DoubleQuadRule(
    tqd0[1,1],
    bqd0[1,2]))

  SS_strategy = SauterSchwabQuadrature.CommonVertex(BEAST._legendre(8,T(0.0),T(1.0)))

  Double_strategy = BEAST.DoubleQuadRule( #doublequadstrat
    tqd0[1,1],
    bqd0[1,2]
  )
  #single layer

  z_cv_se = zeros(complex(T),1,1)
  z_cv_ss = zeros(complex(T),1,1)
  z_cv_double = zeros(complex(T),1,1)

  BEAST.momintegrals!(op1, lag0, lag0, t1, t2, z_cv_se, SE_strategy)
  BEAST.momintegrals!(op1, lag0, lag0, t1, t2, z_cv_ss, SS_strategy)
  BEAST.momintegrals!(op1, lag0, lag0, t1, t2, z_cv_double, Double_strategy)

  @test z_cv_double ≈ z_cv_ss rtol = 1e-4
  @test z_cv_se ≈ z_cv_ss rtol=1e-4
  @test_broken norm(z_cv_ss-z_cv_se)/norm(z_cv_ss) < norm(z_cv_ss-z_cv_double)/norm(z_cv_ss)
  #doublelayer
  z_cv_se = zeros(complex(T),1,1)
  z_cv_ss = zeros(complex(T),1,1)
  z_cv_double = zeros(complex(T),1,1)

  BEAST.momintegrals!(op2, lag0, lag0, t1, t2, z_cv_se, SE_strategy)
  BEAST.momintegrals!(op2, lag0, lag0, t1, t2, z_cv_ss, SS_strategy)
  BEAST.momintegrals!(op2, lag0, lag0, t1, t2, z_cv_double, Double_strategy)

  @test z_cv_double ≈ z_cv_ss rtol = 1e-4
  @test z_cv_se ≈ z_cv_double rtol=1e-4
  @test norm(z_cv_ss-z_cv_se) < norm(z_cv_ss-z_cv_double)

  # doublelayer_transposed
  z_cv_se = zeros(complex(T),1,1)
  z_cv_ss = zeros(complex(T),1,1)
  z_cv_double = zeros(complex(T),1,1)

  BEAST.momintegrals!(op3, lag0, lag0, t1, t2, z_cv_se, SE_strategy)
  BEAST.momintegrals!(op3, lag0, lag0, t1, t2, z_cv_ss, SS_strategy)
  BEAST.momintegrals!(op3, lag0, lag0, t1, t2, z_cv_double, Double_strategy)

  @test z_cv_double ≈ z_cv_ss rtol = 1e-4
  @test z_cv_se ≈ z_cv_ss rtol=1e-4
  @test norm(z_cv_ss-z_cv_se) < norm(z_cv_ss-z_cv_double)

#pyramid pyramid
  # singlelayer
  SE_strategy = BEAST.WiltonSERule(
    tqd1[1,1],
    BEAST.DoubleQuadRule(
    tqd1[1,1],
    bqd1[1,2]))
  Double_strategy = BEAST.DoubleQuadRule( #doublequadstrat
    tqd1[1,1],
    bqd1[1,2]
  )
  z_cv_se = zeros(complex(T),3,3)
  z_cv_ss = zeros(complex(T),3,3)
  z_cv_double = zeros(complex(T),3,3)

  BEAST.momintegrals!(op1, lag1, lag1, t1, t2, z_cv_se, SE_strategy)
  BEAST.momintegrals!(op1, lag1, lag1, t1, t2, z_cv_ss, SS_strategy)
  BEAST.momintegrals!(op1, lag1, lag1, t1, t2, z_cv_double, Double_strategy)

  @test z_cv_double ≈ z_cv_ss rtol = 1e-4
  @test z_cv_se ≈ z_cv_double rtol=1e-4
  @test_broken norm(z_cv_ss-z_cv_se) < norm(z_cv_ss-z_cv_double)

  #doublelayer 
  z_cv_se = zeros(complex(T),3,3)
  z_cv_ss = zeros(complex(T),3,3)
  z_cv_double = zeros(complex(T),3,3)

  BEAST.momintegrals!(op2, lag1, lag1, t1, t2, z_cv_se, SE_strategy)
  BEAST.momintegrals!(op2, lag1, lag1, t1, t2, z_cv_ss, SS_strategy)
  BEAST.momintegrals!(op2, lag1, lag1, t1, t2, z_cv_double, Double_strategy)

  @test z_cv_double ≈ z_cv_ss rtol = 1e-4
  @test z_cv_se ≈ z_cv_ss rtol=1e-4
  @test norm(z_cv_ss-z_cv_se) < norm(z_cv_ss-z_cv_double)

  #doublelayer_transposed
  z_cv_se = zeros(complex(T),3,3)
  z_cv_ss = zeros(complex(T),3,3)
  z_cv_double = zeros(complex(T),3,3)

  BEAST.momintegrals!(op3, lag1, lag1, t1, t2, z_cv_se, SE_strategy)
  BEAST.momintegrals!(op3, lag1, lag1, t1, t2, z_cv_ss, SS_strategy)
  BEAST.momintegrals!(op3, lag1, lag1, t1, t2, z_cv_double, Double_strategy)

  @test z_cv_double ≈ z_cv_ss rtol = 1e-3
  @test z_cv_double ≈ z_cv_se rtol = 1e-5
  @test norm(z_cv_ss-z_cv_se) < norm(z_cv_ss-z_cv_double)

  #hypersingular
  z_cv_se = zeros(complex(T),3,3)
  z_cv_ss = zeros(complex(T),3,3)
  z_cv_double = zeros(complex(T),3,3)

  BEAST.momintegrals!(op4, lag1, lag1, t1, t2, z_cv_se, SE_strategy)
  BEAST.momintegrals!(op4, lag1, lag1, t1, t2, z_cv_ss, SS_strategy)
  BEAST.momintegrals!(op4, lag1, lag1, t1, t2, z_cv_double, Double_strategy)

  # @show z_cv_double
  # @show z_cv_se
  # @show z_cv_ss
  # @show 2 * norm(z_cv_double + z_cv_ss) / norm(z_cv_double - z_cv_ss)
  @test z_cv_double ≈ z_cv_ss rtol = 1e-4
  @test z_cv_se ≈ z_cv_ss rtol=1e-4
  @test_broken norm(z_cv_ss-z_cv_se) < norm(z_cv_ss-z_cv_double)

#pyramid test patch basis
  SE_strategy = BEAST.WiltonSERule(
    tqd1[1,1],
    BEAST.DoubleQuadRule(
    tqd1[1,1],
    bqd0[1,2]))

  SS_strategy = SauterSchwabQuadrature.CommonVertex(BEAST._legendre(8,T(0.0),T(1.0)))

  Double_strategy = BEAST.DoubleQuadRule( #doublequadstrat
    tqd1[1,1],
    bqd0[1,2]
  )
  #singlelayer
  z_cv_se = zeros(complex(T),3,1)
  z_cv_ss = zeros(complex(T),3,1)
  z_cv_double = zeros(complex(T),3,1)

  BEAST.momintegrals!(op1, lag1, lag0, t1, t2, z_cv_se, SE_strategy)
  BEAST.momintegrals!(op1, lag1, lag0, t1, t2, z_cv_ss, SS_strategy)
  BEAST.momintegrals!(op1, lag1, lag0, t1, t2, z_cv_double, Double_strategy)

  @test z_cv_double ≈ z_cv_ss rtol = 1e-4
  @test z_cv_se ≈ z_cv_ss rtol=1e-4
  @test_broken norm(z_cv_ss-z_cv_se) < norm(z_cv_ss-z_cv_double)

  #doublelayer
  z_cv_se = zeros(complex(T),3,1)
  z_cv_ss = zeros(complex(T),3,1)
  z_cv_double = zeros(complex(T),3,1)

  BEAST.momintegrals!(op2, lag1, lag0, t1, t2, z_cv_se, SE_strategy)
  BEAST.momintegrals!(op2, lag1, lag0, t1, t2, z_cv_ss, SS_strategy)
  BEAST.momintegrals!(op2, lag1, lag0, t1, t2, z_cv_double, Double_strategy)

  @test z_cv_double ≈ z_cv_ss rtol = 1e-4
  @test z_cv_se ≈ z_cv_ss rtol=1e-4
  @test norm(z_cv_ss-z_cv_se)/norm(z_cv_ss) < norm(z_cv_ss-z_cv_double)/norm(z_cv_ss)
  #doublelayer_transposed

  z_cv_se = zeros(complex(T),3,1)
  z_cv_ss = zeros(complex(T),3,1)
  z_cv_double = zeros(complex(T),3,1)

  BEAST.momintegrals!(op3, lag1, lag0, t1, t2, z_cv_se, SE_strategy)
  BEAST.momintegrals!(op3, lag1, lag0, t1, t2, z_cv_ss, SS_strategy)
  BEAST.momintegrals!(op3, lag1, lag0, t1, t2, z_cv_double, Double_strategy)

  @test z_cv_double ≈ z_cv_ss rtol = 1e-4
  @test z_cv_se ≈ z_cv_ss rtol=1e-4
  @test norm(z_cv_ss-z_cv_se) < norm(z_cv_ss-z_cv_double)

  #pyramid basis patch test
  SE_strategy = BEAST.WiltonSERule(
    tqd0[1,1],
    BEAST.DoubleQuadRule(
    tqd0[1,1],
    bqd1[1,2]))

  SS_strategy = SauterSchwabQuadrature.CommonVertex(BEAST._legendre(8,T(0.0),T(1.0)))

  Double_strategy = BEAST.DoubleQuadRule( #doublequadstrat
    tqd0[1,1],
    bqd1[1,2]
  )
  #singlelayer
  z_cv_se = zeros(complex(T),1,3)
  z_cv_ss = zeros(complex(T),1,3)
  z_cv_double = zeros(complex(T),1,3)

  BEAST.momintegrals!(op1, lag0, lag1, t1, t2, z_cv_se, SE_strategy)
  BEAST.momintegrals!(op1, lag0, lag1, t1, t2, z_cv_ss, SS_strategy)
  BEAST.momintegrals!(op1, lag0, lag1, t1, t2, z_cv_double, Double_strategy)

  @test z_cv_double ≈ z_cv_ss rtol = 1e-4
  @test z_cv_se ≈ z_cv_ss rtol=1e-4
  @test_broken norm(z_cv_ss-z_cv_se) < norm(z_cv_ss-z_cv_double)

  #double layer
  z_cv_se = zeros(complex(T),1,3)
  z_cv_ss = zeros(complex(T),1,3)
  z_cv_double = zeros(complex(T),1,3)

  BEAST.momintegrals!(op2, lag0, lag1, t1, t2, z_cv_se, SE_strategy)
  BEAST.momintegrals!(op2, lag0, lag1, t1, t2, z_cv_ss, SS_strategy)
  BEAST.momintegrals!(op2, lag0, lag1, t1, t2, z_cv_double, Double_strategy)

  @test z_cv_double ≈ z_cv_ss rtol = 1e-4
  @test z_cv_se ≈ z_cv_ss rtol=1e-4
  @test norm(z_cv_ss-z_cv_se) < norm(z_cv_ss-z_cv_double)

  #double layer transposed
  z_cv_se = zeros(complex(T),1,3)
  z_cv_ss = zeros(complex(T),1,3)
  z_cv_double = zeros(complex(T),1,3)

  BEAST.momintegrals!(op3, lag0, lag1, t1, t2, z_cv_se, SE_strategy)
  BEAST.momintegrals!(op3, lag0, lag1, t1, t2, z_cv_ss, SS_strategy)
  BEAST.momintegrals!(op3, lag0, lag1, t1, t2, z_cv_double, Double_strategy)

  @test z_cv_double ≈ z_cv_ss rtol = 1e-4
  @test z_cv_se ≈ z_cv_ss rtol=1e-4
  @test norm(z_cv_ss-z_cv_se) < norm(z_cv_ss-z_cv_double)

end

@testset "Common edge" begin
  # common edge
  t1 = simplex(
      T.(@SVector [0.180878, -0.941848, -0.283207]),
      T.(@SVector [0.0, -0.980785, -0.19509]),
      T.(@SVector [0.0, -0.92388, -0.382683])
      )
  t2 = simplex(
      T.(@SVector [0.158174, -0.881178, -0.44554]),
      T.(@SVector [0.180878, -0.941848, -0.283207]),
      T.(@SVector [0.0, -0.92388, -0.382683])
      )

  lag0 = BEAST.LagrangeRefSpace{T,0,3,1}() # patch basis
  lag1 = BEAST.LagrangeRefSpace{T,1,3,3}() #pyramid basis

  tqd0 = BEAST.quadpoints(lag0, [t1,t2], (12,)) #test quadrature data
  tqd1 = BEAST.quadpoints(lag1, [t1,t2], (12,))

  bqd0 = BEAST.quadpoints(lag0, [t1,t2], (13,)) #basis quadrature data
  bqd1 = BEAST.quadpoints(lag1, [t1,t2], (13,))

  #patch patch
  SE_strategy = BEAST.WiltonSERule(
    tqd0[1,1],
    BEAST.DoubleQuadRule(
    tqd0[1,1],
    bqd0[1,2]))

  SS_strategy = SauterSchwabQuadrature.CommonEdge(BEAST._legendre(8,T(0.0),T(1.0)))

  Double_strategy = BEAST.DoubleQuadRule( #doublequadstrat
    tqd0[1,1],
    bqd0[1,2]
  )
  #single layer

  z_ce_se = zeros(complex(T),1,1)
  z_ce_ss = zeros(complex(T),1,1)
  z_ce_double = zeros(complex(T),1,1)

  BEAST.momintegrals!(op1, lag0, lag0, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op1, lag0, lag0, t1, t2, z_ce_ss, SS_strategy)
  BEAST.momintegrals!(op1, lag0, lag0, t1, t2, z_ce_double, Double_strategy)

  @test z_ce_double ≈ z_ce_ss rtol = 1e-3
  @test z_ce_se ≈ z_ce_ss rtol=1e-4
  @test norm(z_ce_ss-z_ce_se) < norm(z_ce_ss-z_ce_double)

  #doublelayer
  z_ce_se = zeros(complex(T),1,1)
  z_ce_ss = zeros(complex(T),1,1)
  z_ce_double = zeros(complex(T),1,1)

  BEAST.momintegrals!(op2, lag0, lag0, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op2, lag0, lag0, t1, t2, z_ce_ss, SS_strategy)
  BEAST.momintegrals!(op2, lag0, lag0, t1, t2, z_ce_double, Double_strategy)

  @test z_ce_ss ≈ z_ce_double rtol = 1e-1
  @test z_ce_se ≈ z_ce_ss rtol = 1e-2
  @test norm(z_ce_ss-z_ce_se) < norm(z_ce_ss-z_ce_double)

  # doublelayer_transposed
  z_ce_se = zeros(complex(T),1,1)
  z_ce_ss = zeros(complex(T),1,1)
  z_ce_double = zeros(complex(T),1,1)

  BEAST.momintegrals!(op3, lag0, lag0, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op3, lag0, lag0, t1, t2, z_ce_ss, SS_strategy)
  BEAST.momintegrals!(op3, lag0, lag0, t1, t2, z_ce_double, Double_strategy)

  @test z_ce_double ≈ z_ce_ss rtol = 1e-1
  @test z_ce_se ≈ z_ce_ss rtol=1e-2
  @test norm(z_ce_ss-z_ce_se) < norm(z_ce_ss-z_ce_double)

#pyramid pyramid
  # singlelayer
  SE_strategy = BEAST.WiltonSERule(
    tqd1[1,1],
    BEAST.DoubleQuadRule(
    tqd1[1,1],
    bqd1[1,2]))
  Double_strategy = BEAST.DoubleQuadRule( #doublequadstrat
    tqd1[1,1],
    bqd1[1,2]
  )
  z_ce_se = zeros(complex(T),3,3)
  z_ce_ss = zeros(complex(T),3,3)
  z_ce_double = zeros(complex(T),3,3)

  BEAST.momintegrals!(op1, lag1, lag1, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op1, lag1, lag1, t1, t2, z_ce_ss, SS_strategy)
  BEAST.momintegrals!(op1, lag1, lag1, t1, t2, z_ce_double, Double_strategy)

  @test z_ce_double ≈ z_ce_ss rtol = 1e-3
  @test z_ce_se ≈ z_ce_double rtol=1e-3
  @test norm(z_ce_ss-z_ce_se) < norm(z_ce_ss-z_ce_double)

  #doublelayer 
  z_ce_se = zeros(complex(T),3,3)
  z_ce_ss = zeros(complex(T),3,3)
  z_ce_double = zeros(complex(T),3,3)

  BEAST.momintegrals!(op2, lag1, lag1, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op2, lag1, lag1, t1, t2, z_ce_ss, SS_strategy)
  BEAST.momintegrals!(op2, lag1, lag1, t1, t2, z_ce_double, Double_strategy)

  @test z_ce_double ≈ z_ce_ss rtol = 1e-1
  @test z_ce_se ≈ z_ce_ss rtol=1e-4
  @test norm(z_ce_ss-z_ce_se)/norm(z_ce_ss) < norm(z_ce_ss-z_ce_double)/norm(z_ce_ss)

  #doublelayer_transposed
  z_ce_se = zeros(complex(T),3,3)
  z_ce_ss = zeros(complex(T),3,3)
  z_ce_double = zeros(complex(T),3,3)

  BEAST.momintegrals!(op3, lag1, lag1, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op3, lag1, lag1, t1, t2, z_ce_ss, SS_strategy)
  BEAST.momintegrals!(op3, lag1, lag1, t1, t2, z_ce_double, Double_strategy)

  @test z_ce_double ≈ z_ce_ss rtol = 1e-1
  @test z_ce_ss ≈ z_ce_se rtol = 1e-2
  @test norm(z_ce_ss-z_ce_se) < norm(z_ce_ss-z_ce_double)

  #hypersingular
  z_ce_se = zeros(complex(T),3,3)
  z_ce_ss = zeros(complex(T),3,3)
  z_ce_double = zeros(complex(T),3,3)

  BEAST.momintegrals!(op4, lag1, lag1, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op4, lag1, lag1, t1, t2, z_ce_ss, SS_strategy)
  BEAST.momintegrals!(op4, lag1, lag1, t1, t2, z_ce_double, Double_strategy)

  @test z_ce_double ≈ z_ce_ss rtol = 1e-3
  @test z_ce_se ≈ z_ce_ss rtol=1e-4
  @test norm(z_ce_ss-z_ce_se) < norm(z_ce_ss-z_ce_double)

#pyramid test patch basis
  SE_strategy = BEAST.WiltonSERule(
    tqd1[1,1],
    BEAST.DoubleQuadRule(
    tqd1[1,1],
    bqd0[1,2]))

  SS_strategy = SauterSchwabQuadrature.CommonEdge(BEAST._legendre(8,T(0.0),T(1.0)))

  Double_strategy = BEAST.DoubleQuadRule( #doublequadstrat
    tqd1[1,1],
    bqd0[1,2]
  )
  #singlelayer
  z_ce_se = zeros(complex(T),3,1)
  z_ce_ss = zeros(complex(T),3,1)
  z_ce_double = zeros(complex(T),3,1)

  BEAST.momintegrals!(op1, lag1, lag0, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op1, lag1, lag0, t1, t2, z_ce_ss, SS_strategy)
  BEAST.momintegrals!(op1, lag1, lag0, t1, t2, z_ce_double, Double_strategy)

  @test z_ce_double ≈ z_ce_ss rtol = 1e-3
  @test z_ce_se ≈ z_ce_ss rtol=1e-4
  @test norm(z_ce_ss-z_ce_se) < norm(z_ce_ss-z_ce_double)
  
  #doublelayer
  z_ce_se = zeros(complex(T),3,1)
  z_ce_ss = zeros(complex(T),3,1)
  z_ce_double = zeros(complex(T),3,1)

  BEAST.momintegrals!(op2, lag1, lag0, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op2, lag1, lag0, t1, t2, z_ce_ss, SS_strategy)
  BEAST.momintegrals!(op2, lag1, lag0, t1, t2, z_ce_double, Double_strategy)

  @test z_ce_double ≈ z_ce_ss rtol = 1e-1
  @test z_ce_se ≈ z_ce_ss rtol=1e-4
  @test norm(z_ce_ss-z_ce_se)/norm(z_ce_ss) < norm(z_ce_ss-z_ce_double)/norm(z_ce_ss)
  #doublelayer_transposed
  
  z_ce_se = zeros(complex(T),3,1)
  z_ce_ss = zeros(complex(T),3,1)
  z_ce_double = zeros(complex(T),3,1)

  BEAST.momintegrals!(op3, lag1, lag0, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op3, lag1, lag0, t1, t2, z_ce_ss, SS_strategy)
  BEAST.momintegrals!(op3, lag1, lag0, t1, t2, z_ce_double, Double_strategy)

  @test z_ce_double ≈ z_ce_ss rtol = 1e-1
  @test z_ce_se ≈ z_ce_ss rtol=1e-2
  @test norm(z_ce_ss-z_ce_se)/norm(z_ce_ss) < norm(z_ce_ss-z_ce_double)/norm(z_ce_ss)

  #pyramid basis patch test
  SE_strategy = BEAST.WiltonSERule(
    tqd0[1,1],
    BEAST.DoubleQuadRule(
    tqd0[1,1],
    bqd1[1,2]))

  SS_strategy = SauterSchwabQuadrature.CommonEdge(BEAST._legendre(8,T(0.0),T(1.0)))

  Double_strategy = BEAST.DoubleQuadRule( #doublequadstrat
    tqd0[1,1],
    bqd1[1,2]
  )
  #singlelayer
  z_ce_se = zeros(complex(T),1,3)
  z_ce_ss = zeros(complex(T),1,3)
  z_ce_double = zeros(complex(T),1,3)

  BEAST.momintegrals!(op1, lag0, lag1, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op1, lag0, lag1, t1, t2, z_ce_ss, SS_strategy)
  BEAST.momintegrals!(op1, lag0, lag1, t1, t2, z_ce_double, Double_strategy)

  @test z_ce_double ≈ z_ce_ss rtol = 1e-3
  @test z_ce_se ≈ z_ce_ss rtol=1e-4
  @test norm(z_ce_ss-z_ce_se)/norm(z_ce_ss) < norm(z_ce_ss-z_ce_double)/norm(z_ce_ss)

  #double layer
  z_ce_se = zeros(complex(T),1,3)
  z_ce_ss = zeros(complex(T),1,3)
  z_ce_double = zeros(complex(T),1,3)

  BEAST.momintegrals!(op2, lag0, lag1, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op2, lag0, lag1, t1, t2, z_ce_ss, SS_strategy)
  BEAST.momintegrals!(op2, lag0, lag1, t1, t2, z_ce_double, Double_strategy)

  @test z_ce_double ≈ z_ce_ss rtol = 1e-1
  @test z_ce_se ≈ z_ce_ss rtol=1e-4
  @test norm(z_ce_ss-z_ce_se)/norm(z_ce_ss) < norm(z_ce_ss-z_ce_double)/norm(z_ce_ss)

  #double layer transposed
  z_ce_se = zeros(complex(T),1,3)
  z_ce_ss = zeros(complex(T),1,3)
  z_ce_double = zeros(complex(T),1,3)

  BEAST.momintegrals!(op3, lag0, lag1, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op3, lag0, lag1, t1, t2, z_ce_ss, SS_strategy)
  BEAST.momintegrals!(op3, lag0, lag1, t1, t2, z_ce_double, Double_strategy)

  @test z_ce_double ≈ z_ce_ss rtol = 1e-1
  @test z_ce_se ≈ z_ce_ss rtol=1e-1
  @test norm(z_ce_ss-z_ce_se)/norm(z_ce_ss) < norm(z_ce_ss-z_ce_double)/norm(z_ce_ss)

end

@testset "Seperated triangles" begin
  t1 = simplex(
      T.(@SVector [0.0, -0.980785, -0.19509]),
      T.(@SVector [0.180878, -0.941848, -0.283207]),
      T.(@SVector [0.0, -0.92388, -0.382683]))
  t2 = simplex(
      T.(@SVector [0.0, -0.980785, -0.230102]),
      T.(@SVector [0.180878, -0.941848, -0.383207]),
      T.(@SVector [0.0, -0.92388, -0.411962]))

  lag0 = BEAST.LagrangeRefSpace{T,0,3,1}() # patch basis
  lag1 = BEAST.LagrangeRefSpace{T,1,3,3}() #pyramid basis
  
  tqd0 = BEAST.quadpoints(lag0, [t1,t2], (12,)) #test quadrature data
  tqd1 = BEAST.quadpoints(lag1, [t1,t2], (12,))
  
  bqd0 = BEAST.quadpoints(lag0, [t1,t2], (13,)) #basis quadrature data
  bqd1 = BEAST.quadpoints(lag1, [t1,t2], (13,))
  # singlelayer
  SE_strategy = BEAST.WiltonSERule(
    tqd0[1,1],
    BEAST.DoubleQuadRule(
    tqd0[1,1],
    bqd0[1,2]))
  
  Double_strategy = BEAST.DoubleQuadRule( #doublequadstrat
    tqd0[1,1],
    bqd0[1,2]
  )
  z_ce_se = zeros(complex(T),1,1)
  z_ce_double = zeros(complex(T),1,1)
  
  BEAST.momintegrals!(op1, lag0, lag0, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op1, lag0, lag0, t1, t2, z_ce_double, Double_strategy)
  
  @test z_ce_se ≈ z_ce_double rtol=1e-6

  # calculate a reference value with hcubature to desired accuracy
  #=using HCubature
  r(u,v,t) = (1-u)*t[1]+u*((1-v)*t[2]+v*t[3])
  gamma=1.0+1.0im
  jac(u,v,t) = u * norm(cross(t[1],t[2])+cross(t[2],t[3])+cross(t[3],t[1]))
  function outerf(x)
    r2(x2) = r(x2[1],x2[2],t1) - r(x[1],x[2],t2)
    innerf(x2) = exp(-gamma*norm(r2(x2)))/(4*pi*norm(r2(x2))) * jac(x2[1],x2[2],t1)
    hcubature(innerf, (0.0,0.0), (1.0,1.0), rtol = 1e-8)[1] * jac(x[1],x[2],t2)
  end

  @show refval = hcubature(outerf, (0.0,0.0), (1.0,1.0), rtol = 1e-8)[1] =#
  refval = 0.00033563892481993545 - 2.22592609372769e-5im # can be calculated with above code
  
  @test z_ce_se[1] ≈ refval rtol = 1e-7
  @test z_ce_double[1] ≈ refval rtol = 1e-6


  SE_strategy = BEAST.WiltonSERule(
    tqd1[1,1],
    BEAST.DoubleQuadRule(
    tqd1[1,1],
    bqd1[1,2]))

  Double_strategy = BEAST.DoubleQuadRule( #doublequadstrat
    tqd1[1,1],
    bqd1[1,2]
  )
  z_ce_se = zeros(complex(T),3,3)
  z_ce_double = zeros(complex(T),3,3)

  BEAST.momintegrals!(op1, lag1, lag1, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op1, lag1, lag1, t1, t2, z_ce_double, Double_strategy)

  @test z_ce_se ≈ z_ce_double rtol=1e-4

  # doublelayer
  SE_strategy = BEAST.WiltonSERule(
    tqd0[1,1],
    BEAST.DoubleQuadRule(
    tqd0[1,1],
    bqd1[1,2]))

  Double_strategy = BEAST.DoubleQuadRule(
    tqd0[1,1],
    bqd1[1,2])
  z_ce_se = zeros(complex(T),1,3)
  z_ce_double = zeros(complex(T),1,3)
  BEAST.momintegrals!(op2, lag0, lag1, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op2, lag0, lag1, t1, t2, z_ce_double, Double_strategy)
  
  @test z_ce_se ≈ z_ce_double rtol = 1e-4

  SE_strategy = BEAST.WiltonSERule(
    tqd1[1,1],
    BEAST.DoubleQuadRule(
    tqd1[1,1],
    bqd1[1,2]))

  Double_strategy = BEAST.DoubleQuadRule(
    tqd1[1,1],
    bqd1[1,2])
  z_ce_se = zeros(complex(T),3,3)
  z_ce_double = zeros(complex(T),3,3)
  BEAST.momintegrals!(op2, lag1, lag1, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op2, lag1, lag1, t1, t2, z_ce_double, Double_strategy)

  @test z_ce_se ≈ z_ce_double rtol = 1e-4
  # doublelayer transposed
  SE_strategy = BEAST.WiltonSERule(
      tqd1[1,1],
      BEAST.DoubleQuadRule(
          tqd1[1,1],
          bqd0[1,2]))
  Double_strategy = BEAST.DoubleQuadRule(
    tqd1[1,1],
    bqd0[1,2])
  z_ce_se = zeros(complex(T),3,1)
  z_ce_double = zeros(complex(T),3,1)

  BEAST.momintegrals!(op3, lag1, lag0, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op3, lag1, lag0, t1, t2, z_ce_double, Double_strategy)

  @test z_ce_se ≈ z_ce_double rtol=1e-4

  SE_strategy = BEAST.WiltonSERule(
      tqd1[1,1],
      BEAST.DoubleQuadRule(
          tqd1[1,1],
          bqd1[1,2]))
  Double_strategy = BEAST.DoubleQuadRule(
    tqd1[1,1],
    bqd1[1,2])
  z_ce_se = zeros(complex(T),3,3)
  z_ce_double = zeros(complex(T),3,3)

  BEAST.momintegrals!(op3, lag1, lag1, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op3, lag1, lag1, t1, t2, z_ce_double, Double_strategy)

  @test z_ce_se ≈ z_ce_double rtol=1e-4
  # hypersingular
  SE_strategy = BEAST.WiltonSERule(
      tqd1[1,1],
      BEAST.DoubleQuadRule(
          tqd1[1,1],
          bqd1[1,2]))
  Double_strategy = BEAST.DoubleQuadRule(
    tqd1[1,1],
    bqd1[1,2])

  z_ce_se = zeros(complex(T),3,3)
  z_ce_double = zeros(complex(T),3,3)

  BEAST.momintegrals!(op4, lag1, lag1, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op4, lag1, lag1, t1, t2, z_ce_double, Double_strategy)
  @test z_ce_se ≈ z_ce_double rtol = 1e-5

end
