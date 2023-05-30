using Test
using LinearAlgebra

using BEAST, CompScienceMeshes, SauterSchwabQuadrature, StaticArrays

T=Float64


# operators
op1 = Helmholtz3D.singlelayer(gamma=T(2.0))
op2 = Helmholtz3D.doublelayer(gamma=T(2.0))
op3 = Helmholtz3D.doublelayer_transposed(gamma=T(2.0))
op4 = Helmholtz3D.hypersingular(gamma=T(2.0))
## Common face case:
@testset "Common Face" begin
  # triangles
t1 = simplex(
  T.(@SVector [0.180878, -0.941848, -0.283207]),
  T.(@SVector [0.0, -0.980785, -0.19509]),
  T.(@SVector [0.0, -0.92388, -0.382683]))

  lag0 = BEAST.LagrangeRefSpace{T,0,3,1}() # patch basis
  lag1 = BEAST.LagrangeRefSpace{T,1,3,3}() #pyramid basis

  tqd0 = BEAST.quadpoints(lag0, [t1], (12,)) #test quadrature data
  tqd1 = BEAST.quadpoints(lag1, [t1], (12,))

  bqd0 = BEAST.quadpoints(lag0, [t1], (13,)) #basis quadrature data
  bqd1 = BEAST.quadpoints(lag1, [t1], (13,))

  SE_strategy = BEAST.WiltonSERule( #wilton
    tqd0[1,1],
    BEAST.DoubleQuadRule(
    tqd0[1,1],
    bqd0[1,1],
    ),
  )

  SS_strategy = SauterSchwabQuadrature.CommonFace(BEAST._legendre(8,T(0.0),T(1.0))) #sauter

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
end
# common vertex case
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
  #single layer
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
  z_cv_se = zeros(3,1)
  z_cv_ss = zeros(3,1)
  z_cv_double = zeros(3,1)

  BEAST.momintegrals!(op1, lag0, lag0, t1, t2, z_cv_se, SE_strategy)
  BEAST.momintegrals!(op1, lag0, lag0, t1, t2, z_cv_ss, SS_strategy)
  BEAST.momintegrals!(op1, lag0, lag0, t1, t2, z_cv_double, Double_strategy)

  @test z_cv_double ≈ z_cv_ss rtol = 1e-4
  @test z_cv_se ≈ z_cv_ss rtol=1e-4
  #@test norm(z_cv_ss-z_cv_se)/norm(z_cv_ss) < norm(z_cv_ss-z_cv_double)/norm(z_cv_ss)
  #double layer
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
  z_cv_se = zeros(1,3)
  z_cv_ss = zeros(1,3)
  z_cv_double = zeros(1,3)

  BEAST.momintegrals!(op2, lag0, lag1, t1, t2, z_cv_se, SE_strategy)
  BEAST.momintegrals!(op2, lag0, lag1, t1, t2, z_cv_ss, SS_strategy)
  BEAST.momintegrals!(op2, lag0, lag1, t1, t2, z_cv_double, Double_strategy)

  @test z_cv_double ≈ z_cv_ss rtol = 1e-4
  @test z_cv_se ≈ z_cv_ss rtol=1e-4
  @test norm(z_cv_ss-z_cv_se)/norm(z_cv_ss) < norm(z_cv_ss-z_cv_double)/norm(z_cv_ss)

  #double layer transposed
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
  z_cv_se = zeros(3,1)
  z_cv_ss = zeros(3,1)
  z_cv_double = zeros(3,1)

  BEAST.momintegrals!(op3, lag1, lag0, t1, t2, z_cv_se, SE_strategy)
  BEAST.momintegrals!(op3, lag1, lag0, t1, t2, z_cv_ss, SS_strategy)
  BEAST.momintegrals!(op3, lag1, lag0, t1, t2, z_cv_double, Double_strategy)

  @test z_cv_double ≈ z_cv_ss rtol = 1e-4
  @test z_cv_se ≈ z_cv_ss rtol=1e-4
  @test norm(z_cv_ss-z_cv_se)/norm(z_cv_ss) < norm(z_cv_ss-z_cv_double)/norm(z_cv_ss)

  #hypersingular
  SE_strategy = BEAST.WiltonSERule(
    tqd1[1,1],
    BEAST.DoubleQuadRule(
    tqd1[1,1],
    bqd1[1,2]))

  SS_strategy = SauterSchwabQuadrature.CommonVertex(BEAST._legendre(8,T(0.0),T(1.0)))

  Double_strategy = BEAST.DoubleQuadRule( #doublequadstrat
    tqd1[1,1],
    bqd1[1,2]
  )
  z_cv_se = zeros(3,3)
  z_cv_ss = zeros(3,3)
  z_cv_double = zeros(3,3)

  BEAST.momintegrals!(op4, lag1, lag1, t1, t2, z_cv_se, SE_strategy)
  BEAST.momintegrals!(op4, lag1, lag1, t1, t2, z_cv_ss, SS_strategy)
  BEAST.momintegrals!(op4, lag1, lag1, t1, t2, z_cv_double, Double_strategy)

  @test z_cv_double ≈ z_cv_ss rtol = 1e-4
  @test z_cv_se ≈ z_cv_ss rtol=1e-4
  #@test norm(z_cv_ss-z_cv_se)/norm(z_cv_ss) < norm(z_cv_ss-z_cv_double)/norm(z_cv_ss)

end

@testset "Common edge" begin
  # common edge
  t1 = simplex(
      T.(@SVector [0.180878, -0.941848, -0.283207]),
      T.(@SVector [0.0, -0.980785, -0.19509]),
      T.(@SVector [0.0, -0.92388, -0.382683])
      )
  t2 = simplex(
      T.(@SVector [0.180878, -0.941848, -0.283207]),
      T.(@SVector [0.158174, -0.881178, -0.44554]),
      T.(@SVector [0.0, -0.92388, -0.382683])
      )

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

  SS_strategy = SauterSchwabQuadrature.CommonEdge(BEAST._legendre(8,T(0.0),T(1.0)))

  Double_strategy = BEAST.DoubleQuadRule( #doublequadstrat
    tqd0[1,1],
    bqd0[1,2]
  )
  z_ce_se = zeros(3,1)
  z_ce_ss = zeros(3,1)
  z_ce_double = zeros(3,1)

  BEAST.momintegrals!(op1, lag0, lag0, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op1, lag0, lag0, t1, t2, z_ce_ss, SS_strategy)
  BEAST.momintegrals!(op1, lag0, lag0, t1, t2, z_ce_double, Double_strategy)

  @test z_ce_se ≈ z_ce_ss rtol=1e-4
  @test z_ce_double ≈ z_ce_ss rtol = 1e-3
  @test norm(z_ce_ss-z_ce_se)/norm(z_ce_ss) < norm(z_ce_ss-z_ce_double)/norm(z_ce_ss)

  # doublelayer
  SE_strategy = BEAST.WiltonSERule(
    tqd0[1,1],
    BEAST.DoubleQuadRule(
    tqd0[1,1],
    bqd1[1,2]))
  SS_strategy = SauterSchwabQuadrature.CommonEdge(BEAST._legendre(8,T(0.0),T(1.0)))
  Double_strategy = BEAST.DoubleQuadRule(
    tqd0[1,1],
    bqd1[1,2])
  z_ce_se = zeros(1,3)
  z_ce_ss = zeros(1,3)
  z_ce_double = zeros(1,3)
  BEAST.momintegrals!(op2, lag0, lag1, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op2, lag0, lag1, t1, t2, z_ce_ss, SS_strategy)
  BEAST.momintegrals!(op2, lag0, lag1, t1, t2, z_ce_double, Double_strategy)

  @test z_ce_se ≈ z_ce_ss rtol = 1e-4
  @test z_ce_ss ≈ z_ce_double rtol = 1e-1
  @test norm(z_ce_ss-z_ce_se)/norm(z_ce_ss) < norm(z_ce_ss-z_ce_double)/norm(z_ce_ss)

  # doublelayer transposed
  SE_strategy = BEAST.WiltonSERule(
      tqd1[1,1],
      BEAST.DoubleQuadRule(
          tqd1[1,1],
          bqd0[1,2]))
  SS_strategy = SauterSchwabQuadrature.CommonEdge(BEAST._legendre(12,T(0.0),T(1.0)))
  Double_strategy = BEAST.DoubleQuadRule(
    tqd1[1,1],
    bqd0[1,2])
  z_ce_se = zeros(3,1)
  z_ce_ss = zeros(3,1)
  z_ce_double = zeros(3,1)

  BEAST.momintegrals!(op3, lag1, lag0, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op3, lag1, lag0, t1, t2, z_ce_ss, SS_strategy)
  BEAST.momintegrals!(op3, lag1, lag0, t1, t2, z_ce_double, Double_strategy)

  @test z_ce_se ≈ z_ce_ss rtol=1e-2
  @test z_ce_ss ≈ z_ce_double rtol = 1e-1
  @test norm(z_ce_ss-z_ce_se)/norm(z_ce_ss) < norm(z_ce_ss-z_ce_double)/norm(z_ce_ss)

  # hypersingular

  SE_strategy = BEAST.WiltonSERule(
      tqd1[1,1],
      BEAST.DoubleQuadRule(
          tqd1[1,1],
          bqd1[1,2]))
  SS_strategy = SauterSchwabQuadrature.CommonEdge(BEAST._legendre(12,T(0.0),T(1.0)))
  Double_strategy = BEAST.DoubleQuadRule(
    tqd1[1,1],
    bqd1[1,2])

  z_ce_se = zeros(3,3)
  z_ce_ss = zeros(3,3)
  z_ce_double = zeros(3,3)

  BEAST.momintegrals!(op4, lag1, lag1, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op4, lag1, lag1, t1, t2, z_ce_ss, SS_strategy)
  BEAST.momintegrals!(op4, lag1, lag1, t1, t2, z_ce_double, Double_strategy)
  @test z_ce_se ≈ z_ce_ss rtol = 1e-5
  @test z_ce_double ≈ z_ce_double rtol = 1e-5
  @test norm(z_ce_se-z_ce_ss) < norm(z_ce_double-z_ce_ss)

end

@testset "Seperated triangles" begin
  t1 = simplex(
      T.(@SVector [0.180878, -0.941848, -0.283207]),
      T.(@SVector [0.0, -0.980785, -0.19509]),
      T.(@SVector [0.0, -0.92388, -0.382683]))
  t2 = simplex(
      T.(@SVector [0.373086, -0.881524, -1.289348]),
      T.(@SVector [0.180878, -0.941848, -1.283207]),
      T.(@SVector [0.294908, -0.944921, -1.141962]))

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
  z_ce_se = zeros(3,1)
  z_ce_double = zeros(3,1)
  
  BEAST.momintegrals!(op1, lag0, lag0, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op1, lag0, lag0, t1, t2, z_ce_double, Double_strategy)
  
  @test z_ce_se ≈ z_ce_double rtol=1e-7
  
  # doublelayer
  SE_strategy = BEAST.WiltonSERule(
    tqd0[1,1],
    BEAST.DoubleQuadRule(
    tqd0[1,1],
    bqd1[1,2]))

  Double_strategy = BEAST.DoubleQuadRule(
    tqd0[1,1],
    bqd1[1,2])
  z_ce_se = zeros(1,3)
  z_ce_double = zeros(1,3)
  BEAST.momintegrals!(op2, lag0, lag1, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op2, lag0, lag1, t1, t2, z_ce_double, Double_strategy)
  
  @test z_ce_se ≈ z_ce_double rtol = 1e-8
  
  # doublelayer transposed
  SE_strategy = BEAST.WiltonSERule(
      tqd1[1,1],
      BEAST.DoubleQuadRule(
          tqd1[1,1],
          bqd0[1,2]))
  Double_strategy = BEAST.DoubleQuadRule(
    tqd1[1,1],
    bqd0[1,2])
  z_ce_se = zeros(3,1)
  z_ce_double = zeros(3,1)
  
  BEAST.momintegrals!(op3, lag1, lag0, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op3, lag1, lag0, t1, t2, z_ce_double, Double_strategy)
  
  @test z_ce_se ≈ z_ce_double rtol=1e-7
  
  # hypersingular
  SE_strategy = BEAST.WiltonSERule(
      tqd1[1,1],
      BEAST.DoubleQuadRule(
          tqd1[1,1],
          bqd1[1,2]))
  Double_strategy = BEAST.DoubleQuadRule(
    tqd1[1,1],
    bqd1[1,2])
  
  z_ce_se = zeros(3,3)
  z_ce_double = zeros(3,3)
  
  BEAST.momintegrals!(op4, lag1, lag1, t1, t2, z_ce_se, SE_strategy)
  BEAST.momintegrals!(op4, lag1, lag1, t1, t2, z_ce_double, Double_strategy)
  @test z_ce_se ≈ z_ce_double rtol = 1e-7

end
