using Test
using LinearAlgebra

using BEAST, CompScienceMeshes, StaticArrays


# Float32 not working since hankelh2 returns always F64
for T in [Float64]

#T=Float64

t1 = simplex(
  T.(@SVector [1.0, 0.0]),
  T.(@SVector [6.123233995736766e-17, 1.0]))

#t2 = simplex(
 # T.(@SVector [6.123233995736766e-17, 1.0]),
 # T.(@SVector [-1.0, 1.2246467991473532e-16]))

#t3 = simplex(
 # T.(@SVector [-1.0, 1.2246467991473532e-16]),
 # T.(@SVector [-1.8369701987210297e-16, -1.0]))
  
t4 = simplex(
  T.(@SVector [-1.8369701987210297e-16, -1.0]),
  T.(@SVector [1.0, 0.0]))

t5 = simplex(
  T.(@SVector [1.0, 0.0]),
  T.(@SVector [-1.8369701987210297e-16, -1.0]))  

# Common Edge case:
#L0
Degr = 0
Dim1 = 1
NumF = 1

sp0 = BEAST.LagrangeRefSpace{T,Degr,Dim1,NumF}()

tqd0 = BEAST.quadpoints(sp0, [t1,t1], (1000,))
bqd0 = BEAST.quadpoints(sp0, [t1,t1], (1001,))

op1 = BEAST.Helmholtz2D.singlelayer(; wavenumber=T(1.0)) 
op4 = Helmholtz2D.hypersingular(; wavenumber=T(1.0))

gl_strategy_0 = BEAST.DoubleQuadRule(
	tqd0[1,1],
	bqd0[1,1],
  )

mrw_strategy_0 = BEAST.SauterSchwabQuadrature1D.CommonEdge(BEAST.SauterSchwabQuadrature1D._NRWrules(10,T(0.0),T(1.0)) )

#op1 
z_gl_10 = zeros(Complex{T}, 1, 1)
z_mrw_10 = zeros(Complex{T}, 1, 1)

BEAST.momintegrals!(op1, sp0, sp0, t1, t1, z_gl_10, gl_strategy_0)
BEAST.momintegrals!(op1, sp0, sp0, t1, t1, z_mrw_10, mrw_strategy_0)

@test abs(z_gl_10[1,1] - z_mrw_10[1,1]) / abs(z_gl_10[1,1]) ≈ 0.0 atol=1e-3

#=
@show z_gl_10
@show z_mrw_10
@show abs(z_gl_10[1,1] - z_mrw_10[1,1]) / abs(z_gl_10[1,1])
=#

#op4
z_gl_40 = zeros(Complex{T}, 1, 1)
z_mrw_40 = zeros(Complex{T}, 1, 1)

BEAST.momintegrals!(op4, sp0, sp0, t1, t1, z_gl_40, gl_strategy_0)
BEAST.momintegrals!(op4, sp0, sp0, t1, t1, z_mrw_40, mrw_strategy_0)

@test abs(z_gl_40[1,1] - z_mrw_40[1,1]) / abs(z_gl_40[1,1]) ≈ 0.0 atol=1e-3

#=
@show z_gl_40
@show z_mrw_40
@show abs(z_gl_40[1,1] - z_mrw_40[1,1]) / abs(z_gl_40[1,1])
=#

# L1 
Degr = 1
Dim1 = 2
NumF = 2

sp1 = BEAST.LagrangeRefSpace{T,Degr,Dim1,NumF}()

tqd1 = BEAST.quadpoints(sp1, [t1,t1], (1000,))
bqd1 = BEAST.quadpoints(sp1, [t1,t1], (1001,))

gl_strategy_1 = BEAST.DoubleQuadRule(
	tqd1[1,1],
	bqd1[1,1],
  )

mrw_strategy_1 = BEAST.SauterSchwabQuadrature1D.CommonEdge(BEAST.SauterSchwabQuadrature1D._NRWrules(10,T(0.0),T(1.0)) )

#op1
z_gl_11 = zeros(Complex{T}, 2, 2)
z_mrw_11 = zeros(Complex{T}, 2, 2)
  
BEAST.momintegrals!(op1, sp1, sp1, t1, t1, z_gl_11, gl_strategy_1)
BEAST.momintegrals!(op1, sp1, sp1, t1, t1, z_mrw_11, mrw_strategy_1)
  
@test abs(z_gl_11[1,1] - z_mrw_11[1,1]) / abs(z_gl_11[1,1]) ≈ 0.0 atol=1e-3
@test abs(z_gl_11[1,2] - z_mrw_11[1,2]) / abs(z_gl_11[1,2]) ≈ 0.0 atol=1e-3
@test abs(z_gl_11[2,1] - z_mrw_11[2,1]) / abs(z_gl_11[2,1]) ≈ 0.0 atol=1e-3
@test abs(z_gl_11[2,2] - z_mrw_11[2,2]) / abs(z_gl_11[2,2]) ≈ 0.0 atol=1e-3

#=
@show z_gl_11
@show z_mrw_11
@show abs(z_gl_11[1,1] - z_mrw_11[1,1]) / abs(z_gl_11[1,1])
@show abs(z_gl_11[1,2] - z_mrw_11[1,2]) / abs(z_gl_11[1,2])
@show abs(z_gl_11[1,2] - z_mrw_11[1,2]) / abs(z_gl_11[1,2])
@show abs(z_gl_11[2,1] - z_mrw_11[2,1]) / abs(z_gl_11[2,1])
=#

#op4
z_gl_41 = zeros(Complex{T}, 2, 2)
z_mrw_41 = zeros(Complex{T}, 2, 2)
  
BEAST.momintegrals!(op4, sp1, sp1, t1, t1, z_gl_41, gl_strategy_1)
BEAST.momintegrals!(op4, sp1, sp1, t1, t1, z_mrw_41, mrw_strategy_1)
  
@test abs(z_gl_41[1,1] - z_mrw_41[1,1]) / abs(z_gl_41[1,1]) ≈ 0.0 atol=1e-3
@test abs(z_gl_41[1,2] - z_mrw_41[1,2]) / abs(z_gl_41[1,2]) ≈ 0.0 atol=1e-3
@test abs(z_gl_41[2,1] - z_mrw_41[2,1]) / abs(z_gl_41[2,1]) ≈ 0.0 atol=1e-3
@test abs(z_gl_41[2,2] - z_mrw_41[2,2]) / abs(z_gl_41[2,2]) ≈ 0.0 atol=1e-3

#=
@show z_gl_41
@show z_mrw_41
@show abs(z_gl_41[1,1] - z_mrw_41[1,1]) / abs(z_gl_41[1,1])
@show abs(z_gl_41[1,2] - z_mrw_41[1,2]) / abs(z_gl_41[1,2])
@show abs(z_gl_41[1,2] - z_mrw_41[1,2]) / abs(z_gl_41[1,2])
@show abs(z_gl_41[2,1] - z_mrw_41[2,1]) / abs(z_gl_41[2,1])
=#

# common Vertex case
#L0
Degr = 0
Dim1 = 1
NumF = 1

sp0 = BEAST.LagrangeRefSpace{T,Degr,Dim1,NumF}()

tqd0_14= BEAST.quadpoints(sp0, [t1,t4], (1000,))
bqd0_14 = BEAST.quadpoints(sp0, [t1,t4], (1001,))

κ = T(1.0)
op1 = BEAST.Helmholtz2D.singlelayer(; wavenumber=κ) 
#op2 = Helmholtz2D.doublelayer(; wavenumber=T(1.0)) 
#op3 = Helmholtz2D.doublelayer_transposed(; wavenumber=T(1.0))
op4 = Helmholtz2D.hypersingular(; wavenumber=T(1.0))

gl_strategy_0 = BEAST.DoubleQuadRule(
	tqd0_14[1,1],
	bqd0_14[1,2],
  )


mrw_strategy_0 = BEAST.SauterSchwabQuadrature1D.CommonVertex(BEAST.SauterSchwabQuadrature1D._NRWrules(10,T(0.0),T(1.0)) )

#op1
z_gl_10_14 = zeros(Complex{T}, 1, 1)
z_mrw_10_14 = zeros(Complex{T}, 1, 1)


BEAST.momintegrals!(op1, sp0, sp0, t1, t4, z_gl_10_14, gl_strategy_0)
BEAST.momintegrals!(op1, sp0, sp0, t1, t4, z_mrw_10_14, mrw_strategy_0)

@test abs(z_gl_10_14[1,1] - z_mrw_10_14[1,1]) / abs(z_gl_10_14[1,1]) ≈ 0.0 atol=1e-3

#=
@show z_gl_10_14
@show z_mrw_10_14
@show abs(z_gl_10_14[1,1] - z_mrw_10_14[1,1]) / abs(z_gl_10_14[1,1])
=#

#op4
z_gl_40_14 = zeros(Complex{T}, 1, 1)
z_mrw_40_14 = zeros(Complex{T}, 1, 1)


BEAST.momintegrals!(op4, sp0, sp0, t1, t4, z_gl_40_14, gl_strategy_0)
BEAST.momintegrals!(op4, sp0, sp0, t1, t4, z_mrw_40_14, mrw_strategy_0)

@test abs(z_gl_40_14[1,1] - z_mrw_40_14[1,1]) / abs(z_gl_40_14[1,1]) ≈ 0.0 atol=1e-3

#=
@show z_gl_40_14
@show z_mrw_40_14
@show abs(z_gl_40_14[1,1] - z_mrw_40_14[1,1]) / abs(z_gl_40_14[1,1])
=#

# improper allocation of vertices in [t1,t5], so we need to flip them first

tqd0_15= BEAST.quadpoints(sp0, [t1,t5], (1000,))
bqd0_15 = BEAST.quadpoints(sp0, [t1,t5], (1001,))

gl_strategy_015 = BEAST.DoubleQuadRule(
	tqd0_15[1,1],
	bqd0_15[1,2],
  )

#op1
z_gl_10_15 = zeros(Complex{T}, 1, 1)
z_mrw_10_15 = zeros(Complex{T}, 1, 1)


BEAST.momintegrals!(op1, sp0, sp0, t1, t5, z_gl_10_15, gl_strategy_015)
BEAST.momintegrals!(op1, sp0, sp0, t1, t5, z_mrw_10_15, mrw_strategy_0)

@test abs(z_gl_10_15[1,1] - z_mrw_10_15[1,1]) / abs(z_gl_10_15[1,1]) ≈ 0.0 atol=1e-3

#op4
z_gl_40_15 = zeros(Complex{T}, 1, 1)
z_mrw_40_15 = zeros(Complex{T}, 1, 1)

BEAST.momintegrals!(op4, sp0, sp0, t1, t5, z_gl_40_15, gl_strategy_015)
BEAST.momintegrals!(op4, sp0, sp0, t1, t5, z_mrw_40_15, mrw_strategy_0)

@test abs(z_gl_40_15[1,1] - z_mrw_40_15[1,1]) / abs(z_gl_40_15[1,1]) ≈ 0.0 atol=1e-3

#L1  
Degr = 1
Dim1 = 2
NumF = 2

sp1 = BEAST.LagrangeRefSpace{T,Degr,Dim1,NumF}()

tqd1_14 = BEAST.quadpoints(sp1, [t1,t4], (1000,))
bqd1_14 = BEAST.quadpoints(sp1, [t1,t4], (1001,))

gl_strategy_1 = BEAST.DoubleQuadRule(
	tqd1_14[1,1],
	bqd1_14[1,2],
  )


mrw_strategy_1 = BEAST.SauterSchwabQuadrature1D.CommonVertex(BEAST.SauterSchwabQuadrature1D._NRWrules(10,T(0.0),T(1.0)) )

#op1
z_gl_11_14 = zeros(Complex{T}, 2, 2)
z_mrw_11_14 = zeros(Complex{T}, 2, 2)
  
  
BEAST.momintegrals!(op1, sp1, sp1, t1, t4, z_gl_11_14, gl_strategy_1)
BEAST.momintegrals!(op1, sp1, sp1, t1, t4, z_mrw_11_14, mrw_strategy_1)
  
@test abs(z_gl_11_14[1,1] - z_mrw_11_14[1,1]) / abs(z_gl_11_14[1,1]) ≈ 0.0 atol=1e-3
@test abs(z_gl_11_14[1,2] - z_mrw_11_14[1,2]) / abs(z_gl_11_14[1,2]) ≈ 0.0 atol=1e-3
@test abs(z_gl_11_14[2,1] - z_mrw_11_14[2,1]) / abs(z_gl_11_14[2,1]) ≈ 0.0 atol=1e-3
@test abs(z_gl_11_14[2,2] - z_mrw_11_14[2,2]) / abs(z_gl_11_14[2,2]) ≈ 0.0 atol=1e-3

#=
@show z_gl_11_14
@show z_mrw_11_14
@show abs(z_gl_11_14[1,1] - z_mrw_11_14[1,1]) / abs(z_mrw_11_14[1,1])
@show abs(z_gl_11_14[1,2] - z_mrw_11_14[1,2]) / abs(z_gl_11_14[1,2])
@show abs(z_gl_11_14[1,2] - z_mrw_11_14[1,2]) / abs(z_gl_11_14[1,2])
@show abs(z_gl_11_14[2,1] - z_mrw_11_14[2,1]) / abs(z_gl_11_14[2,1])
=#

 #op4
 z_gl_41_14 = zeros(Complex{T}, 2, 2)
 z_mrw_41_14 = zeros(Complex{T}, 2, 2)
  
 BEAST.momintegrals!(op4, sp1, sp1, t1, t4, z_gl_41_14, gl_strategy_1)
 BEAST.momintegrals!(op4, sp1, sp1, t1, t4, z_mrw_41_14, mrw_strategy_1)
 
 @test abs(z_gl_41_14[1,1] - z_mrw_41_14[1,1]) / abs(z_gl_41_14[1,1]) ≈ 0.0 atol=1e-3
 @test abs(z_gl_41_14[1,2] - z_mrw_41_14[1,2]) / abs(z_gl_41_14[1,2]) ≈ 0.0 atol=1e-3
 @test abs(z_gl_41_14[2,1] - z_mrw_41_14[2,1]) / abs(z_gl_41_14[2,1]) ≈ 0.0 atol=1e-3
 @test abs(z_gl_41_14[2,2] - z_mrw_41_14[2,2]) / abs(z_gl_41_14[2,2]) ≈ 0.0 atol=1e-3
 
 #=
 @show z_gl_41_14
 @show z_mrw_41_14
 @show abs(z_gl_41_14[1,1] - z_mrw_41_14[1,1]) / abs(z_gl_41_14[1,1])
 @show abs(z_gl_41_14[1,2] - z_mrw_41_14[1,2]) / abs(z_gl_41_14[1,2])
 @show abs(z_gl_41_14[1,2] - z_mrw_41_14[1,2]) / abs(z_gl_41_14[1,2])
 @show abs(z_gl_41_14[2,1] - z_mrw_41_14[2,1]) / abs(z_gl_41_14[2,1])
=#

# improper allocation of vertices in [t1,t5], so we need to flip them first
#op1
tqd1_15 = BEAST.quadpoints(sp1, [t1,t5], (1000,))
bqd1_15 = BEAST.quadpoints(sp1, [t1,t5], (1001,))

gl_strategy_115 = BEAST.DoubleQuadRule(
	tqd1_15[1,1],
	bqd1_15[1,2],
  )

z_gl_11_15 = zeros(Complex{T}, 2, 2)
z_mrw_11_15 = zeros(Complex{T}, 2, 2)
    
BEAST.momintegrals!(op1, sp1, sp1, t1, t5, z_gl_11_15, gl_strategy_115)
BEAST.momintegrals!(op1, sp1, sp1, t1, t5, z_mrw_11_15, mrw_strategy_1)
  
@test abs(z_gl_11_15[1,1] - z_mrw_11_15[1,1]) / abs(z_gl_11_15[1,1]) ≈ 0.0 atol=1e-3
@test abs(z_gl_11_15[1,2] - z_mrw_11_15[1,2]) / abs(z_gl_11_15[1,2]) ≈ 0.0 atol=1e-3
@test abs(z_gl_11_15[2,1] - z_mrw_11_15[2,1]) / abs(z_gl_11_15[2,1]) ≈ 0.0 atol=1e-3
@test abs(z_gl_11_15[2,2] - z_mrw_11_15[2,2]) / abs(z_gl_11_15[2,2]) ≈ 0.0 atol=1e-3

#op4
z_gl_41_15 = zeros(Complex{T}, 2, 2)
z_mrw_41_15 = zeros(Complex{T}, 2, 2)
 
BEAST.momintegrals!(op4, sp1, sp1, t1, t4, z_gl_41_15, gl_strategy_115)
BEAST.momintegrals!(op4, sp1, sp1, t1, t4, z_mrw_41_15, mrw_strategy_1)

@test abs(z_gl_41_15[1,1] - z_mrw_41_15[1,1]) / abs(z_gl_41_15[1,1]) ≈ 0.0 atol=1e-3
@test abs(z_gl_41_15[1,2] - z_mrw_41_15[1,2]) / abs(z_gl_41_15[1,2]) ≈ 0.0 atol=1e-3
@test abs(z_gl_41_15[2,1] - z_mrw_41_15[2,1]) / abs(z_gl_41_15[2,1]) ≈ 0.0 atol=1e-3
@test abs(z_gl_41_15[2,2] - z_mrw_41_15[2,2]) / abs(z_gl_41_15[2,2]) ≈ 0.0 atol=1e-3

end