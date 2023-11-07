using CompScienceMeshes, BEAST, LinearAlgebra
r = 1.0
h = 0.1
Γ = CompScienceMeshes.meshsphere(radius=r, h=h)
X = raviartthomas(Γ)
X1 = BEAST.nedelec(Γ)
Δt0 = 0.07376
reduce_tstep = 1
Δt = Δt0/reduce_tstep
Nt = Int(floor(100*reduce_tstep))
T = timebasisshiftedlagrange(Δt, Nt, 3)
U = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = X ⊗ U
V1 = X1 ⊗ T
W1 = X1 ⊗ U

duration = 20 * Δt0 * 2
delay = 1.0 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)
direction, polarisation = ŷ, x̂
E = planewave(polarisation, direction, derive(gaussian), 1.0)
chr = BEAST.Polynomial(0.681,-1.620,1.200)
SL = TDMaxwell3D.singlelayer(speedoflight=1.0, numdiffs=1)
idST = Identity()⊗Identity()
σ = BEAST.conductivityfunc(chr,numdiffs=1)
@hilbertspace j
@hilbertspace k
@hilbertspace l
@hilbertspace m

tdefie = @discretise 1.0SL[k,j] == -1.0E[k] j∈V k∈W
ohmslaw = @discretise 1.0idST[m,l] == 1.0σ[m] l∈V1 m∈W1

xj, xe, xj_all, xe_all = BEAST.td_solve(tdefie, ohmslaw)

#post processing

fcr, geo = facecurrents(xj, V)
fef, geo = facecurrents(xe, V1)

xaxs = 0:0.01:1

using Plots
Plots.plotly()
p = Plots.plot(xaxs, xaxs.*chr.(xaxs))

face = 1300
display(Plots.plot!(norm.(fef[face, :]), norm.(fcr[face, :]), marker = '.', markersize=1))