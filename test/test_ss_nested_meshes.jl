using Test
# Write a test that compares the momintegrals_nested approach with applying the
# Wilton singularity extraction to a non-conforming mesh.

# Conclusion: WiltonSERule should not be applied to deal with overlapping
# charts in the case of geometrically non-conforming meshes; the result of the
# inner integral is not analytic, crippling the accuracy of the outer integration,
# which is done with triangular Gauss-Legendre points.

using CompScienceMeshes, BEAST
using LinearAlgebra

o = point(0,0,0)
c = point(1/3,1/3,0)
e = point(1/2,0,0)
d = point(1/3,-1/3,0)

# Relative weights chosen to make the two terms equally important
𝒜 = Maxwell3D.singlelayer(wavenumber=0.0, alpha=22.0, beta=1.0)

coarse = Mesh([o,x̂,ŷ], [CompScienceMeshes.index(1,2,3)])
fine = barycentric_refinement(coarse)

X = raviartthomas(coarse, skeleton(coarse,1))
Y = raviartthomas(fine, skeleton(fine,1))

𝒳 = refspace(X)
T = scalartype(𝒜, 𝒳, 𝒳)
@test 𝒳 == refspace(Y)
@test CompScienceMeshes.parent(geometry(Y)) == geometry(X)

p = 6
trial_chart = chart(coarse,1)
test_chart = chart(fine,p)
# @show test_chart.vertices

# test_chart_old = simplex(e,x̂,c)
# trial_chart_old = simplex(o,x̂,ŷ)

test_quadpoints  = BEAST.quadpoints(𝒳, [test_chart],  (12,))[1,1]
trial_quadpoints = BEAST.quadpoints(𝒳, [trial_chart], (13,))[1,1]

function BEAST.quaddata(op::BEAST.MaxwellOperator3D,
    test_local_space::BEAST.RefSpace, trial_local_space::BEAST.RefSpace,
    test_charts, trial_charts)

    a, b = 0.0, 1.0
    tqd = quadpoints(test_local_space, test_charts, (12,12))
    bqd = quadpoints(trial_local_space, trial_charts, (13,13))
    leg = (BEAST._legendre(10,a,b), BEAST._legendre(10,a,b), BEAST._legendre(10,a,b),)
    return (tpoints=tqd, bpoints=bqd, gausslegendre=leg)
end

qs_strat = BEAST.DoubleNumWiltonSauterQStrat(1,1,6,7,10,10,10,10)

# sauterschwab = BEAST.SauterSchwabQuadrature.CommonFace(BEAST._legendre(10,0.0,1.0))
out_ss = zeros(T, numfunctions(𝒳, domain(test_chart)), numfunctions(𝒳, domain(trial_chart)))
BEAST.momintegrals!(out_ss, 𝒜,
    Y, p, test_chart,
    X, 1, trial_chart,
    BEAST.TestRefinesTrialQRule(qs_strat))


wiltonsingext = BEAST.WiltonSERule(test_quadpoints, BEAST.DoubleQuadRule(test_quadpoints, trial_quadpoints))
out_dw = zeros(T, numfunctions(𝒳, domain(test_chart)), numfunctions(𝒳, domain(trial_chart)))
BEAST.momintegrals!(out_dw, 𝒜,
    Y, p, test_chart,
    X, 1, trial_chart,
    wiltonsingext)

@show norm(out_ss-out_dw) / norm(out_dw)

test_chart = simplex(e,x̂,c)
trial_chart = simplex(e,x̂,c)

test_quadpoints  = BEAST.quadpoints(𝒳, [test_chart],  (12,))[1,1]
trial_quadpoints = BEAST.quadpoints(𝒳, [trial_chart], (13,))[1,1]

sauterschwab = BEAST.SauterSchwabQuadrature.CommonFace(BEAST._legendre(10,0.0,1.0))
out_ss1 = zeros(T, numfunctions(𝒳, domain(test_chart)), numfunctions(𝒳, domain(trial_chart)))
BEAST.momintegrals!(𝒜,𝒳,𝒳,test_chart,trial_chart,out_ss1,sauterschwab)

wiltonsingext = BEAST.WiltonSERule(test_quadpoints, BEAST.DoubleQuadRule(test_quadpoints, trial_quadpoints))
out_dw1 = zeros(T, numfunctions(𝒳, domain(test_chart)), numfunctions(𝒳, domain(trial_chart)))
BEAST.momintegrals!(𝒜,𝒳,𝒳,test_chart,trial_chart,out_dw1,wiltonsingext)

@test out_ss1 ≈ out_dw1 rtol=3e-6

# Test accuracy CommonEdge
test_chart = simplex(e,x̂,c)
trial_chart = simplex(c,x̂,(x̂+ŷ)/2)

test_quadpoints  = BEAST.quadpoints(𝒳, [test_chart],  (12,))[1,1]
trial_quadpoints = BEAST.quadpoints(𝒳, [trial_chart], (13,))[1,1]

sauterschwab = BEAST.SauterSchwabQuadrature.CommonEdge(BEAST._legendre(10,0.0,1.0))
out_ss2 = zeros(T, numfunctions(𝒳, domain(test_chart)), numfunctions(𝒳, domain(trial_chart)))
BEAST.momintegrals!(𝒜,𝒳,𝒳,test_chart,trial_chart,out_ss2,sauterschwab)

wiltonsingext = BEAST.WiltonSERule(test_quadpoints, BEAST.DoubleQuadRule(test_quadpoints, trial_quadpoints))
out_dw2 = zeros(T, numfunctions(𝒳, domain(test_chart)), numfunctions(𝒳, domain(trial_chart)))
BEAST.momintegrals!(𝒜,𝒳,𝒳,test_chart,trial_chart,out_dw2,wiltonsingext)

@test out_ss2 ≈ out_dw2 rtol=4e-6


test_chart = simplex(e,x̂,c)
trial_chart = simplex(o,e,c)

test_quadpoints  = BEAST.quadpoints(𝒳, [test_chart],  (12,))[1,1]
trial_quadpoints = BEAST.quadpoints(𝒳, [trial_chart], (13,))[1,1]

sauterschwab = BEAST.SauterSchwabQuadrature.CommonEdge(BEAST._legendre(10,0.0,1.0))
out_ss3 = zeros(T, numfunctions(𝒳, domain(test_chart)), numfunctions(𝒳, domain(trial_chart)))
BEAST.momintegrals!(𝒜,𝒳,𝒳,test_chart,trial_chart,out_ss3,sauterschwab)

wiltonsingext = BEAST.WiltonSERule(test_quadpoints, BEAST.DoubleQuadRule(test_quadpoints, trial_quadpoints))
out_dw3 = zeros(T, numfunctions(𝒳, domain(test_chart)), numfunctions(𝒳, domain(trial_chart)))
BEAST.momintegrals!(𝒜,𝒳,𝒳,test_chart,trial_chart,out_dw3,wiltonsingext)

@test out_ss3 ≈ out_dw3 rtol=3e-6

# Test the CommonVertex cases

test_chart = simplex(e,x̂,c)
trial_chart = simplex(c,(x̂+ŷ)/2,ŷ)

test_quadpoints  = BEAST.quadpoints(𝒳, [test_chart],  (12,))[1,1]
trial_quadpoints = BEAST.quadpoints(𝒳, [trial_chart], (13,))[1,1]

sauterschwab = BEAST.SauterSchwabQuadrature.CommonVertex(BEAST._legendre(10,0.0,1.0))
out_ss4 = zeros(T, numfunctions(𝒳, domain(test_chart)), numfunctions(𝒳, domain(trial_chart)))
BEAST.momintegrals!(𝒜,𝒳,𝒳,test_chart,trial_chart,out_ss4,sauterschwab)

wiltonsingext = BEAST.WiltonSERule(test_quadpoints, BEAST.DoubleQuadRule(test_quadpoints, trial_quadpoints))
out_dw4 = zeros(T, numfunctions(𝒳, domain(test_chart)), numfunctions(𝒳, domain(trial_chart)))
BEAST.momintegrals!(𝒜,𝒳,𝒳,test_chart,trial_chart,out_dw4,wiltonsingext)

@test out_ss4 ≈ out_dw4 rtol=2e-9


test_chart = simplex(e,x̂,c)
trial_chart = simplex(c,ŷ,ŷ/2)

test_quadpoints  = BEAST.quadpoints(𝒳, [test_chart],  (12,))[1,1]
trial_quadpoints = BEAST.quadpoints(𝒳, [trial_chart], (13,))[1,1]

sauterschwab = BEAST.SauterSchwabQuadrature.CommonVertex(BEAST._legendre(10,0.0,1.0))
out_ss5 = zeros(T, numfunctions(𝒳, domain(test_chart)), numfunctions(𝒳, domain(trial_chart)))
BEAST.momintegrals!(𝒜,𝒳,𝒳,test_chart,trial_chart,out_ss5,sauterschwab)

wiltonsingext = BEAST.WiltonSERule(test_quadpoints, BEAST.DoubleQuadRule(test_quadpoints, trial_quadpoints))
out_dw5 = zeros(T, numfunctions(𝒳, domain(test_chart)), numfunctions(𝒳, domain(trial_chart)))
BEAST.momintegrals!(𝒜,𝒳,𝒳,test_chart,trial_chart,out_dw5,wiltonsingext)

@test out_ss4 ≈ out_dw4 rtol=3e-8


test_chart = simplex(e,x̂,c)
trial_chart = simplex(c,ŷ/2,o)

test_quadpoints  = BEAST.quadpoints(𝒳, [test_chart],  (12,))[1,1]
trial_quadpoints = BEAST.quadpoints(𝒳, [trial_chart], (13,))[1,1]

sauterschwab = BEAST.SauterSchwabQuadrature.CommonVertex(BEAST._legendre(10,0.0,1.0))
out_ss6 = zeros(T, numfunctions(𝒳, domain(test_chart)), numfunctions(𝒳, domain(trial_chart)))
BEAST.momintegrals!(𝒜,𝒳,𝒳,test_chart,trial_chart,out_ss6,sauterschwab)

wiltonsingext = BEAST.WiltonSERule(test_quadpoints, BEAST.DoubleQuadRule(test_quadpoints, trial_quadpoints))
out_dw6 = zeros(T, numfunctions(𝒳, domain(test_chart)), numfunctions(𝒳, domain(trial_chart)))
BEAST.momintegrals!(𝒜,𝒳,𝒳,test_chart,trial_chart,out_dw6,wiltonsingext)

@test out_ss4 ≈ out_dw4 rtol=6e-8

out_ss = out_ss1 + out_ss2 + out_ss3 + out_ss4 + out_ss5 + out_ss6
out_dw = out_dw1 + out_dw2 + out_dw3 + out_dw4 + out_dw5 + out_dw6

@test out_ss ≈ out_dw rtol=4e-7
