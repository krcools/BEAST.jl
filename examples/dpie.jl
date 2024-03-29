using CompScienceMeshes
using BEAST

Γ = readmesh(joinpath(dirname(pathof(BEAST)),"../examples/sphere2.in"))
X = raviartthomas(Γ)

d = BEAST.CurlSingleLayerDP3D(1.0im, 1.0)
s = Helmholtz3D.singlelayer(wavenumber=1.0)

X = raviartthomas(Γ,BEAST.Continuity{:none})
Y = lagrangec0d1(Γ)

x = refspace(X)
y = refspace(Y)
charts = [chart(Γ,c) for c in Γ]

qs = BEAST.defaultquadstrat(s, x, x)
qd = quaddata(s, x, x, charts, charts, qs)
m1,m2 = 10,11
ql = quadrule(s, x, x, m1, charts[m1], m2, charts[m2], qd, qs)
@show @which quadrule(s, x, x, m1, charts[m1], m2, charts[m2], qd, qs)
@show typeof(ql)
BEAST.momintegrals!(s, x, x, charts[10], charts[20],
    zeros(ComplexF64,3,3), ql)

# sop = BEAST.singularpart(s)
ctr = CompScienceMeshes.center(charts[m1])
# BEAST.momintegrals!(sop, ctr, x, x, charts[m1], charts[m1], zeros(3,3), ql, 1.0)

dvg = divergence

@hilbertspace j p
@hilbertspace k q

κ = 1.0
curlA = BEAST.ScalarTrace{ComplexF64}(((x,y,z),)->-exp(-im*κ*z)*ŷ)
ndotA = BEAST.NDotTrace{ComplexF64}(((x,y,z),)->-x*exp(-im*κ*z)*ẑ)
a = @varform s[k,j] + s[dvg(k),p] +
    s[q,dvg(j)] - κ^2*s[q,p] == curlA[k]+ndotA[q]

Eq = @discretise a k∈X q∈Y j∈X p∈Y

u = solve(Eq)

# lform = Eq.equation.rhs
# @which scalartype(lform.terms[2].functional.field)
# @enter assemble(Eq.equation.rhs, Eq.test_space_dict)
# assemble(Eq.equation.lhs, Eq.test_space_dict, Eq.trial_space_dict)

# uj = u[1:numfunctions(X)]
# up = u[numfunctions(X)+1:end]

uj = u[j]
up = u[p]

xrange = range(-2.0, stop=2.0, length=40)
# xrange = [0.0]
# yrange = [0.0]
zrange = range(-2.0, stop=2.0, length=40)
pts = [point(x,0,z) for x ∈ xrange, z ∈ zrange]

A1 = BEAST.DPVectorialTerm(1.0, 1.0*im)
A2 = BEAST.DPScalarTerm(1.0, 1.0*im)

p1 = potential(A1,pts,uj,X)
p2 = potential(A2,pts,up,Y)

p = p1 + p2
