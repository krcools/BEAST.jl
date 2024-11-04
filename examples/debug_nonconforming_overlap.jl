using BEAST, CompScienceMeshes
import Plotly

mesh_h = 0.2/5
h_ref = 0.2/10
m = meshcuboid(1.0,1.0,0.25,mesh_h)
mref = meshcuboid(1.0,1.0,0.25,h_ref)

X = raviartthomas(m)

m1 = CompScienceMeshes.Mesh(m.vertices, [m.faces[120]])
m2 = CompScienceMeshes.Mesh(mref.vertices, [mref.faces[452]])

wfm1 = wireframe(m1;color="rgb(100,0,0)")
wfm2 = wireframe(m2;color="rgb(0,0,100)")
Plotly.plot([wfm1, wfm2])

basis_chart = CompScienceMeshes.simplex(m1.vertices[m1.faces[1]])
test_chart = CompScienceMeshes.simplex(m2.vertices[m2.faces[1]])

test_charts, tclps = CompScienceMeshes.intersection_keep_clippings(test_chart, basis_chart)
_, bclps = CompScienceMeshes.intersection_keep_clippings(basis_chart, test_chart)
bsis_charts = copy(test_charts)

for tclp in tclps append!(test_charts, tclp) end
for bclp in bclps append!(bsis_charts, bclp) end

T = coordtype(test_chart)
h = max(volume(test_chart), volume(basis_chart))
test_charts = [ch for ch in test_charts if volume(ch) .> 1e6 * eps(T) * h]
bsis_charts = [ch for ch in bsis_charts if volume(ch) .> 1e6 * eps(T) * h]

test_charts = [ch for ch in test_charts if volume(ch) .> 1e6 * eps(T)]
bsis_charts = [ch for ch in bsis_charts if volume(ch) .> 1e6 * eps(T)]

zlocal = zeros(ComplexF64,3,3)
dqstrat = BEAST.defaultquadstrat(t, X, X)
qstrat = BEAST.NonConformingIntegralOpQStrat(dqstrat)
innerqstrat = BEAST.CommonFaceOverlappingEdgeQStrat(dqstrat)
qdata = BEAST.quaddata(t, refspace(X), refspace(X), test_charts, bsis_charts, innerqstrat)
qrule = BEAST.quadrule(t, refspace(X), refspace(X), 5, test_charts[5], 4, bsis_charts[4], qdata, innerqstrat)
BEAST.momintegrals!(t, refspace(X), refspace(X), test_charts[3], bsis_charts[3], zlocal, qrule, true)

wft = wireframe(test_charts;color="rgb(100,100,0)")
wfb = Plotly.GenericTrace[]
for bsis in bsis_charts2 
    push!(wfb, wireframe([bsis];color="rgb(0,100,100)",width=3)) 
end

push!(wfb, wfm1)
push!(wfb, wfm2)
Plotly.plot(wfb)

trtest = wireframe([test_charts[3]], color=:red, width=3)
trbsis = wireframe([bsis_charts2[3]], color=:blue, width=3)
Plotly.plot([trtest, trbsis, wfm1, wfm2])