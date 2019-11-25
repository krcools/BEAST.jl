using CompScienceMeshes
using BEAST

m = meshrectangle(1.0, 1.0, 0.2, 3)
X = raviartthomas(m, skeleton(m,1))

# n = boundary(m)
# x = ntrace(X,n)
#
# V₀ = 1.0
# f = ScalarTrace(p -> V₀)
#
# assemble(f, x)

els, ad = BEAST.assemblydata(X)

Ad = [ maximum(m for (m,a) in ad[c,1]) for c in 1:numcells(m)]
