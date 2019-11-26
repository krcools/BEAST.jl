using CompScienceMeshes
using BEAST

m = meshsphere(1.0, 0.45)
# m = CompScienceMeshes.rotate(m, rand(3))
X = raviartthomas(m)

Y = buffachristiansen(m)
# Y = X
Y = BEAST.subset(Y,1:135)
geo = geometry(Y)

els, ad = BEAST.assemblydata(Y)
Base.length(it::BEAST.ADIterator) = something(findfirst(x->x[1]<1, it.ad.data[:,it.r,it.c]), size(it.ad.data,1)+1)-1

Ad = [Int[] for c in 1:length(els)]
for c in 1:length(els)
    for r in 1:3
        length(ad[c,r]) == 0 && continue
        push!(Ad[c], maximum(m for (m,a) in ad[c,r]))
    end
end
top = maximum.(Ad)
bot = minimum.(Ad)

using Plots

plot()
plot!(top)
plot!(bot)

import PlotlyJS

t1 = patch(Mesh(geo.mesh.vertices, geo.mesh.faces[1:150]), rand(150))
P = [p[i] for p in Y.pos, i in 1:3]
t2 = PlotlyJS.scatter3d(x=P[1:38,1], y=P[1:38,2], z=P[1:38,3])
t3 = PlotlyJS.scatter3d(x=P[264:264,1], y=P[264:264,2], z=P[264:264,3])
t4 = PlotlyJS.scatter3d(x=P[1:180,1], y=P[1:180,2], z=P[1:180,3])
PlotlyJS.plot(t4)

PlotlyJS.plot([t1,t2,t3])
