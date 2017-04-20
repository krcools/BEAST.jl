using Base.Test

import CompScienceMeshes; CM = CompScienceMeshes
import BEAST; BE = BEAST

ω = 1.0
ui, wi = legendre(4, 0.0, 1.0); ui = transpose(ui)
uo, wo = legendre(3, 0.0, 1.0); uo = transpose(uo)
universe = Universe(ω, ui, wi, uo, wo)

fn = Pkg.dir("BEAST", "test", "sphere.in")
mesh = CM.meshfromfile(fn)

rwg = BE.raviartthomas(m);
bc = BE.buffachristiansen(m);

T = BE.MWSingleLayer3D(universe);

trr = BE.assemble(rwg, rwg, t);
tbb = BE.assemble(bc, bc, t);
