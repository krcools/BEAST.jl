using CompScienceMeshes
using BEAST

Γ1 = readmesh(joinpath(dirname(pathof(BEAST)),"../examples/sphere2.in"))
Γ2 = CompScienceMeshes.translate(Γ1, [10,0,0])

SL = Maxwell3D.singlelayer(wavenumber=1.0)
X1 = raviartthomas(Γ1)
X2 = raviartthomas(Γ2)

x1 = refspace(X1)
x2 = refspace(X2)

function BEAST.quaddata(op::typeof(SL), tref::typeof(x1), bref::typeof(x2),
    tels, bels, qs::BEAST.DoubleNumQStrat)

    qs = BEAST.DoubleNumWiltonSauterQStrat(qs.outer_rule, qs.inner_rule, 1, 1, 1, 1, 1, 1)
    BEAST.quaddata(op, tref, bref, tels, bels, qs)
end

function BEAST.quadrule(op::typeof(SL), tref::typeof(x1), bref::typeof(x2),
    i ,τ, j, σ, qd, qs::BEAST.DoubleNumQStrat)

    return BEAST.DoubleQuadRule(
        qd.tpoints[1,i],
        qd.bpoints[1,j])
end

BEAST.quadinfo(SL,X1,X2)
BEAST.quadinfo(SL,X1,X2, quadstrat=BEAST.DoubleNumQStrat(2,3))

@time Z1 = assemble(SL,X1,X2);
@time Z2 = assemble(SL,X1,X2,quadstrat=BEAST.DoubleNumQStrat(2,3));

ba1 = blockassembler(SL,X1,X2)
ba2 = blockassembler(SL,X1,X2,quadstrat=BEAST.DoubleNumQStrat(2,3))

T = scalartype(SL,X1,X2)
n1 = numfunctions(X1)
n2 = numfunctions(X2)

W1 = zeros(T,n1,n2);
W2 = zeros(T,n1,n2);

idcs1 = collect(1:n1)
idcs2 = collect(1:n2)
@time ba1(idcs1,idcs2, (v,m,n)->(W1[m,n]+=v))
@time ba2(idcs1,idcs2, (v,m,n)->(W2[m,n]+=v));