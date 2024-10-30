
mutable struct ND2RefSpace{T} <: RefSpace{T} end

function (ϕ::ND2RefSpace)(nbd)

    u, v = parametric(nbd)
    n = normal(nbd)
    j = jacobian(nbd)

    tu = tangents(nbd,1)
    tv = tangents(nbd,2)

    d = 2/j

    inv_j = 1/j

    return SVector((
        (value= n × ((8*u^2+8*u*v-12*u-6*v+4)*tu + (2*v*(4*u+4*v-3))*tv)*inv_j, curl=(24*u+24*v-18)*inv_j),
        (value= n × ((-8*u*v+2*u+6*v-2)*tu + (4*v*(-2*v+1))*tv)*inv_j, curl=(6-24*v)*inv_j),
        (value= n × ((4*u*(1-2*u))*tu + (-8*u*v+6*u+2*v-2)*tv)*inv_j, curl=(6-24*u)*inv_j),
        (value= n × ((2*u*(4*u+4*v-3))*tu + (8*u*v-6*u+8*v^2-12*v+4)*tv)*inv_j, curl=(24*u+24*v-18)*inv_j),
        (value= n × ((2u*(1-4v))*tu + (4v*(1-2v))*tv)*inv_j, curl=(6-24*v)*inv_j),
        (value= n × ((4u*(1-2u))*tu + (2v*(1-4u))*tv)*inv_j, curl=(6-24*u)*inv_j),
        (value= n × ((8*u*(-2*u-v+2))*tu + (8*v*(-2*u-v+1))*tv)*inv_j, curl=(-48*u-24*v+24)*inv_j),
        (value= n × ((8*u*(-u-2*v+1))*tu + (8*v*(-u-2*v+2))*tv)*inv_j, curl=(-24*u-48*v+24)*inv_j)
    ))

end

numfunctions(x::ND2RefSpace, dom::CompScienceMeshes.ReferenceSimplex{2}) = 8