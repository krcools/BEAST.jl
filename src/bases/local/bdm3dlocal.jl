struct BDM3DRefSpace{T} <: RefSpace{T} end

numfunctions(f::BDM3DRefSpace, ch::CompScienceMeshes.ReferenceSimplex{3}) = 12

function (f::BDM3DRefSpace)(p)

    u,v,w = parametric(p)

    tu = tangents(p,1)
    tv = tangents(p,2)
    tw = tangents(p,3)

    j = jacobian(p)
    d=1/j

    return SVector((

      
        (value= 2*(-v*tu+v*tv)/j , divergence= 2*d),
        (value= 2*(-w*tu+w*tw)/j , divergence= 2*d),
        (value= 2*((u+v+w-1)*tu)/j , divergence= 2*d),

        (value= 2*(-w*tv+w*tw)/j , divergence= 2*d),
        (value= 2*((u+v+w-1)*tv)/j , divergence= 2*d),
        (value= 2*(u*tu-u*tv)/j , divergence= 2*d),


        (value= 2*((u+v+w-1)*tw)/j , divergence= 2*d),
        (value= 2*(u*tu-u*tw)/j , divergence= 2*d),
        (value= 2*(v*tv-v*tw)/j  , divergence= 2*d),

        (value= 2*(u*tu)/j      , divergence= 2*d),      
        (value= 2*(v*tv)/j       , divergence= 2*d),       
        (value= 2*(w*tw)/j       , divergence= 2*d)

    ))
end