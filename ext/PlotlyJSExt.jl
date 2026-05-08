module PlotlyJSExt

using BEAST
using PlotlyJS
using CompScienceMeshes
using LinearAlgebra

function PlotlyJS.mesh3d(u::AbstractVector, U::BEAST.AbstractSpace; kwargs...)

    fcr, geo = BEAST.facecurrents(u, U)

    v = CompScienceMeshes.vertexarray(geo)
    c = CompScienceMeshes.cellarray(geo)

    x = v[:,1];    y = v[:,2];    z = v[:,3]
    i = c[:,1].-1; j = c[:,2].-1; k = c[:,3].-1

    return PlotlyJS.mesh3d(;
        x=x, y=y, z=z,
        i=i, j=j, k=k,
        intensitymode="cell",
        intensity=norm.(fcr),
        kwargs...)
end


PlotlyJS.mesh3d(u::BEAST.FEMFunction; kwargs...) = PlotlyJS.mesh3d(u.coeffs, u.space; kwargs...)

end