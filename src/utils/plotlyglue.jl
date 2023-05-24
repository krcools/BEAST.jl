using Requires

export functionvals

function __init__()
    @require PlotlyJS="f0f68f2c-4968-5e81-91da-67840de0976a" begin     
        @eval function PlotlyJS.cone(xyz::Vector, uvw::Vector;kwargs...)
            PlotlyJS.cone(;
                x=getindex.(xyz,1),
                y=getindex.(xyz,2),
                z=getindex.(xyz,3),
                u=getindex.(uvw,1),
                v=getindex.(uvw,2),
                w=getindex.(uvw,3),
                kwargs...)
        end
    end
end