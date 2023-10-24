import Base.length

struct NormalChart{U,D,C,N,T}
    chart::CompScienceMeshes.Simplex{U,D,C,N,T}
    normals::SVector{C,SVector{U,T}}
end

normals(chart::NormalChart) = chart.normals
chart(chart::NormalChart) = chart.chart
dimtype(chart::NormalChart) = dimtype(chart(chart))
coordtype(chart::NormalChart) = coordtype(chart(chart))
volume(chart::NormalChart) = volume(chart(chart))
dimension(chart::NormalChart) = dimension(chart(chart))
Base.length(chart::NormalChart) = length(chart(chart))
universedimension(chart::NormalChart) = universedimension(chart(chart))
getindex(p::NormalChart,I::Union{Number,SVector,Array}) = chart(p)[I]


function normalchart(chart::Simplex)
    normals = normals(chart)
    NormalChart(chart,normals)
end

normalchart(chart::Simplex,normals) = NormalChart(chart,normals)