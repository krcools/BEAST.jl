# import Base.length
# import CompScienceMeshes: chart,dimtype,coordtype,volume,dimension,universedimension
# struct NormalChart{U,V}
#     chart::U
#     normals::V
# end

# normals(chart::NormalChart) = chart.normals
# chart(c::NormalChart) = c.chart
# dimtype(c::NormalChart) = dimtype(chart(c))
# coordtype(c::NormalChart) = coordtype(chart(chart))

# volume(c::NormalChart) = volume(chart(chart))
# dimension(chart::NormalChart) = dimension(chart(chart))
# Base.length(chart::NormalChart) = length(chart(chart))
# universedimension(chart::NormalChart) = universedimension(chart(chart))
# getindex(p::NormalChart,I::Union{Number,SVector,Array}) = chart(p)[I]


# function normalchart(chart::CompScienceMeshes.Simplex)
#     normals = chart.normals
#     NormalChart(chart,normals)
# end

# normalchart(chart::CompScienceMeshes.Simplex,normals) = NormalChart(chart,normals)

#CompScienceMeshes.normal(t::CompScienceMeshes.Simplex{3,2,1,3,<:Number}) = t.normals[1]


# function permute_barycentric(perm,bary::Tuple{T,T}) where {T}
#     last_coef = 1-sum(bary)
#     total_bary = SVector{3,T}([bary[1],bary[2],last_coef])
#     return Tuple{T,T}(total_bary[perm][1:end-1])
# end
# function permute_barycentric(perm,bary::Tuple{T,T,T}) where {T}
#     last_coef = 1-sum(bary)
#     total_bary = SVector{4,T}([bary[1],bary[2],bary[3],last_coef])
#     return Tuple{T,T,T}(total_bary[perm][1:end-1])
# end