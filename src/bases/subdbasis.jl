import CompScienceMeshes.chart,CompScienceMeshes.subd_chart,CompScienceMeshes.GSubdMesh,CompScienceMeshes.subdMesh


mutable struct subdBasis{T,M,P} <: Space{T}
  geo::M
  fns::Vector{Vector{Shape{T}}}
  pos::Vector{P}
  data::Tuple{Array{subd_chart},Array{Array{Array{Tuple}}}}
end

function subset(s::S,I) where {S<:subdBasis}
    @warn "No parallel assembly supported when using subdivision bases!"
    return s
end


mutable struct subReferenceSpace{T,D} <: RefSpace{T,D} end

refspace(s::subdBasis{T,M,P}) where {T,M,P} = subReferenceSpace{T,12}()
# 12 only for regular case

function (rs::subReferenceSpace)(nbhd)
    # u = parametric(nbhd)
    # x = cartesian(nbhd)
    s = shapefuns(nbhd)
    scurl = get_shape_curl(nbhd)
    # TODO: this needs fixing, stat!
    return s #,scurl
end

function subdsurface(mesh)
    subd_mesh = GSubdMesh(mesh)
    vertices = mesh.vertices
    subd_elements = subd_mesh.elements
    nvertices = length(vertices)
    nelem = length(subd_elements)
    funs=Array{Vector{BEAST.Shape{Float64}}}(undef,nvertices)
    for i = 1 : nvertices funs[i] = [] end
    for ie = 1:nelem
        inodes = subd_elements[ie].RingNodes
        N = length(inodes)
        for ib = 1:N
            nodeid = inodes[ib]
            coeff = 1.0
            ifun = BEAST.Shape{Float64}(ie,ib,coeff)
            push!(funs[nodeid], ifun)
        end
    end
    data = assembly(subd_mesh)
    return subdBasis(mesh,funs,vertices,data)
end

function assemblydata(subdB::subdBasis)
    return subdB.data
end

function assembly(subdG::subdMesh)
    nelem = length(subdG.elements)
    assemblydata = Array{Array{Array{Tuple}}}(undef,nelem)
    ElementCharts = Array{subd_chart}(undef,nelem)
    for e = 1:nelem
        cha = chart(subdG,e)
        ringnodes = cha.RingNodes
        nnodes = length(ringnodes)
        line=Array{Array{Tuple}}(undef,nnodes)
        for inode = 1:nnodes
            Node = ringnodes[inode]
            entries=[]
            entry=tuple(Node,1.0)
            push!(entries,entry)
            line[inode] = entries
        end
        assemblydata[e] = line
        ElementCharts[e] = cha
    end
    return(ElementCharts,assemblydata)
end
