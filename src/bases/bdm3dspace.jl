mutable struct BDM3DBasis{T,M,P} <: VectorVolumeSpace{T}
    geo::M
    fns::Vector{Vector{Shape{T}}}
    pos::Vector{P}
end

refspace(s::BDM3DBasis{T}) where {T} = BDM3DRefSpace{T}()


function brezzidouglasmarini3d(mesh, faces)
    T = coordtype(mesh)
    P = vertextype(mesh)
    num_faces = numcells(faces)

    C = connectivity(faces, mesh, identity)

    #println(C)
    rows = rowvals(C)
    vals = nonzeros(C)

    #println(rows)

    fns = Vector{Vector{Shape{T}}}(undef, 3*num_faces)
    pos = Vector{P}(undef, 3*num_faces)

    #println(cells(faces))


    for (i,face) in enumerate(faces)

        fns[3*(i-1)+1] = Vector{Shape{T}}()
        fns[3*(i-1)+2] = Vector{Shape{T}}()
        fns[3*(i-1)+3] = Vector{Shape{T}}()

        pos[3*(i-1)+1] = cartesian(center(chart(faces,face)))
        pos[3*(i-1)+2] = cartesian(center(chart(faces,face)))
        pos[3*(i-1)+3] = cartesian(center(chart(faces,face)))

        #println( cartesian(center(chart(faces,face))))
        #println(face)

        # add shape function for each cell to the basis 
        for k in nzrange(C,i)

            #Cell index
            j = rows[k]
            #face index (sign gives local orientation)
            s = vals[k]

            cell = cells(mesh)[j]

            #local face index
            lf = abs(s);
            #local vertex index
            lv1 = mod1(lf+1,4)
            lv2 = mod1(lf+2,4)
            lv3 = mod1(lf+3,4)

            #println("\t",lv1," ",lv2," ",lv3)
            #println("\t=>",cell[lv1]," ",cell[lv2]," ",cell[lv3])
            #global vertex index
            gv1 = cell[lv1]
            gv2 = cell[lv2]
            gv3 = cell[lv3]

            if (gv1 < gv2) && (gv2 < gv3)
                push!(fns[3*(i-1)+1], Shape{T}(j, 3*(abs(s)-1)+1, sign(s)))
                push!(fns[3*(i-1)+2], Shape{T}(j, 3*(abs(s)-1)+2, sign(s)))
                push!(fns[3*(i-1)+3], Shape{T}(j, 3*(abs(s)-1)+3, sign(s)))
                #println("Case 1")
            elseif (gv1 < gv3) && (gv3 < gv2)
                push!(fns[3*(i-1)+1], Shape{T}(j, 3*(abs(s)-1)+1, sign(s)))
                push!(fns[3*(i-1)+2], Shape{T}(j, 3*(abs(s)-1)+3, sign(s)))
                push!(fns[3*(i-1)+3], Shape{T}(j, 3*(abs(s)-1)+2, sign(s))) 
                #println("Case 2")
            elseif (gv2 < gv1) && (gv1 < gv3)
                push!(fns[3*(i-1)+1], Shape{T}(j, 3*(abs(s)-1)+2, sign(s)))
                push!(fns[3*(i-1)+2], Shape{T}(j, 3*(abs(s)-1)+1, sign(s)))
                push!(fns[3*(i-1)+3], Shape{T}(j, 3*(abs(s)-1)+3, sign(s)))
                #println("Case 3")
            elseif (gv2 < gv3) && (gv3 < gv1)
                push!(fns[3*(i-1)+1], Shape{T}(j, 3*(abs(s)-1)+2, sign(s)))
                push!(fns[3*(i-1)+2], Shape{T}(j, 3*(abs(s)-1)+3, sign(s)))
                push!(fns[3*(i-1)+3], Shape{T}(j, 3*(abs(s)-1)+1, sign(s)))
                #println("Case 4")
            elseif (gv3 < gv2) && (gv2 < gv1)
                push!(fns[3*(i-1)+1], Shape{T}(j, 3*(abs(s)-1)+3, sign(s)))
                push!(fns[3*(i-1)+2], Shape{T}(j, 3*(abs(s)-1)+2, sign(s)))
                push!(fns[3*(i-1)+3], Shape{T}(j, 3*(abs(s)-1)+1, sign(s)))
                #println("Case 5")
            elseif (gv3 < gv1) && (gv1 < gv2)
                push!(fns[3*(i-1)+1], Shape{T}(j, 3*(abs(s)-1)+3, sign(s)))
                push!(fns[3*(i-1)+2], Shape{T}(j, 3*(abs(s)-1)+1, sign(s)))
                push!(fns[3*(i-1)+3], Shape{T}(j, 3*(abs(s)-1)+2, sign(s)))
                #println("Case 6")
            end

            
            #println(s)

            #push!(fns[3*(i-1)+1], Shape{T}(j, 3*(abs(s)-1)+1, sign(s)))
            #push!(fns[3*(i-1)+2], Shape{T}(j, 3*(abs(s)-1)+2, sign(s)))
            #push!(fns[3*(i-1)+3], Shape{T}(j, 3*(abs(s)-1)+3, sign(s)))
        end

    end

    BDM3DBasis(mesh, fns, pos)
end 

function brezzidouglasmarini3d(mesh)
    faces = skeleton(mesh, 2)
    brezzidouglasmarini3d(mesh, faces)
end
