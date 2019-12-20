struct NDLCCRefSpace{T} <: RefSpace{T,6} end

function valuetype(ref::NDLCCRefSpace{T}, charttype::Type) where {T}
    SVector{universedimension(charttype),T}
end

function (ϕ::NDLCCRefSpace)(ndlc)

    u,v,w = parametric(ndlc)
    j = jacobian(ndlc)

    tu = tangents(ndlc,1)
    tv = tangents(ndlc,2)
    tw = tangents(ndlc,3)

    #B = [tv-tu tw-tu -tu]
    B = [tu tv tw]
    BMT = transpose(inv(B))

    #Choose 1-(0,0,0), 2-(1,0,0), 3-(0,1,0), 4-(0,0,1)
    #Now it is listed as:
    #(1,2) [1-v-w,u,u]  (direction: from 1 to 2)
    #(1,3) [v,1-u-w,v]
    #(1,4) [w,w,1-u,v]
    #(2,3) [-v,u,0]
    #(2,4) [-w,0,u]
    #(3,4) [0,-w,v]

    return SVector((
        (value=(BMT*[-v,u,0])       ,curl=(B*[0,0,2])/j) ,
        (value=(BMT*[-w,0,u])       ,curl=(B*[0,-2,0])/j),
        (value=(BMT*[(w+v-1),-u,-u]),curl=(B*[0,2,-2])/j),
        (value=(BMT*[0,-w,v])       ,curl=(B*[2,0,0])/j) ,
        (value=(BMT*[-v,(u+w-1),-v]),curl=(B*[-2,0,2])/j),
        (value=(BMT*[-w,-w,(u+v-1)]),curl=(B*[2,-2,0])/j)
    ))

##    	return SVector((
##        (value=(BMT*[(w+v-1),-u,-u]),curl=(B*[0,2,-2])/j),
##        (value=(BMT*[-v,(u+w-1),-v])  ,curl=(B*[-2,0,2])/j),
##        (value=(BMT*[-w,-w,(u+v-1)])  ,curl=(B*[2,-2,0])/j),
##        (value=(BMT*[-v,u,0])       ,curl=(B*[0,0,2])/j) ,
##        (value=(BMT*[-w,0,u])       ,curl=(B*[0,-2,0])/j),
##        (value=(BMT*[0,-w,v])       ,curl=(B*[2,0,0])/j)
##        ))
end

#check orientation
function curl(ref::NDLCCRefSpace, sh, el)
    a = [4,2,3,4,1,2]##[2,1,1,1,2,4]#[2,1,4,4,3,2]#[4,2,3,4,1,2]
    b = [3,4,2,1,3,1]##[3,2,4,3,4,3]#[4,3,3,1,2,1]#[3,4,2,1,3,1]
    sh1 = Shape(sh.cellid, b[sh.refid], -sh.coeff)
    sh2 = Shape(sh.cellid, a[sh.refid], sh.coeff)
    return [sh1,sh2]
end
