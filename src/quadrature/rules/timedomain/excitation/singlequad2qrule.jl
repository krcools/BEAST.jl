# TODO: consolidate with the existing definition of SingleQuadStrategy
mutable struct SingleQuadStrategy2{P} #<: NumQuadStrategy
    quad_points::P
end



function timeintegrals!(z, exc::TDFunctional,
    testrefs, timerefs,
    testpoint, timeelement, dx, qr::SingleQuadStrategy2, f)

    num_tshapes = numfunctions(testrefs, domain(chart(testpoint)))
    for p in qr.quad_points
        t = p.point
        w = p.weight
        U = p.value
        dt = w #* jacobian(t) # * volume(timeelement)

        for i in 1 : num_tshapes
            for k in 1 : numfunctions(timerefs)
                z[i,k] += dot(f[i][1]*U[k], exc(testpoint,t)) * dt * dx
            end
        end
    end
end


function timeintegrals!(z, exc::TDFunctional,
        spacerefs, timerefs::DiracBoundary,
        testpoint, timeelement,
        dx, qr::Nothing, testvals)

        num_tshapes = numfunctions(spacerefs, domain(chart(testpoint)))
        # since timeelement uses barycentric coordinates,
        # the first/left vertex has coords u = 1.0!
        testtime = neighborhood(timeelement, point(0.0))
        @assert cartesian(testtime)[1] ≈ timeelement.vertices[2][1]

        for i in 1 : num_tshapes
            z[i,1] += dot(testvals[i][1], exc(testpoint, testtime)) * dx
        end
end