mutable struct MultiQuadStrategy{P,R} #<: NumQuadStrategy
    quad_points::P
    inner_rule::R
end


timequadrule(qr::MultiQuadStrategy, p) = qr.inner_rule

function momintegrals!(z, exc::TDFunctional, testrefs, timerefs, τ, ρ, qr::MultiQuadStrategy)

    for p in qr.quad_points
        x = p.point
        w = p.weight
        f = p.value
        dx = w

        tqr = timequadrule(qr,p)
        timeintegrals!(z, exc, testrefs, timerefs, x, ρ, dx, tqr, f)

    end

end