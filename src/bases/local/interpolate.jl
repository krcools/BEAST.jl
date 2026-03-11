"""
    interpolate(interpolant::RefSpace, chart1, interpolee::RefSpace, chart2)

Computes by interpolation approximations of the local shape functions for
`interpolee` on `chart2` in terms of the local shape functions for `interpolant`
on `chart1`. The returned value is a matrix `Q` such that

```math
\\phi_i \\approx \\sum_j Q_{ij} \\psi_j
```

with ``\\phi_i`` the i-th local shape function for `interpolee` and ``\\psi_j`` the
j-th local shape function for `interpolant`.
"""
function interpolate(interpolant::RefSpace, chart1, interpolee::RefSpace, chart2)
    function fields(p)
        x = cartesian(p)
        v = carttobary(chart2, x)
        r = neighborhood(chart2, v)
        fieldvals = [f.value for f in interpolee(r)]
    end

    interpolate(fields, interpolant, chart1)
end

function interpolate!(out, interpolant::RefSpace, chart1, interpolee::RefSpace, chart2)
    function fields(p)
        x = cartesian(p)
        v = carttobary(chart2, x)
        r = neighborhood(chart2, v)
        fieldvals = [f.value for f in interpolee(r)]
    end

    interpolate!(out, fields, interpolant, chart1)
end


function interpolate(interpolant::RefSpace, chart1, interpolee::RefSpace, chart2, ch1toch2)
    function fields(p1)
        u1 = parametric(p1)
        u2 = cartesian(ch1toch2, u1)
        p2 = neighborhood(chart2, u2)
        fieldvals = [f.value for f in interpolee(p2)]
    end

    interpolate(fields, interpolant, chart1)
end

function interpolate!(out, interpolant::RefSpace, chart1, interpolee::RefSpace, chart2, ch1toch2)
    function fields(p1)
        u1 = parametric(p1)
        u2 = cartesian(ch1toch2, u1)
        p2 = neighborhood(chart2, u2)
        fieldvals = [f.value for f in interpolee(p2)]
    end

    interpolate!(out, fields, interpolant, chart1)
end


function restrict(ϕ::RefSpace, dom1, dom2)
    interpolate(ϕ, dom2, ϕ, dom1)
end

function restrict!(out, ϕ::RefSpace, dom1, dom2)
    interpolate!(out, ϕ, dom2, ϕ, dom1)
end

function restrict(ϕ::RefSpace, dom1, dom2, dom2todom1)
    interpolate(ϕ, dom2, ϕ, dom1, dom2todom1)
end

function restrict!(out, ϕ::RefSpace, dom1, dom2, dom2todom1)
    interpolate!(out, ϕ, dom2, ϕ, dom1, dom2todom1)
end