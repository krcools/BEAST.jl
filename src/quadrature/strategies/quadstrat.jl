abstract type AbstractQuadStrat end

function (qs::AbstractQuadStrat)(a, X, Y)
    qs
end

function (qs::AbstractQuadStrat)(l, X)
    qs
end