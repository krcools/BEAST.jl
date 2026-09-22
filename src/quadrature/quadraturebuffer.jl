quadraturebuffer(quadstrat, test_space, trial_space) = (;)

const SauterSchwabBuffer = NamedTuple{(:I, :J, :K, :L)}
const SauterSchwabBufferSet = NamedTuple{(:edge, :triangle, :quadrilateral)}

function quadraturebuffer(qs::Union{
    TestRefinesTrialQStrat,
    TrialRefinesTestQStrat,
    NonConformingIntegralOpQStrat,
    CommonFaceOverlappingEdgeQStrat,
}, test_space, trial_space)
    return quadraturebuffer(qs.conforming_qstrat, test_space, trial_space)
end

function quadraturebuffer(
    qs::NonConfTestBaryRefOfTrialQStrat, test_space, trial_space
)
    return quadraturebuffer(qs.conforming_qstrat, test_space, trial_space)
end

function _sauterschwab_buffer(n)
    return (;
        I=Vector{Int}(undef, n),
        J=Vector{Int}(undef, n),
        K=Vector{Int}(undef, n),
        L=Vector{Int}(undef, n),
    )
end

function quadraturebuffer(::Union{
    DoubleNumSauterQstrat,
    DoubleNumWiltonSauterQStrat,
    SelfSauterOtherwiseDNumQStrat,
    CommonFaceVertexSauterCommonEdgeWiltonPostitiveDistanceNumQStrat,
}, test_space, trial_space)
    return (;
        edge=_sauterschwab_buffer(2),
        triangle=_sauterschwab_buffer(3),
        quadrilateral=_sauterschwab_buffer(4),
    )
end

function quadraturebuffer(::Union{
    SauterSchwabQuadrature.CommonVertex,
    SauterSchwabQuadrature.CommonEdge,
    SauterSchwabQuadrature.CommonFace,
}, test_space, trial_space)
    return _sauterschwab_buffer(3)
end

function quadraturebuffer(::Union{
    SauterSchwabQuadrature.CommonVertexQuad,
    SauterSchwabQuadrature.CommonEdgeQuad,
    SauterSchwabQuadrature.CommonFaceQuad,
}, test_space, trial_space)
    return _sauterschwab_buffer(4)
end
