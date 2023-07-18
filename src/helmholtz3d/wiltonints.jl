# single layer
function (igd::Integrand{<:HH3DSingleLayerReg})(x,y,f,g)
    α = igd.operator.alpha
    γ = igd.operator.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1 / R
    green = (expm1(-γ*R) #= + γ*R =# - 0.5*γ^2*R^2) / (4pi*R)
    αG = α * green

    _integrands(f,g) do fi, gi
        dot(gi.value, αG*fi.value)
    end
end

function innerintegrals!(op::HH3DSingleLayerSng, test_neighborhood,
    test_refspace::LagrangeRefSpace{T,0} where {T},
    trial_refspace::LagrangeRefSpace{T,0} where {T},
    test_elements, trial_element, zlocal, quadrature_rule::WiltonSERule, dx)

γ = op.gamma
α = op.alpha

s1, s2, s3 = trial_element.vertices

x = cartesian(test_neighborhood)
n = normalize((s1-s3)×(s2-s3))
ρ = x - dot(x - s1, n) * n

scal, vec = WiltonInts84.wiltonints(s1, s2, s3, x, Val{1})
∫G = (scal[2]#=  - γ*scal[3] =# + 0.5*γ^2*scal[4]) / (4π)

zlocal[1,1] += α * ∫G * dx
return nothing
end



# double layer transposed
function (igd::Integrand{<:HH3DDoubleLayerTransposedReg})(x,y,f,g)
    γ = igd.operator.gamma
    
    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    γR = γ*R
    expo = exp(-γR)

    gradgreen = ( -(γR + 1)*expo + (1 - 0.5*γR^2) ) * (i4pi*iR^3) * r

    n = normal(x)
    fvalue = getvalue(f)
    gvalue = getvalue(g)

    return _krondot(fvalue,gvalue) * dot(n, gradgreen)
end

function innerintegrals!(op::HH3DDoubleLayerTransposedSng, test_neighborhood,
    test_refspace::LagrangeRefSpace{T,1,3,3} where {T},
    trial_refspace::LagrangeRefSpace{T,0,3,1} where {T},
    test_elements, trial_element, zlocal, quadrature_rule::WiltonSERule, dx)
    γ = op.gamma
    α = op.alpha

    s1, s2, s3 = trial_element.vertices
    t1, t2, t3 = test_elements.vertices
    x = cartesian(test_neighborhood)
    n = normalize((t1-t3)×(t2-t3))
    ρ = x - dot(x - s1, n) * n

    scal, vec, grad = WiltonInts84.wiltonints(s1, s2, s3, x, Val{1})

    ∫∇G = -(grad[1]+0.5*γ^2*grad[3])/(4π)
    ∫n∇G = dot(n,∫∇G)
    Atot = 1/2*norm(cross(t3-t1,t3-t2))
    for i in 1:numfunctions(test_refspace)
        Ai = 1/2*norm(cross(test_elements.vertices[mod1(i-1,3)]-x,test_elements.vertices[mod1(i+1,3)]-x))
        g = Ai/Atot
        for j in 1:numfunctions(trial_refspace)
            zlocal[i,j] += α * ∫n∇G * g * dx
        end
    end

    return nothing
end

# double layer
function (igd::Integrand{<:HH3DDoubleLayerReg})(x,y,f,g)
    γ = igd.operator.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    γR = γ*R
    expo = exp(-γR)

    gradgreen = ( -(γR + 1)*expo + (1 - 0.5*γR^2) ) * (i4pi*iR^3) * r

    n = normal(y)
    fvalue = getvalue(f)
    gvalue = getvalue(g)

    return _krondot(fvalue,gvalue) * dot(n, -gradgreen)

end

function innerintegrals!(op::HH3DDoubleLayerSng, p,
    g::LagrangeRefSpace{T,0} where {T},
    f::LagrangeRefSpace{T,1} where {T},
    t, s, z, quadrature_rule::WiltonSERule, dx)
    γ = op.gamma
    α = op.alpha

    s1, s2, s3 = s.vertices

    x = cartesian(p)
    n = normalize((s1-s3)×(s2-s3))
    ρ = x - dot(x - s1, n) * n

    _, _, _, grad = WiltonInts84.higherorder(s1,s2,s3,x,3)

    ∫∇G = (-grad[1] - 0.5*γ^2*grad[2]) / (4π)

  
for i in 1:numfunctions(g)
    for j in 1:numfunctions(f)     
       z[i,j] += α * dot(n,∫∇G[j]) * dx
    end
end

    return nothing
end


function (igd::Integrand{<:HH3DHyperSingularReg})(x,y,f,g)
    α = igd.operator.alpha
    β = igd.operator.beta
    γ = igd.operator.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1 / R
    green = (expm1(-γ*R) - 0.5*γ^2*R^2) / (4pi*R)
    nx = normal(x)
    ny = normal(y)

    _integrands(f,g) do fi, gi
        α*dot(nx,ny)*gi.value*fi.value*green + β*dot(gi.curl,fi.curl)*green
    end
end

function innerintegrals!(op::HH3DHyperSingularSng, p,
    g::LagrangeRefSpace{T,1} where {T},
    f::LagrangeRefSpace{T,1} where {T},
    t, s, z, quadrature_rule::WiltonSERule, dx)
    α = op.alpha
    β = op.beta
    γ = op.gamma
    s1, s2, s3 = s.vertices
    t1, t2, t3 = t.vertices
    x = cartesian(p)
    nx = normalize((s1-s3)×(s2-s3))
    ny = normalize((t1-t3)×(t2-t3))

    ∫Rⁿ, ∫RⁿN = WiltonInts84.higherorder(s1,s2,s3,x,3)
    greenconst = (∫Rⁿ[2] + 0.5*γ^2*∫Rⁿ[3]) / (4π)
    greenlinear = (∫RⁿN[2] + 0.5*γ^2*∫RⁿN[3] ) / (4π)

    jt = volume(t) * factorial(dimension(t))
    js = volume(s) * factorial(dimension(s))
    curlt = [(t3-t2)/jt,(t1-t3)/jt,(t2-t1)/jt]
    curls = [(s3-s2)/js,(s1-s3)/js,(s2-s1)/js]
    Atot = 1/2*norm(cross(t3-t1,t3-t2))
    for i in 1:numfunctions(g)
        Ai = 1/2*norm(cross(t.vertices[mod1(i-1,3)]-x,t.vertices[mod1(i+1,3)]-x))
        g = Ai/Atot
        for j in 1:numfunctions(f)
           z[i,j] += β * dot(curlt[i],curls[j])*greenconst*dx + α*dot(nx,ny) * greenlinear[j]*g*dx
        end
    end
end
