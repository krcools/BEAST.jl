# single layer
function (igd::Integrand{<:HH3DSingleLayerReg})(x,y,f,g)  
    # TODO: It is better to dispatch on γ
    # to handle the static case more efficiently
    γ = gamma(igd.operator)
    α = igd.operator.alpha
    
    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1 / R
    green = (expm1(-γ*R) - 0.5*γ^2*R^2) / (4pi*R)
    αG = α * green

    _integrands(f,g) do fi, gi
        dot(gi.value, αG*fi.value)
    end
end

function innerintegrals!(op::HH3DSingleLayerSng, test_neighborhood,
    test_refspace::LagrangeRefSpace{T,0} where {T},
    trial_refspace::LagrangeRefSpace{T,0} where {T},
    test_elements, trial_element, zlocal, quadrature_rule::WiltonSERule, dx)
    # TODO: It is better to dispatch on γ
    # to handle the static case more efficiently
    γ = gamma(op)
    α = op.alpha
    
    s1, s2, s3 = trial_element.vertices

    x = cartesian(test_neighborhood)
    n = normalize((s1-s3)×(s2-s3))
    ρ = x - dot(x - s1, n) * n

    scal, vec = WiltonInts84.wiltonints(s1, s2, s3, x, Val{1})
    ∫G = (scal[2] + 0.5*γ^2*scal[4]) / (4π)

    zlocal[1,1] += α * ∫G * dx

    return nothing
end

# single layer with patch basis and pyramid testing
function innerintegrals!(op::HH3DSingleLayerSng, test_neighborhood,
    test_refspace::LagrangeRefSpace{T,1} where {T},
    trial_refspace::LagrangeRefSpace{T,0} where {T},
    test_elements, trial_element, zlocal, quadrature_rule::WiltonSERule, dx)
    # TODO: It is better to dispatch on γ
    # to handle the static case more efficiently
    γ = gamma(op)
    α = op.alpha

    s1, s2, s3 = trial_element.vertices
    t1, t2, t3 = test_elements.vertices

    num_tshapes = numfunctions(test_refspace, domain(test_elements))
    num_bshapes = numfunctions(trial_refspace, domain(trial_element))

    x = cartesian(test_neighborhood)
    n = normalize((s1-s3)×(s2-s3))
    ρ = x - dot(x - s1, n) * n

    scal, vec = WiltonInts84.wiltonints(s1, s2, s3, x, Val{1})
    ∫G = (scal[2] + 0.5*γ^2*scal[4]) / (4π)
    Atot = 1/2*norm(cross(t3-t1,t3-t2))
    for i in 1:num_tshapes
        Ai = 1/2*norm(cross(test_elements.vertices[mod1(i-1,3)]-x,test_elements.vertices[mod1(i+1,3)]-x))
        g = Ai/Atot
        for j in 1:num_bshapes
            zlocal[i,j] += α * ∫G * g * dx
        end
    end

    return nothing
end

# single layer with pyramid basis and patch testing
function innerintegrals!(op::HH3DSingleLayerSng, test_neighborhood,
    test_refspace::LagrangeRefSpace{T,0} where {T},
    trial_refspace::LagrangeRefSpace{T,1} where {T},
    test_elements, trial_element, zlocal, quadrature_rule::WiltonSERule, dx)
    # TODO: It is better to dispatch on γ
    # to handle the static case more efficiently
    γ = gamma(op)
    α = op.alpha

    num_tshapes = numfunctions(test_refspace, domain(test_elements))
    num_bshapes = numfunctions(trial_refspace, domain(trial_element))

    s1, s2, s3 = trial_element.vertices

    x = cartesian(test_neighborhood)
    n = normalize((s1-s3)×(s2-s3))
    ρ = x - dot(x - s1, n) * n

    ∫Rⁿ, ∫RⁿN = WiltonInts84.higherorder(s1,s2,s3,x,3)
    ∫G = (∫RⁿN[2] + 0.5*γ^2*∫RⁿN[3]) / (4π)

    for i in 1:num_bshapes
    zlocal[1,i] += α * ∫G[i] * dx
    end

    return nothing
end

#single layer with pyramid basis and pyramid testing
function innerintegrals!(op::HH3DSingleLayerSng, test_neighborhood,
    test_refspace::LagrangeRefSpace{T,1} where {T},
    trial_refspace::LagrangeRefSpace{T,1} where {T},
    test_elements, trial_element, zlocal, quadrature_rule::WiltonSERule, dx)
    # TODO: It is better to dispatch on γ
    # to handle the static case more efficiently    
    γ = gamma(op)
    α = op.alpha

    s1, s2, s3 = trial_element.vertices
    t1, t2, t3 = test_elements.vertices

    num_tshapes = numfunctions(test_refspace, domain(test_elements))
    num_bshapes = numfunctions(trial_refspace, domain(trial_element))

    x = cartesian(test_neighborhood)
    n = normalize((s1-s3)×(s2-s3))
    ρ = x - dot(x - s1, n) * n

    ∫Rⁿ, ∫RⁿN = WiltonInts84.higherorder(s1,s2,s3,x,3)
    ∫G = (∫RⁿN[2] + 0.5*γ^2*∫RⁿN[3]) / (4π)

    Atot = 1/2*norm(cross(t3-t1,t3-t2))
    for i in 1:num_tshapes
        Ai = 1/2*norm(cross(test_elements.vertices[mod1(i-1,3)]-x,test_elements.vertices[mod1(i+1,3)]-x))
        g = Ai/Atot
        for j in 1:num_bshapes
            zlocal[i,j] += α * ∫G[j] * g * dx
        end
    end

    return nothing
end

# double layer transposed
function (igd::Integrand{<:HH3DDoubleLayerTransposedReg})(x,y,f,g)
    # TODO: It is better to dispatch on γ
    # to handle the static case more efficiently
    γ = gamma(igd.operator)
    
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

#double layer transposed with patch basis and pyramid testing
function innerintegrals!(op::HH3DDoubleLayerTransposedSng, test_neighborhood,
    test_refspace::LagrangeRefSpace{T,0,3,1} where {T},
    trial_refspace::LagrangeRefSpace{T,0,3,1} where {T},
    test_elements, trial_element, zlocal, quadrature_rule::WiltonSERule, dx)
    # TODO: It is better to dispatch on γ
    # to handle the static case more efficiently    
    γ = gamma(op)
    α = op.alpha

    s1, s2, s3 = trial_element.vertices
    t1, t2, t3 = test_elements.vertices
    x = cartesian(test_neighborhood)
    n = normalize((t1-t3)×(t2-t3))
    ρ = x - dot(x - s1, n) * n

    num_tshapes = numfunctions(test_refspace, domain(test_elements))
    num_bshapes = numfunctions(trial_refspace, domain(trial_element))

    scal, vec, grad = WiltonInts84.wiltonints(s1, s2, s3, x, Val{1})

    ∫∇G = -(grad[1]+0.5*γ^2*grad[3])/(4π)
    ∫n∇G = dot(n,∫∇G)
    for i in 1:num_tshapes
        for j in 1:num_bshapes
            zlocal[i,j] += α * ∫n∇G * dx
        end
    end

    return nothing
end
#double layer transposed with patch basis and pyramid testing
function innerintegrals!(op::HH3DDoubleLayerTransposedSng, test_neighborhood,
    test_refspace::LagrangeRefSpace{T,1,3,3} where {T},
    trial_refspace::LagrangeRefSpace{T,0,3,1} where {T},
    test_elements, trial_element, zlocal, quadrature_rule::WiltonSERule, dx)
    # TODO: It is better to dispatch on γ
    # to handle the static case more efficiently    
    γ = gamma(op)
    α = op.alpha

    s1, s2, s3 = trial_element.vertices
    t1, t2, t3 = test_elements.vertices
    x = cartesian(test_neighborhood)
    n = normalize((t1-t3)×(t2-t3))
    ρ = x - dot(x - s1, n) * n

    num_tshapes = numfunctions(test_refspace, domain(test_elements))
    num_bshapes = numfunctions(trial_refspace, domain(trial_element))

    scal, vec, grad = WiltonInts84.wiltonints(s1, s2, s3, x, Val{1})

    ∫∇G = -(grad[1]+0.5*γ^2*grad[3])/(4π)
    ∫n∇G = dot(n,∫∇G)
    Atot = 1/2*norm(cross(t3-t1,t3-t2))
    for i in 1:num_tshapes
        Ai = 1/2*norm(cross(test_elements.vertices[mod1(i-1,3)]-x,test_elements.vertices[mod1(i+1,3)]-x))
        g = Ai/Atot
        for j in 1:num_bshapes
            zlocal[i,j] += α * ∫n∇G * g * dx
        end
    end

    return nothing
end
#double layer transposed with pyramid basis and patch testing
function innerintegrals!(op::HH3DDoubleLayerTransposedSng, test_neighborhood,
    test_refspace::LagrangeRefSpace{T,0,3,1} where {T},
    trial_refspace::LagrangeRefSpace{T,1,3,3} where {T},
    test_elements, trial_element, zlocal, quadrature_rule::WiltonSERule, dx)
    # TODO: It is better to dispatch on γ
    # to handle the static case more efficiently    
    γ = gamma(op)
    α = op.alpha

    s1, s2, s3 = trial_element.vertices
    t1, t2, t3 = test_elements.vertices
    x = cartesian(test_neighborhood)
    n = normalize((t1-t3)×(t2-t3))
    ρ = x - dot(x - s1, n) * n

    num_tshapes = numfunctions(test_refspace, domain(test_elements))
    num_bshapes = numfunctions(trial_refspace, domain(trial_element))

    _, _, _, grad = WiltonInts84.higherorder(s1,s2,s3,x,3)

    ∫∇G = -(grad[1] + 0.5*γ^2*grad[2]) / (4π)
    for i in 1:num_tshapes
        for j in 1:num_bshapes
            ∫n∇G = dot(n,∫∇G[j])
            zlocal[i,j] += α * ∫n∇G  * dx
        end
    end

    return nothing
end
#double layer transposed with pyramid basis and pyramid testing
function innerintegrals!(op::HH3DDoubleLayerTransposedSng, test_neighborhood,
    test_refspace::LagrangeRefSpace{T,1,3,3} where {T},
    trial_refspace::LagrangeRefSpace{T,1,3,3} where {T},
    test_elements, trial_element, zlocal, quadrature_rule::WiltonSERule, dx)
    # TODO: It is better to dispatch on γ
    # to handle the static case more efficiently    
    γ = gamma(op)
    α = op.alpha

    s1, s2, s3 = trial_element.vertices
    t1, t2, t3 = test_elements.vertices
    x = cartesian(test_neighborhood)
    n = normalize((t1-t3)×(t2-t3))
    ρ = x - dot(x - s1, n) * n

    num_tshapes = numfunctions(test_refspace, domain(test_elements))
    num_bshapes = numfunctions(trial_refspace, domain(trial_element))

    _, _, _, grad = WiltonInts84.higherorder(s1,s2,s3,x,3)

    ∫∇G = -(grad[1] + 0.5*γ^2*grad[2]) / (4π)

    Atot = 1/2*norm(cross(t3-t1,t3-t2))
    for i in 1:num_tshapes
        Ai = 1/2*norm(cross(test_elements.vertices[mod1(i-1,3)]-x,test_elements.vertices[mod1(i+1,3)]-x))
        g = Ai/Atot
        for j in 1:num_bshapes
            ∫n∇G = dot(n,∫∇G[j])
            zlocal[i,j] += α * ∫n∇G * g * dx
        end
    end

    return nothing
end
# double layer
function (igd::Integrand{<:HH3DDoubleLayerReg})(x,y,f,g)
    # TODO: It is better to dispatch on γ
    # to handle the static case more efficiently    
    γ = gamma(igd.operator)

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
#double layer with pyramid basis and pyramid testing
function innerintegrals!(op::HH3DDoubleLayerSng, p,
    g::LagrangeRefSpace{T,1} where {T},
    f::LagrangeRefSpace{T,1} where {T},
    t, s, z, quadrature_rule::WiltonSERule, dx)
    # TODO: It is better to dispatch on γ
    # to handle the static case more efficiently    
    γ = gamma(op)
    α = op.alpha

    s1, s2, s3 = s.vertices
    t1, t2, t3 = t.vertices

    num_tshapes = numfunctions(g, domain(t))
    num_bshapes = numfunctions(f, domain(s))

    x = cartesian(p)
    n = normalize((s1-s3)×(s2-s3))
    ρ = x - dot(x - s1, n) * n

    _, _, _, grad = WiltonInts84.higherorder(s1,s2,s3,x,3)

    ∫∇G = -(grad[1] + 0.5*γ^2*grad[2]) / (4π)

    Atot = 1/2*norm(cross(t3-t1,t3-t2))
    for i in 1:num_tshapes
        Ai = 1/2*norm(cross(t.vertices[mod1(i-1,3)]-x,t.vertices[mod1(i+1,3)]-x))
        g = Ai/Atot
        for j in 1:num_bshapes     
        z[i,j] += α * dot(n,-∫∇G[j]) * g * dx
        end
    end

    return nothing
end
#double layer with patch basis and pyramid testing
function innerintegrals!(op::HH3DDoubleLayerSng, p,
    g::LagrangeRefSpace{T,0} where {T},
    f::LagrangeRefSpace{T,0} where {T},
    t, s, z, quadrature_rule::WiltonSERule, dx)
    # TODO: It is better to dispatch on γ
    # to handle the static case more efficiently    
    γ = gamma(op)
    α = op.alpha

    num_tshapes = numfunctions(g, domain(t))
    num_bshapes = numfunctions(f, domain(s))

    s1, s2, s3 = s.vertices
    t1, t2, t3 = t.vertices

    x = cartesian(p)
    n = normalize((s1-s3)×(s2-s3))
    ρ = x - dot(x - s1, n) * n

    scal, vec, grad = WiltonInts84.wiltonints(s1, s2, s3, x, Val{1})

    ∫∇G = -(grad[1]+0.5*γ^2*grad[3])/(4π)

    for i in 1:num_tshapes
        for j in 1:num_bshapes
            z[i,j] += α * dot(n,-∫∇G) * dx #why the minus?
        end
    end
    return nothing
end
#double layer with patch basis and pyramid testing
function innerintegrals!(op::HH3DDoubleLayerSng, p,
    g::LagrangeRefSpace{T,1} where {T},
    f::LagrangeRefSpace{T,0} where {T},
    t, s, z, quadrature_rule::WiltonSERule, dx)
    # TODO: It is better to dispatch on γ
    # to handle the static case more efficiently    
    γ = gamma(op)
    α = op.alpha

    s1, s2, s3 = s.vertices
    t1, t2, t3 = t.vertices

    num_tshapes = numfunctions(g, domain(t))
    num_bshapes = numfunctions(f, domain(s))

    x = cartesian(p)
    n = normalize((s1-s3)×(s2-s3))
    ρ = x - dot(x - s1, n) * n

    scal, vec, grad = WiltonInts84.wiltonints(s1, s2, s3, x, Val{1})

    ∫∇G = -(grad[1]+0.5*γ^2*grad[3])/(4π)

    Atot = 1/2*norm(cross(t3-t1,t3-t2))
    for i in 1:num_tshapes
        Ai = 1/2*norm(cross(t.vertices[mod1(i-1,3)]-x,t.vertices[mod1(i+1,3)]-x))
        g = Ai/Atot
        for j in 1:num_bshapes
        z[i,j] += α * dot(n,-∫∇G) * g * dx
        end
    end
    return nothing
end

#double layer with pyramid basis and patch testing
function innerintegrals!(op::HH3DDoubleLayerSng, p,
    g::LagrangeRefSpace{T,0} where {T},
    f::LagrangeRefSpace{T,1} where {T},
    t, s, z, quadrature_rule::WiltonSERule, dx)
    # TODO: It is better to dispatch on γ
    # to handle the static case more efficiently    
    γ = gamma(op)
    α = op.alpha

    s1, s2, s3 = s.vertices

    x = cartesian(p)
    n = normalize((s1-s3)×(s2-s3))
    ρ = x - dot(x - s1, n) * n

    _, _, _, grad = WiltonInts84.higherorder(s1,s2,s3,x,3)

    ∫∇G = -(grad[1] + 0.5*γ^2*grad[2]) / (4π)

    num_tshapes = numfunctions(g, domain(t))
    num_bshapes = numfunctions(f, domain(s))
  
    for i in 1:num_tshapes
        for j in 1:num_bshapes
        z[i,j] += α * dot(n,-∫∇G[j]) * dx
        end
    end

    return nothing
end


function (igd::Integrand{<:HH3DHyperSingularReg})(x,y,f,g)
    α = igd.operator.alpha
    β = igd.operator.beta
    # TODO: It is better to dispatch on γ
    # to handle the static case more efficiently    
    γ = gamma(igd.operator)

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
    # TODO: It is better to dispatch on γ
    # to handle the static case more efficiently    
    γ = gamma(op)
    s1, s2, s3 = s.vertices
    t1, t2, t3 = t.vertices
    x = cartesian(p)
    nx = normalize((s1-s3)×(s2-s3))
    ny = normalize((t1-t3)×(t2-t3))

    ∫Rⁿ, ∫RⁿN = WiltonInts84.higherorder(s1,s2,s3,x,3)
    greenconst = (∫Rⁿ[2] + 0.5*γ^2*∫Rⁿ[3]) / (4π)
    greenlinear = (∫RⁿN[2] + 0.5*γ^2*∫RⁿN[3] ) / (4π)

    num_tshapes = numfunctions(g, domain(t))
    num_bshapes = numfunctions(f, domain(s))

    jt = volume(t) * factorial(dimension(t))
    js = volume(s) * factorial(dimension(s))
    curlt = [(t3-t2)/jt,(t1-t3)/jt,(t2-t1)/jt]
    curls = [(s3-s2)/js,(s1-s3)/js,(s2-s1)/js]
    Atot = 1/2*norm(cross(t3-t1,t3-t2))
    for i in 1:num_tshapes
        Ai = 1/2*norm(cross(t.vertices[mod1(i-1,3)]-x,t.vertices[mod1(i+1,3)]-x))
        g = Ai/Atot
        for j in 1:num_bshapes
           z[i,j] += β * dot(curlt[i],curls[j])*greenconst*dx + α*dot(nx,ny) * greenlinear[j]*g*dx
        end
    end
end
