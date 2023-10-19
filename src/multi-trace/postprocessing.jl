import Base.Threads: @spawn
function nearfield_A(dom::Domain{HomogeneousDomain},pts,strat::VectorStrat)
    p = physicalconstants(dom)
    parent_const = physicalconstants(dom.parent)
    k = sqrt(p.ϵ*p.μ)*p.ω

    nxb_out1,σ_out1,a_out1,γ_out1 = dom.results
    nxb1,σ1,a1,γ1 = dom.data.trialbasises

    vector_potential = p.μ/parent_const.μ* potential(BEAST.HHHGreenField(wavenumber=k),pts,a_out1,a1) .- 
    potential(BEAST.HHHGradGreenCrossField(wavenumber=k),pts,nxb_out1,nxb1) .-
    parent_const.ϵ/p.ϵ*potential(BEAST.HHHGradGreenField(wavenumber=k),pts,γ_out1,γ1) .-#compensatie omdat normal naar buiten wijst
    potential(BEAST.HHHGreenField(k*im,BEAST.HHHBasisNtimesField(BEAST.HHHIdentityField())),pts,σ_out1,σ1);
    
    for Ω in dom.children
    nxb_out,σ_out,a_out,γ_out = Ω.results
    nxb,σ,a,γ = Ω.data.trialbasises
    vector_potential = vector_potential .+ -1*potential(BEAST.HHHGreenField(wavenumber=k),pts,a_out,a) .+ 
    potential(BEAST.HHHGradGreenCrossField(wavenumber=k),pts,nxb_out,nxb) .+
    potential(BEAST.HHHGradGreenField(wavenumber=k),pts,γ_out,γ) .+#compensatie omdat normal naar buiten wijst
    potential(BEAST.HHHGreenField(k*im,BEAST.HHHBasisNtimesField(BEAST.HHHIdentityField())),pts,σ_out,σ);
    end
    return vector_potential
end
function nearfield_A(dom::Domain{BackgroundDomain},pts,strat::VectorStrat)
    p = physicalconstants(dom)
    k = sqrt(p.ϵ*p.μ)*p.ω
    vector_potential = zeros(SVector{3,Complex},size(pts))
    for Ω in dom.children
    nxb_out,σ_out,a_out,γ_out = Ω.results
    nxb,σ,a,γ = Ω.data.trialbasises
    vector_potential = vector_potential + -1*potential(BEAST.HHHGreenField(wavenumber=k),pts,a_out,a) + 
    potential(BEAST.HHHGradGreenCrossField(wavenumber=k),pts,nxb_out,nxb) +
    potential(BEAST.HHHGradGreenField(wavenumber=k),pts,γ_out,γ) +#compensatie omdat normal naar buiten wijst
    potential(BEAST.HHHGreenField(k*im,BEAST.HHHBasisNtimesField(BEAST.HHHIdentityField())),pts,σ_out,σ);
    end
    return vector_potential
end
function nearfield_A(config::Configuration,points,strat::VectorStrat)
    tasks = []
    for (id,dom) in config.domains
        push!(tasks,@spawn nearfield_A(dom,points,strat))
    end
    fields = fetch.(tasks)
    out = zeros(SVector{3,Complex},size(points))
    for field in fields
        out = out .+ field
    end
    return out
end
function nearfield_B(dom::Domain{HomogeneousDomain},pts,strat::VectorStrat)
    p = physicalconstants(dom)
    parent_const = physicalconstants(dom.parent)
    k = sqrt(p.ϵ*p.μ)*p.ω

    nxb_out1,σ_out1,a_out1,γ_out1 = dom.results
    nxb1,σ1,a1,γ1 = dom.data.trialbasises

    BField = p.μ/parent_const.μ* potential(BEAST.HHHGradGreenCrossField(wavenumber=k),pts,a_out1,a1) .- 
    (potential(BEAST.HHHGradGreenField(k*im,BEAST.HHHDivergenceField()),pts,nxb_out1,nxb1)+k^2*potential(BEAST.HHHGreenField(wavenumber=k),pts,nxb_out1,nxb1)) .-
    potential(BEAST.HHHGradGreenCrossField(k*im,BEAST.HHHBasisNtimesField(BEAST.HHHIdentityField())),pts,σ_out1,σ1);
    
    for Ω in dom.children
    nxb_out,σ_out,a_out,γ_out = Ω.results
    nxb,σ,a,γ = Ω.data.trialbasises
    BField = BField .+ -1*potential(BEAST.HHHGradGreenCrossField(wavenumber=k),pts,a_out,a) .+ 
    (potential(BEAST.HHHGradGreenField(k*im,BEAST.HHHDivergenceField()),pts,nxb_out,nxb)+k^2*potential(BEAST.HHHGreenField(wavenumber=k),pts,nxb_out,nxb)) .+
    potential(BEAST.HHHGradGreenCrossField(k*im,BEAST.HHHBasisNtimesField(BEAST.HHHIdentityField())),pts,σ_out,σ);
    end
    return BField
end

function nearfield_H(dom::Domain{HomogeneousDomain},points,strat::VectorStrat)
    nearfield_B(dom,points,strat)./physicalconstants(dom).μ
end


function nearfield_B(dom::Domain{BackgroundDomain},pts,strat::VectorStrat)
    p = physicalconstants(dom)
    k = sqrt(p.ϵ*p.μ)*p.ω
    BField = zeros(SVector{3,Complex},size(pts))
    for Ω in dom.children
        nxb_out,σ_out,a_out,γ_out = Ω.results
        nxb,σ,a,γ = Ω.data.trialbasises
        BField = BField + -1*potential(BEAST.HHHGradGreenCrossField(wavenumber=k),pts,a_out,a) + 
        (potential(BEAST.HHHGradGreenField(k*im,BEAST.HHHDivergenceField()),pts,nxb_out,nxb)+k^2*potential(BEAST.HHHGreenField(wavenumber=k),pts,nxb_out,nxb))  +
        potential(BEAST.HHHGradGreenCrossField(k*im,BEAST.HHHBasisNtimesField(BEAST.HHHIdentityField())),pts,σ_out,σ)
    end
    return BField
end
function nearfield_B(config::Configuration,points,strat::VectorStrat)
    tasks = []
    for (id,dom) in config.domains
        push!(tasks,@spawn nearfield_B(dom,points,strat))
    end
    fields = fetch.(tasks)
    out = zeros(SVector{3,Complex},size(points))
    for field in fields
        out = out .+ field
    end
    return out
end
function nearfield_H(dom::Domain{BackgroundDomain},points,strat::VectorStrat)
    nearfield_B(dom,points,strat)./physicalconstants(dom).μ
end

function nearfield_H(config::Configuration,points,strat::VectorStrat)
    tasks = []
    for (id,dom) in config.domains
        push!(tasks,nearfield_H(dom,points,strat))
    end
    #fields = fetch.(tasks)
    out = zeros(SVector{3,Complex},size(points))
    for field in tasks
        out = out .+ field
    end
    return out
end