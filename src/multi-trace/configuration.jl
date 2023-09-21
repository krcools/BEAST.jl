abstract type Domain end

abstract type PhysicalInformation end #depends on the strategy
abstract type DomainData end

mutable struct HomogeneousDomain <: DomainData
ϵr
μr
testbasises
trialbasises
testindex = []
trialindex = []
end

struct BackgroundDomain <: DomainData
    ϵ0
    μ0
end

mutable struct SubDomain{T<:DomainData} <: Domain
    id::Int
    children::Vector{Domain}
    parent::Domain
    data::T
end
mutable struct RootDomain{T<:DomainData} <: Domain
    id::Int
    children::Vector{Domain}
    data::T
end

mutable struct Configuration #TODO add dict contining all subbasis info af touching objects
    domains::Dict{Int,Domain}
    root::RootDomain
    touching::Dict{Tuple{Int,Int},<:Vector}
    testdirectproductspace 
    testcounter 
    trialcounter 
    trialdirectproductspace 
end
function Configuration(dom,root,touching)
Configuration(dom,root,touching,nothing,0,0,nothing)
end

function _adddomain(config::Configuration, newdom::Domain)
    id = newdom.id
    @assert !(id in keys(config.domains))
    dom = newdom.parent
    config.domains[id] = newdom
    for space in newdom.data.testbasises
        config.testcounter += 1
        push!(newdom.data.testindex,config.testcounter)
        config.testdirectproductspace = space × config.test_direct_productspace
    end
    for space in newdom.data.trialbasises
        config.trialcounter += 1
        push!(newdom.data.trialindex,config.trialcounter)
        config.trialdirectproductspace = space × config.trial_direct_productspace
    end
    push!(dom.children,newdom)
end

function _createdomain(config::Configuration,id::Int,parentid::Int,values)
    dom = SubDomain(id,Domain[],config.domains[parentid],values::DomainData)

   # push!(config.domains[parentid].children,dom)
    return dom
end
function _createdomain(config::Configuration,id::Int,parentdom::Domain,values)
    _createdomain(config,id,parentdom.id,values)
end

function createdomain(config::Configuration,id::Int,parent::Union{Int,Domain},values::DomainData)
    _adddomain(config,_createdomain(config,id,parent,values))

end
"""
a child of b
"""
function is_child_of(a,b)
    return a∈b
end

function brothers(a,b)
return is_child_of(a,b.parent)
end
function createconfiguration(data::BackgroundDomain)
    r = RootDomain(0,Domain[],data)
    Configuration(Dict{Int,Domain}(0=>r),r,Dict{Tuple{Int,Int},Any}()) # typle of the indices, returns the spaces: [test1,trial1,test2,trial2]
end

function createdomains(config::Configuration,confdict::Dict) end#the confdict is the dictionary containing the information about the mesh, physical parameters,...

"""
all objects in a domain touch each other, not with subdomains, well with the boundary

"""
function alltouching(config::Configuration)
    for i in keys(config)
        #config.touching[(i,i)] = [config.domains[i].data.testbasises,config.domains[i].data.trialbasises,config.domains[i].data.testbasises,config.domains[i].data.trialbasises]
        for child in config.domains[i].children
            j = child.index
            if typeof(config.domains[i]) <: Subdomain
            config.touching[(i,j)] = [config.domains[i].data.testbasises,config.domains[i].data.trialbasises,child.data.testbasises,child.data.trialbasises]
            config.touching[(j,i)] = [config.touching[(i,j)][3],config.touching[(i,j)][4],config.touching[(i,j)][1],config.touching[(i,j)][2]]
            end
            for child2 in config.domains[i].children
                k = child2.index
                config.touching[(j,k)] = [child.data.testbasises,child.data.trialbasises,child2.data.testbasises,child2.data.trialbasises]
                config.touching[(k,j)] = [child2.data.testbasises,child2.data.trialbasises,child.data.testbasises,child.data.trialbasises]
            end
        end 
    end
end

function (dom::Domain)(n::NormalVector,dom2::Domain)
    if dom==dom2
        return 1.0
    elseif dom2 ∈ dom.children
        return -1.0
    else
        @error "domain is not in children of parent domain"
    end

end

function generate_configuration(typelist,id_of_parentlist,background)
    conf = createconfiguration(background)
    l = length(typelist)
    for (index,t,parent_id) in zip(index,typelist,id_of_parentlist)
        createdomain(conf,index,parent_id,t)
    end
    return conf
end
function _create_bilform(operator_matrix,test_direct_productspace,trial_direct_productspace)
#TODO create the bilform 

end

function generate_problem(config::Configuration)
    @assert length(config.testdirectproductspace.factors)==length(config.trialdirectproductspace.factors)
    N = length(config.testdirectproductspace.factors)
    OperatorMatrix = fill!(Array{AbstractOperator}(undef,N,N),ZeroOperator())
    for (id,Ω) in config.domains
        if id != 0
            inter = Interaction(config,Ω,Ω,Ω)
            OperatorMatrix[Ω.data.testindex[1]:last(Ω.data.testindex),Ω.data.trialindex[1]:last(Ω.data.trialindex)] += inter()
            for child in Ω.children
                inter = Interaction(config,Ω,child,Ω)
                OperatorMatrix[Ω.data.testindex[1]:last(Ω.data.testindex),child.data.trialindex[1]:last(child.data.trialindex)] += inter()
                inter = Interaction(config,child,Ω,Ω)
                OperatorMatrix[child.data.testindex[1]:last(child.data.testindex),Ω.data.trialindex[1]:last(Ω.data.trialindex)] += inter()
            end

        end
        for Ω1 in Ω.children
            for Ω2 in Ω.children
                if Ω1===Ω2
                    inter = Interaction(config,Ω1,Ω2,Ω)
                    OperatorMatrix[Ω1.data.testindex[1]:last(Ω1.data.testindex),Ω2.data.trialindex[1]:last(Ω2.data.trialindex)] += inter()
                else
                    inter = Interaction(config,Ω1,Ω2,Ω)
                    OperatorMatrix[Ω1.data.testindex[1]:last(Ω1.data.testindex),Ω2.data.trialindex[1]:last(Ω2.data.trialindex)] += inter()
                    inter = Interaction(config,Ω2,Ω1,Ω)
                    OperatorMatrix[Ω2.data.testindex[1]:last(Ω2.data.testindex),Ω1.data.trialindex[1]:last(Ω1.data.trialindex)] += inter()

                end

            end
        end
    end

#TODO the right hand side

end