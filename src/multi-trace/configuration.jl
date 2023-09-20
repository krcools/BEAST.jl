abstract type Domain end

abstract type PhysicalInformation end #depends on the strategy
abstract type DomainData end

struct HomogeneousDomain <: DomainData
ϵr
μr
testbasises
trialbasises
end

struct BackgroundDomain <: DomainData
    ϵ0
    μ0
end

struct SubDomain <: Domain
    id::Int
    children::Vector{Domain}
    parent::Domain
    data::DomainData
end
struct RootDomain <: Domain
    id::Int
    children::Vector{Domain}
    data::DomainData
end

struct Configuration #TODO add dict contining all subbasis info af touching objects
    domains::Dict{Int,Domain}
    root::RootDomain
    touching::Dict{Tuple{Int,Int},<:Vector}
end

function _adddomain(config::Configuration, newdom::Domain)
    id = newdom.id
    @assert !(id in keys(config.domains))
    dom = newdom.parent
    config.domains[id] = newdom
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


