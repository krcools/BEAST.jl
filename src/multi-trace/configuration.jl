abstract type DomainData end
abstract type Domain{T<:DomainData} end
using PrettyTables
"""
A is the type of the test structure, B is the type of the trial structure, C the type of ths structure embedding
"""
struct Interaction{A,B,C}# wordt gecrieerd bij de itteratie over alle volumes.
    config
    testvol::A
    trialvol::B
    embedvol::C
end
# abstract type PhysicalInformation end #depends on the strategy
"""
subtypes of the NumericalStrategy type descirbe the followed strategy. for example a preconditioning strategy for a vector potential based problem,
 a vector potential based problem, a fields based problem,...
"""
abstract type NumericalStrategy end

export NumericalStrategy
mutable struct HomogeneousDomain <: DomainData
physicalinformation
testbasises
trialbasises
coeff
testindex 
trialindex 
HomogeneousDomain(a,b,c,d) = new(a,b,c,d,[],[])
end

struct BackgroundDomain <: DomainData
    physicalinformation
    coeff
end
physicalconstants(d::DomainData) = d.physicalinformation
physicalconstants(d::Domain) = physicalconstants(d.data)
testspace(d::DomainData) = d.testbasises
trialspace(d::DomainData) = d.trialbasises

mutable struct SubDomain{T} <: Domain{T}
    id::Int
    children::Vector{Domain}
    parent::Domain
    data::T
    excitation
    results
end
mutable struct RootDomain{T} <: Domain{T}
    id::Int
    children::Vector{Domain}
    data::T
    excitation
end

mutable struct Configuration #TODO add dict contining all subbasis info af touching objects
    domains::Dict{Int,Domain}
    root::RootDomain
    touching::Dict{Tuple{Int,Int},Vector}
    testdirectproductspace 
    trialdirectproductspace 
end
Configuration(dom,root,touching) = Configuration(dom,root,touching,nothing,nothing)


function _adddomain(config::Configuration, newdom::Domain)
    id = newdom.id
    @assert !(id in keys(config.domains))
    dom = newdom.parent
    config.domains[id] = newdom
    for space in newdom.data.testbasises
        config.testdirectproductspace = config.testdirectproductspace × space
        push!(newdom.data.testindex,length(config.testdirectproductspace.factors))
    end
    for space in newdom.data.trialbasises
        config.trialdirectproductspace = config.trialdirectproductspace × space
        push!(newdom.data.trialindex,length(config.trialdirectproductspace.factors))
    end
    push!(dom.children,newdom)
end

function _createdomain(config::Configuration,id::Int,parentid::Int,values;excitation)
    dom = SubDomain(id,Domain[],config.domains[parentid],values::DomainData,excitation,[])
    return dom
end
function _createdomain(config::Configuration,id::Int,parentdom::Domain,values;excitation)
    _createdomain(config,id,parentdom.id,values;excitation)
end

function createdomain(config::Configuration,id::Int,parent::Union{Int,Domain},values::DomainData;excitation=Dict())
    _adddomain(config,_createdomain(config,id,parent,values;excitation))
end
"""
a child of b
"""
function is_child_of(a,b)
    return a∈b.children
end

function brothers(a,b)
return is_child_of(a,b.parent)
end
function createconfiguration(data::BackgroundDomain,excitation=Dict())
    r = RootDomain(0,Domain[],data,excitation)
    Configuration(Dict{Int,Domain}(0=>r),r,Dict{Tuple{Int,Int},Any}()) # typle of the indices, returns the spaces: [test1,trial1,test2,trial2]
end

function createdomains(config::Configuration,confdict::Dict) end#the confdict is the dictionary containing the information about the mesh, physical parameters,...

"""
all objects in a domain touch each other, not with subdomains, well with the boundary

"""
function alltouching(config::Configuration)
    for i in keys(config.domains)
        #config.touching[(i,i)] = [config.domains[i].data.testbasises,config.domains[i].data.trialbasises,config.domains[i].data.testbasises,config.domains[i].data.trialbasises]
        for child in config.domains[i].children
            j = child.id
            if typeof(config.domains[i]) <: SubDomain
            config.touching[(i,j)] = [config.domains[i].data.testbasises,config.domains[i].data.trialbasises,child.data.testbasises,child.data.trialbasises]
            config.touching[(j,i)] = [config.touching[(i,j)][3],config.touching[(i,j)][4],config.touching[(i,j)][1],config.touching[(i,j)][2]]
            end
            for child2 in config.domains[i].children
                k = child2.id
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

function generate_configuration(typelist,id_of_parentlist,background,backgroundexcitation=Dict())
    conf = createconfiguration(background,backgroundexcitation)
    l = length(typelist)
    for (ind,t,parent_id) in zip(1:l,typelist,id_of_parentlist)
        createdomain(conf,ind,parent_id,t)
    end
    return conf
end
# function _create_bilform(operator_matrix,test_direct_productspace,trial_direct_productspace)
# #TODO create the bilform 

# end
#function convert_inside_to_outside_basis(Ωchild,Ωparent,strat) end
function interaction_matrix(config::Configuration,id,strat::NumericalStrategy)
    @assert length(config.testdirectproductspace.factors)==length(config.trialdirectproductspace.factors)
    N = length(config.testdirectproductspace.factors)
    OperatorMatrix = fill!(Array{AbstractOperator}(undef,N,N),ZeroOperator())
    indexen = []
    Ω = config.domains[id]
    if id != 0
        indexen = [indexen; Ω.data.testindex]
        inter = Interaction(config,Ω,Ω,Ω)
        OperatorMatrix[Ω.data.testindex[1]:last(Ω.data.testindex),Ω.data.trialindex[1]:last(Ω.data.trialindex)] += inv(convert_outside_to_inside_basis(Ω,Ω.parent,strat))*inter(strat)*convert_outside_to_inside_basis(Ω,Ω.parent,strat)
        for child in Ω.children
            inter = Interaction(config,Ω,child,Ω)
            OperatorMatrix[Ω.data.testindex[1]:last(Ω.data.testindex),child.data.trialindex[1]:last(child.data.trialindex)] += inv(convert_outside_to_inside_basis(Ω,Ω.parent,strat))*inter(strat)
            inter = Interaction(config,child,Ω,Ω)
            OperatorMatrix[child.data.testindex[1]:last(child.data.testindex),Ω.data.trialindex[1]:last(Ω.data.trialindex)] += inter(strat)*convert_outside_to_inside_basis(Ω,Ω.parent,strat)
        end

    end
    for Ω1 in Ω.children
        indexen = [indexen; Ω1.data.testindex]
        inter = Interaction(config,Ω1,Ω1,Ω)
        OperatorMatrix[Ω1.data.testindex[1]:last(Ω1.data.testindex),Ω1.data.trialindex[1]:last(Ω1.data.trialindex)] += inter(strat)

        for Ω2 in Ω.children
            if Ω1!==Ω2
                inter = Interaction(config,Ω1,Ω2,Ω)
                OperatorMatrix[Ω1.data.testindex[1]:last(Ω1.data.testindex),Ω2.data.trialindex[1]:last(Ω2.data.trialindex)] += inter(strat)
            end
        end
    end

    
    id = fill!(Array{AbstractOperator}(undef,N,N),ZeroOperator())
    for i in 1:N
        i∈indexen && (id[i,i] = Identity())
    end
    pretty_table(id-OperatorMatrix, noheader=true, backend=Val(:latex))
    return assemble(dot(config.testdirectproductspace,id-OperatorMatrix,config.trialdirectproductspace),config.testdirectproductspace,config.trialdirectproductspace)
end

# function generate_problem_lhs(config::Configuration,strat::NumericalStrategy)#Make more correct by assambling for eacht domain an add with coeffeicient, add idetity at that place
#     @assert length(config.testdirectproductspace.factors)==length(config.trialdirectproductspace.factors)
#     N = length(config.testdirectproductspace.factors)
#     OperatorMatrix = fill!(Array{AbstractOperator}(undef,N,N),ZeroOperator())
#     for (id,Ω) in config.domains
#         c = Ω.data.coeff
#         if id != 0
#             inter = Interaction(config,Ω,Ω,Ω)
#             OperatorMatrix[Ω.data.testindex[1]:last(Ω.data.testindex),Ω.data.trialindex[1]:last(Ω.data.trialindex)] += c*inv(convert_outside_to_inside_basis(Ω,Ω.parent,strat))*inter(strat)*convert_outside_to_inside_basis(Ω,Ω.parent,strat)
#             for child in Ω.children
#                 inter = Interaction(config,Ω,child,Ω)
#                 OperatorMatrix[Ω.data.testindex[1]:last(Ω.data.testindex),child.data.trialindex[1]:last(child.data.trialindex)] += c*inv(convert_outside_to_inside_basis(Ω,Ω.parent,strat))*inter(strat)
#                 inter = Interaction(config,child,Ω,Ω)
#                 OperatorMatrix[child.data.testindex[1]:last(child.data.testindex),Ω.data.trialindex[1]:last(Ω.data.trialindex)] += c*inter(strat)*convert_outside_to_inside_basis(Ω,Ω.parent,strat)
#             end

#         end
#         for Ω1 in Ω.children

#             inter = Interaction(config,Ω1,Ω1,Ω)
#             OperatorMatrix[Ω1.data.testindex[1]:last(Ω1.data.testindex),Ω1.data.trialindex[1]:last(Ω1.data.trialindex)] += c*inter(strat)

#             for Ω2 in Ω.children
#                 if Ω1!==Ω2
#                     inter = Interaction(config,Ω1,Ω2,Ω)
#                     OperatorMatrix[Ω1.data.testindex[1]:last(Ω1.data.testindex),Ω2.data.trialindex[1]:last(Ω2.data.trialindex)] += c*inter(strat)
#                 end
#             end
#         end

#     end
#     id = fill!(Array{AbstractOperator}(undef,N,N),ZeroOperator())
#     for i in 1:N
#         id[i,i] = 2*Identity()
#     end
#     pretty_table(id-OperatorMatrix, noheader=true, backend=Val(:latex))
#     return assemble(dot(config.testdirectproductspace,id-OperatorMatrix,config.trialdirectproductspace),config.testdirectproductspace,config.trialdirectproductspace)

# end
function generate_problem_lhs(config::Configuration,strat::NumericalStrategy)
    out = interaction_matrix(config,0,strat)
    for (id,Ω) in config.domains
        id==0 && continue
        out += interaction_matrix(config,id,strat)
    end
    return out
end

function generate_problem_rhs(config::Configuration)
    linterms = Vector{LinTerm}()
    @warn "only excitations in root domain implemented, check first if same aproach for others is appropriate"
    for (id,Ω) in config.domains end

    for (key,func) in config.root.excitation
        for child in config.root.children


            push!(linterms,LinTerm(child.data.testindex[key],[],1,func))
        end
    end

    return assemble(LinForm([],linterms),config.testdirectproductspace)
end

function map_solution_to_volumes!(config,solution)
    for (id,dom) in config.domains
        id==0 && continue
        out = []
        for i in dom.data.trialindex
            push!(out,solution[Block(i)])
        end
        dom.results = out
    end
end




