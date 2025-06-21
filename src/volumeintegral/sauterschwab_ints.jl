

# 6D integral: ∫∫∫_Ω ∫∫∫_Ω, 5D integral: ∫∫∫_Ω ∫∫_Γ
# function (f::PulledBackIntegrand)(u,v) <--- see quadrature/sauterschwabints.jl
#     ...
# end


# 5D integral: ∫∫_Γ ∫∫∫_Ω
function (f::PulledBackIntegrand{<:Integrand{<:BoundaryOperator}, 
    <:CompScienceMeshes.Simplex{<:Any,2,<:Any,3}, <:CompScienceMeshes.Simplex{<:Any,3,<:Any,4} })(v,u) 
    # "In general I think a Jacobian determinant needs to be included. For Simplical and
    # Quadrilateral charts this is not needed because they are 1."
    f.igd(cartesian(f.chart1,u), cartesian(f.chart2,v))
end



function reorder(test_chart::CompScienceMeshes.Simplex{<:Any,2,<:Any,3}, trial_chart::CompScienceMeshes.Simplex{<:Any,3,<:Any,4}, 
    strat::SauterSchwab3DStrategy) 
    J, I = SauterSchwab3D.reorder(strat.sing) # 5D integral: ∫∫_Γ ∫∫∫_Ω
    return I, J
end

function reorder(test_chart::CompScienceMeshes.Simplex, trial_chart::CompScienceMeshes.Simplex, 
    strat::SauterSchwab3DStrategy)
    I, J = SauterSchwab3D.reorder(strat.sing) # 6D integral: ∫∫∫_Ω ∫∫∫_Ω, 5D integral: ∫∫∫_Ω ∫∫_Γ
    return I, J
end



#TODO: use trial_ptr to get the cell material using material_array[trial_ptr] ...

function momintegrals!(out, op::VIEOperator,
    test_functions::Space, test_ptr, test_chart,
    trial_functions::Space, trial_ptr, trial_chart,
    strat::SauterSchwab3DStrategy)

    test_local_space = refspace(test_functions)
    trial_local_space = refspace(trial_functions)

    num_tshapes = numfunctions(test_local_space, domain(test_chart))
    num_bshapes = numfunctions(trial_local_space, domain(trial_chart))

    I, J = reorder(test_chart, trial_chart, strat)
        
    igd = Integrand(op, test_local_space, trial_local_space, test_chart, trial_chart)
    igdp = pulledback_integrand(igd, I, test_chart, J, trial_chart)
    Q = SauterSchwab3D.sauterschwab_parameterized(igdp, strat)


    for j in 1:num_bshapes
        for i in 1:num_tshapes
            out[i,j] += Q[i,j]
        end
    end

    return nothing
end



#TODO: use trial_ptr to get the cell material using material_array[trial_ptr] ...

# function momintegrals!(z, op::VIEOperator,
#     test_functions::Space, test_cellptr, test_chart,
#     trial_functions::Space, trial_cellptr, trial_chart,
#     strat::DoubleQuadRule)

#     tshs = refspace(test_functions)
#     bshs = refspace(trial_functions) 

#     igd = Integrand(op, tshs, bshs, test_chart, trial_chart) 

#     womps = strat.outer_quad_points
#     wimps = strat.inner_quad_points
    
#     for womp in womps
#         tgeo = womp.point
#         tvals = womp.value
#         M = length(tvals)
#         jx = womp.weight
        
#         for wimp in wimps
#             bgeo = wimp.point
#             bvals = wimp.value
#             N = length(bvals)
#             jy = wimp.weight

#             j = jx * jy

#             z1 = j * igd(tgeo, bgeo, tvals, bvals)
#             for n in 1:N
#                 for m in 1:M
#                     z[m,n] += z1[m,n]
#             end end
#         end
#     end

#     return z
# end