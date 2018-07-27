# function momintegrals!(op::MWSingleLayer3D, g::RTRefSpace, f::RTRefSpace,
#                             t, s, z, strat::sauterschwabstrategy)
#     accuracy = 12
#     accuracy_pd = 8
#
#     ĝ = op.gamma
#     b = op.β
#     a = op.α
#
#
#     function G(x,y)
#         exp(-ĝ*norm(x-y))/(4pi*norm(x-y))
#     end
#
#
#     for i = 1:3
#         for j = 1:3
#
#             function Integrand(x,y)
#                 a*(((x-t.vertices[i])'*(y-s.vertices[j]))/(2*(volume(t))*2*volume(s)))*G(x,y)+
#                 b*((2*2)/(2*volume(t)*2*volume(s)))*G(x,y)
#             end
#
#             z[i,j] = sauterschwabintegral(s, t, Integrand, accuracy, accuracy_pd)
#
#         end
#     end
# end

function momintegrals!(op::MWSingleLayer3D,
    test_local_space::RTRefSpace, trial_local_space::RTRefSpace,
    test_triangular_element, trial_triangular_element, out, strat::SauterSchwabStrategy)

    I, J, K, L = SauterSchwabQuadrature.reorder(
        test_triangular_element.vertices,
        trial_triangular_element.vertices)

    test_triangular_element  = simplex(test_triangular_element.vertices[I])
    trial_triangular_element = simplex(trial_triangular_element.vertices[I])

    α = op.α
    β = op.β
    γ = op.gamma

    function igd(u,v)

        x = neighborhood(test_triangular_element,u)
        y = neighborhood(trial_triangular_element,v)

        r = cartesian(x) - cartesian(y)
        R = norm(r)
        G = exp(-γ*R)/(4π*R)

        f = test_local_space(x)
        g = trial_local_space(y)

        jx = jacobian(x)
        jy = jacobian(y)

        jx*jy*SMatrix{3,3}(
            α*dot(f[1][1], G*g[1][1]) + β*dot(f[1][2], G*g[1][2]),
            α*dot(f[2][1], G*g[1][1]) + β*dot(f[2][2], G*g[1][2]),
            α*dot(f[3][1], G*g[1][1]) + β*dot(f[3][2], G*g[1][2]),
            α*dot(f[1][1], G*g[2][1]) + β*dot(f[1][2], G*g[2][2]),
            α*dot(f[2][1], G*g[2][1]) + β*dot(f[2][2], G*g[2][2]),
            α*dot(f[3][1], G*g[2][1]) + β*dot(f[3][2], G*g[2][2]),
            α*dot(f[1][1], G*g[3][1]) + β*dot(f[1][2], G*g[3][2]),
            α*dot(f[2][1], G*g[3][1]) + β*dot(f[2][2], G*g[3][2]),
            α*dot(f[3][1], G*g[3][1]) + β*dot(f[3][2], G*g[3][2]),);

    end


    ssq = sauterschwab_parameterized(igd, strat)
    out .+= ssq[K,L]

    nothing
end



# function momintegrals!(op::MWDoubleLayer3D, g::RTRefSpace, f::RTRefSpace,
#                         t, s, z, strat::sauterschwabstrategy)
#     accuracy = 12
#     accuracy_pd = 8
#
#     ĝ = op.gamma
#
#
#     function G(x,y)
#         exp(-ĝ*norm(x-y))/(4pi*norm(x-y))
#     end
#
#     function u(x,y)
#         -G(x,y)*(1/norm(x-y))*((1/norm(x-y))+ĝ)*(x-y)
#     end
#
#
#     for i = 1:3
#         for j = 1:3
#
#             function Integrand(x,y)
#                 ((x-t.vertices[i])/(2*volume(t)))'*cross(u(x,y) , (y-s.vertices[j])/(2*volume(s)))
#             end
#
#             z[i,j] = sauterschwabintegral(s, t, Integrand, accuracy, accuracy_pd)
#
#         end
#     end
# end


function momintegrals!(op::MWDoubleLayer3D,
    test_local_space::RTRefSpace, trial_local_space::RTRefSpace,
    test_triangular_element, trial_triangular_element, out, strat::SauterSchwabStrategy)

    I = indexin(test_triangular_element, trial_triangular_element)
    @assert I == [1,2,3]

    γ = op.gamma

    function igd(u,v)

        x = test_triangular_element(u)
        y = trial_triangular_element(v)

        r = cartesian(x) - cartesian(y)
        R = norm(r)

        expn = exp(-γ*R)
        G = expn / (4pi*R)
        GG = -(γ + 1/R) * G / R * r
        T = @SMatrix [
            0 -GG[3] GG[2]
            GG[3] 0 -GG[1]
            -GG[2] GG[1] 0 ]

        f = test_local_space(x)
        g = trial_local_space(y)

        jx = jacobian(x)
        jy = jacobian(y)

        out .+= jx*jy*SMatrix(
            (
                α*dot(f[1][1], T*g[1][1]),
                α*dot(f[2][1], T*g[1][1]),
                α*dot(f[3][1], T*g[1][1]),),
            (
                α*dot(f[1][1], T*g[2][1]),
                α*dot(f[2][1], T*g[2][1]),
                α*dot(f[3][1], T*g[2][1]),),
            (
                α*dot(f[1][1], T*g[3][1]),
                α*dot(f[2][1], T*g[3][1]),
                α*dot(f[3][1], T*g[3][1]),),);
    end

    sauterschwab_parameterized(igd, strat)

    nothing
end
