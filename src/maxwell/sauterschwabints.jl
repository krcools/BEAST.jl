function momintegrals!(op::MWSingleLayer3D, g::RTRefSpace, f::RTRefSpace,
                            t, s, z, strat::sauterschwabstrategy)
    accuracy = 12
    accuracy_pd = 8

    ĝ = op.gamma
    b = op.β
    a = op.α


    function G(x,y)
        exp(-ĝ*norm(x-y))/(4pi*norm(x-y))
    end


    for i = 1:3
        for j = 1:3

            function Integrand(x,y)
                a*(((x-t.vertices[i])'*(y-s.vertices[j]))/(2*(volume(t))*2*volume(s)))*G(x,y)+
                b*((2*2)/(2*volume(t)*2*volume(s)))*G(x,y)
            end

            z[i,j] = sauterschwabintegral(s, t, Integrand, accuracy, accuracy_pd)

        end
    end
end




function momintegrals!(op::MWDoubleLayer3D, g::RTRefSpace, f::RTRefSpace,
                        t, s, z, strat::sauterschwabstrategy)
    accuracy = 12
    accuracy_pd = 8

    ĝ = op.gamma


    function G(x,y)
        exp(-ĝ*norm(x-y))/(4pi*norm(x-y))
    end

    function u(x,y)
        -G(x,y)*(1/norm(x-y))*((1/norm(x-y))+ĝ)*(x-y)
    end


    for i = 1:3
        for j = 1:3

            function Integrand(x,y)
                ((x-t.vertices[i])/(2*volume(t)))'*cross(u(x,y) , (y-s.vertices[j])/(2*volume(s)))
            end

            z[i,j] = sauterschwabintegral(s, t, Integrand, accuracy, accuracy_pd)

        end
    end
end
