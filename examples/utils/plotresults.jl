@isdefined(plotresults) || (plotresults = false)

if plotresults

    postproc || error("Cannot plot results without going through post-processing.")
    @eval begin
        import Plots
        using LinearAlgebra

        p1 = Plots.scatter(Θ, real.(norm.(ffd)))
        p2 = Plots.heatmap(clamp.(real.(norm.(nfd)), 0.0, 2.0), clims=(0,2))
        p3 = Plots.contour(clamp.(real.(norm.(nfd)), 0.0, 2.0), clims=(0,2))
        Plots.plot(p1,p2,p3,layout=(3,1))
    end
end
