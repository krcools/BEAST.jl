@isdefined(plotresults) || (plotresults = false)

if plotresults

    postproc || error("Cannot plot results without going through post-processing.")
    @eval begin
        using Plots
        using LinearAlgebra

        p1 = scatter(Θ, real.(norm.(ffd)))
        p2 = heatmap(clamp.(real.(norm.(nfd)), 0.0, 2.0))
        p3 = contour(clamp.(real.(norm.(nfd)), 0.0, 2.0))
        plot(p1,p2,p3,layout=(3,1))
    end
end
