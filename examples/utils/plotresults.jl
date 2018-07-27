@isdefined(plotresults) || (plotresults = false)

if plotresults

    postproc || error("Cannot plot results without going through post-processing.")
    @eval begin
        using PlotlyJS
        #include(Pkg.dir("CompScienceMeshes","examples","plotlyjs_patches.jl"))

        t1 = scatter(x=Θ, y=real.(norm.(ffd)))
        t2 = patch(geo, real.(norm.(fcr)))
        t3 = heatmap(x = xs, y = zs, z = clamp.(real.(norm.(nfd)), 0.0, 2.0))
    end
end
