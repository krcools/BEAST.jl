using BEAST
using CompScienceMeshes

o, x, y, z = euclidianbasis(3)

Γ = readmesh(Pkg.dir("BEAST","examples","sphere2.in"))
RT = raviartthomas(Γ)

κ = 1.0
t = MWSingleLayer3D(im*κ)
E = planewavemw3d(direction=z, polarization=x, wavenumber=κ)
e = (n × E) × n

@hilbertspace j
@hilbertspace k

efie = @discretise t[k,j]==e[k]  j∈RT k∈RT

u = solve(efie)
println("Solution computed")

isdefined(:plotresults) || (plotresults = false)
plotresults && (postproc = true)
isdefined(:postproc) || (postproc = false)

if postproc

    Φ, Θ = [0.0], linspace(0,π,100)
    pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for ϕ in [0.0] for θ in Θ]
    ffd = potential(MWFarField3D(im*κ), pts, u, RT)
    println("far field computed")

    fcr, geo = facecurrents(u, RT)
    println("Face currents computed")

    nx, nz = 50, 100
    X, Z = linspace(-2,2,nx), linspace(-4,4,nz)
    grid = [point(x,0,z) for x in X, z in Z]
    nfd = potential(MWSingleLayerField3D(κ), grid, u, RT)
    nfd = reshape(nfd, (nx,nz))
    nfd .-= E.(grid)
    println("near field computed.")

end

if plotresults
    @eval begin
        using PlotlyJS
        include(Pkg.dir("CompScienceMeshes","examples","plotlyjs_patches.jl"))
        t1 = scatter(x=Θ, y=real.(norm.(ffd)))
        t2 = patch(geo, real.(norm.(fcr)))
        t3 = heatmap(x = X, y = Z, z = clamp.(real.(norm.(nfd)), 0.0, 2.0))
    end
end

#plot(t1) # uncomment to plot the far field (line2d)
#plot(t2) # uncomment to plot the induced current (mesh3d)
#plot(t3) # uncomment to plot the near field (heatmap)

nothing;
