@isdefined(postproc) || (postproc = false)
@isdefined(plotresults) || (plotresults = false)
plotresults && (postproc = true)

if postproc

    Φ, Θ = [0.0], range(0,stop=π,length=100)
    pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for ϕ in Φ for θ in Θ]
    ffd = potential(MWFarField3D(wavenumber=κ), pts, u, X)

    fcr, geo = facecurrents(u, X)

    nx, nz = 50, 100
    x = 0.0
    ys, zs = range(-2,stop=2,length=nx), range(-4,stop=4,length=nz)
    gridpoints = [point(x,y,z) for y in ys, z in zs]
    nfd = potential(MWSingleLayerField3D(wavenumber = κ), gridpoints, u, X)
    nfd = reshape(nfd, (nx,nz))
    nfd .-= E.(gridpoints)
end
