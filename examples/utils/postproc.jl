isdefined(:postproc) || (postproc = false)

if postproc

    Φ, Θ = [0.0], linspace(0,π,100)
    pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for ϕ in Φ for θ in Θ]
    ffd = potential(MWFarField3D(im*κ), pts, u, X)

    fcr, geo = facecurrents(u, X)

    nx, nz = 50, 100
    xs, zs = linspace(-2,2,nx), linspace(-4,4,nz)
    grid = [point(x,0,z) for x in xs, z in zs]
    nfd = potential(MWSingleLayerField3D(κ), grid, u, X)
    nfd = reshape(nfd, (nx,nz))
    nfd .-= E.(grid)
end
