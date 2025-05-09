# # EFIE
#
# This example demonstrates the modelling of scattering by a perfectly conducting
# sphere using the Electric Field Integral Equation (EFIE).
#

using LinearAlgebra
using CompScienceMeshes
using BEAST
import PlotlyDocumenter # hide

# We call upon GMsh to create a flat-faceted triangular mesh for the sphere. 
# Subordinate to this mesh the space of lowest order Raviart-Thomas elements is
# constructed. For reproducibility of this example, a pre-constructed mesh is
# loaded from disk.

Γ = readmesh(joinpath(dirname(pathof(BEAST)),"../examples/sphere2.in"))
X = raviartthomas(Γ);

# The integral equation can be wrtiten down in terms of the single layer operaotr
# and a plane wave excitation.

κ, η = 1.0, 1.0;
t = Maxwell3D.singlelayer(wavenumber=κ);
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ);
e = (n × E) × n;

# To write down the final variational formulation, placeholders for the test and
# trial functions are required.

@hilbertspace j;
@hilbertspace k;
efie = @discretise t[k,j]==e[k]  j∈X k∈X;

# The system is assembled and solved by GMRES. When this processs terminated, the
# solution and the convergence history are returned.

u, ch = BEAST.gmres_ch(efie; restart=1500);
@show ch

# With the solution in hand, various secondary quatities can be computed. Here we
# query the far field pattern, the near field, and the norm of the solution at the 
# barycenters of the triangles in the mesh.

Φ, Θ = [0.0], range(0,stop=π,length=100);
pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for ϕ in Φ for θ in Θ];
ffd = potential(MWFarField3D(wavenumber=κ), pts, u, X);

fcr, geo = facecurrents(u, X);

ys = range(-2,stop=2,length=50);
zs = range(-4,stop=4,length=100);
gridpoints = [point(0,y,z) for y in ys, z in zs];
Esc = potential(MWSingleLayerField3D(wavenumber = κ), gridpoints, u, X);
Ein = E.(gridpoints);
"" #hide

# The results can be exported or visualised using e.g. Plotly
import PlotlyBase as Plt
tr1 = Plt.scatter(x=Θ, y=norm.(ffd))
pl1 = Plt.Plot(tr1)

tr2 = Plt.heatmap(x=ys, y=zs, z=norm.(Esc-Ein),
    colorscale="Viridis", zmin=0, zmax=2, showscale=false)
pl2 = Plt.Plot(tr2)

tr3 = CompScienceMeshes.patch(geo, norm.(fcr);
    caxis=(0,2), colorscale="Viridis")
pl3 = Plt.Plot(tr3)
pl = [[pl1; pl2] pl3]
PlotlyDocumenter.to_documenter(pl) #hide