using BEAST
using CompScienceMeshes
using Base.Test

m = readmesh(joinpath(dirname(@__FILE__),"sphere.in"))
X = raviartthomas(m)
κ = ω = c = 1.0

fn = Pkg.dir("BEAST","test","efie_solution.txt")
u = map(x->eval(parse(x)), readcsv(fn))

fcr = facecurrents(u,X)

T, P = linspace(0.0,π,7), 0.0
points = vec([point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in T, ϕ in P])
utheta = vec([point(cos(ϕ)*cos(θ), sin(ϕ)*cos(θ), -sin(θ)) for θ in T, ϕ in P])
uphi   = vec([point(-sin(ϕ)*sin(θ), cos(ϕ)*sin(θ), 0) for θ in T, ϕ in P])
for i in eachindex(uphi)
    uphi[i] /= norm(uphi[i])
end

FF = BEAST.MWFarField3D(1.0im)
ff = potential(FF, points, u, X)
ff *= im*ω/(4π*c)

a = Float64[real(norm(f)) for f in ff]
b = Float64[sqrt(abs(dot(ff[i],utheta[i]))^2 + abs(dot(ff[i],uphi[i]))^2) for i in eachindex(ff)]
c = [dot(uphi[i], uphi[i]) for i in eachindex(ff)]

tol = eps() * 10000
@test a[1] ≈ 0.6426645162880508  atol=tol
@test a[2] ≈ 0.5222644345134998  atol=tol
@test a[3] ≈ 0.28287430186444407 atol=tol
@test a[4] ≈ 0.3873439919314161  atol=tol
@test a[5] ≈ 0.6788129041674158  atol=tol
@test a[6] ≈ 0.8800967217249486  atol=tol
@test a[7] ≈ 0.9487895695366841  atol=tol

# NF = BEAST.MWSingleLayerField3D(κ)
# xs = linspace(-8.0, 8.0, 100)
# ys = linspace(-8.0, 8.0, 100)
# points = [point(x,y,0) for x in xs, y in ys]
# nf = BEAST.potential(NF, points, u, X)
#
# pw = PlaneWaveMW(point(0,0,1),point(1,0,0),1.0,complex(1.0))
# for i in eachindex(points); nf[i] -= pw(points[i]); end
# a = Float64[abs(norm(f[1])) for f in nf]
# b = reshape(a,size(points)...)
