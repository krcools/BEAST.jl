using CompScienceMeshes
using BEAST
using LinearAlgebra

trias = CompScienceMeshes.meshrectangle(1.0,1.0,0.05)
trias = CompScienceMeshes.meshcircle(1.0,0.5)
tetrs = CompScienceMeshes.tetmeshcuboid(1.0,1.0,1.0,0.2;generator = :gmsh)

f = BEAST.ScalarTrace{Float64}(x -> sin(pi*x[1]*x[2])*sin(pi*x[2]))
f2 = BEAST.ScalarTrace{Float64}(x -> (-x[2]*sin(2*pi*(x[1]^2+x[2]^2))*sin(atan(x[2],x[1])), x[1]*sin(2*pi*(x[1]^2+x[2]^2))*sin(atan(x[2],x[1])), 0) )
g = BEAST.ScalarTrace{Float64}(x -> (sin(pi*x[1])*cos(pi*x[2]), -cos(pi*x[1])*sin(pi*x[2]), 0 ))
g = BEAST.ScalarTrace{Float64}(x -> (sin(pi*x[1])*cos(pi*x[2]), sin(pi*x[2])*cos(pi*x[2]), sin(pi*x[3])*cos(pi*x[2])))
g2 = BEAST.ScalarTrace{Float64}(x -> ( sin(pi*x[1])*cos(pi*x[3]),  0, -cos(pi*x[1])*sin(pi*x[3])))


N = 4

h = Vector(undef, N)
errorL = Vector(undef, N)
errorLnew = Vector(undef, N)
errorBDM = Vector(undef, N)
errorRT = Vector(undef, N)

for n in 1:N

    h[n] = 0.5^(n)

    trias = CompScienceMeshes.meshrectangle(1.0,1.0, h[n])
    
    # trias = CompScienceMeshes.meshcircle(1.0, h[n], 3)

    Z = BEAST.lagrangec0d1(trias, dirichlet=false)
    Y = BEAST.brezzidouglasmarini(trias)
    X = BEAST.raviartthomas(trias)

    u_exact = DofInterpolate(Z,f)
    b_exact = DofInterpolate(Y,f2)
    r_exact = DofInterpolate(X,f2)

    identity = BEAST.Identity()
    
    M = assemble(identity, Z, Z)
    b = assemble(f,Z)
    u_n = M \ b

    errorL[n] = sqrt((u_n-u_exact)'*M*(u_n-u_exact))

 
    M = assemble(identity, Y, Y)
    b = assemble(f2,Y)
    b_n = M \ b

    errorBDM[n] = sqrt((b_n-b_exact)'*M*(b_n-b_exact))


    M = assemble(identity, X, X)
    b = assemble(f2,X)
    r_n = M \ b

    errorRT[n] = sqrt((r_n-r_exact)'*M*(r_n-r_exact))

end

N = 15

delta = Vector(undef, N)
error = Vector(undef, N)


for i in 1:N

    delta[i] = 0.5-0.025*i

    tetrs =  tetmeshcuboid(1.0,1.0,1.0, delta[i])

    bndry = boundary(tetrs)
    faces = skeleton(tetrs, 2)

    # onbnd1 = CompScienceMeshes.in(bnd_edges)
    # onbnd0 = CompScienceMeshes.in(bnd_verts)

    # interior_edges = submesh(!onbnd1, all_edges)
    # interior_verts = submesh(!onbnd0, all_verts)

    # bndry_faces = [sort(c) for c in cells(skeleton(bndry, 2))]
    # function is_interior(faces)
    #     !(sort(faces) in bndry_faces)
    # end

    onbnd = CompScienceMeshes.in(bndry)
    interior_faces = submesh(!onbnd, faces)

    # interior_faces = submesh(is_interior, faces)

    W = BEAST.brezzidouglasmarini3d(tetrs, interior_faces)

    u_exact = DofInterpolate(W,g2)

    identity = BEAST.Identity()

    M = assemble(identity, W, W)

    b = assemble(g2,W)
    u_n = M \ b

    error[i] = real(sqrt((u_n-u_exact)'*M*(u_n-u_exact)))
end


N=15

delta2 = Vector(undef, N)
error2 = Vector(undef, N)

for i in 1:N

    delta2[i] =0.5-0.025*i

    tetrs =  tetmeshcuboid(1.0,1.0,1.0, delta2[i];generator = :gmsh)

    bndry = boundary(tetrs)
    faces = skeleton(tetrs, 2)

    # bndry_faces = [sort(c) for c in cells(skeleton(bndry, 2))]
    # function is_interior(faces)
    #     !(sort(faces) in bndry_faces)
    # end

    # interior_faces = submesh(is_interior, faces)

    onbnd = CompScienceMeshes.in(bndry)
    interior_faces = submesh(!onbnd, faces)

    V = BEAST.nedelecd3d(tetrs, interior_faces)
    

  
    q_exact = DofInterpolate(V,g)

    identity = BEAST.Identity()

    M2 = assemble(identity, V, V)

   
    b2 = assemble(g,V)
 
    q_n = M2 \ b2

   
    error2[i] = real(sqrt((q_n-q_exact)'*M2*(q_n-q_exact)))
end

using Plots

p2 = Plots.plot(delta,error, label="BDM3D", title="L2-Error on [0,1]^3",xlabel="h", ylabel="Error", legend=:bottomright, xaxis=:log, yaxis=:log)
Plots.plot!(p2, delta2, error2, label="Nedelec Div", xlims=[0.075,0.55])

p1 = Plots.plot(h,real(errorL), label="Lagrange", title="L2-Error on [0,1]^2",xlabel="h", ylabel="Error", legend=:bottomright, xaxis=:log, yaxis=:log, xlims=[0.025,0.5])
Plots.plot!(p1,h,real(errorBDM),label="BDM")
Plots.plot!(p1,h,real(errorRT),label="Raviart-Thomas")

Plots.plot(p1,p2, layout=2)

