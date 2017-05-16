export facecurrents
export potential
export get_scatter_parameters

"""
	fcr, geo = facecurrents(coeffs, basis)

Compute the value of the function with the given collection of coeffient in the provided basis in all the centroids of the mesh underlying the basis. The mesh is returned together with the currents.
"""
function facecurrents(coeffs, basis)

	T = eltype(coeffs)
	RT = real(T)

	refs = refspace(basis)
	numrefs = numfunctions(refs)

	cells, tad = assemblydata(basis)

	mesh = geometry(basis)
	D = dimension(mesh)
	U = D+1

	# TODO: express relative to input types
	PT = SVector{U, T}
	fcr = zeros(PT, numcells(mesh))

	for (t,cell) in enumerate(cells)

		mps = neighborhood(cell, ones(RT,D)/(D+1))
        vals = refs(mps)

		# assemble in the right hand side vector
		for i in 1 : numrefs
			fx = vals[i][1]
            for (m,a) in tad[t,i]
				fcr[t] += coeffs[m] * a * fx
			end
		end
	end

	return fcr, geometry(basis)
end


function facecurrents(u, X1, Xs...)

    offset = 0

    n = numfunctions(X1)
    fcrs = [ facecurrents(u[offset + (1:n)], X1)[1] ]
    offset += n

    for i in 1:length(Xs)
        n = numfunctions(Xs[i])
        push!(fcrs, facecurrents(u[offset + (1:n)], Xs[i])[1])
        offset += n
    end

    fcrs
end


function facecurrents(u, X::DirectProductSpace)

	if length(X.factors) == 1
		return facecurrents(u, X.factors[1]), geometry(X.factors[1])[1]
	end

	fcrs = facecurrents(u, X.factors...)
	fcr = append!(fcrs...)
	m = weld([geometry(x) for x in X.factors]...)

	fcr, m
end

function potential(op, points, coeffs, basis)
	T = SVector{3,Complex128}
	ff = zeros(T, length(points))
	store(v,m,n) = (ff[m] += v*coeffs[n])
	potential!(store, op, points, basis)
	return ff
end

function potential(op, points, coeffs, space::DirectProductSpace)
	T = SVector{3,Complex128}
	ff = zeros(T, length(points))

	@assert length(coeffs) == numfunctions(space)

	offset = 0
	for fct in space.factors
		store(v,m,n) = (ff[m] += v*coeffs[offset+n])
		potential!(store, op, points, fct)
		offset += numfunctions(fct)
	end

	ff
end

function potential!(store, op, points, basis)

	T = SVector{3,Complex128}
	z = zeros(T,length(points))

	els, ad = assemblydata(basis)
	rs = refspace(basis)

	zlocal = Array{T}(numfunctions(rs))
	qdata = quaddata(op,rs,els)

	print("Computing nearfield.")
	print("dots out of 10: ")

	for (p,y) in enumerate(points)
		for (q,el) in enumerate(els)

			fill!(zlocal,zero(T))
			qr = quadrule(op,rs,p,y,q,el,qdata)
			farfieldlocal!(zlocal,op,rs,y,el,qr)

			# assemble from local contributions
			for (r,z) in enumerate(zlocal)
                for (n,b) in ad[q,r]
					store(z*b,p,n)
				end
			end

			done += 1
			new_pctg = round(Int, done / todo * 100)
			if new_pctg > pctg + 9
					print(".")
					pctg = new_pctg
			end
		end
	end
end

function farfieldlocal!(zlocal,op,refspace,y,el,qr)

    for q in qr
        x = q.point
        F = q.value
        dx = q.weight

        krn = kernelvals(op, y, x)
        for r in 1 : length(zlocal)
            zlocal[r] += integrand(op,krn,y,F[r],x) * dx
        end

    end

end
