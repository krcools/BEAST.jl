

"""
	fcr, geo = facecurrents(coeffs, basis)

Compute the value of the function with the given collection of coeffient in the provided basis in all the centroids of the mesh underlying the basis. The mesh is returned together with the currents.
"""
function facecurrents(coeffs, basis)

	T = eltype(coeffs)
	RT = real(T)

	refs = refspace(basis)
	numrefs = numfunctions(refs)

	cells, tad, a2g = assemblydata(basis)

	mesh = geometry(basis)
	D = dimension(mesh)
	# U = D+1
	U = 3

	# TODO: remove ugliness
	vals = refs(center(first(cells)))
	PT = typeof(first(coeffs)*vals[1][1])

	fcr = zeros(PT, numcells(mesh))

	for (t,cell) in enumerate(cells)

		mps = neighborhood(cell, ones(RT,D)/(D+1))
        vals = refs(mps)

		# assemble in the right hand side vector
		for i in 1 : numrefs
			fx = vals[i][1]
            for (m,a) in tad[t,i]
				# fcr[t] += coeffs[m] * a * fx
				fcr[a2g[t]] += coeffs[m] * a * fx
			end
		end
	end

	return fcr, geometry(basis)
end

function facecurrents(coeffs, basis::SpaceTimeBasis)

	space_basis = basis.space
	time_basis = basis.time

	Nt = numfunctions(time_basis)
	Δt = timestep(time_basis)

	refs = refspace(space_basis)
	trefs = refspace(time_basis)
	numrefs = numfunctions(refs)
	tnumrefs = numfunctions(trefs)

	cells, ad = assemblydata(space_basis)
	tcells, tad = assemblydata(time_basis)

	mesh = geometry(space_basis)
	T = eltype(coeffs)
	D = dimension(mesh)
	U = D+1

	# TODO: express relative to input types
	PT = SVector{U, T}
	fcr = zeros(PT, numcells(mesh), Nt)

	for (k,tcell) in enumerate(tcells)
		tmps = neighborhood(tcell,1)
		tvals = trefs(tmps)
		for (p,cell) in enumerate(cells)
			mps = center(cell)
	        vals = refs(mps)

			# assemble
			for i in 1:numrefs
				fx = vals[i][1]
				for j in 1:tnumrefs
					tfx = tvals[j]
	    			for (m,a) in ad[p,i]
						for (n,b) in tad[k,j]
							fcr[p,k] += (coeffs[m,n] * tfx * b) * a * fx
						end
					end
				end
			end

		end
	end

	return fcr, geometry(space_basis)
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
	#fcr = append!(fcrs...)
	fcr = vcat(fcrs...)
	m = weld([geometry(x) for x in X.factors]...)

	fcr, m
end

function potential(op, points, coeffs, basis)
	T = SVector{3,ComplexF64}
	ff = zeros(T, size(points))
	store(v,m,n) = (ff[m] += v*coeffs[n])
	potential!(store, op, points, basis)
	return ff
end

function potential(op, points,coeffs, basis::SpaceTimeBasis)
	T = SVector{3,ComplexF64}
	ff = zeros(T, length(points), size(coeffs)[2])
	store(v,m,n,k,o) = (ff[m,k] += v*coeffs[n,o])
	potential!(store, op, points, basis)
	return ff
end

function potential(op, points, coeffs, space::DirectProductSpace)
	T = SVector{3,ComplexF64}
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

	T = SVector{3,ComplexF64}
	z = zeros(T,length(points))

	els, ad = assemblydata(basis)
	rs = refspace(basis)

	zlocal = Array{T}(undef,numfunctions(rs))
	qdata = quaddata(op,rs,els)

	#println("Computing nearfield.")
	print("dots out of 10: ")

	todo, done, pctg = length(points), 0, 0

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
		end

		done += 1
		new_pctg = round(Int, done / todo * 100)
		if new_pctg > pctg + 9
				#println(todo," ",done," ",new_pctg)
				print(".")
				pctg = new_pctg
		end
	end

	println("")

end

function potential!(store, op, points, basis::SpaceTimeBasis)

    space_basis = basis.space
	tb = time_basis = basis.time

	Nt = numfunctions(time_basis)
	Δt = timestep(time_basis)

	T = SVector{3,ComplexF64}

	refs = refspace(space_basis)
	trefs = refspace(time_basis)

	numrefs = numfunctions(refs)
	tnumrefs = numfunctions(trefs)

	els, ad = assemblydata(space_basis)
	tels, tad = assemblydata(time_basis)


	zlocal = Array{T}(numrefs,Nt)
	qdata = quaddata(op,refs,trefs,els,tels)

	#println("Computing nearfield.")
	print("dots out of 10: ")

	todo, done, pctg = Nt, 0, 0
    for (k,tel) in enumerate(tels)
    	for (p,y) in enumerate(points)
    		for (q,el) in enumerate(els)

    			fill!(zlocal,zero(T))
    			qr = quadrule(op,refs,trefs,p,y,q,el,k,tel,qdata)
	    		farfieldlocal!(zlocal,op,refs,trefs,p,y,q,el,k,tel,tb,qr)

    			# assemble from local contributions
    			 for r in 1:numrefs
					 for s in 1:Nt
	                    for (n,b) in ad[q,r]
	        				store(zlocal[r,s]*b,p,n,k,s)
							# end
	    				end
					end
    			end
            end
    	end

		done += 1
		new_pctg = round(Int, done / todo * 100)
		if new_pctg > pctg + 9
				#println(todo," ",done," ",new_pctg)
				print(".")
				pctg = new_pctg
		end
    end

	println("")

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

function farfieldlocal!(zlocal,op,trialrefs, timerefs,
        p, testel, q, trialel, k, timeel, tb,quadrule)

	M,N = size(zlocal)
	Δt = tb.timestep

    for qr in quadrule[1]
        x = qr.point
        F = qr.value
        dx = qr.weight

        krn = kernelvals(op, testel, x)

        for r = 1:M
			for l = 1:N
				t = (k - l)*Δt
            	zlocal[r,l] += integrand(op,krn,testel,F[r], t, tb) * dx
			end
		end

    end
end
