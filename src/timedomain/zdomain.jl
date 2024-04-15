

"""
    laplace_to_z(rho, n, N, dt, A, b)

Returns the complex matrix valued Laplace variable s that correspond to the
variable z = rho*exp(2*im*pi*n/N) for a given Butcher tableau (A,b,c) and a time step dt.
"""
function laplace_to_z(rho, n, N, dt, A, b)
	z = rho * exp(2*im*pi*n/N)
	s = inv(dt * (A + ones(size(b)) * b' / (z-1)))
	return s
end

"""
    inverse_z_transform(k, rho, N, X)

Returns the k-th term of the inverse z-transform. X is an array of the z-transform
evaluated in the points z=rho*exp(2*im*pi*n/N) for n in 0:(N-1).
"""
function inverse_z_transform(k, rho, N, X::AbstractArray{T,1}) where T
	return ((rho^k) / N) * sum(n -> X[n+1] * exp(2*im*pi*k*n/N), 0:(N-1))
end

"""
    real_inverse_z_transform(k, rho, N, X)

Returns the k-th term of the inverse z-transform.
It is assumed that X[n+1] = conj(X[N-n]) for each n in 1:(N-1)
so that Nmax = N/2+1 or (N+1)/2 (resp. if N%2==0 or N%2==1) terms are used in X
X is an array of the z-transform
evaluated in the points z=rho*exp(2*im*pi*n/N) for n in 0:(Nmax-1).
"""
function real_inverse_z_transform(k, rho, N, X::AbstractArray{T,1}) where T
	Nmax = (N+1)>>1
	realTerms = (N%2==0) ? real(X[1]) + (-1)^k * real(X[Nmax+1]) : real(X[1])
	return ((rho^k) / N) * (realTerms + 2*sum(n -> real(X[n+1] * exp(2*im*pi*k*n/N)), 1:(Nmax-1)))
end
