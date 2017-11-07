export laplace_to_z, inverse_z_transform;

"""
    laplace_to_z(rho, n, N, dt, A, b)

Returns the complex matrix valued Laplace variable s that correspond to the
variable z = rho*exp(2*im*pi*n/N) for a given Butcher tableau (A,b,c) and a time step dt.
"""
function laplace_to_z(rho, n, N, dt, A, b)
	z = rho * exp(2*im*pi*n/N);
	s = inv(dt * (A + ones(b) * b' / (z-1)));
	return s;
end

"""
    inverse_z_transform(k, rho, N, X)

Returns the k-th term of the inverse z-transform. X is an array of the z-transform
evaluated in the points z=rho*exp(2*im*pi*n/N) for n in 0:(N-1).
"""
function inverse_z_transform{T}(k, rho, N, X::AbstractArray{T,1})
	return ((rho^k) / N) * sum(n -> X[n+1] * exp(2*im*pi*k*n/N), 0:(N-1));
end
