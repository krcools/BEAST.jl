BE = BEAST

T = Float64
tol = eps(T) * 10^3

# Test on cubic polynomial in two variables
c = T[7,-12,3,2,-1,10,-3,2,8,-7];
p = BE.Poly{2,T}(c);
x = [2.0, -3.0];
y = BE.polyval(p,x);
u = x[1]; v = x[2];
y2 = 7 - 12u + 3v + 2(u^2) -u*v +10(v^2) -3(u^3) +2(u^2*v) + 8(u*v^2) -7(v^3)

@test abs(y2-y) < tol
