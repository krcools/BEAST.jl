using RungeKutta, LinearAlgebra

tableRadau2a = TableauRadauIIA(3)
A = tableRadau2a.a
b = tableRadau2a.b
inv_A = inv(A)
function s(z, inv_A=inv_A, b=b)
    btinv_A = b'*inv_A
    p1 = ones(size(inv_A,1),1)
    M = inv_A - inv_A*p1*btinv_A./(z-1.0 .+btinv_A*p1)
    return M
    Meig = eigen(M)
    Λ = Meig.values.*I(size(inv_A,1))
    Q = Meig.vectors
    inv_Q = inv(Q)
    return Λ, Q, inv_Q
end