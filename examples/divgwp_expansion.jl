using BEAST
using CompScienceMeshes
const CSM = CompScienceMeshes
using  StaticArrays

T = Float64
D = 4
NF = binomial(2+D, 2)
gwp = BEAST.GWPDivRefSpace{T,D}()
lgx = BEAST.LagrangeRefSpace{T,D,3,10}()

function fields(p)
    map(gwp(p)) do x
        x.divergence
    end
end

ch = CSM.simplex(
    point(1,0,0),
    point(0,1,0),
    point(0,0,0))

# p = CSM.center(ch)
# v = fields(p)

coeffs = BEAST.interpolate(fields, lgx, ch)
Q = rationalize.(coeffs; tol=sqrt(eps(T)))

using LinearAlgebra
norm(Q - coeffs)

Q3 = Rational{Int64}[89//3 223//81 -7//81 -2 223//81 -7//81 -2 -7//81 -2 -2; -59//2 236//27 -13//54 13 -83//9 -4//27 8 -1//18 3 -2; 13 -13//54 236//27 -59//2 8 -4//27 -83//9 3 -1//18 -2; -2 -7//81 223//81 89//3 -2 -7//81 223//81 -2 -7//81 -2; 89//3 223//81 -7//81 -2 223//81 -7//81 -2 -7//81 -2 -2; -59//2 -83//9 -1//18 -2 236//27 -4//27 3 -13//54 8 13; 13 8 3 -2 -13//54 -4//27 -1//18 236//27 -83//9 -59//2; -2 -2 -2 -2 -7//81 -7//81 -7//81 223//81 223//81 89//3; -2 -7//81 223//81 89//3 -2 -7//81 223//81 -2 -7//81 -2; -2 -1//18 -83//9 -59//2 3 -4//27 236//27 8 -13//54 13; -2 3 8 13 -1//18 -4//27 -13//54 -83//9 236//27 -59//2; -2 -2 -2 -2 -7//81 -7//81 -7//81 223//81 223//81 89//3; -50 890//81 55//81 35//3 -565//81 125//162 20//3 70//81 5//3 -10//3; 50 565//81 -70//81 10//3 -890//81 -125//162 -5//3 -55//81 -20//3 -35//3; 30 -355//27 355//27 -30 80//9 0 -80//9 -10//9 10//9 0; -40 440//27 85//27 -10 5//2 -125//27 205//18 35//3 95//9 -25//2; -35//3 -55//81 -890//81 50 -20//3 -125//162 565//81 -5//3 -70//81 10//3; 15 -5//27 485//27 0 5 0 -485//27 -5 5//27 -15; 40 -5//2 -35//3 25//2 -440//27 125//27 -95//9 -85//27 -205//18 10; -30 -80//9 10//9 0 355//27 0 -10//9 -355//27 80//9 30; -25//2 35//3 5//2 -40 95//9 -125//27 440//27 205//18 85//27 -10; 25//2 -95//9 -205//18 10 -35//3 125//27 -85//27 -5//2 -440//27 40; -15 -5 5 15 5//27 0 -5//27 -485//27 485//27 0; 35//3 20//3 5//3 -10//3 55//81 125//162 70//81 890//81 -565//81 -50]
function BEAST.divergence(localspace::BEAST.GWPDivRefSpace, sh, ch)
    BEAST.divergence(localspace, sh, ch, BEAST.dimtype(localspace, domain(ch)))
end

function BEAST.divergence(localspace::BEAST.GWPDivRefSpace{T,D}, cellid, ch, ::Type{Val{N}}) where {N,D}
    function fields(p)
        map(localspace(p)) do x
            x.divergence
        end
    end
    T = coordtype(ch)
    Dim = 2
    NFout = div((Dim+1)*(Dim+2), 2)
    lag = BEAST.LagrangeRefSpace{T,D,Dim+1,NFout}()
    coeffs = BEAST.interpolate(fields, lag, ch)
    S = BEAST.Shape{T}
    A = Vector{Vector{S}}(undef, size(coeffs,1))
    for r in axes(coeffs,1)
        A[r] = collect(S(cellid, c, coeffs[r,c]) for c in axes(coeffs,2))
    end
    return SVector{N}(A)
end

divfns = BEAST.divergence(gwp, 1, ch)

p = neighborhood(ch, (0.2, 0.6))
gwp_vals = gwp(p)
val1 = gwp_vals[1].divergence

lgx_vals = lgx(p)
val2 = zero(T)
for sh in divfns[1]
    val2 += sh.coeff * lgx_vals[sh.refid].value
end