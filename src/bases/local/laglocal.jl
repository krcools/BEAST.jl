# T: coeff type
# Degree: degree
# Dim1: dimension of the support + 1
struct LagrangeRefSpace{T,Degree,Dim1,NF} <: RefSpace{T} end

numfunctions(s::LagrangeRefSpace{T,D,2}, ch::CompScienceMeshes.ReferenceSimplex{1}) where {T,D} = D+1
numfunctions(s::LagrangeRefSpace{T,0,3}, ch::CompScienceMeshes.ReferenceSimplex{2}) where {T} = 1
numfunctions(s::LagrangeRefSpace{T,1,3}, ch::CompScienceMeshes.ReferenceSimplex{2}) where {T} = 3
numfunctions(s::LagrangeRefSpace{T,2,3}, ch::CompScienceMeshes.ReferenceSimplex{2}) where {T} = 6
numfunctions(s::LagrangeRefSpace{T,Dg}, ch::CompScienceMeshes.ReferenceSimplex{D}) where {T,Dg,D} = binomial(D+Dg,Dg)

# valuetype(ref::LagrangeRefSpace{T}, charttype) where {T} =
#         SVector{numfunctions(ref), Tuple{T,T}}
valuetype(ref::LagrangeRefSpace{T}, charttype) where {T} = T

# Evaluate constant lagrange elements on anything
(ϕ::LagrangeRefSpace{T,0})(tp) where {T} = SVector(((value=one(T), derivative=zero(T)),))
(ϕ::LagrangeRefSpace{T,0,3})(tp) where {T} = SVector(((value=one(T), derivative=zero(T)),))


# Evaluate linear Lagrange elements on a segment
# The derivative denotes the tangential derivative
function (f::LagrangeRefSpace{T,1,2})(mp) where T
    u = mp.bary[1]
    j = jacobian(mp)
    SVector(
        (value=  u, derivative=-1/j),
        (value=1-u, derivative= 1/j))
end

function (f::LagrangeRefSpace{T,2,2})(mp) where T
    u = mp.bary[1]
    j = jacobian(mp)
    SVector(
        (value=-u*(1 - 2*u), derivative=-(4*u - 1)/j),
        (value=(1 - 2*u)*(1 - u), derivative=-(4*u - 3)/j),
        (value=4*u*(1 - u), derivative=-(4 - 8*u)/j)
    )
end

function (f::LagrangeRefSpace{T,3,2})(mp) where T
    u = mp.bary[1]
    j = jacobian(mp)
    SVector(
        (value=u*(1 - 3*u)*(2 - 3*u)/2, derivative=-(27*u^2/2 - 9*u + 1)/j),
        (value=(1 - 3*u)*(1 - u)*(2 - 3*u)/2, derivative=(27*u^2/2 - 18*u + 11/2)/j),
        (value=-9*u*(1 - 3*u)*(1 - u)/2, derivative=(81*u^2/2 - 36*u + 9/2)/j),
        (value=9*u*(1 - u)*(2 - 3*u)/2, derivative=-(81*u^2/2 - 45*u + 9)/j)
    )
end

function (f::LagrangeRefSpace{T,4,2})(mp) where T
    u = mp.bary[1]
    j = jacobian(mp)
    SVector(
        (value=-u*(1 - 4*u)*(1 - 2*u)*(3 - 4*u)/3,
        derivative=-(128*u^3/3 - 48*u^2 + 44*u/3 - 1)/j),
        (value=(1 - 4*u)*(1 - 2*u)*(1 - u)*(3 - 4*u)/3,
        derivative=-(128*u^3/3 - 80*u^2 + 140*u/3 - 25/3)/j),
        (value=16*u*(1 - 4*u)*(1 - 2*u)*(1 - u)/3,
        derivative=-( -512*u^3/3 + 224*u^2 - 224*u/3 + 16/3)/j),
        (value=-4*u*(1 - 4*u)*(1 - u)*(3 - 4*u),
        derivative=-(256*u^3 - 384*u^2 + 152*u - 12)/j),
        (value=16*u*(1 - 2*u)*(1 - u)*(3 - 4*u)/3,
        derivative=-( -512*u^3/3 + 288*u^2 - 416*u/3 + 16)/j)
    )
end

function (f::LagrangeRefSpace{T,5,2})(mp) where T
    u = mp.bary[1]
    j = jacobian(mp)
    SVector(
        (value=u*(1 - 5*u)*(2 - 5*u)*(3 - 5*u)*(4 - 5*u)/24,
        derivative=-(3125*u^4/24 - 625*u^3/3 + 875*u^2/8 - 125*u/6 + 1)/j),
        (value=(1 - 5*u)*(1 - u)*(2 - 5*u)*(3 - 5*u)*(4 - 5*u)/24,
        derivative=(3125*u^4/24 - 625*u^3/2 + 2125*u^2/8 - 375*u/4 + 137/12)/j),
        (value=-25*u*(1 - 5*u)*(1 - u)*(2 - 5*u)*(3 - 5*u)/24,
        derivative=(15625*u^4/24 - 6875*u^3/6 + 5125*u^2/8 - 1525*u/12 + 25/4)/j),
        (value=25*u*(1 - 5*u)*(1 - u)*(2 - 5*u)*(4 - 5*u)/12,
        derivative=-(15625*u^4/12 - 2500*u^3 + 6125*u^2/4 - 325*u + 50/3)/j),
        (value=-25*u*(1 - 5*u)*(1 - u)*(3 - 5*u)*(4 - 5*u)/12,
        derivative=(15625*u^4/12 - 8125*u^3/3 + 7375*u^2/4 - 2675*u/6 + 25)/j),
        (value=25*u*(1 - u)*(2 - 5*u)*(3 - 5*u)*(4 - 5*u)/24,
        derivative=-(15625*u^4/24 - 4375*u^3/3 + 8875*u^2/8 - 1925*u/6 + 25)/j)
    )
end

function (f::LagrangeRefSpace{T,6,2})(mp) where T
    u = mp.bary[1]
    j = jacobian(mp)
    SVector(
        (value=-u*(1 - 6*u)*(1 - 3*u)*(1 - 2*u)*(2 - 3*u)*(5 - 6*u)/10,
        derivative=-(1944*u^5/5 - 810*u^4 + 612*u^3 - 405*u^2/2 + 137*u/5 - 1)/j),
        (value=(1 - 6*u)*(1 - 3*u)*(1 - 2*u)*(1 - u)*(2 - 3*u)*(5 - 6*u)/10,
        derivative=(1944*u^5/5 - 1134*u^4 + 1260*u^3 - 1323*u^2/2 + 812*u/5 - 147/10)/j),
        (value=18*u*(1 - 6*u)*(1 - 3*u)*(1 - 2*u)*(1 - u)*(2 - 3*u)/5,
        derivative=-( -11664*u^5/5 + 5184*u^4 - 4104*u^3 + 1404*u^2 - 972*u/5 + 36/5)/j),
        (value=-9*u*(1 - 6*u)*(1 - 3*u)*(1 - 2*u)*(1 - u)*(5 - 6*u)/2,
        derivative=-(5832*u^5 - 13770*u^4 + 11556*u^3 - 8289*u^2/2 + 594*u - 45/2)/j),
        (value=4*u*(1 - 6*u)*(1 - 3*u)*(1 - u)*(2 - 3*u)*(5 - 6*u),
        derivative=-( -7776*u^5 + 19440*u^4 - 17424*u^3 + 6696*u^2 - 1016*u + 40)/j),
        (value=-9*u*(1 - 6*u)*(1 - 2*u)*(1 - u)*(2 - 3*u)*(5 - 6*u)/2,
        derivative=-(5832*u^5 - 15390*u^4 + 14796*u^3 - 12447*u^2/2 + 1053*u - 45)/j),
        (value=18*u*(1 - 3*u)*(1 - 2*u)*(1 - u)*(2 - 3*u)*(5 - 6*u)/5,
        derivative=-( -11664*u^5/5 + 6480*u^4 - 6696*u^3 + 3132*u^2 - 3132*u/5 + 36)/j)
    )
end

function (f::LagrangeRefSpace{T,7,2})(mp) where T
    u = mp.bary[1]
    j = jacobian(mp)
    SVector(
        (value=u*(7*u - 6)*(7*u - 5)*(7*u - 4)*(7*u - 3)*(7*u - 2)*(7*u - 1)/720,
        derivative=-(823543*u^6/720 - 117649*u^5/40 + 420175*u^4/144 - 16807*u^3/12 + 9947*u^2/30 - 343*u/10 + 1)/j),
        (value=-(u - 1)*(7*u - 6)*(7*u - 5)*(7*u - 4)*(7*u - 3)*(7*u - 2)*(7*u - 1)/720,
        derivative=(823543*u^6/720 - 117649*u^5/30 + 386561*u^4/72 - 33614*u^3/9 + 331681*u^2/240 - 22981*u/90 + 363/20)/j),
        (value=-49*u*(u - 1)*(7*u - 5)*(7*u - 4)*(7*u - 3)*(7*u - 2)*(7*u - 1)/720,
        derivative=(5764801*u^6/720 - 1294139*u^5/60 + 1596665*u^4/72 - 98441*u^3/9 + 634207*u^2/240 - 49931*u/180 + 49/6)/j),
        (value=49*u*(u - 1)*(7*u - 6)*(7*u - 4)*(7*u - 3)*(7*u - 2)*(7*u - 1)/240,
        derivative=-(5764801*u^6/240 - 2705927*u^5/40 + 1159683*u^4/16 - 444185*u^3/12 + 45962*u^2/5 - 9849*u/10 + 147/5)/j),
        (value=-49*u*(u - 1)*(7*u - 6)*(7*u - 5)*(7*u - 3)*(7*u - 2)*(7*u - 1)/144,
        derivative=(5764801*u^6/144 - 117649*u^5 + 9495955*u^4/72 - 211288*u^3/3 + 872935*u^2/48 - 2009*u + 245/4)/j),
        (value=49*u*(u - 1)*(7*u - 6)*(7*u - 5)*(7*u - 4)*(7*u - 2)*(7*u - 1)/144,
        derivative=-(5764801*u^6/144 - 2941225*u^5/24 + 20756645*u^4/144 - 2926819*u^3/36 + 133427*u^2/6 - 46501*u/18 + 245/3)/j),
        (value=-49*u*(u - 1)*(7*u - 6)*(7*u - 5)*(7*u - 4)*(7*u - 3)*(7*u - 1)/240,
        derivative=(5764801*u^6/240 - 1529437*u^5/20 + 756315*u^4/8 - 170471*u^3/3 + 1347647*u^2/80 - 43071*u/20 + 147/2)/j),
        (value=49*u*(u - 1)*(7*u - 6)*(7*u - 5)*(7*u - 4)*(7*u - 3)*(7*u - 2)/720,
        derivative=-(5764801*u^6/720 - 1058841*u^5/40 + 4958065*u^4/144 - 88837*u^3/4 + 109417*u^2/15 - 10927*u/10 + 49)/j)
    )
end

function (f::LagrangeRefSpace{T,8,2})(mp) where T
    u = mp.bary[1]
    j = jacobian(mp)
    SVector(
        (value=u*(2*u - 1)*(4*u - 3)*(4*u - 1)*(8*u - 7)*(8*u - 5)*(8*u - 3)*(8*u - 1)/315,
        derivative=-(1048576*u^7/315 - 458752*u^6/45 + 188416*u^5/15 - 71680*u^4/9 + 123776*u^3/45 - 7504*u^2/15 + 1452*u/35 - 1)/j),
        (value=(u - 1)*(2*u - 1)*(4*u - 3)*(4*u - 1)*(8*u - 7)*(8*u - 5)*(8*u - 3)*(8*u - 1)/315,
        derivative=-(1048576*u^7/315 - 65536*u^6/5 + 106496*u^5/5 - 18432*u^4 + 136832*u^3/15 - 12816*u^2/5 + 118124*u/315 - 761/35)/j),
        (value=-64*u*(u - 1)*(2*u - 1)*(4*u - 3)*(4*u - 1)*(8*u - 5)*(8*u - 3)*(8*u - 1)/315,
        derivative=(8388608*u^7/315 - 3801088*u^6/45 + 1605632*u^5/15 - 624640*u^4/9 + 1097728*u^3/45 - 67456*u^2/15 + 13184*u/35 - 64/7)/j),
        (value=16*u*(u - 1)*(2*u - 1)*(4*u - 1)*(8*u - 7)*(8*u - 5)*(8*u - 3)*(8*u - 1)/45,
        derivative=-(4194304*u^7/45 - 917504*u^6/3 + 1998848*u^5/5 - 266240*u^4 + 1435136*u^3/15 - 17952*u^2 + 68576*u/45 - 112/3)/j),
        (value=-64*u*(u - 1)*(2*u - 1)*(4*u - 3)*(4*u - 1)*(8*u - 7)*(8*u - 3)*(8*u - 1)/45,
        derivative=(8388608*u^7/45 - 28442624*u^6/45 + 12812288*u^5/15 - 5285888*u^4/9 + 9773056*u^3/45 - 626048*u^2/15 + 18048*u/5 - 448/5)/j),
        (value=4*u*(u - 1)*(4*u - 3)*(4*u - 1)*(8*u - 7)*(8*u - 5)*(8*u - 3)*(8*u - 1)/9,
        derivative=-(2097152*u^7/9 - 7340032*u^6/9 + 3424256*u^5/3 - 7331840*u^4/9 + 2814208*u^3/9 - 186496*u^2/3 + 5528*u - 140)/j),
        (value=-64*u*(u - 1)*(2*u - 1)*(4*u - 3)*(4*u - 1)*(8*u - 7)*(8*u - 5)*(8*u - 1)/45,
        derivative=(8388608*u^7/45 - 10092544*u^6/15 + 4882432*u^5/5 - 727040*u^4 + 4390912*u^3/15 - 306048*u^2/5 + 256384*u/45 - 448/3)/j),
        (value=16*u*(u - 1)*(2*u - 1)*(4*u - 3)*(8*u - 7)*(8*u - 5)*(8*u - 3)*(8*u - 1)/45,
        derivative=-(4194304*u^7/45 - 15597568*u^6/45 + 7831552*u^5/15 - 3665920*u^4/9 + 7827968*u^3/45 - 587296*u^2/15 + 19872*u/5 - 112)/j),
        (value=-64*u*(u - 1)*(2*u - 1)*(4*u - 3)*(4*u - 1)*(8*u - 7)*(8*u - 5)*(8*u - 3)/315,
        derivative=(8388608*u^7/315 - 917504*u^6/9 + 2392064*u^5/15 - 1177600*u^4/9 + 2695168*u^3/45 - 44672*u^2/3 + 61568*u/35 - 64)/j)
    )
end

function (f::LagrangeRefSpace{T,9,2})(mp) where T
    u = mp.bary[1]
    j = jacobian(mp)
    SVector(
        (value=u*(3*u - 2)*(3*u - 1)*(9*u - 8)*(9*u - 7)*(9*u - 5)*(9*u - 4)*(9*u - 2)*(9*u - 1)/4480,
        derivative=-(43046721*u^8/4480 - 4782969*u^7/140 + 16120377*u^6/320 - 1594323*u^5/40 + 2337903*u^4/128 - 194643*u^3/40 + 797337*u^2/1120 - 6849*u/140 + 1)/j),
        (value=-(u - 1)*(3*u - 2)*(3*u - 1)*(9*u - 8)*(9*u - 7)*(9*u - 5)*(9*u - 4)*(9*u - 2)*(9*u - 1)/4480,
        derivative=(43046721*u^8/4480 - 4782969*u^7/112 + 5137263*u^6/64 - 2657205*u^5/32 + 6589431*u^4/128 - 623295*u^3/32 + 122121*u^2/28 - 58635*u/112 + 7129/280)/j),
        (value=-81*u*(u - 1)*(3*u - 2)*(3*u - 1)*(9*u - 7)*(9*u - 5)*(9*u - 4)*(9*u - 2)*(9*u - 1)/4480,
        derivative=(387420489*u^8/4480 - 176969853*u^7/560 + 152523567*u^6/320 - 61470009*u^5/160 + 22878207*u^4/128 - 7712091*u^3/160 + 3986901*u^2/560 - 275967*u/560 + 81/8)/j),
        (value=81*u*(u - 1)*(3*u - 2)*(3*u - 1)*(9*u - 8)*(9*u - 5)*(9*u - 4)*(9*u - 2)*(9*u - 1)/1120,
        derivative=-(387420489*u^8/1120 - 90876411*u^7/70 + 80247591*u^6/40 - 66075831*u^5/40 + 25043337*u^4/32 - 2142531*u^3/10 + 8968887*u^2/280 - 78327*u/35 + 324/7)/j),
        (value=-9*u*(u - 1)*(3*u - 1)*(9*u - 8)*(9*u - 7)*(9*u - 5)*(9*u - 4)*(9*u - 2)*(9*u - 1)/160,
        derivative=(129140163*u^8/160 - 62178597*u^7/20 + 197164611*u^6/40 - 166341033*u^5/40 + 64448703*u^4/32 - 22480173*u^3/40 + 6828867*u^2/80 - 60381*u/10 + 126)/j),
        (value=81*u*(u - 1)*(3*u - 2)*(3*u - 1)*(9*u - 8)*(9*u - 7)*(9*u - 4)*(9*u - 2)*(9*u - 1)/320,
        derivative=-(387420489*u^8/320 - 4782969*u^7 + 249245829*u^6/32 - 54029835*u^5/8 + 215023653*u^4/64 - 3844017*u^3/4 + 2386017*u^2/16 - 21465*u/2 + 1134/5)/j),
        (value=-81*u*(u - 1)*(3*u - 2)*(3*u - 1)*(9*u - 8)*(9*u - 7)*(9*u - 5)*(9*u - 2)*(9*u - 1)/320,
        derivative=(387420489*u^8/320 - 196101729*u^7/40 + 1313190711*u^6/160 - 586888011*u^5/80 + 241241409*u^4/64 - 89119521*u^3/80 + 14257053*u^2/80 - 526419*u/40 + 567/2)/j),
        (value=9*u*(u - 1)*(3*u - 2)*(9*u - 8)*(9*u - 7)*(9*u - 5)*(9*u - 4)*(9*u - 2)*(9*u - 1)/160,
        derivative=-(129140163*u^8/160 - 33480783*u^7/10 + 115322697*u^6/20 - 213107841*u^5/40 + 91020753*u^4/32 - 8776431*u^3/10 + 5878089*u^2/40 - 56601*u/5 + 252)/j),
        (value=-81*u*(u - 1)*(3*u - 2)*(3*u - 1)*(9*u - 8)*(9*u - 7)*(9*u - 5)*(9*u - 4)*(9*u - 1)/1120,
        derivative=(387420489*u^8/1120 - 205667667*u^7/140 + 26040609*u^6/10 - 99733761*u^5/40 + 44529507*u^4/32 - 18152829*u^3/40 + 45570519*u^2/560 - 475389*u/70 + 162)/j),
        (value=81*u*(u - 1)*(3*u - 2)*(3*u - 1)*(9*u - 8)*(9*u - 7)*(9*u - 5)*(9*u - 4)*(9*u - 2)/4480,
        derivative=-(387420489*u^8/4480 - 52612659*u^7/140 + 219485133*u^6/320 - 13640319*u^5/20 + 51221727*u^4/128 - 5589243*u^3/40 + 30921993*u^2/1120 - 373329*u/140 + 81)/j)
    )
end

function (f::LagrangeRefSpace{T,10,2})(mp) where T
    u = mp.bary[1]
    j = jacobian(mp)
    SVector(
        (value=u*(2*u - 1)*(5*u - 4)*(5*u - 3)*(5*u - 2)*(5*u - 1)*(10*u - 9)*(10*u - 7)*(10*u - 3)*(10*u - 1)/4536,
        derivative=-(15625000*u^9/567 - 781250*u^8/7 + 36250000*u^7/189 - 546875*u^6/3 + 1883125*u^5/18 - 296875*u^4/8 + 4523000*u^3/567 - 162875*u^2/168 + 7129*u/126 - 1)/j),
        (value=(u - 1)*(2*u - 1)*(5*u - 4)*(5*u - 3)*(5*u - 2)*(5*u - 1)*(10*u - 9)*(10*u - 7)*(10*u - 3)*(10*u - 1)/4536,
        derivative=-(15625000*u^9/567 - 8593750*u^8/63 + 55000000*u^7/189 - 9453125*u^6/27 + 4695625*u^5/18 - 26846875*u^4/216 + 42711625*u^3/1134 - 10511875*u^2/1512 + 177133*u/252 - 7381/252)/j),
        (value=-25*u*(u - 1)*(2*u - 1)*(5*u - 4)*(5*u - 3)*(5*u - 2)*(5*u - 1)*(10*u - 7)*(10*u - 3)*(10*u - 1)/1134,
        derivative=(156250000*u^9/567 - 71875000*u^8/63 + 377500000*u^7/189 - 52062500*u^6/27 + 10090625*u^5/9 - 21709375*u^4/54 + 49435250*u^3/567 - 4033825*u^2/378 + 13150*u/21 - 100/9)/j),
        (value=25*u*(u - 1)*(2*u - 1)*(5*u - 3)*(5*u - 2)*(5*u - 1)*(10*u - 9)*(10*u - 7)*(10*u - 3)*(10*u - 1)/504,
        derivative=-(78125000*u^9/63 - 36718750*u^8/7 + 590000000*u^7/63 - 82796875*u^6/9 + 32584375*u^5/6 - 142028125*u^4/72 + 54486625*u^3/126 - 2990025*u^2/56 + 88325*u/28 - 225/4)/j),
        (value=-50*u*(u - 1)*(2*u - 1)*(5*u - 4)*(5*u - 3)*(5*u - 2)*(5*u - 1)*(10*u - 9)*(10*u - 3)*(10*u - 1)/189,
        derivative=(625000000*u^9/189 - 100000000*u^8/7 + 1640000000*u^7/63 - 234625000*u^6/9 + 15662500*u^5 - 52006250*u^4/9 + 242639000*u^3/189 - 1121950*u^2/7 + 200600*u/21 - 1200/7)/j),
        (value=25*u*(u - 1)*(2*u - 1)*(5*u - 4)*(5*u - 2)*(5*u - 1)*(10*u - 9)*(10*u - 7)*(10*u - 3)*(10*u - 1)/108,
        derivative=-(156250000*u^9/27 - 76562500*u^8/3 + 47500000*u^7 - 437281250*u^6/9 + 89384375*u^5/3 - 403334375*u^4/36 + 68357750*u^3/27 - 11544725*u^2/36 + 174025*u/9 - 350)/j),
        (value=-u*(u - 1)*(5*u - 4)*(5*u - 3)*(5*u - 2)*(5*u - 1)*(10*u - 9)*(10*u - 7)*(10*u - 3)*(10*u - 1)/9,
        derivative=(62500000*u^9/9 - 31250000*u^8 + 535000000*u^7/9 - 560000000*u^6/9 + 117216250*u^5/3 - 135371875*u^4/9 + 31274500*u^3/9 - 448875*u^2 + 27508*u - 504)/j),
        (value=25*u*(u - 1)*(2*u - 1)*(5*u - 4)*(5*u - 3)*(5*u - 1)*(10*u - 9)*(10*u - 7)*(10*u - 3)*(10*u - 1)/108,
        derivative=-(156250000*u^9/27 - 26562500*u^8 + 155000000*u^7/3 - 498968750*u^6/9 + 107321875*u^5/3 - 510353125*u^4/36 + 91073375*u^3/27 - 1792225*u^2/4 + 168775*u/6 - 525)/j),
        (value=-50*u*(u - 1)*(2*u - 1)*(5*u - 4)*(5*u - 3)*(5*u - 2)*(5*u - 1)*(10*u - 9)*(10*u - 7)*(10*u - 1)/189,
        derivative=(625000000*u^9/189 - 325000000*u^8/21 + 1940000000*u^7/63 - 305375000*u^6/9 + 67737500*u^5/3 - 83431250*u^4/9 + 433739000*u^3/189 - 20028950*u^2/63 + 1308200*u/63 - 400)/j),
        (value=25*u*(u - 1)*(2*u - 1)*(5*u - 4)*(5*u - 3)*(5*u - 2)*(10*u - 9)*(10*u - 7)*(10*u - 3)*(10*u - 1)/504,
        derivative=-(78125000*u^9/63 - 41406250*u^8/7 + 758750000*u^7/63 - 122828125*u^6/9 + 56396875*u^5/6 - 289909375*u^4/72 + 66191750*u^3/63 - 8694225*u^2/56 + 153025*u/14 - 225)/j),
        (value=-25*u*(u - 1)*(2*u - 1)*(5*u - 4)*(5*u - 3)*(5*u - 2)*(5*u - 1)*(10*u - 9)*(10*u - 7)*(10*u - 3)/1134,
        derivative=(156250000*u^9/567 - 9375000*u^8/7 + 527500000*u^7/189 - 29312500*u^6/9 + 20965625*u^5/9 - 18878125*u^4/18 + 165985250*u^3/567 - 1997825*u^2/42 + 243050*u/63 - 100)/j),
    )
end

# Evaluete linear lagrange elements on a triangle
function (f::LagrangeRefSpace{T,1,3})(t) where T
    u,v,w, = barycentric(t)

    j = jacobian(t)
    p = t.patch
    σ = sign(dot(normal(t), cross(p[1]-p[3],p[2]-p[3])))
    SVector(
        (value=u, curl=σ*(p[3]-p[2])/j),
        (value=v, curl=σ*(p[1]-p[3])/j),
        (value=w, curl=σ*(p[2]-p[1])/j))
end



# Evaluate constant Lagrange elements on a triangle, with their curls
function (f::LagrangeRefSpace{T,0,3})(t, ::Type{Val{:withcurl}}) where T
    i = one(T)
    z = zero(cartesian(t))
    SVector(((value=i, curl=z,),))
end

#= Replaced by generalized curl using GWP functions
function curl(ref::LagrangeRefSpace{T,1,3} where {T}, sh, el)
    sh1 = Shape(sh.cellid, mod1(sh.refid+1,3), -sh.coeff)
    sh2 = Shape(sh.cellid, mod1(sh.refid+2,3), +sh.coeff)
    return [sh1, sh2]
end

function curl(ref::LagrangeRefSpace{T,2,3}, sh, el) where T

    j = 1.0 #volume(el) * factorial(dimension(el))

    if sh.refid < 4
        sh1 = Shape(sh.cellid, mod1(2*sh.refid+1,6), -sh.coeff*j)
        sh2 = Shape(sh.cellid, mod1(2*sh.refid+2,6), 3*sh.coeff*j)
        sh3 = Shape(sh.cellid, mod1(2*sh.refid+3,6), -3*sh.coeff*j)
        sh4 = Shape(sh.cellid, mod1(2*sh.refid+4,6), sh.coeff*j)

        return [sh1, sh2, sh3, sh4]
    else 
        sh1 = Shape(sh.cellid, 2*mod1(sh.refid,3)-1, 4*sh.coeff*j)
        sh2 = Shape(sh.cellid, 2*mod1(sh.refid,3), -4*sh.coeff*j)

        return [sh1, sh2]
    end 

end
=#

function gradient(ref::LagrangeRefSpace{T,1,4}, sh, tet) where {T}
    this_vert = tet.vertices[sh.refid]
    # other_verts = deleteat(tet.vertices, sh.refid)
    # opp_face = simplex(other_verts...)
    opp_face = faces(tet)[sh.refid]
    ctr_opp_face = center(opp_face)
    n = normal(ctr_opp_face)
    h = -dot(this_vert - cartesian(ctr_opp_face), n)
    @assert h > 0
    gradval = -(1/h)*n
    output = Vector{Shape{T}}()
    for (i,edge) in enumerate(CompScienceMeshes.edges(tet))
        ctr_edge = center(edge)
        tgt = tangents(ctr_edge,1)
        tgt = normalize(tgt)
        lgt = volume(edge)
        cff = -lgt * dot(tgt, gradval)
        isapprox(cff, 0, atol=sqrt(eps(T))) && continue
        push!(output, Shape(sh.cellid, i, sh.coeff * cff))
    end
    return output
end

function gradient(ref::LagrangeRefSpace{T,1,3} where {T}, sh, el)
    sh1 = Shape(sh.cellid, mod1(sh.refid+1,3), +sh.coeff)
    sh2 = Shape(sh.cellid, mod1(sh.refid+2,3), -sh.coeff)
    return [sh1, sh2]
end

function gradient(ref::LagrangeRefSpace{T,1,2}, sh, seg) where {T}

    sh.refid == 1 && return [Shape(sh.cellid, 1, +sh.coeff/volume(seg))]
    @assert sh.refid == 2
    return [Shape(sh.cellid, 1, -sh.coeff/volume(seg))]

end



function strace(x::LagrangeRefSpace, cell, localid, face)

    Q = zeros(scalartype(x),2,3)

    p1 = neighborhood(face, 1)
    p2 = neighborhood(face, 0)

    u1 = carttobary(cell, cartesian(p1))
    u2 = carttobary(cell, cartesian(p2))

    P1 = neighborhood(cell, u1)
    P2 = neighborhood(cell, u2)

    vals1 = x(P1)
    vals2 = x(P2)

    num_shapes = numfunctions(x, domain(cell))
    for j in 1:num_shapes
        Q[1,j] = vals1[j].value
        Q[2,j] = vals2[j].value
    end

    Q
end

function strace(x::LagrangeRefSpace{T, 1, 4, 4}, cell, localid, face) where {T}

    #T = scalartype(x)
    t = zeros(T, 3, 4)
    for (k,fvert) in enumerate(face.vertices)
        for (l,cvert) in enumerate(cell.vertices)
            nrm = norm(fvert - cvert)
            if isapprox(nrm, 0, atol=sqrt(eps(T)))
                t[k,l] = T(1.0)
                break
            end
        end
    end

    return t
end


function restrict(refs::LagrangeRefSpace{T,0}, dom1, dom2) where T
    n = numfunctions(refs, domain(dom1))
    Q = Matrix{T}(I, n, n)
end

function restrict(f::LagrangeRefSpace{T,1}, dom1, dom2) where T

    D = numfunctions(f, domain(dom1))
    Q = zeros(T, D, D)

    # for each point of the new domain
    for i in 1:D
        v = dom2.vertices[i]

        # find the barycentric coordinates in dom1
        uvn = carttobary(dom1, v)

        # evaluate the shape functions in this point
        x = neighborhood(dom1, uvn)
        fx = f(x)

        for j in 1:D
            Q[j,i] = fx[j][1]
        end
    end

    return Q
end



## Quadratic Lagrange element on a triangle
function (f::LagrangeRefSpace{T,2,3})(t) where T
    u,v,w, = barycentric(t)

    j = jacobian(t)
    p = t.patch


    σ = sign(dot(normal(t), cross(p[1]-p[3],p[2]-p[3])))
     SVector(
        (value=u*(2*u-1), curl=σ*(p[3]-p[2])*(4u-1)/j),
        (value=4*u*v, curl=4*σ*(u*(p[1]-p[3])+v*(p[3]-p[2]))/j),
        (value=v*(2*v-1), curl=σ*(p[1]-p[3])*(4v-1)/j),
        (value=4*w*u, curl=4*σ*(w*(p[3]-p[2])+u*(p[2]-p[1]))/j),
        (value=4*v*w, curl=4*σ*(w*(p[1]-p[3])+v*(p[2]-p[1]))/j),
        (value=w*(2*w-1), curl=σ*(p[2]-p[1])*(4w-1)/j),
    )
end

#= Not sure if that is correct anymore
function curl(ref::LagrangeRefSpace{T,2,3} where {T}, sh, el)
    #curl of lagc0d2 as combination of bdm functions 
    z=zero(typeof(sh.coeff))
    if sh.refid < 4
        sh1 = Shape(sh.cellid, mod1(2*sh.refid+1,6), +sh.coeff)
        sh2 = Shape(sh.cellid, mod1(2*sh.refid+2,6), -3*sh.coeff)
        sh3 = Shape(sh.cellid, mod1(2*sh.refid+3,6), +3*sh.coeff)
        sh4 = Shape(sh.cellid, mod1(2*sh.refid+4,6), -sh.coeff)
    else
        sh1 = Shape(sh.cellid, mod1(2*sh.refid+4,6), z*sh.coeff)
        sh2 = Shape(sh.cellid, mod1(2*sh.refid+5,6), -4*sh.coeff)
        sh3 = Shape(sh.cellid, mod1(2*sh.refid+6,6), +4*sh.coeff)
        sh4 = Shape(sh.cellid, mod1(2*sh.refid+7,6), z*sh.coeff)
    end
    return [sh1, sh2, sh3, sh4]
end
=#


function curl(localspace::LagrangeRefSpace, sh, ch)
    fns = curl(localspace, sh.cellid, ch)
    α = sh.coeff
    S = typeof(sh)
    return S[S(s.cellid, s.refid, α*s.coeff) for s in fns[sh.refid]]
end


function curl(localspace::LagrangeRefSpace{T,D,3,N}, cellid::Int, ch) where {T,D,N}
    function fields(p)
        map(localspace(p)) do x
            x.curl
        end
    end
    atol = sqrt(eps(T))
    gwp = GWPDivRefSpace{T,D-1}()
    coeffs = interpolate(fields, gwp, ch)
    S = BEAST.Shape{T}
    A = Vector{Vector{S}}(undef, size(coeffs,1))
    for r in axes(coeffs,1)
        A[r] = collect(S(cellid, c, coeffs[r,c]) for c in axes(coeffs,2) if abs(coeffs[r,c]) > atol)
    end
    return SVector{N}(A)
end



@testitem "curl - chartwise" begin
    using LinearAlgebra
    using CompScienceMeshes
    const CSM = CompScienceMeshes

    T = Float64
    D = 1
    NF = binomial(2+D, 2)
    gwp = BEAST.GWPDivRefSpace{T,D-1}()
    lgc = BEAST.LagrangeRefSpace{T,D,3,NF}()

    ch = CSM.simplex(
        point(1,0,0),
        point(0,1,0),
        point(0,0,0))

    curlfns = BEAST.curl(lgc, 1, ch)

    p = neighborhood(ch, (0.2341, 0.4321))
    gwp_vals = gwp(p)
    lgc_vals = lgc(p)

    err = similar(Vector{T}, axes(gwp_vals))
    for i in eachindex(lgc_vals)
        val1 = lgc_vals[i].curl
        val2 = zero(val1)
        for sh in curlfns[i]
            val2 += sh.coeff * gwp_vals[sh.refid].value
        end
        err[i] = norm(val1 - val2)
    end
    atol = sqrt(eps(T))
    @test all(err .< atol)
end


function restrict(f::LagrangeRefSpace{T,2}, dom1, dom2) where T

    D = numfunctions(f)
    Q = zeros(T, D, D)

    # for each point of the new domain
    for i in 1:3

        #vertices
        v = dom2.vertices[i]

        # find the barycentric coordinates in dom1
        uvn = carttobary(dom1, v)

        # evaluate the shape functions in this point
        x = neighborhood(dom1, uvn)
        fx = f(x)

        for j in 1:D
            Q[j,i] = fx[j][1]
        end
        
            
        #edges
        # find the center of edge i of dom2
        a = dom2.vertices[mod1(i+1,3)]
        b = dom2.vertices[mod1(i+2,3)]
        v = (a + b) / 2

        # find the barycentric coordinates in dom1
        uvn = carttobary(dom1, v)

        # evaluate the shape functions in this point
        x = neighborhood(dom1, uvn)
        fx = f(x)
  
        for j in 4:D
            Q[j,i+3] = fx[j][1]
        end
    end

    return Q
end


const _vert_perms_lag = [
    (1,2,3),
    (2,3,1),
    (3,1,2),
    (2,1,3),
    (1,3,2),
    (3,2,1),
]

const _dof_perms_lag0 = [
    (1),
    (1),
    (1),
    (1),
    (1),
    (1),
]
const _dof_perms_lag1 = [
    (1,2,3),
    (3,1,2),
    (2,3,1),
    (2,1,3),
    (1,3,2),
    (3,2,1),
]

function dof_permutation(::LagrangeRefSpace{<:Any,0}, vert_permutation)
    i = findfirst(==(tuple(vert_permutation...)), _vert_perms_lag)
    return _dof_perms_lag0[i]
end

function dof_permutation(::LagrangeRefSpace{<:Any,1}, vert_permutation)
    i = findfirst(==(tuple(vert_permutation...)), _vert_perms_lag)
    return _dof_perms_lag1[i]
end

function dof_perm_matrix(::LagrangeRefSpace{<:Any,0}, vert_permutation)
    i = findfirst(==(tuple(vert_permutation...)), _vert_perms_rt)
    @assert i != nothing
    return _dof_lag0perm_matrix[i]
end

function dof_perm_matrix(::LagrangeRefSpace{<:Any,1}, vert_permutation)
    i = findfirst(==(tuple(vert_permutation...)), _vert_perms_rt)
    @assert i != nothing
    return _dof_rtperm_matrix[i]
end

const _dof_lag0perm_matrix = [
    @SMatrix[1],         # 1. {1,2,3}
    @SMatrix[1],         # 2. {2,3,1}
    @SMatrix[1],         # 3. {3,1,2}
    @SMatrix[1],         # 4. {2,1,3}
    @SMatrix[1],         # 5. {1,3,2}
    @SMatrix[1]         # 6. {3,2,1}
]

function (ϕ::LagrangeRefSpace{T, 1, 4, 4})(lag) where T

    u, v, w = parametric(lag)

    tu = tangents(lag, 1)
    tv = tangents(lag, 2)
    tw = tangents(lag, 3)

    B = [tu tv tw]
    A = inv(transpose(B))

    # gradient in u,v,w (unit tetrahedron)
    gr1=SVector{3, T}(1.0, 0.0, 0.0)
    gr2=SVector{3, T}(0.0, 1.0, 0.0)
    gr3=SVector{3, T}(0.0, 0.0, 1.0)
    gr4=SVector{3, T}(-1.0, -1.0, -1.0)

    return SVector((
        (value = u, gradient = A*gr1),
        (value = v, gradient = A*gr2),
        (value = w, gradient = A*gr3),
        (value = T(1.0)-u-v-w, gradient = A*gr4)
    ))
end



# Evaluate higher order Lagrange elements on triangles
# TODO: Optimise using code generation
function (ϕ::LagrangeRefSpace{T,Degree,3})(p) where {T,Degree}

    u, v = parametric(p)
    w = 1 - u - v
    idx = 0

    suppdim = 2
    localdim = binomial(suppdim+Degree, suppdim)
    vals = T[]
    diffus = T[]
    diffvs = T[]

    D1 = Degree + 1
    s = range(zero(T), one(T), length=D1)
    for i in 0:Degree
        ui = i/Degree
        for j in 0:Degree
            vj = j/Degree
            for k in 0:Degree
                wk = k/Degree
                i + j + k == Degree || continue

                prod_p = one(T)
                for p in 0:i-1
                    up = p / Degree
                    prod_p *= (u-up) / (ui-up)
                end
                prod_q = one(T)
                for q in 0:j-1
                    vq = q / Degree
                    prod_q *= (v-vq) / (vj-vq)
                end
                prod_r = one(T)
                for r in 0:k-1
                    wr = r / Degree
                    prod_r *= (w-wr) / (wk-wr)
                end
                push!(vals, prod_p * prod_q * prod_r)

                diffu = zero(T)
                diffv = zero(T)
                for l in 0:i-1
                    ul = l/Degree
                    prod_pl = one(T)
                    for p in 0:i-1
                        p == l && continue
                        up = p/Degree
                        prod_pl *= (u-up) / (ui-up)
                    end
                    diffu += prod_pl * prod_q * prod_r / (ui-ul)
                end
                for m in 0:j-1
                    vm = m/Degree
                    prod_qm = one(T)
                    for q in 0:j-1
                        q == m && continue
                        vq = q/Degree
                        prod_qm *= (v-vq) / (vj-vq)
                    end
                    diffv += prod_p * prod_qm * prod_r / (vj-vm)
                end
                for n in 0:k-1
                    wn = n/Degree
                    prod_rn = one(T)
                    for r in 0:k-1
                        r == n && continue
                        wr = r/Degree
                        prod_rn *= (w-wr) / (wk-wr)
                    end
                    diffu -= prod_p * prod_q * prod_rn / (wk-wn)
                    diffv -= prod_p * prod_q * prod_rn / (wk-wn)
                end

                push!(diffus, diffu)
                push!(diffvs, diffv)

                idx += 1
    end end end
 
    tu = tangents(p,1)
    tv = tangents(p,2)
    j = jacobian(p)
    NF = length(vals)
    SVector{NF}([(value=f, curl=(-dv*tu+du*tv)/j) for (f,du,dv) in zip(vals, diffus, diffvs)])
end

# fields[i] ≈ sum(Q[j,i] * interpolant[j].value for j in 1:numfunctions(interpolant))
function interpolate(fields, interpolant::LagrangeRefSpace{T,Degree,3}, chart) where {T,Degree}

    dim = binomial(2+Degree, Degree)
    vals = Vector{Vector{T}}()
    if Degree > 0
        I = 0:Degree
        s = range(0,1,length=Degree+1)
        Is = zip(I,s)
        for (i,ui) in Is
            for (j,vj) in Is
                for (k,wk) in Is
                    i + j + k == Degree || continue
                    @assert ui + vj + wk ≈ 1
                    p = neighborhood(chart, (ui,vj))
                    push!(vals, fields(p))
        end end end
    else
        p = center(chart)
        push!(vals, fields(p))
    end
    # Q = hcat(vals...)
    Q = Matrix{T}(undef, length(vals[1]), length(vals))
    for i in eachindex(vals)
        Q[:,i] .= vals[i]
    end
    return Q
end

function interpolate(fields, interpolant::LagrangeRefSpace{T,1,3}, chart) where {T}
    vals = Vector{Vector{T}}()


    p1 = neighborhood(chart, (1,0))
    p2 = neighborhood(chart, (0,1))
    p3 = neighborhood(chart, (0,0))
    push!(vals, fields(p1))
    push!(vals, fields(p2))
    push!(vals, fields(p3))
    
   

    # Q = hcat(vals...)
    Q = Matrix{T}(undef, length(vals[1]), length(vals))
    for i in eachindex(vals)
        Q[:,i] .= vals[i]
    end
    return Q
end


function interpolate!(out, fields, interpolant::LagrangeRefSpace{T,Degree,3}, chart) where {T,Degree}
    Is = zip((0:Degree), range(0,1,length=Degree+1))
    idx = 0
    for (i,ui) in Is
        for (j,vj) in Is
            for (k,wk) in Is
                i + j + k == Degree || continue; idx += 1
                @assert ui + vj + wk ≈ 1
                p = neighborhood(chart, (ui,vj))
                vals = fields(p)
                for (g, val) in zip(axes(out, 1), vals)
                    out[g,idx] = val
end end end end end

function interpolate!(out, fields, interpolant::LagrangeRefSpace{T,0,3}, chart) where {T}
    p = center(chart)
    vals = fields(p)
    for (g, val) in zip(axes(out, 1), vals)
        out[g,1] = val
end end