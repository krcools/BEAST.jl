using Test
using BEAST

const B3D = BEAST.SparseND.Banded3D

k0 = [2 2; 3 2]
k1 = [2 3; 3 3]
M = 2
N = 2
bandwidth = maximum(k1 - k0 .+ 1)
data = collect(begin
    k = k0[i,j] + l - 1
    v = (k > k1[i,j]) ? 0.0 : Float64((i-1) + (j-1)*M + (k-1)*M*N)
end for l in 1:bandwidth, i in axes(k0,1), j in axes(k0,2))
A = B3D(k0,data,maximum(k1))

@test size(A) == (2,2,3)
@test eltype(A) == Float64
@test A[1,1,1] == 0
@test A[1,1,2] == M*N
@test A[1,1,2] == M*N
@test A[1,1,3] == 0

slice = A[:,:,2]
@test slice == [4 6; 0 7]

F = Array(A)
for m in axes(F,1)
    for n in axes(F,2)
        for l in axes(F,3)
            @test F[m,n,l] ≈ A[m,n,l]
        end
    end
end

st_data = rand(2,5)
cv1 = BEAST.convolve(F, st_data, 4, 2)
cv2 = BEAST.convolve(A, st_data, 4, 2)
@test cv1 ≈ cv2

cv1 = BEAST.convolve(F, st_data, 4, 1)
cv2 = BEAST.convolve(A, st_data, 4, 1)
@test cv1 ≈ cv2

cv1 = BEAST.convolve(F, st_data, 4, 3)
cv2 = BEAST.convolve(A, st_data, 4, 3)
@test cv1 ≈ cv2
