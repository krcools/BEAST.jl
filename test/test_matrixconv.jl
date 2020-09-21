using Test

using BEAST

A11 = BEAST.MatrixConvolution(rand(2,3,10))
A21 = BEAST.MatrixConvolution(rand(4,3,8))
A12 = BEAST.MatrixConvolution(rand(2,5,12))
A22 = BEAST.MatrixConvolution(rand(4,5,14))

A = BEAST.hvcat((2,2), A11, A21, A12, A22)

@test A[3,5,6] == A22[1,2,6]