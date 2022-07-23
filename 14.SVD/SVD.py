import numpy as np
from numpy.linalg import svd

A = np.loadtxt('B:/Nguyễn Quang Minh - Lớp 133584 - Bài tập nhóm GTS/14.SVD/SVD.txt', dtype = float)
print("Ma trận A : \n", A)
U, S, VT = svd(A, full_matrices = True)
print(U.shape, S.shape, VT.shape)
print('Ma tran U la: ')
print(U)
print('Duong cheo tren ma tran S la: ')
print(S)
print('Ma tran chuyen vi cua V la ')
print(VT.T)