import numpy as np
import math

# Gói kiểm tra tính chéo trội của ma trận
def dominantMatrix(matrix):
    diag_A = np.diag(np.abs(matrix)) 

    # tong theo hang
    sum_of_row_except_diag_A = np.sum(np.abs(matrix), axis = 1) - diag_A 

    # tong theo cot
    sum_of_col_except_diag_A = np.sum(np.abs(matrix), axis = 0) - diag_A 

    # kiem tra ma tran cheo troi
    if np.all(diag_A > sum_of_row_except_diag_A) :
        return 0

    if np.all(diag_A > sum_of_col_except_diag_A):
        return 1
    
    return -1

# Gói tính và trả về giá trị chuẩn vô cùng
def normInf(matrix):
    return np.linalg.norm(matrix, np.inf)

# Gói tính và trả về giá trị chuẩn 1
def norm1(matrix):
    return np.linalg.norm(matrix, 1)

def normMatrix(matrix):
    A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi.txt", dtype = float)
    if dominantMatrix(A) == 0:
        return normInf(matrix)
    if dominantMatrix(A) == 1:
        return norm1(matrix)
    if dominantMatrix(A) == -1:
        return 0
    return 0 

# Gói xác định B, B1, d

def defintion_B(A):
    n = A.shape[0]
    B = np.arange(n ** 2)
    B = B.reshape((n,n))
    B = np.zeros_like(B, dtype = float)
    
    for i in range(n):
        for j in range(n):
            if i == j:
                B[i,j] = 0
            else:
                B[i,j] = -A[i,j] / A[i,i]

    return B

def defintion_B1(A):
    n = A.shape[0]
    B1 = np.arange(n ** 2)
    B1 = B1.reshape((n,n))
    B1 = np.zeros_like(B1, dtype = float)
    
    for i in range(n):
        for j in range(n):
            if i == j:
                B1[i,j] = 0
            else:
                B1[j,i] = -A[j,i] / A[i,i]

    return B1

def definition_d(A, b):
    n = A.shape[0]
    d = np.zeros(n)
    for i in range(n):
        d[i] = b[i] / A[i,i]

    return d

# Gói lặp đơn
def singleLoop(B,d,eps):
    n = B.shape[0]
    X = np.zeros(n)
    X = d   
    X_k = np.zeros(n)
    X_k = B @ X + d

    i = 0
    while normMatrix(X - X_k) > eps:
        print("Lần lặp thứ ", i)
        print(f'X{i} = ')
        print(X_k)
        X = X_k
        X_k = B @ X + d
        i += 1

    print("số lần lặp là: ", i)
    print("Nghiệm của hệ PT là: ")
    print("X = ")
    
    return X_k
    




# =========================================================================
# chương trình chính

A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/single loop.txt", dtype = float)
b = np.array([1,2,3,4])

if dominantMatrix(A) == -1:
    print("Ma trận A không chéo trội")

print("Nhập sai số epxilon: ")
eps = float(input())

d = definition_d(A, b)
B = defintion_B(A)
B1 = defintion_B1(A)

if dominantMatrix(A) == 0:
    print("A là ma trận chéo trội hàng")
    norm_B = normMatrix(B)
    print("q = ", norm_B)
    eps_0 = eps * (1 - norm_B) / norm_B
    print(singleLoop(B,d,eps_0))

if dominantMatrix(A) == 1:
    norm_B1 = normMatrix(B1)
    print("A là ma trận chéo trội cột")
    print("q = ", norm_B1)

    norm_T = 0
    n = A.shape[0]
    for i in range(n):
        if 1 / math.fabs(A[i,i]) > norm_T:
            norm_T = 1 / math.fabs(A[i,i])

    norm_D = 0
    for i in range(n):
        if math.fabs(A[i,i]) > norm_D:
            norm_D = math.fabs(A[i,i])

    eps0 = eps * (1 - norm_B1) / (norm_B1 * norm_D * norm_T)

    print(singleLoop(B,d,eps0))