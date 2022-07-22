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

def getMatrix():
    # A = np.loadtxt("B:/Nguyễn Quang Minh - Lớp 133584 - Bài tập nhóm GTS/Giải gần đúng hệ đại số tuyến tính/5. Lặp gauss-seidel/gauss_seidel.txt", dtype = float)
    # A = np.loadtxt("B:/Nguyễn Quang Minh - Lớp 133584 - Bài tập nhóm GTS/Giải gần đúng hệ đại số tuyến tính/5. Lặp gauss-seidel/gauss_seidel_cheo_troi_hang.txt", dtype = float)
    # A = np.loadtxt("B:/Nguyễn Quang Minh - Lớp 133584 - Bài tập nhóm GTS/Giải gần đúng hệ đại số tuyến tính/5. Lặp gauss-seidel/gauss_seidel_ko_cheo_troi.txt", dtype = float)
    A = np.loadtxt("B:/Nguyễn Quang Minh - Lớp 133584 - Bài tập nhóm GTS/Giải gần đúng hệ đại số tuyến tính/5. Lặp gauss-seidel/gauss_seidel_cheo_troi_cot.txt", dtype = float)
    return A

# gói tính chuẩn
def normMatrix(matrix):
    A = getMatrix()
    if dominantMatrix(A) == 0:
        return normInf(matrix)
    if dominantMatrix(A) == 1:
        return norm1(matrix)
    if dominantMatrix(A) == -1:
        return 0
    return 0 
    
# Gói biến đổi từ AX + B = 0 về dạng X = CX + D hoặc từ AX + B = 0 về dạng Y = CY + D
def setAtoC(A):
    n = A.shape[0]
    C = np.arange(n ** 2)
    C = C.reshape((n,n))
    C = np.zeros_like(C, dtype = float)

    if dominantMatrix(A) == 0:
        for i in range(n):
            for j in range(n):
                if i != j:
                    C[i,j] = -A[i,j] / A[i,i]
                else:
                    C[i,j] = 0

    if dominantMatrix(A) == 1:
        for i in range(n):
            for j in range(n):
                if i != j:
                    C[i,j] = -A[j,i] / A[i,i]
                else:
                    C[i,j] = 0

    return C

def setBtoD(A, b):
    n = A.shape[1]
    D = np.zeros(n)
    for i in range(n):
        D[i] = b[i] / A[i,i]

    return D 

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

# Thiết lập công thức sai số, hệ số co
def getqCoeff(A, C):
    n = C.shape[1]
    a = np.zeros(n)
    b = np.zeros(n)
    sum = q = 0
    if dominantMatrix(A) == 0:
        for i in range(n):
            sum = 0
            for j in range(i,n):
                sum += math.fabs(C[i,j])
            a[i] = sum 
            sum = 0
            for k in range(i-1):
                sum += math.fabs(C[i,k])
            b[i] = 1 - sum 
            if a[i] / b[i] > q:
                q = a[i] / b[i]
    print(q)

    if dominantMatrix(A) == 1:
        for i in range(n):
            sum = 0
            for j in range(1,i):
                sum += math.fabs(C[j,i])
            a[i] = sum 
            sum = 0
            for k in range(i+1, n):
                sum += math.fabs(C[k,i])
            b[i] = 1 - sum 
            if a[i] / b[i] > q:
                q = a[i] / b[i]
    print(q)
    
    return q 

def getSCoeff(A,C):
    n = A.shape[0]
    S = sum = 0
    if dominantMatrix(A) == 0:
        S = 0
    if dominantMatrix(A) == 1:
        for i in range(n):
            sum = 0
            for j in range(i+1,n):
                sum += math.fabs(C[j,i])
            if sum > S:
                S = sum

    return S 

# Gói lặp Gauss Seidel
def seidelLoop(A, D, eps):
    n = A.shape[0]
    X = np.zeros(n)
    X_k = np.zeros(n)
    Col = np.arange(n**2)
    Col = Col.reshape((n,n))
    Col = np.zeros_like(Col, dtype = float)

    C_ = setAtoC(A)
    C = defintion_B(A)

    for i in range(n):
        Col[i,i] = 1 / A[i,i]

    
    S = getSCoeff(A, C_)
    q = getqCoeff(A, C_)

    C1 = np.arange(n**2)
    C1 = C1.reshape((n,n))
    C1 = np.zeros_like(C1, dtype = float)

    C2 = np.arange(n**2)
    C2 = C2.reshape((n,n))
    C2 = np.zeros_like(C2, dtype = float)

    for i in range(n):
        for j in range(0, i):
            C1[i,j] = C[i,j]
        for k in range(i, n):
            C2[i,k] = C[i,k]
    
    i = 0
            
    X = C1 @ X + C2 @ X + D

    while q * normMatrix(X - X_k) > eps * (1 - q) * (1 - S):
        print(X)
        for i in range(n):
            X_k[i] = X[i]
        X = C1 @ X + C2 @ X + D
        i += 1
        
    print("\nNghiệm của hệ PT là: ")
    print("X = ")
    print(X)

    return X

# =======================================================================
# chương trình chính

A = getMatrix()
print("A = \n", A)

# b = np.array([1,2,3,4])
b = np.array([1,2,3,4,5,6])
if dominantMatrix(A) == -1:
    print("Ma trận A ko chéo trội\n")
else:
    print("Nhập sai số epxilon: ")
    eps = float(input())

    D = setBtoD(A, b)

    if dominantMatrix(A) == 0:
        print("A là ma trận chéo trội hàng \n")
        seidelLoop(A, D, eps)

    if dominantMatrix(A) == 1:
        print("A là ma trận chéo trội cột \n")
        seidelLoop(A, D, eps)

