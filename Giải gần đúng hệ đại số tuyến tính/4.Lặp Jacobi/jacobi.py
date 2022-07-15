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

def get_B(A):
    A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi.txt", dtype = float)
    n = A.shape[0]
    I = np.identity(n) 
    T = np.linalg.inv(np.diag(np.diag(A)))
    B = I - np.dot(T, A)
    
    return B

def get_d(A,b):
    A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi.txt", dtype = float)
    b = np.array([1,2,3,4], dtype = float)
    T = np.linalg.inv(np.diag(np.diag(A)))
    d = np.dot(T, b)
    
    return d
    
    

def jacobiDominantRow(A,b,eps):
    A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi.txt", dtype = float)
    
    B = get_B(A)
    d = get_d(A,b)
    norm_inf = normMatrix(B)
    if norm_inf < 1:
        print('thoa man dieu kien hoi tu')
        eps0 = ((1 - norm_inf) * eps) / norm_inf
        x0 = np.zeros_like(b)
        x1 = np.dot(B, x0) + d

        i = 0

        while normMatrix(x1 - x0) > eps0 :
            x0 = x1
            x1 = np.dot(B, x0) + d
            i += 1
            print(x1)
        print("so lan lap", i)    


def jacobiDominantCol(A,b,eps):
    A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi.txt", dtype = float)
    
    B = get_B(A)
    d = get_d(A,b)
    norm_1 = normMatrix(B)
    if norm_1 < 1 :
        print('thoa man dieu kien hoi tu')

        lst = np.diag(A)

        lamda = max(abs(lst)) / min(abs(lst))

        eps0 = ((1 - norm_1) * eps) / ( lamda * norm_1)
        x0 = np.zeros_like(b)
        x1 = np.dot(B, x0) + d

        i = 0

        while normMatrix(x1 - x0) > eps0:
            x0 = x1
            x1 = np.dot(B, x0) + d
            i += 1
            print(x1)
        print("so lan lap", i)

def jacobi(A,b,eps):
    if dominantMatrix(A) == 0 :
        jacobiDominantRow(A,b,eps)
    elif dominantMatrix(A) == 1:
        jacobiDominantCol(A,b,eps)
    else:
        print("Ma trận không chéo trội")

# main
A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi.txt", dtype = float)
b = np.array([1,2,3,4], dtype = float)
print("Nhập eps:")
eps = float(input())
jacobi(A,b,eps)
