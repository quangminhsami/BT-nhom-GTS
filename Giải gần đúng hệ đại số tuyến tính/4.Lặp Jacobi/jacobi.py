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

    elif np.all(diag_A > sum_of_col_except_diag_A):
        return 1
    else:
        return -1

# Gói tính và trả về giá trị chuẩn vô cùng
def normInf(matrix):
    return np.linalg.norm(matrix, np.inf)

# Gói tính và trả về giá trị chuẩn 1
def norm1(matrix):
    return np.linalg.norm(matrix, 1)

# gói tính chuẩn
def normMatrix(matrix):
    # A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi.txt", dtype = float)
    # A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi_cheo_troi_hang.txt", dtype = float)
    A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi_cheo_troi_cot.txt", dtype = float)
    # A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi_ko_cheo_troi.txt", dtype = float)
    if dominantMatrix(A) == 0:
        return normInf(matrix)
    elif dominantMatrix(A) == 1:
        return norm1(matrix)
    else:
        return 0 

def get_B(A):
    # A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi.txt", dtype = float)
    # A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi_cheo_troi_hang.txt", dtype = float)
    A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi_cheo_troi_cot.txt", dtype = float)
    # A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi_ko_cheo_troi.txt", dtype = float)
    n = A.shape[0]
    I = np.identity(n) 
    T = np.linalg.inv(np.diag(np.diag(A)))
    B = I - np.dot(T, A)
    
    return B

def get_d(A,b):
    # A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi.txt", dtype = float)
    # A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi_cheo_troi_hang.txt", dtype = float)
    A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi_cheo_troi_cot.txt", dtype = float)
    # A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi_ko_cheo_troi.txt", dtype = float)
    b = np.array([1,2,3,4,5,6], dtype = float)
    T = np.linalg.inv(np.diag(np.diag(A)))
    d = np.dot(T, b)
    
    return d
    
    
# gói lặp jacobi TH chéo trội hàng
def jacobiDominantRow(A,b,eps):
    # A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi.txt", dtype = float)
    # A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi_cheo_troi_hang.txt", dtype = float)
    A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi_cheo_troi_cot.txt", dtype = float)
    # A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi_ko_cheo_troi.txt", dtype = float)
    
    B = get_B(A)
    d = get_d(A,b)
    
    norm_inf = normMatrix(B)
    print("q = ", norm_inf)
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

# gói lặp jacobi TH chéo trội cột
def jacobiDominantCol(A,b,eps):
    # A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi.txt", dtype = float)
    # A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi_cheo_troi_hang.txt", dtype = float)
    A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi_cheo_troi_cot.txt", dtype = float)
    # A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi_ko_cheo_troi.txt", dtype = float)
    
    B = get_B(A)
    d = get_d(A,b)
    norm_1 = normMatrix(B)
    print("q = ", norm_1)
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


# main

# A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi.txt", dtype = float)
# A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi_cheo_troi_hang.txt", dtype = float)
A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi_cheo_troi_cot.txt", dtype = float)
# A = np.loadtxt("B:/BT nhóm GTS/Giải gần đúng hệ đại số tuyến tính/4.Lặp Jacobi/jacobi_ko_cheo_troi.txt", dtype = float)
b = np.array([1,2,3,4,5,6], dtype = float)

if dominantMatrix(A) == -1:
    print("Ma trận A ko chéo trội")
else:
    print("Nhập sai số epxilon: ")
    eps = float(input())

    if dominantMatrix(A) == 0:
        print("A là ma trận chéo trội hàng")
        jacobiDominantRow(A,b,eps)

    if dominantMatrix(A) == 1:
        print("A là ma trận chéo trội cột")
        jacobiDominantCol(A,b,eps)
    
