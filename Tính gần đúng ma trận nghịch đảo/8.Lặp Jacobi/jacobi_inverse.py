# khai báo thư viện
import numpy as np
from math import *

class Jacobi_Loop_inverse:

    def __init__(self, A, n, eps):  # khởi tạo
        self.A = np.reshape(np.array(A),(n,n))
        self.n = n
        self.eps = eps
        self.iters = 0

    def norm(self, __A, __norm_type = 2):   # chuẩn ma trận
        return np.linalg.norm(__A, __norm_type)

    # Kiểm tra ma trận chéo trội
    def checkDominant(self):

        # đường chéo của ma trận A
        diag_A = np.diag(np.abs(A)) 

        # tổng theo hàng trừ phần tử đường chéo
        sum_of_row_except_diag_A = np.sum(np.abs(A), axis = 1) - diag_A 

        # tổng theo cột trừ phần tử đường chéo
        sum_of_col_except_diag_A = np.sum(np.abs(A), axis = 0) - diag_A 

        if np.all(diag_A > sum_of_row_except_diag_A) and np.all(diag_A > sum_of_col_except_diag_A):
            return True
        else:
            return False

    def Jacobi_Loop(self):

        A = self.A
        n = self.n
        iters = self.iters

        # ma trận đơn vị của A
        I = np.identity(n) 

        # ma trận nghịch đảo đường chéo của A
        T = np.linalg.inv(np.diag(np.diag(A)))

        B = I - np.dot(T, A)

        norm_1 = self.norm(B, 1) # chuẩn 1 của ma trận B
        norm_inf = self.norm(B, np.inf) # chuẩn vô cùng của ma trận B

        # chéo trội cột
        if norm_1 < 1 :
            print('thảo mãn điều kiện hội tụ')

            lst = np.diag(A)

            lamda = max(abs(lst)) / min(abs(lst))

            eps0 = ((1 - norm_1) * eps) / ( lamda * norm_1)

            x0 = np.zeros(n)

            x1 = np.dot(B, x0) + T

            while self.norm((x1 - x0), 1) > eps0:
                x0 = x1
                x1 = np.dot(B, x0) + T
                iters += 1
                print("lần lặp thứ "+ str(iters))
                print(x1)
        
        # chéo trội hàng
        elif norm_inf < 1:
            eps0 = ((1 - norm_inf) * eps) / norm_inf

            x0 = np.zeros(n)

            x1 = np.dot(B, x0) + T

            while self.norm((x1 - x0), np.inf) > eps0 :
                x0 = x1
                x1 = np.dot(B, x0) + T
                iters += 1
                print("lần lặp thứ "+ str(iters))
                print(x1)     
        else:
            print("Không thỏa mãn điều kiện hội tụ")

    def Jacobi_Solve(self):  
        if(np.linalg.det(self.A) == 0):
            print("A không khả nghịch nên không đưa ra được ma trận nghịch đảo")
            return np.full((self.n, self.n), float("NaN"))
        
        elif self.checkDominant() == False:
            print("A không là ma trận chéo trội")
        else:
            return self.Jacobi_Loop()

# chương trình chính
eps = float(input("Nhập sai số eps : "))

# nhập ma trận A từ file
# A = np.loadtxt('B:/BT nhóm GTS/Tính gần đúng ma trận nghịch đảo/8.Lặp Jacobi/jacobi inverse.txt', dtype = float)
A = np.loadtxt('B:/BT nhóm GTS/Tính gần đúng ma trận nghịch đảo/8.Lặp Jacobi/cau 2(2018.2).txt', dtype = float)
n,_ = np.shape(A)

test = Jacobi_Loop_inverse(A,n,eps)

A1 = test.Jacobi_Loop()

print("================================================================")

#sử dụng numpy.linalg để kiểm tra kết quả

print("Nghich dao cua ma tran A theo numpy:")
print(np.linalg.inv(A))
