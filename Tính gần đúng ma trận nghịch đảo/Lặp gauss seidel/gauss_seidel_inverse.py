import numpy as np
from math import *

# Phần thuật toán chính
class gaussseidel_mat_inversion:
    def __init__(self, A, n, eps):      # khởi tạo
        self.A = np.reshape(np.array(A), (n, n))
        self.eps = eps
        self.n = n

    def __norm(self, __A, __norm_type = 2):        # tính chuẩn của ma trận
        return np.linalg.norm(__A, __norm_type)

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

    def __getNorm(self, __A, __domination_status):  # Tính chuẩn
        if(__domination_status == 1): return self.__norm(__A, inf)
        return self.__norm(__A, 1)
    
    def __get_S_coeff(self, __A, __domination_status):  # Tính S
        if(__domination_status == 1): return 0
        
        S = 0
        n = int(__A.shape[0])
        for i in range(n): 
            tmp = 0
            for j in range(i, n): tmp += abs(__A[j, i])
            S = max(S, tmp)
        return S

    def __get_q_coeff(self, __A, __domination_status): # Tính q

        q = 0
        n = int(__A.shape[0])
        for i in range(n):

            Q1 = Q2 = 0
            for j in range(n):
                if(__domination_status == 1):
                    if(j >= i-1): Q1 += abs(__A[i, j])
                    else: Q2 += abs(__A[i, j])
                else:
                    if(j <= i-1): Q1 += abs(__A[j, i])
                    else: Q2 += abs(__A[j, i])
            q = max(q, Q1 / (1 - Q2))
        return q

    # Tiến trình lặp
    def _________next_iteration(self, old_X, B, T, relax_factor):
        n = int(B.shape[0])
        new_X = np.zeros_like(B)
        for i in range(n):
            for j in range(i): 
                new_X[i] += B[i, j] * new_X[j]
            for j in range(i, n):
                 new_X[i] += B[i, j] * old_X[j]
            new_X[i] += T[i]

        return (1 - relax_factor) * old_X + relax_factor * new_X

    def __predecessor_iteration(self, X_0, B, T, S, q, p, relax_factor):    #SD công thức sai số tiên nghiệm
        X_1 = self._________next_iteration(X_0, B, T, relax_factor)
        predecessor_norm = self.__getNorm(X_1 - X_0, p)
        X = X_0
        qk = 1
        nr_iteration = 0

        while(qk * predecessor_norm > self.eps * (1 - q) * (1 - S)):
            nr_iteration += 1
            X   = self._________next_iteration(X, B, T, relax_factor)
            qk *= q
            print("Lần lặp thứ {0}".format(nr_iteration))
            print(X)

        print(f"Phương pháp Gauss-Seidel đánh giá tiên nghiệm kết thúc sau {nr_iteration} bước lặp")
        return X

    def __successor_iteration(self, X_0, B, T, S, q, p, relax_factor):      #SD công thức sai số hậu nghiệm
        new_X = self._________next_iteration(X_0, B, T, relax_factor)
        old_X = X_0

        nr_iteration = 0

        while(q * self.__getNorm(new_X - old_X, p) > self.eps * (1 - q) * (1 - S)):
            nr_iteration += 1
            old_X = new_X
            new_X = self._________next_iteration(old_X, B, T, relax_factor)

        print(f"Phương pháp Gauss-Seidel đánh giá hậu nghiệm kết thúc sau {nr_iteration} bước lặp")
        return new_X
 
    # PP lặp Gauss-Seidel với chế độ đánh giá tiên nghiệm hoặc hậu nghiệm và cải tiến hệ số điều chỉnh
    def gauss_seidel_iteration(self, mode = 1, relax_factor = 1):
        # Gán các biến cơ bản
        A = self.A
        n = self.n
        E = np.eye(self.n)

        # Kiểm tra chéo trội và khả nghịch
        p = self.checkDominant()
        if(np.linalg.det(self.A) == 0):
            print("A không khả nghịch nên không đưa ra được ma trận nghịch đảo")
            return np.full((self.n, self.n), float("NaN"))
        elif(p == False): 
            print("A không chéo trội nên không đưa ra được ma trận nghịch đảo")
            return np.full((self.n, self.n), float("NaN"))

        else:
            # Tính T, B, q, S
            T = np.diag(1 / np.diag(A))
            B = E - T @ A
            S = self.__get_S_coeff(B, p)
            q = self.__get_q_coeff(B, p)
            
            # Đưa ra ma trận cuối cùng
            if(mode == 1): return self.__predecessor_iteration(A, B, T, S, q, p, relax_factor)
            return self.__successor_iteration(A, B, T, S, q, p, relax_factor)
   
# chuong trinh chinh
# Nhập sai số eps
eps = float(input("Nhập sai số eps: "))

A = np.loadtxt('seidel_inverse.txt', dtype = float) 

n,_ = np.shape(A)

# In ra ma trận A
print("--------------------------------------------------------------------")
print("Ma trận A:")
print(A)
print("--------------------------------------------------------------------")

# Nghịch đảo ma trận
uu = gaussseidel_mat_inversion(A, n, eps)

# In ra ma trận nghịch đảo
print("--------------------------------------------------------------------")
A1 = uu.gauss_seidel_iteration(1)
print("Ma trận nghịch đảo của A theo PP Gauss-Seidel tiên nghiệm:")
print(A1)

# In ra ma trận nghịch đảo
print("--------------------------------------------------------------------")
A2 = uu.gauss_seidel_iteration(2)
print("Ma trận nghịch đảo của A theo PP Gauss-Seidel hậu nghiệm:");
print(A2)