# khai báo thư viện 
import numpy as np
from math import *

class newton_matrix_inversion:
    max_attempt = 2

    def __init__(self, A, n, eps):    # Khởi tạo
        self.A = np.reshape(np.array(A), (n, n))
        self.n = n
        self.eps = eps
        self.iters= 0

    def norm(self, __A, __norm_type = 2):   # Chuẩn ma trận
        return np.linalg.norm(__A, __norm_type)

    def initial_approx(self):  # Tìm xấp xỉ đầu:

        E  = np.identity(self.n)    # Ma trận đơn vị của A
        A  = self.A
        X  = self.A
        
        # PP xấp xỉ đầu 
        t1 = self.norm(X, 1)
        t2 = self.norm(X, inf)
        X  = (X / (t1 * t2)).T

        # Hiệu chỉnh lại giá trị q của ma trận xấp xỉ đầu
        attempt = 0
        while(attempt <= newton_matrix_inversion.max_attempt):
            X = X @ (2 * E - A @ X)
            if(self.norm(E - A @ X) < 1):
                attempt += 1
            self.iters += 1 * bool(attempt == 0)

        # Trả về ma trận xấp xỉ đầu
        return X

    def __pure_newton(self, X_0):  # Lặp Newton nguyên bản

        norm_X0 = self.norm(X_0)   
        E       = np.identity(self.n)
        A       = self.A
        eps     = self.eps
        
        q2k = q = self.norm(E - A @ X_0)
        print(q)
        X = X_0

        # Kiểm tra điều kiện hội tụ
        if(q >= 1):
            print("Xấp xỉ đầu không thỏa mãn nên không đưa ra được ma trận nghịch đảo.")
            return np.full((self.n, self.n), float("NaN"))
        
        # Lặp
        while(norm_X0 * q2k >= self.eps * (1 - q)):
            self.iters += 1
            X = X @ (2 * E - A @ X)
            q2k = q2k ** 2
            print("X{0} = {1}".format(self.iters, X))

        # Đưa ra ma trận cuối cùng
        print(f"Phương pháp Newton kết thúc sau {self.iters} bước lặp")
        print(X)
        return X
        

    def improved_newton(self):  # PP Newton cải tiến với xấp xỉ đầu
        if(np.linalg.det(self.A) == 0): # định thức của ma trận A
            print("A không khả nghịch nên không đưa ra được ma trận nghịch đảo")
            return np.full((self.n, self.n), float("NaN"))
        return self.__pure_newton(self.initial_approx())

# chương trình chính
# Nhập sai số eps
eps = float(input("Nhập sai số eps: "));

# nhập ma trận A từ file

# A = np.loadtxt('B:/Nguyễn Quang Minh - Lớp 133584 - Bài tập nhóm GTS/Tính gần đúng ma trận nghịch đảo/10.Lặp Newton/newton_inverse.txt', dtype = float)
# A = np.loadtxt('B:/Nguyễn Quang Minh - Lớp 133584 - Bài tập nhóm GTS/Tính gần đúng ma trận nghịch đảo/10.Lặp Newton/cau 2(2018.2).txt', dtype = float)
A = np.loadtxt('B:/Nguyễn Quang Minh - Lớp 133584 - Bài tập nhóm GTS/Tính gần đúng ma trận nghịch đảo/10.Lặp Newton/matrix det get equal 0.txt', dtype = float)

n,_ = np.shape(A)

# In ra ma trận A
print("--------------------------------------------------------------------")
print("Ma trận A:")
print(A)
print("--------------------------------------------------------------------")

# Nghịch đảo ma trận
uu = newton_matrix_inversion(A, n, eps);

# In ra ma trận nghịch đảo
print("--------------------------------------------------------------------")
B = uu.improved_newton()
# lưu ma trận nghịch đảo vào file
np.savetxt('B:/Nguyễn Quang Minh - Lớp 133584 - Bài tập nhóm GTS/Tính gần đúng ma trận nghịch đảo/10.Lặp Newton/inverse_matrix_A.txt', B, delimiter = ' ')
# np.savetxt('B:/Nguyễn Quang Minh - Lớp 133584 - Bài tập nhóm GTS/Tính gần đúng ma trận nghịch đảo/10.Lặp Newton/dap an cau 2(2018.2).txt', B, delimiter = ' ')
    
