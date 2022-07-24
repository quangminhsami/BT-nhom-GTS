import numpy as np
np.set_printoptions(suppress = True, linewidth = np.inf, precision = 12) # sai so

print("Chương Trình Tìm nghịch đảo của ma trận bằng phương pháp viền quanh")

# nhap ma tran A tu file
# A = np.loadtxt("B:/Nguyễn Quang Minh - Lớp 133584 - Bài tập nhóm GTS/Tính gần đúng ma trận nghịch đảo/9.Viền quanh/bordering.txt", dtype = 'float', delimiter = ' ')
A = np.loadtxt("B:/Nguyễn Quang Minh - Lớp 133584 - Bài tập nhóm GTS/Tính gần đúng ma trận nghịch đảo/9.Viền quanh/matrix det get equal 0.txt", dtype = float)
# kich thuoc cua ma tran A
n = len(A)
# chuyen vi cua ma tran A
b = np.transpose(A)
M = b.dot(A)
print("Ma trận vừa nhập tu file là:")
print(A)

def bordering(A, n):    # Nhap ma tran A va kich thuoc cua ma tran
    if n == 2: 
        inv = np.zeros([2,2])   # Khởi tạo ma trận nghịch đảo alpha^-1
        if (A[0, 0] * A[1, 1] - A[1, 0] * A[0, 1]) == 0:    # dinh thuc cap 2
            print("Ma tran ko kha nghich")
            quit()
        else:   # tinh ma tran nghich dao alpha^(-1)
            inv[0, 0] = A[1, 1] / (A[0, 0] * A[1, 1] - A[1, 0] * A[0, 1])
            inv[0, 1] = -A[0, 1] / (A[0, 0] * A[1, 1] - A[1, 0] * A[0, 1])
            inv[1, 0] = -A[1, 0] / (A[0, 0] * A[1, 1] - A[1, 0] * A[0, 1])
            inv[1, 1] = A[0, 0] / (A[0, 0] * A[1, 1] - A[1, 0] * A[0, 1])
            return inv
    A_11 = np.zeros([n - 1, n - 1])  # Khởi tạo các ma trận con
    A_12 = np.zeros([n - 1, 1])
    A_21 = np.zeros([1, n - 1])
    A_22 = A[n - 1, n - 1]
    for i in range(0, n - 1):
        for j in range(0, n - 1):
            A_11[i, j] = A[i, j]
    for i in range(0, n - 1):
        A_12[i, 0] = A[i, n - 1]
        A_21[0, i] = A[n - 1, i]
    if n > 2: # Nếu chưa gặp trường hợp chặn
        a_111 = bordering(A_11, n-1)
        X = a_111.dot(A_12)
        Y = A_21.dot(a_111)
        theta = A_22 - Y.dot(A_12)
        if theta == 0:
            print('Ma Tran Khong Kha Nghich')
            quit()
        x = X.dot(Y)
        for i in range(0, n - 1):
            for j in range(0, n - 1):
                A[i, j] = a_111[i, j] + x[i, j] / theta
                A[n - 1, i] = Y[0, i] / (- theta)
                A[i, n - 1] = - X[i, 0] / theta
        A[n - 1, n - 1] = 1 / theta
    return A

s = bordering(M, n)

print("Ma trận nghịch đảo là:")
print(s.dot(b))
# luu ket qua vao file
np.savetxt("B:/Nguyễn Quang Minh - Lớp 133584 - Bài tập nhóm GTS/Tính gần đúng ma trận nghịch đảo/9.Viền quanh/test_output_bordering.txt", s.dot(b), fmt= '%-.8f', delimiter = ' ')

