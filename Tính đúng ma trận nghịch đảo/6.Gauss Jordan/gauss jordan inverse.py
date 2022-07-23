import numpy as np
import sys

# nhập ma trận b
b = np.loadtxt('B:/Nguyễn Quang Minh - Lớp 133584 - Bài tập nhóm GTS/Tính đúng ma trận nghịch đảo/6.Gauss Jordan/gauss jordan inverse.txt', dtype = float)
#b = np.loadtxt('B:/Nguyễn Quang Minh - Lớp 133584 - Bài tập nhóm GTS/Tính đúng ma trận nghịch đảo/6.Gauss Jordan/gauss_jordan_inverse_det_equal0.txt', dtype = float)

# Kiểm tra điều kiện ma trận khả nghịch
if np.linalg.det(b) != 0:
    n = b.shape[0]

    # tìm ma trận bổ sung 
    a = np.zeros((n,2*n))

    for i in range(n):        
        for j in range(n):
            if i == j:
                a[i][j+n] = 1

    for i in range(n):        
        for j in range(n):
            a[i][j] = b[i][j]

    #Biến đổi Gauss Jordan
    print("Bien doi Gauss Jordan tim ma tran nghich dao: \n")
    for i in range(n):
        if a[i][i] == 0.0:
            sys.exit('ERROR')
        for j in range(n):
            if i != j:
                ratio = a[j][i] / a[i][i]

                for k in range(2*n):
                    a[j][k] = a[j][k] - ratio * a[i][k]
        print("\n", a)
    #Cho vế trái thành mt đơn vị
    for i in range(n):
        temp = a[i][i]
        for j in range(2*n):
            a[i][j] = a[i][j] / temp
        print("\n", a)
    print("Vay ta thu duoc ma tran: \n",a)
    print("==================================================================== \n")

    print("ma trận nghịch đảo:")
    for i in range(n):
        for j in range(n, 2*n):
            print(a[i][j], end='\t')
        print()

    print("==================================================================== \n")
else:
    print("Ma trận ko khả nghịch, ko tồn tại ma trận nghịch đảo!")