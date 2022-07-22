import numpy as np
import sys

print("PP Gauss Jordan")
print("=========================================")
# Nhập ma trận bổ sung/ liên kết a   
 
# TH có nghiệm duy nhát
# a = np.loadtxt('B:/BT nhóm GTS/Giải đúng hệ đại số tuyến tính/2.Gauss Jordan/gauss_jordan.txt', dtype = float)

# TH vo so nghiem
# a = np.loadtxt('B:/BT nhóm GTS/Giải đúng hệ đại số tuyến tính/2.Gauss Jordan/gauss_jordan_vo_so_nghiem.txt', dtype = float)

# TH vô nghiệm
a = np.loadtxt('B:/BT nhóm GTS/Giải đúng hệ đại số tuyến tính/2.Gauss Jordan/gauss_jordan_vo_nghiem.txt', dtype = float)

# Hạng của ma trận a
rank_a = np.linalg.matrix_rank(a)
            
print("1.INPUT: \n")
print("Ma trận a: \n", a)

n = a.shape[0]

# Ma trận hệ số B
b = np.zeros((n,n),float)

for i in range(n):        
    for j in range(n):
        b[i][j]=a[i][j]

# Hạng của ma trận b
rank_b = np.linalg.matrix_rank(b)

#Tạo ma trận 0
x = np.zeros(n,float)

# Biến đổi Gauss Jordan
def gauss_jordan():
    for i in range(n):
        if a[i][i] == 0.0:
            sys.exit('ERROR')
            
        for j in range(n):
            if i != j:
                ratio = a[j][i] / a[i][i]

                for k in range(n+1):
                    a[j][k] = a[j][k] - ratio * a[i][k]
        print("Lan ", i+1)   
        print(a)
    print("Vay ta thu duoc ma tran \n",a)

    # Tìm nghiệm
    for i in range(n):
        x[i] = a[i][n]/a[i][i]
    print("\nKet qua: ",x)
    for i in range(n):
        print('X%d = %0.2f' %(i,x[i]), end = '\t')

# Đánh giá
if rank_a == rank_b:
    if rank_a == n:
        print("2.OUTPUT: Co nghiem duy nhat\n")
    else:
        print("2.OUTPUT: Co vo so nghiem")
else:
    print("2.OUTPUT: Vo nghiem")

gauss_jordan()