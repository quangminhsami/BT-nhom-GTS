import numpy as np
import sys

b = np.loadtxt('B:/BT nhóm GTS/Tính đúng ma trận nghịch đảo/6.Gauss Jordan/gauss jordan inverse.txt', dtype = float)

n = b.shape[0]

# tìm ma trận bổ sung 
a = np.zeros((n,2*n))

for i in range(n):        
    for j in range(n):
        if i == j:
            a[i][j+n] = 1

for i in range(n):        
    for j in range(n):
        a[i][j]=b[i][j]


#Biến đổi Gauss Jordan
print("Bien doi Gauss Jordan tim ma tran nghich dao: \n")
for i in range(n):
    if a[i][i] == 0.0:
        sys.exit('ERROR')
    for j in range(n):
        if i != j:
            ratio = a[j][i]/a[i][i]

            for k in range(2*n):
                a[j][k] = a[j][k] - ratio * a[i][k]
    print("\n", a)
#Cho vế trái thành mt đơn vị
for i in range(n):
    temp = a[i][i]
    for j in range(2*n):
        a[i][j] = a[i][j]/temp
    print("\n", a)
print("Vay ta thu duoc ma tran: \n",a)
print("==================================================================== \n")

print("ma trận nghịch đảo:")
for i in range(n):
    for j in range(n, 2*n):
        print(a[i][j], end='\t')
    print()

print("==================================================================== \n")