import numpy as np
import sys

print("PP Gauss Jordan tim ma tran nghich dao\n")
print("=========================================")
#Nhap b
b = np.array([[4, 12, -16],
             [12, 37, -43],
             [-16, -43, 98]],float)

print("1.INPUT: \n")
print(b)
#Tinh số hàng
n,_=np.shape(b)
print("n=",n,"\n")

# Ma trận hợp
a = np.zeros((n,2*n))
##Tạo đơn vị
for i in range(n):        
    for j in range(n):
        if i == j:
            a[i][j+n] = 1
##Chép ma trận đề bài
for i in range(n):        
    for j in range(n):
        a[i][j]=b[i][j]
print(a)
print("MT nghich dao: \n", np.linalg.inv(b))
print("==================================================================== \n")

#Biến đổi Gauss Jordan
print("2.Bien doi Gauss Jordan tim ma tran nghich dao: \n")
for i in range(n):
    if a[i][i] == 0.0:
        sys.exit('ERROR')
    for j in range(n):
        if i != j:
            ratio = a[j][i]/a[i][i]

            for k in range(2*n):
                a[j][k] = a[j][k] - ratio * a[i][k]
    print(a)
#Cho vế tría thành mt đơn vị
for i in range(n):
    temp = a[i][i]
    for j in range(2*n):
        a[i][j] = a[i][j]/temp
    print(a)
print("Vay ta thu duoc ma tran: \n",a)
print("==================================================================== \n")

#Kết quả
print("3.KET QUA:\n")
print("Chua lam tron")
for i in range(n):
    for j in range(n, 2*n):
        print((a[i][j]), end='\t')
    print()

print("\nLam tron")
for i in range(n):
    for j in range(n, 2*n):
        print(round(a[i][j]), end='\t')
    print()

print("==================================================================== \n")