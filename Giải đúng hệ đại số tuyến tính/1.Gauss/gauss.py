import numpy as np
import sys

print("PP Gauss")
print("=========================================")
# Nhập ma trận bổ sung/ liên kết A    

# Th có nghiệm duy nhất
# a = np.loadtxt('B:/BT nhóm GTS/Giải đúng hệ đại số tuyến tính/1.Gauss/gauss.txt', dtype = float)

# th vô số nghiệm
# a = np.loadtxt('B:/BT nhóm GTS/Giải đúng hệ đại số tuyến tính/1.Gauss/gauss_vo_so_nghiem.txt', dtype = float)

# th vô nghiệm
a = np.loadtxt('B:/BT nhóm GTS/Giải đúng hệ đại số tuyến tính/1.Gauss/gauss_vo_nghiem.txt', dtype = float)
rank_a=np.linalg.matrix_rank(a)
            
print("1.INPUT: \n")
print(a)
# TÍNH n a
n,_=np.shape(a)
print("n=",n,"\n")

#Ma trận hệ số B
b = np.zeros((n,n),float)

for i in range(n):        
    for j in range(n):
        b[i][j]=a[i][j]

rank_b=np.linalg.matrix_rank(b)

#Tạo ma trận nghiệm 0
x = np.zeros(n,float)


#Biến đổi Gauss 
def gauss():
    for i in range(n):
        if a[i][i] == 0.0:
            sys.exit('ERROR')
         
        for j in range(i+1, n):
            ratio = a[j][i]/a[i][i]
         
            for k in range(n+1):
             a[j][k] = a[j][k] - ratio * a[i][k]
        print("Lan ",i+1)
        print(a,"\n")
    print("Vay ta thu duoc ma tran \n",a)

    #Tính X
    x[n-1] = a[n-1][n]/a[n-1][n-1]
    for i in range(n-2,-1,-1):
        x[i] = a[i][n]
     
        for j in range(i+1,n):
            x[i] = x[i] - a[i][j]*x[j]
     
        x[i] = x[i]/a[i][i]

    print("\nKet qua:")
    for i in range(n):
        print('X%d = %0.7f' %(i,x[i]), end = '\t')
 
#Ket qua

if rank_a == rank_b:
    if rank_a == n:
        print("2.OUTPUT: Co nghiem duy nhat\n")
    else:
        print("2.OUTPUT: Co vo so nghiem")
else:
    print("2.OUTPUT: Vo nghiem")

gauss()

