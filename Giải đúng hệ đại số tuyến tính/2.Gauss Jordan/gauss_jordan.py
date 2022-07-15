import numpy as np
import sys

def gauss_jordan_forward(a):
    # Biến đổi Gauss Jordan
    for i in range(n):
        if a[i][i] == 0.0:
            sys.exit('no solutions')
            
        for j in range(n):
            if i != j:
                ratio = a[j][i] / a[i][i]

                for k in range(n+1):
                    a[j][k] = a[j][k] - ratio * a[i][k]
            print(a)

def solve_x(a):
    # Tính X
    for i in range(n):
        x[i] = a[i][n] / a[i][i]

    # KET QUA
    print("\nKet qua: ",x)
    for i in range(n):
        print('X%d = %0.2f' %(i,x[i]), end = '\t')


a = np.loadtxt('B:/BT nhóm GTS/Giải đúng hệ đại số tuyến tính/2.Gauss Jordan/gauss_jordan.txt', dtype = float)
n = a.shape[0]
x = np.zeros(n,float)
gauss_jordan_forward(a)
solve_x(a)
