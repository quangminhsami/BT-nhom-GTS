from turtle import back
import numpy as np

def forwardElim(a):
    n = a.shape[0]
    for i in range(n):

        for j in range(i+1, n):
            ratio = a[j][i] / a[i][i]
            
            for k in range(n+1):
                a[j][k] = a[j][k] - ratio * a[i][k]

        print(a)
        
    return -1
 
def backSub(a):
    x[n-1] = a[n-1][n] / a[n-1][n-1]
    
    for i in range(n-2,-1,-1):
        x[i] = a[i][n]
        
        for j in range(i+1,n):
            x[i] = x[i] - a[i][j] * x[j]
        
        x[i] = x[i] / a[i][i]
    
    print("\nKet qua:",x)
    for i in range(n):
        print('X%d = %0.7f' %(i,x[i]), end = '\t')

def gaussElimination(a):
    n = a.shape[0]
    singular_flag = forwardElim(a)
    
    if singular_flag != -1:
        print("Singular Matrix \n")

        if a[singular_flag][n] == True:
            print("No solutions")
        else:
            print("Infinitive solutions")   

# main

# Nhập ma trận bổ sung/ liên kết A    
a = np.loadtxt('B:/BT nhóm GTS/Giải đúng hệ đại số tuyến tính/1.Gauss/gauss.txt', dtype = float)

n = a.shape[0]
x = np.zeros(n,float)
forwardElim(a)
backSub(a)
# gaussElimination(a)
