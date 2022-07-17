import numpy as np
from numpy import poly1d

A = np.loadtxt('B:/BT nhóm GTS/Tìm trị riêng và vector riêng/Danielevsky/matrixA.txt', dtype = float)
print("A1 : \n", A)
n = A.shape[0] 
I = np.eye(n)

M = np.eye(n)
i = 1

def coeff(M):
    return M @ coeff(M)

while A[n-1][n-2] != 0:
    k = A[n-1]
    M[n-2] = k
    print("M{0}: \n {1}".format(i, M))
    invM = np.linalg.inv(M)
    print("invM{0}: \n {1}".format(i, invM))
    A = M @ A @ invM
    print("A{0}: \n {1}".format(i+1, A))
    
    n -= 1
    i += 1

    q = I[n-1]
    M[n-1] = q

    if n == 1 :
        break

P = A 
print("Ma trận dạng Frobenius : \n", P)

equation = [1]
for i in range(P.shape[0]):
    equation.append(-P[0][i])
print(equation)

equation = poly1d(equation)
print("phương trình đặc trưng: \n{0} = 0\n".format(equation))

solve = np.roots(equation)
print("Giá trị riêng : \n ", solve)

M1 = np.loadtxt('B:/BT nhóm GTS/Tìm trị riêng và vector riêng/Danielevsky/matrixM.txt', dtype = float)
M2 = np.loadtxt('B:/BT nhóm GTS/Tìm trị riêng và vector riêng/Danielevsky/matrixM2.txt', dtype = float)
M3 = np.loadtxt('B:/BT nhóm GTS/Tìm trị riêng và vector riêng/Danielevsky/matrixM3.txt', dtype = float)

coeff = np.linalg.inv(M3 @ M2 @ M1)

for i in range(P.shape[0]):
    Y = np.array([solve[i] ** 3, solve[i] ** 2, solve[i], 1]).T

    X = coeff @ Y
    print("vector riêng ứng với trị riêng {0} là : {1}".format(solve[i], X))
