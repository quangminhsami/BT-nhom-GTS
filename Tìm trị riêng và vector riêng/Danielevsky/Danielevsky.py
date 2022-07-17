import numpy as np
from numpy import poly1d

A1 = A = np.loadtxt('B:/Giai tich so/Danielevsky/matrixA.txt', dtype = float)
M1 = M = np.loadtxt('B:/Giai tich so/Danielevsky/matrixM.txt', dtype = float)

detM1 = np.linalg.det(M1)
print("det M1 : \n", detM1)

invM1 = np.linalg.inv(M1)
print("inv M1 : \n", invM1)

A2 = M1 @ A1 @ invM1
print("A2 : \n", A2)

M2 = np.loadtxt('B:/Giai tich so/Danielevsky/matrixM2.txt', dtype = float)
detM2 = np.linalg.det(M2)
print("det M2 : \n", detM2)

invM2 = np.linalg.inv(M2)
print("inv M2 : \n", invM2)

A3 = M2 @ A2 @ invM2
print("A3 : \n", A3)

M3 = np.loadtxt('B:/Giai tich so/Danielevsky/matrixM3.txt', dtype = float)
detM3 = np.linalg.det(M3)
print("det M3 : \n", detM3)

invM3 = np.linalg.inv(M3)
print("inv M3 : \n", invM3)

A4 = M3 @ A3 @ invM3
print("A4 : \n", A4)

P = (M3 @ M2 @ M1) @ A @ (np.linalg.inv(M3 @ M2 @ M1))
print("Ma trận dạng Frobenius : \n", P)

equation = [1, - A4[0][0], - A4[0][1],- A4[0][2], - A4[0][3]]

solve = poly1d(equation)
print("phương trình đặc trưng: \n", solve)

solve = np.roots(equation)
print("Giá trị riêng : \n ", solve)

coeff = np.linalg.inv(M3 @ M2 @ M1)
print(coeff)

Y1 = np.array([solve[0] ** 3, solve[0] ** 2, solve[0], 1]).T

X1 = coeff @ Y1
print("vector riêng ứng với trị riêng {0} là : {1}".format(solve[0], X1))
