import numpy as np

# nhập ma trận A
A = np.loadtxt('B:/Giai tich so/Danielevsky/matrixA.txt', dtype = float)
print("A : \n", A)

print("Eigenvalues :", np.linalg.eigvals(A))

eigenvalue, eigenvector = np.linalg.eig(A)

vector_riêng = eigenvector.T
# giá trị riêng
print("giá trị riêng : \n ", eigenvalue)

# vector riêng
for i in range(eigenvector.shape[0]):
    print("vector riêng ứng với trị riêng {0} là : {1}".format(eigenvalue[i], vector_riêng[i]))