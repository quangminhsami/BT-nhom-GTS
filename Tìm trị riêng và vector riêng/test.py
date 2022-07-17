import numpy as np

# nhập ma trận A
A = np.mat("3 2; 1 0 ")
print("A : \n", A)

print("Eigenvalues :", np.linalg.eigvals(A))

eigenvalue, eigenvector = np.linalg.eig(A)

# 2 giá trị riêng
print("First tuple of eig:", eigenvalue)

# 2 vector riêng
print("Second tuple of eig", eigenvector)